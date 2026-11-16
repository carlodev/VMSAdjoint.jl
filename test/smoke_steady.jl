# Coarse-mesh end-to-end smoke: mesh generation -> steady primal -> objective ->
# steady adjoint -> shape sensitivity. Checks the pipeline runs and yields finite
# numbers; NOT a discretization-accuracy test.

using SegregatedVMSSolver.ParametersDef
using VMSAdjoint.IteratorTools: compute_sensitivity, generate_model, obj_fun, dJobj_fun,
    perturbed_airfoil_model

ap0 = get_airfoil_coordinates(joinpath(@__DIR__, "..", "examples", "n0012.csv"))

control_px = collect(LinRange(0.1, 0.9, 4))
rbfd = RBFDesign(RBFGeometry(ControlPoints(control_px, control_px),
        RBFFunctionLocalSupport(RBF_CP4, 0.25)), ap0)

AoA = 2.5
sprob = StabilizedProblem(VMS(1))
physicalp = PhysicalParameters(Re=1000, u_in=[1.0, 0.0])
timep = TimeParameters(dt=0.05, tF=1.0, time_window=(0.5, 1.0))
meshinfo = AirfoilMesh{Structured}(AoA=AoA, MS=MeshSize(meshref=0.5))
meshp = MeshParameters((1, 1), 2, meshinfo)
exportp = ExportParameters(printinitial=false, printmodel=false,
    name_tags=["airfoil"], fieldexport=[["uh", "ph"]])
solverp = SolverParameters(θ=1.0, M=100, matrix_freq_update=2)
airfoil_case = Airfoil(meshp, SimulationParameters(timep, physicalp, solverp, exportp), sprob)

J(CDCL) = 8 * CDCL[1]^2 + (CDCL[2] - 1.0)^2

workdir = mktempdir()
meshdir = joinpath(workdir, "MeshFiles")

model = generate_model(rbfd, 0, 0.0, meshinfo, physicalp, meshdir)
am = AirfoilModel(model, airfoil_case)

cd(workdir) do
    uh, ph = solve_inc_primal(am, airfoil_case, :steady)
    fval, CDCL = obj_fun(am, airfoil_case, uh, ph, ThickPenalty(), J)
    @test isfinite(fval)
    @test all(isfinite, CDCL)

    adj_bc = -dJobj_fun(J, CDCL)
    uhadj, phadj = solve_inc_adj(am, airfoil_case, adj_bc, :steady, uh, ph)
    @test all(isfinite, uhadj.free_values)

    # shape sensitivity of one suction-side design parameter, via remeshing
    δ = 1e-4
    model1 = generate_model(rbfd, 2, δ, meshinfo, physicalp, meshdir)
    am1 = AirfoilModel(model1, airfoil_case)
    Ju, Jt = compute_sensitivity(am, am1, rbfd, 2, δ, airfoil_case, ThickPenalty(), uh, uhadj)
    @test isfinite(Ju)
    @test isfinite(Jt)

    # ... and via mesh deformation, which must reach the same gradient. The two
    # routes share only the adjoint solution: the mesh-motion route builds δβ
    # from the displacement field it actually applied, the remeshing route
    # differences nodal positions between two independently generated meshes.
    # They are not expected to agree to many digits (the meshes differ), but a
    # sign flip or an order-of-magnitude gap would mean one of them is wrong.
    solver_mm = AdjSolver(δ=δ, mesh_motion=true)
    am1_mm, dh = perturbed_airfoil_model(rbfd, 2, δ, am, airfoil_case, solver_mm, meshdir)
    @test !isnothing(dh)
    Ju_mm, Jt_mm = compute_sensitivity(am, am1_mm, rbfd, 2, δ, airfoil_case,
        ThickPenalty(), uh, uhadj; mesh_velocity=dh)
    @test isfinite(Ju_mm)
    @test isfinite(Jt_mm)
    @test sign(Ju_mm) == sign(Ju)
    @test isapprox(Ju_mm, Ju; rtol=0.5)

    # the deformed mesh is a topological clone of the baseline (this is what
    # makes an FD check against it a clean shape-derivative comparison)
    @test num_cells(am1_mm.model) == num_cells(am.model)
    @test Gridap.Geometry.get_cell_node_ids(Gridap.Geometry.get_grid(am1_mm.model)) ==
          Gridap.Geometry.get_cell_node_ids(Gridap.Geometry.get_grid(am.model))

    # the mesh-motion operator is cached on the baseline model, so a full
    # Ndes sweep pays one assembly + one factorization
    @test haskey(am.params, :mesh_motion)

    # unsteady path: exercises the chunked frozen-operator driver in the real
    # solver (transient multifield spaces, chunk-boundary state handoff) for
    # both the primal and the backward adjoint
    # tF must leave room for the 10*dt adjoint start-up transient the averaging
    # window skips: time_window[2] < tF - 10*dt
    timep_u = TimeParameters(dt=0.05, t0=0.0, tF=1.0, time_window=(0.2, 0.4))
    solverp_u = SolverParameters(θ=0.5, M=100, matrix_freq_update=3)
    case_u = Airfoil(meshp, SimulationParameters(timep_u, physicalp, solverp_u, exportp), sprob)

    am_u = AirfoilModel(generate_model(rbfd, 0, 0.0, meshinfo, physicalp, meshdir), case_u)
    uh_u, ph_u = solve_inc_primal(am_u, case_u, "unsteady-smoke", :unsteady)
    @test all(isfinite, uh_u.free_values)
    # one snapshot per step, plus the initial condition
    @test length(am_u.params[:UH]) == round(Int, (case_u.simulationp.timep.tF -
                                                  case_u.simulationp.timep.t0) /
                                                 case_u.simulationp.timep.dt) + 1

    _, CDCL_u = obj_fun(am_u, case_u, uh_u, ph_u, ThickPenalty(), J)
    uhadj_u, _ = solve_inc_adj(am_u, case_u, -dJobj_fun(J, CDCL_u), :unsteady, uh_u, ph_u)
    @test all(isfinite, uhadj_u.free_values)
    # the unsteady adjoint stores the time-averaged wall-shear correlation that
    # the gradient consumes instead of a product of time averages
    @test haskey(am_u.params, :wall_shear_corr)
end
