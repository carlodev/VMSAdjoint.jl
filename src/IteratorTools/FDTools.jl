"""
    finite_difference_analysis(adjoint_airfoil_problem::AdjointProblem)
Iterator for Finite Difference Analysis
"""
function finite_difference_analysis(adjoint_airfoil_problem::AdjointProblem; idxs::Vector{Int64}=Int64[])
    @unpack timesol, adesign, J, vbcase, solver= adjoint_airfoil_problem
    @unpack thick_penalty = solver
    @sunpack order = vbcase

    tsol = timesol[1]

    fddir = "FD"
    mkpath(fddir)

    airfoil_case = vbcase
    w = get_DesignParameters(adesign)

    meshinfo = airfoil_case.meshp.meshinfo
    physicalp =airfoil_case.simulationp.physicalp

    Ndes = length(w) #number of design parameters

    if isempty(idxs)
        idxs = collect(1:1:Ndes)
    else
        @assert maximum(idxs) <= Ndes "The IDX of Finite Difference has to be in the range of tha maximum number of design parameters; maximum(idxs) = $(maximum(idxs)), Number of design parameters: $Ndes"
    end

    
    δ = solver.δ #0.001
    shift = CSTweights(Int(Ndes/2), δ)
    shiftv =   vcat(shift)
    
    #create the new airfoil model from the weights w
    adesign = create_AirfoilDesign(adesign,w)

    model = generate_model(adesign, 0, 0.0, meshinfo, physicalp, "MeshFiles")

     #modelname = create_msh(meshinfo,adesign, physicalp ; iter = 0)
     #model = GmshDiscreteModel(modelname)
     
    writevtk(model, joinpath(fddir,"model_FD_0"))
    am =  AirfoilModel(model, airfoil_case)
    filename = joinpath(fddir, "FD_0")

   

    uh0,ph0 = solve_inc_primal(am, airfoil_case, filename, tsol; linear_backend=solver.linear_backend)
    fval0, CLCD0 = obj_fun(am, airfoil_case, uh0,ph0,thick_penalty, J)
    fval_fd, CLCD_fd = iterate_fd(shiftv,idxs, adesign, am, airfoil_case,tsol,thick_penalty, J, solver)

    fval_grad = (fval_fd .- fval0)./shiftv
    CLCD_grad= map(CLCDi-> CLCDi .- CLCD0, CLCD_fd) ./shiftv
    
    sol_fd = Dict(:CLCD_grad=>CLCD_grad, :fval_grad=>fval_grad, :idxs=>idxs)
    jldsave( joinpath(fddir, "FD_SOL.jld2"); sol_fd)


    return fval_grad, CLCD_grad

end

"""
    iterate_fd(shift, idxs, adesign, am0, airfoil_case, tsol, thick_penalty, J, solver)

Re-solve the primal flow on each perturbed geometry, for the finite-difference
gradient.

The perturbed geometries come from [`perturbed_airfoil_model`](@ref), i.e. from
the *same* code path the adjoint sensitivity uses. With
`AdjSolver.mesh_motion = true` that means both sides of an FD-vs-adjoint
comparison see meshes with identical topology, differing only by the node
displacement being differentiated. Remeshing each perturbation instead
(`mesh_motion = false`) injects a mesh-dependent difference into `ΔJ` that is
unrelated to the shape derivative and that, at the small `δ` an FD check needs,
can be the same order as the signal.
"""
function iterate_fd(shift::Vector{Float64}, idxs::Vector{Int64}, adesign::AirfoilDesign, am0::AirfoilModel, airfoil_case::Airfoil,tsol::Symbol,thick_penalty, J::Function, solver::AdjSolver)
    fddir = "FD"
    mkpath(fddir)

    # sized by the TOTAL number of design parameters: the loop indexes by the
    # design index i ∈ idxs, which may be a non-contiguous subset
    Ndes = length(shift)
    fval_fd = zeros(Ndes)
    CLCD_fd = [zeros(2) for _ in 1:Ndes]


    for i in idxs

        ss = shift[i]

        println("Perturbation Domain $i")

        am_tmp, _ = perturbed_airfoil_model(adesign, i, ss, am0, airfoil_case, solver, "MeshPerturb")

        filename = joinpath(fddir, "FD_$(i)")


        uh_tmp,ph_tmp = solve_inc_primal(am_tmp, airfoil_case, filename, tsol; linear_backend=solver.linear_backend)

        fval_fd[i], CLCD_fd[i] = obj_fun(am_tmp, airfoil_case, uh_tmp,ph_tmp,thick_penalty, J)

    end

    return fval_fd, CLCD_fd
end

