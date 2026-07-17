using AirfoilTools
using VMSAdjoint
using Gridap, GridapGmsh
using SegregatedVMSSolver.ParametersDef

# =============================================================================
# Adjoint shape sensitivity of a NACA0012 (AoA = 2.5 deg, Re = 1000)
# computed on TWO mesh configurations:
#   1. structured C-type mesh
#   2. unstructured mesh with boundary layer
#
# Element types: structured meshes may be :TRI or :QUAD (transfinite, optionally
# recombined); unstructured meshes are :TRI only (GridapGmsh requires a single
# 2D cell type). CST designs additionally require a structured mesh.
# =============================================================================

# ---------------------------------------------------------------------------
# Geometry & design parametrization (RBF control points on both surfaces)
# ---------------------------------------------------------------------------
fname = "n0012.csv"
AoA   = 2.5
ap0   = get_airfoil_coordinates(joinpath(@__DIR__, fname))

ap0.yu[5] = ap0.yu[5] * 1

control_px     = collect(LinRange(0.05, 0.95, 20))     # 20 + 20 = 40 control points
control_points = ControlPoints(control_px, control_px)

R      = 0.15                                           # RBF support radius
RBFfun = RBFFunctionLocalSupport(RBF_CP4, R)
rbfg   = RBFGeometry(control_points, RBFfun)
rbfd   = RBFDesign(rbfg, ap0)

# ---------------------------------------------------------------------------
# Physical / numerical parameters (shared by all four cases)
# ---------------------------------------------------------------------------
sprob     = StabilizedProblem(VMS(2))
physicalp = PhysicalParameters(Re = 1000, u_in = [1.0, 0.0])
timep     = TimeParameters(dt = 0.05, tF = 10.0, time_window = (8.0, 10.0))
solverp   = SolverParameters(θ = 1.0, M = 1, matrix_freq_update = 4)
exportp   = ExportParameters(printinitial = true, printmodel = true,
                             name_tags = ["airfoil"], fieldexport = [["uh", "ph", "friction"]])

# Solve mode: (primal, adjoint).
#   (:steady,   :steady)   -> fast; exact for this attached, steady flow.
#   (:unsteady, :unsteady) -> time-accurate gradient. The unsteady adjoint builds the
#                             time-AVERAGED wall-shear correlation ⟨(∂v_t/∂n)(∂ψ_t/∂n)⟩,
#                             keeping the primal–adjoint covariance term (Srinath &
#                             Mittal, JCP 2010) that a product of time-averages drops.
timesol = (:steady, :steady)

# Objective: minimize drag. Adjoint BC on the airfoil = -dJ/d[CD,CL].
# (Try `8*(CD)^2 + (CL - 1.0)^2` for a lift-target inverse design instead.)
J(CDCL) = CDCL[1]

thick_penalty = ThickPenalty()   # inactive by default (valid=false) -> no penalty contribution
δ             = 1e-4             # design-parameter perturbation


meshinfo = AirfoilMesh{Unstructured}(elements=:TRI, AoA = AoA, MS = MeshSize(BL_fl = 1e-4, BL_tt = 0.01, meshref = 1.0))
tag = "tri_unstructured_BL"

# ---------------------------------------------------------------------------
# Full adjoint-sensitivity pass for one mesh configuration
# ---------------------------------------------------------------------------



meshp        = MeshParameters((1, 1), 2, meshinfo)
simparams    = SimulationParameters(timep, physicalp, solverp, exportp)
airfoil_case = Airfoil(meshp, simparams, sprob)

# baseline mesh + model
model = generate_model(rbfd, 0, 0.0, meshinfo, physicalp, "MeshFiles_$tag")
am    = AirfoilModel(model, airfoil_case)




# primal
uh, ph = solve_inc_primal(am, airfoil_case, "primal_$tag", timesol[1])

    # objective value + adjoint boundary condition
    fval, CDCL = obj_fun(am, airfoil_case, uh, ph, thick_penalty, J)
    adj_bc     = -dJobj_fun(J, CDCL)

    # adjoint (for :unsteady this also stores the wall-shear correlation in am.params)
    uhadj, phadj = solve_inc_adj(am, airfoil_case, adj_bc, "adj_$tag", timesol[2], uh, ph)

    # sensitivity for every design parameter
    Ndes   = length(get_DesignParameters(rbfd))
    shiftv = vcat(CSTweights(Int(Ndes ÷ 2), δ))   # [+δ … , -δ …] -> outward on both surfaces
    grad   = zeros(Ndes)
    for (i, ss) in enumerate(shiftv)
        model_tmp = generate_model(rbfd, i, ss, meshinfo, physicalp, "MeshPerturb_$tag")
        am_tmp    = AirfoilModel(model_tmp, airfoil_case)
        Jsens, _  = compute_sensitivity(am, am_tmp, rbfd, i, ss, airfoil_case, thick_penalty, uh, uhadj)
        grad[i]   = Jsens
    end

using Plots

shift = CSTweights(Int(Ndes/2), δ)
δv=   vcat(shift)
fd = ([ 0.1216404544975305, 0.12164022538139854, 0.121639265678663, 0.12163864470430423, 0.12163829871723085, 0.12163808109652317, 0.12163789670091155, 0.12163769948157467, 0.12163746739003616, 0.1216371925430035, 0.12163687744654388, 0.12163653066004938, 0.12163616545553438, 0.12163579709642698, 0.12163544348892595, 0.12163511230041991, 0.1216348114898485, 0.1216345602199478, 0.12163435681023815, 0.12163427508648679, 0.12163706730354933, 0.12163706106603764, 0.12163663022813836, 0.1216362323852443, 0.12163596143757201, 0.12163582081275652, 0.12163578439796892, 0.12163582287820633, 0.12163591042708378, 0.12163602986232984, 0.12163617245015604, 0.12163633750628997, 0.12163652509682873, 0.12163673532585907, 0.12163697580637903, 0.12163722805357227, 0.12163748254041373, 0.12163779866991345, 0.12163817398616644, 0.12163895833965868] .- 0.12163386560377128) ./ δv

plot(grad)
scatter!(fd)