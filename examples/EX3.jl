using Revise
using AirfoilTools
using VMSAdjoint
using Gridap, GridapGmsh
using SegregatedVMSSolver.ParametersDef


fname = "n0012.csv" #airfoil coordinates to load
AoA = 2.5

ap0 = get_airfoil_coordinates(joinpath(@__DIR__, fname))

control_px = collect([0.05, 0.1, 0.25, 0.4, 0.55, 0.7, 0.85, 0.92])

control_points = ControlPoints(control_px,control_px)


R = 0.25 #support radius. If is bigger, the deformation radius increases, try increasing it
RBFfun = RBFFunctionLocalSupport(RBF_CP4, R)

rbfg = RBFGeometry(control_points,RBFfun)
rbfd = RBFDesign(rbfg, ap0)


sprob = StabilizedProblem(VMS(2))

physicalp = PhysicalParameters(Re=1000, u_in=[1.0,0.0])
timep = TimeParameters(dt=0.05, tF=10.0)

meshinfo = AirfoilMesh(AoA= AoA, meshref=1)
meshp = MeshParameters((1,1), 2, meshinfo)
exportp = ExportParameters(printinitial=true,printmodel=true,name_tags=["airfoil"], fieldexport=[["uh","ph","friction"]])

solverp = SolverParameters(θ=1.0)

simparams = SimulationParameters(timep,physicalp,solverp,exportp)

airfoil_case = Airfoil(meshp,simparams,sprob)

#Define the Objective Function, first argument [CD,CL], then you can define whatever you want

function J(CDCL; CLtarget=0.35)
    CD,CL=CDCL
    return 0.5 * (CL - CLtarget)^2
end  

adjoint_airfoil_problem = AdjointProblem( rbfd,airfoil_case,:steady, J)
solve_adjoint_optimization(adjoint_airfoil_problem)




