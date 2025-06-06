using Revise
using AirfoilTools
using VMSAdjoint
using Gridap, GridapGmsh
using SegregatedVMSSolver.ParametersDef


AoA = 0.0
meshinfo = AirfoilMesh(AoA=AoA, meshref=2)
physicalp = PhysicalParameters(Re=1000, u_in=[1.0,0.0])

xx = collect(0.0:0.01:1.0)
yy = (0.5^2 .-(xx.-0.5).^2).^0.5
ap1 = AirfoilPoints(reverse(xx),xx,reverse(yy),-yy)

control_px = collect(LinRange(0.05,0.95,20))

control_points = ControlPoints(control_px,control_px)

R = 0.15
RBFfun = RBFFunctionLocalSupport(RBF_CP4, R)


rbfg = RBFGeometry(control_points,RBFfun)
rbfd = RBFDesign(rbfg, ap1)


sprob = StabilizedProblem(VMS(2))

physicalp = PhysicalParameters(Re=1000, u_in=[1.0,0.0])
timep = TimeParameters(dt=0.25, tF=10.0, time_window=(8.0, 10.0))

meshinfo = AirfoilMesh(AoA= AoA, meshref=2)
meshp = MeshParameters((1,1), 2, meshinfo)
exportp = ExportParameters(printinitial=true,printmodel=true,name_tags=["airfoil"], fieldexport=[["uh","ph","friction"]])

solverp = SolverParameters(θ=1.0)

simparams = SimulationParameters(timep,physicalp,solverp,exportp)

airfoil_case = Airfoil(meshp,simparams,sprob)

#Define the Objective Function, first argument [CD,CL], then you can define whatever you want

function J(CDCL; CLtarget=0.75)
    CD,CL=CDCL
    return CD #0.5 * (CL - CLtarget)^2
end  



adj_solver = AdjSolver(δ=0.0001, αg = 2.0)

adjoint_airfoil_problem = AdjointProblem( rbfd,airfoil_case,adj_solver,:unsteady, J)
solve_adjoint_optimization(adjoint_airfoil_problem)