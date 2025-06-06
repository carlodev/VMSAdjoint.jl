using Revise
using AirfoilTools
using VMSAdjoint
using Gridap, GridapGmsh
using SegregatedVMSSolver.ParametersDef
"""
Example taken from:
Sorgiovanni, G., Quadrio, M., Ponzini, R., 2016. A robust open-source adjoint optimization method for external aerodynamics. Politecnico di Milano, Milan.

Starting from a NACA0012 finding the derivaties with finite differences
"""

fname = "n0012.csv" #airfoil coordinates to load
AoA = 2.5

ap0 = get_airfoil_coordinates(joinpath(@__DIR__, fname))

control_px = collect(LinRange(0.05,0.95,20))

control_points = ControlPoints(control_px,control_px)


R = 0.15 #support radius. If is bigger, the deformation radius increases, try increasing it
RBFfun = RBFFunctionLocalSupport(RBF_CP4, R)

rbfg = RBFGeometry(control_points,RBFfun)
rbfd = RBFDesign(rbfg, ap0)


sprob = StabilizedProblem(VMS(2))

physicalp = PhysicalParameters(Re=1000, u_in=[1.0,0.0])
timep = TimeParameters(dt=0.05, tF=1.0, time_window=(0.8, 1.0))

meshinfo = AirfoilMesh(AoA= AoA, meshref=2)
meshp = MeshParameters((1,1), 2, meshinfo)
exportp = ExportParameters(printinitial=true,printmodel=true,name_tags=["airfoil"], fieldexport=[["uh","ph","friction"]])

solverp = SolverParameters(θ=1.0)

simparams = SimulationParameters(timep,physicalp,solverp,exportp)

airfoil_case = Airfoil(meshp,simparams,sprob)

#Define the Objective Function, first argument [CD,CL], then you can define whatever you want

function J(CDCL; CLtarget=0.75)
    CD,CL=CDCL
    return -CL #0.5 * (CL - CLtarget)^2
end  



adj_solver = AdjSolver(δ=0.0001)

adjoint_airfoil_problem = AdjointProblem( rbfd,airfoil_case,adj_solver,:unsteady, J)
finite_difference_analysis(adjoint_airfoil_problem)




