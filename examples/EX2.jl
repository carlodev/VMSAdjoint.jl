using Revise
using AirfoilTools
using VMSAdjoint
using Gridap, GridapGmsh
using SegregatedVMSSolver.ParametersDef


AoA = 15.0 #2.5
Re = 15_000
NW0 = 15 #number of w0 to parametrize top and bottom

cst0 = CST_NACA0012(;N=NW0, t=0.0)

cstd0 = AirfoilCSTDesign(cst0,200)

sprob = StabilizedProblem(VMS(1))

physicalp = PhysicalParameters(Re=Re, u_in=[1.0,0.0])
timep = TimeParameters(dt=0.05, tF=15.0, time_window=(10.0, 15.0))

meshinfo = AirfoilMesh(AoA= AoA,meshref=2)
meshp = MeshParameters((1,1), 2, meshinfo)
exportp = ExportParameters(printinitial=true,printmodel=true,name_tags=["airfoil"], fieldexport=[["uh","ph","friction"]])

solverp = SolverParameters(θ=1.0)

simparams = SimulationParameters(timep,physicalp,solverp,exportp)

airfoil_case = Airfoil(meshp,simparams,sprob)

#Define the Objective Function, first argument [CD,CL], then you can define whatever you want

function J(CDCL; CLtarget=1.25)
    CD,CL=CDCL
    return 0.5 * (CL - CLtarget)^2
end  

adjoint_airfoil_problem = AdjointProblem( cstd0,airfoil_case,:unsteady,J)


solve_adjoint_optimization(adjoint_airfoil_problem)





