using Revise
using AirfoilTools
using VMSAdjoint
using Gridap, GridapGmsh
using SegregatedVMSSolver.ParametersDef

"""
    Adjoint 2D optimization starting from a NACA0012
    Primal and Adjoint solved as steady
"""

AoA = 2.5
NW0 = 20 #number of w0 to parametrize top and bottom

cst0 = CST_NACA0012(;N=NW0, t=0.0)

cstd0 = AirfoilCSTDesign(cst0,100)

sprob = StabilizedProblem(VMS(1))

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

adjoint_airfoil_problem = AdjointProblem( cstd0,airfoil_case,:steady, J)
solve_adjoint_optimization(adjoint_airfoil_problem)




