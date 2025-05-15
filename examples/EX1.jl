using Revise
using AirfoilTools
using VMSAdjoint
using Gridap, GridapGmsh
using SegregatedVMSSolver.ParametersDef


AoA = 2.5
NW0 = 15 #number of w0 to parametrize top and bottom

cst0 = CST_NACA0012(;N=20, t=0.0)
cstd0 = AirfoilCSTDesign(cst0,200)


sprob = StabilizedProblem(VMS(1))

physicalp = PhysicalParameters(Re=1000, u_in=[1.0,0.0])
timep = TimeParameters(t0=0.0, dt=0.01, tF = 5.0)
meshinfo = AirfoilMesh(AoA= AoA, cstdesign=cstd0)
meshp = MeshParameters((1,1), 2, meshinfo)
exportp = ExportParameters(printinitial=true,printmodel=true,name_tags=["airfoil"], fieldexport=[["uh","ph","friction"]])

solverp = SolverParameters(θ=1.0)

simparams = SimulationParameters(timep,physicalp,solverp,exportp)

airfoil_case = Airfoil(meshp,simparams,sprob)

adjoint_airfoil_problem = AdjointProblem(AdjBC("airfoil", [-1.0,0.0]), airfoil_case)


modelname =create_msh(meshinfo, physicalp ; iter = 0)
model = GmshDiscreteModel(modelname)
writevtk(model, "model_original")


am = AirfoilModel(model, airfoil_case)

as = AdjSolution(adjoint_airfoil_problem,am)


uh,ph = solve_inc_primal(am, airfoil_case, :steady)

# solve_inc_primal(am, airfoil_case, :unsteady; uh0=uh,ph0=ph)



using ForwardDiff

#Define the Objective Function, first argument [CD,CL], then you can define whatever you want
function J(CDCL; CLtarget=0.75)
    CD,CL=CDCL
    return 0.5 * (CL - CLtarget)^2
end  


fval, CDCL = obj_fun(am, adjoint_airfoil_problem, uh,ph, J)
dJobj_fun(J, CDCL)


