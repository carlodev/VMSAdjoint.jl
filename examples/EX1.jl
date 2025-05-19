using Revise
using AirfoilTools
using VMSAdjoint
using Gridap, GridapGmsh
using SegregatedVMSSolver.ParametersDef


AoA = 2.5
NW0 = 15 #number of w0 to parametrize top and bottom

cst0 = CST_NACA0012(;N=20, t=0.0)

cstd0 = AirfoilCSTDesign(cst0,100)

sprob = StabilizedProblem(VMS(1))

physicalp = PhysicalParameters(Re=1000, u_in=[1.0,0.0])
timep = TimeParameters(dt=0.05, tF=1.0)

meshinfo = AirfoilMesh(AoA= AoA)
meshp = MeshParameters((1,1), 2, meshinfo)
exportp = ExportParameters(printinitial=true,printmodel=true,name_tags=["airfoil"], fieldexport=[["uh","ph","friction"]])

solverp = SolverParameters()

simparams = SimulationParameters(timep,physicalp,solverp,exportp)

airfoil_case = Airfoil(meshp,simparams,sprob)

#Define the Objective Function, first argument [CD,CL], then you can define whatever you want

function J(CDCL; CLtarget=0.35)
    CD,CL=CDCL
    return 0.5 * (CL - CLtarget)^2
end  

adjoint_airfoil_problem = AdjointProblem( cstd0,airfoil_case,J)


solve_adjoint_optimization(adjoint_airfoil_problem)




modelname =create_msh(meshinfo,cstd0, physicalp ; iter = 0)
model = GmshDiscreteModel(modelname)
writevtk(model, "model_original")


am = AirfoilModel(model, airfoil_case)



uh,ph = solve_inc_primal(am, airfoil_case, :steady)
 
PressureCoefficient = get_aerodynamic_features(am,uh,ph)


fval, CDCL = obj_fun(am, airfoil_case, uh,ph, J)
d_bc = dJobj_fun(J, CDCL)
uhadj,phadj = solve_inc_adj(am, airfoil_case, d_bc, :steady, uh, ph)







