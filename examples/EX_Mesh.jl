using Revise
using AirfoilTools
using VMSAdjoint
using Gridap, GridapGmsh
using SegregatedVMSSolver.ParametersDef


fname = "n0012.csv" #airfoil coordinates to load
AoA = 2.5

ap0 = get_airfoil_coordinates(joinpath(@__DIR__, fname))

#40 control points in total, 20 on the suction side, 20 on the pressure side
control_px = collect(LinRange(0.05,0.95,20))
control_points = ControlPoints(control_px,control_px)

#support radius:
R = 0.15 
RBFfun = RBFFunctionLocalSupport(RBF_CP4, R)

rbfg = RBFGeometry(control_points,RBFfun)
rbfd = RBFDesign(rbfg, ap0)


sprob = StabilizedProblem(VMS(2))

physicalp = PhysicalParameters(Re=1000, u_in=[1.0,0.0])
timep = TimeParameters(dt=0.05, tF=10.0, time_window=(8.0, 10.0)) #the time-window define the time-span for time-averaging

msh_size = MeshSize(BL_fl=1e-4,BL_tt=0.01 )

meshinfo = AirfoilMesh(AoA= AoA, MS=msh_size)
meshp = MeshParameters((1,1), 2, meshinfo)
exportp = ExportParameters(printinitial=true,printmodel=true,name_tags=["airfoil"], fieldexport=[["uh","ph","friction"]])

solverp = SolverParameters(θ=1.0,  M=1, matrix_freq_update=4) #keep θ=1.0 for stability

simparams = SimulationParameters(timep,physicalp,solverp,exportp)

airfoil_case = Airfoil(meshp,simparams,sprob)


#Define the Objective Function, it only takes one argument [CD,CL], then you can define all the keywords you want
#The boundary conditions are defined as -dJ/dCDCL
function J(CDCL; CLtarget=0.75)
    CD,CL=CDCL
    return CD/CL  #0.5 * (CL - CLtarget)^2
end  


#Perturbation of your design parameters. The one on the pressure side are pertubed by -δ; so the deformation is always outward
adj_solver = AdjSolver(δ=0.0001)

#here you can choose :steady or :unsteady for the resolution of the primal flow. The unsteady solution is always initialized from a steady one.
adjoint_airfoil_problem = AdjointProblem( rbfd,airfoil_case,adj_solver,(:unsteady,:steady), J )

#solve_adjoint_optimization(adjoint_airfoil_problem)
using Gridap, VMSAdjoint
modelname =create_msh(meshinfo,rbfd, physicalp ; iter = 0)
model = GmshDiscreteModel(modelname)
writevtk(model, "model_0n")



Γ = BoundaryTriangulation(model; tags="airfoil")
dΓ = Measure(Γ, 2*2)
nΓ = -get_normal_vector(Γ) #beacuse they point inward 

writevtk(Γ, "nΓ", cellfields=["nΓ"=>nΓ])


