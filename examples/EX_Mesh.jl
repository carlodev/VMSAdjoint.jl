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
# modelname =create_msh(meshinfo,rbfd, physicalp ; iter = 0)
# model = GmshDiscreteModel(modelname)
# writevtk(model, "model_0n")

w1 = [0.037170259602030145, 0.04300665266635907, 0.05308859181629484, 0.05627205461260461, 0.056742942526510305, 0.05661162468594478, 0.05557170742655077, 0.053792533523437804, 0.051424724170843004, 0.04857465955769425, 0.04683297123598739, 0.04333555659019143, 0.039505887781925146, 0.03541140840420925, 0.029239701060206503, 0.02505659812017946, 0.02099978380709525, 0.01525513060884453, 0.010254059336041524, 0.003911595940572393, -0.005671327805635765, -0.024784692307439866, -0.029007725436627767, -0.03461022305930837, -0.03741093243149193, -0.039543405092036674, -0.04026019915947145, -0.040115948579654864, -0.03926776401278267, -0.03779087194103642, -0.03780489171003431, -0.035835110298814705, -0.03376356382144574, -0.031666661168062354, -0.02720091543454356, -0.02568116265331385, -0.024887537268359402, -0.02592145877139583, -0.05159846443455528, -0.06199643509034411]
rbfd01 = create_AirfoilDesign(rbfd,w1)
for  i = 1:1:40
rbfd1 = perturb_DesignParameter(rbfd01, i, 0.0001)

modelname =create_msh(meshinfo,rbfd1, physicalp ; iter = 1)

model = GmshDiscreteModel(modelname)
writevtk(model, "model_1n")
end





plot(w2)

rbfd2 = create_AirfoilDesign(rbfd,w2)
modelname =create_msh(meshinfo,rbfd2, physicalp ; iter = 1)
model = GmshDiscreteModel(modelname)
writevtk(model, "model_2n")


Γ = BoundaryTriangulation(model; tags="airfoil")
dΓ = Measure(Γ, 2*2)
nΓ = -get_normal_vector(Γ) #beacuse they point inward 
writevtk(Γ, "nΓ", cellfields=["nΓ"=>nΓ])


