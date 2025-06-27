using Revise
using AirfoilTools
using VMSAdjoint
using Gridap, GridapGmsh
using SegregatedVMSSolver.ParametersDef

using Plots

"""
Example taken from:
Sorgiovanni, G., Quadrio, M., Ponzini, R., 2016. A robust open-source adjoint optimization method for external aerodynamics. Politecnico di Milano, Milan.

Starting from a NACA0012 finding the derivaties with the adjoint method
"""

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
timep = TimeParameters(dt=0.05, tF=1.0, time_window=(0.8, 1.0)) #the time-window define the time-span for time-averaging


meshinfo = AirfoilMesh(AoA= AoA, meshref=1)
meshp = MeshParameters((1,1), 2, meshinfo)
exportp = ExportParameters(printinitial=true,printmodel=true,name_tags=["airfoil"], fieldexport=[["uh","ph","friction"]])

solverp = SolverParameters(θ=1.0,  M=1, matrix_freq_update=4) #keep θ=1.0 for stability

simparams = SimulationParameters(timep,physicalp,solverp,exportp)

airfoil_case = Airfoil(meshp,simparams,sprob)



modelname =create_msh(meshinfo,rbfd, physicalp ; iter = 0)
model = GmshDiscreteModel(modelname)

using JLD2, Plots, VMSAdjoint
sol = jldopen(joinpath(@__DIR__,  "ADJ_SOL2.jld2"))["adj_sol"]

rbfd2= create_AirfoilDesign(rbfd, sol.βv)
modelname2 =create_msh(meshinfo,rbfd2, physicalp,"MeshPerturb"; iter = 2+100)
model2 = GmshDiscreteModel(modelname2)
am2 =  AirfoilModel(model2, airfoil_case)

get_aerodynamic_features(am2, uh,ph;tag="airfoil")


plot(am2.ap.xu, am2.ap.yu, aspect_ratio=:equal)
plot!(am2.ap.xl, am2.ap.yl)