using Revise
using AirfoilTools
using VMSAdjoint
using Gridap, GridapGmsh
using SegregatedVMSSolver.ParametersDef

"""
Example taken from:
Sorgiovanni, G., Quadrio, M., Ponzini, R., 2016. A robust open-source adjoint optimization method for external aerodynamics. Politecnico di Milano, Milan.

Starting from a NACA0012 finding the derivaties with the adjoint method
"""



fname = "../n0012.csv" #airfoil coordinates to load
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
timep = TimeParameters(dt=0.05, tF=10.0, time_window=(8.0, 10.0))


msh_size = MeshSize(BL_fl=1e-4,BL_tt=0.01, meshref=2  )

meshinfo = AirfoilMesh{Unstructured}(AoA= AoA, MS=msh_size)

meshp = MeshParameters((1,1), 2, meshinfo)
exportp = ExportParameters(printinitial=true,printmodel=true,name_tags=["airfoil"], fieldexport=[["uh","ph","friction"]])

solverp = SolverParameters(θ=1.0,  M=1, matrix_freq_update=4) #keep θ=1.0 for stability

simparams = SimulationParameters(timep,physicalp,solverp,exportp)

airfoil_case = Airfoil(meshp,simparams,sprob)



#Define the Objective Function, it only takes one argument [CD,CL], then you can define all the keywords you want
#The boundary conditions are defined as -dJ/dCDCL
function J(CDCL)
    CD,CL=CDCL
    CLtarget = 1.0
    CDtarget=0.0
    return 8* (CD-CDtarget)^2 + (CL - CLtarget)^2
end 



#Perturbation of your design parameters. The one on the pressure side are pertubed by -δ; so the deformation is always outward
adj_solver = AdjSolver(δ=0.0001)
#here you can choose :steady or :unsteady for the resolution of the primal flow. The unsteady solution is always initialized from a steady one.
adjoint_airfoil_problem = AdjointProblem( rbfd,airfoil_case,adj_solver,(:unsteady,:steady), J )



thick_penalty = ThickPenalty()

modelname = create_msh(meshinfo, rbfd, physicalp, "MeshPerturb")
model = GmshDiscreteModel(modelname)
am0 = AirfoilModel(model, airfoil_case)


δ = adj_solver.δ
ss = δ*1

model1 = generate_regularized_model(rbfd, 2, ss, meshinfo, physicalp, "MeshPerturb")
am1 =  AirfoilModel(model1, airfoil_case)

using VMSAdjoint.IteratorTools
using VMSAdjoint.IteratorTools: compute_∇penalty
shift = CSTweights(20, δ)
shiftv =   vcat(shift) #[δ,δ,δ,δ,δ...., -δ,-δ,-δ,-δ,.....]

Jthickness = zeros(40)
for (i,ss) in enumerate(shiftv)
    @info "Perturbation Domain $i"
    model_tmp = generate_regularized_model(rbfd, i, ss, meshinfo, physicalp, "MeshPerturb")
    am_tmp =  AirfoilModel(model_tmp, airfoil_case)
    Jthickness[i] = compute_∇penalty(am0,am_tmp,  δ, thick_penalty) #-> fix bug here
    println(Jthickness[i])
end


[println(t) for t in Jthickness]