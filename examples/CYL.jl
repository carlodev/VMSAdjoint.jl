using AirfoilTools
using VMSAdjoint
using SegregatedVMSSolver.ParametersDef


AoA = 2.5

a,b = 0.5, 0.25
NW0 = 15 #number of w0 to parametrize top and bottom
N1 = 0.5
N2 = 1.0

cst0 = CSTGeometry(cstw=CSTweights(NW0,0.1), dz=0.0, N1=N1,N2=N2)

using Plots
function refined_symmetric_distribution(n::Int)
    # Generate half of the points with clustering near 0 using cosine transform
    θ = range(0, π/2; length=n)
    x_left = 0.5 * (1 .- cos.(θ))  # from 0 to 0.5, refined at 0

    # Mirror around x = 0.5
    x_right = 1 .- x_left[end:-1:1]

    # Combine (excluding the midpoint duplicate if n is odd)
    return vcat(x_left, x_right[n % 2 == 0 ? 1 : 2:end])
end

xx = refined_symmetric_distribution(51)

# xx = collect(0.01:0.005:0.99)

# yy = (0.5^2 .-(xx.-0.5).^2).^0.5
yy = circle(xx)
# @. yy = b * sqrt(1 - ( (xx-0.5) /a)^2)

ap0 = AirfoilPoints([reverse(xx); 0.0], xx[1:end], [reverse(yy);0.0],-yy[1:end])
cst01 = CSTGeometry(ap0, cst0)
cstd0 = AirfoilCSTDesign(cst01, ap0)


sprob = StabilizedProblem(VMS(2))
physicalp = PhysicalParameters(Re=1000, u_in=[1.0,0.0])
timep = TimeParameters(dt=0.25, tF=30.0, time_window=(20.0, 30.0))

meshinfo = AirfoilMesh(AoA= AoA, meshref=2, H = 5.0, Lback = 10)
meshp = MeshParameters((1,1), 2, meshinfo)
exportp = ExportParameters(printinitial=true,printmodel=true,name_tags=["airfoil"], fieldexport=[["uh","ph","friction"]])

solverp = SolverParameters(θ=1.0, M=10, matrix_freq_update=4) #M = 1 export every time-step

simparams = SimulationParameters(timep,physicalp,solverp,exportp)

airfoil_case = Airfoil(meshp,simparams,sprob)



#Define the Objective Function, first argument [CD,CL], then you can define whatever you want
function J(CDCL; CLtarget=0.75)
    CD,CL=CDCL
    return CD/CL #0.5 * (CL - CLtarget)^2
end  

function Jf(CDCL; CLtarget=0.75)
    CD,CL=CDCL
    return 2.0

end  


adj_solver = AdjSolver(δ=0.01, αg = 0.15, scaled=true, bounds = DesignBounds(upper=1.0, lower=-1.0))

adjoint_airfoil_problem = AdjointProblem( cstd0,airfoil_case,adj_solver,:unsteady, (J,Jf))
# solve_adjoint_optimization(adjoint_airfoil_problem)


using Gridap, GridapGmsh
modelname =create_msh(meshinfo,cstd0, physicalp )
model = GmshDiscreteModel(modelname)

writevtk(model, "CylinderCST")

w0 = vcat(cstd0.cstg.cstw)
w1 = [0.9641190385287466, 0.9999999999999986, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.9999999999999987, -0.964119037985676, -0.9999999999999983, -1.0, -0.9999999999999993, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -0.9999999999999999, -1.0, -1.0]

adesign = create_AirfoilDesign(cstd0, w0)

norm(cstd0.ap.xu - adesign.ap.xu)
norm(cstd0.ap.xl - adesign.ap.xl)
norm(cstd0.ap.yu - adesign.ap.yu)
norm(cstd0.ap.yl - adesign.ap.yl)


modelname =create_msh(meshinfo,adesign, physicalp )
model = GmshDiscreteModel(modelname)
writevtk(model, "CylinderCST_T")

adesignp = perturb_DesignParameter(adesign, 2, 0.01 )

modelname =create_msh(meshinfo,adesignp, physicalp )
model = GmshDiscreteModel(modelname)
writevtk(model, "CylinderCST_Tp")