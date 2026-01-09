using AirfoilTools
using VMSAdjoint
using SegregatedVMSSolver.ParametersDef
using KissSmoothing

AoA = 2.5
R = 0.15
a,b = 0.5, 0.15


xx0 = collect(0.01:0.005:0.99)
xx =  tanh.(((xx0.-0.5))*7)./2 .+0.5
yy = (0.5^2 .-(xx.-0.5).^2).^0.5
@. yy = b * sqrt(1 - ( (xx-0.5) /a)^2)
ap1 = AirfoilPoints([reverse(xx); 0.0], xx[1:end], [reverse(yy);0.0],-yy[1:end])




control_px = collect(LinRange(0.05,0.95,20))
control_points = ControlPoints(control_px,control_px)

RBFfun = RBFFunctionLocalSupport(RBF_CP4, R)


rbfg = RBFGeometry(control_points,RBFfun)
rbfd = RBFDesign(rbfg, ap1)

sprob = StabilizedProblem(VMS(2))
physicalp = PhysicalParameters(Re=1000, u_in=[1.0,0.0])
timep = TimeParameters(dt=0.25, t_endramp=5.0, tF=40.0, time_window=(30.0, 40.0))

msh_size = MeshSize(BL_fl=1e-4, BL_tt=0.01, H = 5.0, Lback = 10)
meshinfo = AirfoilMesh{Unstructured}(AoA= AoA, MS=msh_size)



meshp = MeshParameters((1,1), 2, meshinfo)
exportp = ExportParameters(printinitial=true,printmodel=true,name_tags=["airfoil"], fieldexport=[["uh","ph","friction"]])

solverp = SolverParameters(θ=1.0, M=10, matrix_freq_update=2) #M = 1 export every time-step

simparams = SimulationParameters(timep,physicalp,solverp,exportp)

airfoil_case = Airfoil(meshp,simparams,sprob)

function regf(x0,y0)
    fact= 0.25 #0.05
    y1, _ = denoise(y0; factor=0.05 )
    return y1
end
reg = Regularization(active=true, iter_reg=2, fun=regf)

#Define the Objective Function, first argument [CD,CL], then you can define whatever you want
function J(CDCL)
    CD,CL=CDCL
    CLtarget=0.25
    CDtarget=0.12
    return 8 * (CD-CDtarget)^2 + 0.5 * (CL - CLtarget)^2
end  


function Jpenalty(pp)
    xx, yu,yl = pp

    N = length(xx)
    @assert length(yu) == N
    @assert length(yl) == N

    yc = yu .- yl

    Δx = diff(xx)

    λ1 = 1e-4
    grad_penalty = sum((diff(yc).^2) ./ (Δx .^2))
    # curv_penalty = sum((yc[i+1] - 2*yc[i] + yc[i-1])^2 /  Δx[i]^4 for i in 2:N-1)


    return λ1 * grad_penalty 
end 

TP = ThickPenalty(valid = true, thickness_penalty = Jpenalty)
adj_solver = AdjSolver(δ=0.0001,bounds = DesignBounds(upper=0.25, lower=-0.25),regularization=reg, thick_penalty = TP)

adjoint_airfoil_problem = AdjointProblem( rbfd,airfoil_case,adj_solver,(:unsteady, :steady), J)



solve_adjoint_optimization(adjoint_airfoil_problem)