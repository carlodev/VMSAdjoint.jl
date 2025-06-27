using AirfoilTools
using VMSAdjoint
using SegregatedVMSSolver.ParametersDef


AoA = 0.0

R = 0.15

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


yy = circle(xx)

ap1 = AirfoilPoints([reverse(xx); 0.0], xx[1:end], [reverse(yy);0.0],-yy[1:end])

control_px = collect(LinRange(0.05,0.95,20))
control_points = ControlPoints(control_px,control_px)
RBFfun = RBFFunctionLocalSupport(RBF_CP4, R)

rbfg = RBFGeometry(control_points,RBFfun)
rbfd = RBFDesign(rbfg, ap1)


sprob = StabilizedProblem(VMS(2))
physicalp = PhysicalParameters(Re=1000, u_in=[1.0,0.0])
timep = TimeParameters(dt=0.10, tF=50.0, t_endramp=1.0, time_window=(20.0, 50.0))

meshinfo = AirfoilMesh(AoA= AoA, meshref=1, H = 5.0, Lback = 10)
meshp = MeshParameters((1,1), 2, meshinfo)
exportp = ExportParameters(printinitial=true,printmodel=true,name_tags=["airfoil"], fieldexport=[["uh","ph","friction"]])

solverp = SolverParameters(θ=1.0, M=1, matrix_freq_update=2) #M = 1 export every time-step

simparams = SimulationParameters(timep,physicalp,solverp,exportp)

airfoil_case = Airfoil(meshp,simparams,sprob)



#Define the Objective Function, first argument [CD,CL], then you can define whatever you want
function J(CDCL; CLtarget=0.75)
    CD,CL=CDCL
    return CD #0.5 * (CL - CLtarget)^2
end  


adj_solver = AdjSolver(δ=0.001, αg = 0.15, scaled=true, bounds = DesignBounds(upper=1.0, lower=-1.0))

adjoint_airfoil_problem = AdjointProblem( rbfd,airfoil_case,adj_solver,(:unsteady,:unsteady), J)
# solve_adjoint_optimization(adjoint_airfoil_problem)

using JLD2, Statistics

am = jldopen("UnsteadyPrimalFields.jld2")["am"]

adjf = jldopen("UnsteadyAdjointFields.jld2")

am.params[:uh].free_values .= mean(am.params[:UH][350:500])

am.params[:ph].free_values .= mean(am.params[:PH][350:500])


adj_bc = [-1.0, 0.0]

uh,ph = am.params[:uh], am.params[:ph]

uhadj,phadj = solve_inc_adj(am, airfoil_case, adj_bc, "inc-adj-steady-0", :steady, am.params[:uh],am.params[:ph])

uhadj.free_values .= mean(adjf["UH_ADJ"][5:150])
phadj.free_values .= mean(adjf["PH_ADJ"][5:150])
# uhadj,phadj = solve_inc_adj(am, airfoil_case, adj_bc, "inc-adj-unsteady-0", :unsteady, 0.0,0.0)

order,D = 2, 2


function rotation(n::VectorValue{2,Float64})
    n1, n2 = [n...] ./ norm(n)
    VectorValue(-n2, n1)
end

δv = 0.0001 .* [ones(20);-ones(20)]
Ω = Triangulation(am.model)
dΩ = Measure(Ω, order*2)
Γ = BoundaryTriangulation(am.model; tags="airfoil")
dΓ = Measure(Γ, order*2)
nΓ = -get_normal_vector(Γ) #beacuse they point inward 
tΓ = rotation ∘ nΓ

ν = 1/1000
reffe = ReferenceFE(lagrangian, VectorValue{D,Float64}, order)
VV0 = FESpace(am.model, reffe; conformity=:H1) 
function fx(x)
    return x
end

m0 = get_free_dof_values(interpolate_everywhere(fx, VV0))#space where to interpolate the variation


Jnew = zeros(40)
Jbound = zeros(40)

for (IIDX, δ) in enumerate(δv)
    
    
    rbfd1 = perturb_DesignParameter(rbfd, [IIDX], [δ])
    modelname_tmp = create_msh(meshinfo, rbfd1, physicalp, "MeshPerturb"; iter=IIDX + 100)
    model_tmp = GmshDiscreteModel(modelname_tmp)
    am_tmp = AirfoilModel(model_tmp, airfoil_case; am=am)


    VV1 = FESpace(am_tmp.model, reffe; conformity=:H1)
    m1 = get_free_dof_values(interpolate_everywhere(fx, VV1))

    vi = (m1 - m0) ./ δ
    v_field = FEFunction(VV0, vi)

    #You want to visualize it?
    #writevtk(Ω, "vfield", nsubcells=order, cellfields=["vfield"=>v_field])


    #-----------------------------------------------------
    #Now computing the Sensitivity
    #-----------------------------------------------------

    Jbound[IIDX] =     sum(-ν* ∫( ((∇(uh) ⋅ tΓ) ⋅ nΓ) * ((∇(uhadj) ⋅ tΓ) ⋅ nΓ) ⋅ (v_field ⋅nΓ) )dΓ)
    Tj_visc =  -ν * ((∇(uh) ⋅ tΓ) ⋅ nΓ) * ((∇(uhadj) ⋅ tΓ) ⋅ nΓ) #* nΓ -ν*(∇(uh)⊙∇(uhadj)) ⋅ nΓ
    Tj =  -ν*(∇(uh)⊙∇(uhadj)) ⋅ nΓ - ph * (uhadj⋅ nΓ) + phadj * (uh ⋅ nΓ) 
    Jnew[IIDX] = sum(∫(Tj⋅v_field)dΓ )
end


using Plots, LaTeXStrings

grad_c = -2 .* [Jbound[1:20]; -Jbound[21:end]]
grad_n = -2 .* [Jnew[1:20]; -Jnew[21:end]]


# grad_c_st = copy(grad_c)

Plots.default(linewidth=2)

plot(grad_c_st[1:20], label = "Steady Adjoint", linecolor=:black)
plot!(grad_c[1:20], label = "Unsteady Time Averaged Adjoint", linecolor=:blue)
# plot!(grad_n, label = "New Adjoint", linecolor=:blue)
plot!(xlabel="Design variable β", ylabel = L"C_D"* " gradient (outward deformation >0)")



