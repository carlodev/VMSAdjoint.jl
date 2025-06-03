using Revise
using AirfoilTools
using VMSAdjoint
using Gridap, GridapGmsh
using SegregatedVMSSolver.ParametersDef
using Plots, JLD2


fname = "n0012.csv" #airfoil coordinates to load
AoA = 4.0

ap0 = get_airfoil_coordinates(joinpath(@__DIR__, fname))

control_px = collect(LinRange(0.05,0.95,20))

control_points = ControlPoints(control_px,control_px)


R = 0.15 #support radius. If is bigger, the deformation radius increases, try increasing it
RBFfun = RBFFunctionLocalSupport(RBF_CP4, R)

rbfg = RBFGeometry(control_points,RBFfun)
rbfd = RBFDesign(rbfg, ap0)


sprob = StabilizedProblem(VMS(2))

physicalp = PhysicalParameters(Re=1000, u_in=[1.0,0.0])
timep = TimeParameters(dt=0.05, tF=15.0, time_window=(10.0, 15.0))

meshinfo = AirfoilMesh(AoA= AoA, meshref=2)
meshp = MeshParameters((1,1), 2, meshinfo)
exportp = ExportParameters(printinitial=true,printmodel=true,name_tags=["airfoil"], fieldexport=[["uh","ph","friction"]])

solverp = SolverParameters(θ=1.0)

simparams = SimulationParameters(timep,physicalp,solverp,exportp)

airfoil_case = Airfoil(meshp,simparams,sprob)

#Define the Objective Function, first argument [CD,CL], then you can define whatever you want

function J(CDCL; CLtarget=0.75)
    CD,CL=CDCL
    return -CL #0.5 * (CL - CLtarget)^2
end  



adj_solver = AdjSolver(δ=0.0001, opt_alg=:LD_LBFGS)

# adjoint_airfoil_problem = AdjointProblem( rbfd,airfoil_case,adj_solver,:steady, J)
# solve_adjoint_optimization(adjoint_airfoil_problem)

δv = 0.0001 .* [ones(20); -ones(20)]

fd1 = load("FD/FD_RES.jld2")["fd"]



ADJ_grad = zeros(40)
res0 = load("results/SOLUHPH2.jld2")



adjsol0 = load("results/ADJ_SOL0.jld2")["adj_sol"]
adjsol2 = load("results/ADJ_SOL2.jld2")["adj_sol"]




plot(adjsol0.airfoil_model.ap.xu,adjsol0.pressure_distribution.uval, linewidth = 2, label="Sol0", linecolor=:red)
plot!(adjsol0.airfoil_model.ap.xl,adjsol0.pressure_distribution.lval, linewidth = 2, label=false, linecolor=:red)

plot!(adjsol2.airfoil_model.ap.xu,adjsol2.pressure_distribution.uval, linewidth = 2, label="Sol2", linecolor=:black)
plot!(adjsol2.airfoil_model.ap.xl,adjsol2.pressure_distribution.lval, linewidth = 2, label=false, linecolor=:black)
yflip!()

# Jtot0= load("results/GJ2.jld2")["JtotD"][:J2s]
# plot(Jtot0)



for (IIDX, aa) in enumerate(ADJ_grad)


rbfd0 = rbfd
modelname0 =create_msh(meshinfo,rbfd0, physicalp ; iter = 0)
model0 = GmshDiscreteModel(modelname0)

rbfd1 = perturb_DesignParameter(rbfd0, [IIDX], [δv[IIDX]])
modelname1 =create_msh(meshinfo,rbfd1, physicalp ; iter = 0)
model1 = GmshDiscreteModel(modelname1)

# writevtk(model, "model_$iter")
Ω0 = Triangulation(model0)
Ω1 = Triangulation(model1)



order = 2
D = 2 
dΩ = Measure(Ω0, order*2)
ν=0.001


function fx(x)
  return x
end


reffe = ReferenceFE(lagrangian, VectorValue{D,Float64}, order)
VV0 = FESpace(model0, reffe; conformity=:H1)
m0 = interpolate_everywhere(fx,VV0)
VV1 = FESpace(model1, reffe; conformity=:H1)
m1 = interpolate_everywhere(fx,VV1)

v_field = interpolate_everywhere(VectorValue(0.0,0.0),VV0)
vi = (m1.free_values-m0.free_values)./ [δv[IIDX]]

v_field.free_values .= [x for u in vi for x in u]




reffeᵤ = ReferenceFE(lagrangian, VectorValue{D,Float64}, order )
V = TestFESpace(model0, reffeᵤ, conformity=:H1, dirichlet_tags=["inlet", "limits", "airfoil"],  dirichlet_masks=[(true,true), (false,true),(true,true) ])
u0 = VectorValue(1.0, 0.0)
u_walls = VectorValue(zeros(D)...) 
U = TrialFESpace(V, [u0, u0, u_walls])
reffeₚ = ReferenceFE(lagrangian, Float64, order)
Q = TestFESpace(model0, reffeₚ, conformity=:H1, dirichlet_tags=["outlet"])
P = TrialFESpace(Q, [0.0])


reffe_u_adj = ReferenceFE(lagrangian, VectorValue{D,Float64}, order)
V_adj = TestFESpace(model0, reffe_u_adj, conformity=:H1, dirichlet_tags=["airfoil", "outlet","limits"], dirichlet_masks=[(true,true), (true,true),(false,true) ])
reffe_p_adj = ReferenceFE(lagrangian, Float64, order)
Q_adj = TestFESpace(model0, reffe_p_adj, conformity=:H1, dirichlet_tags="inlet")

d_boundary =  [0.0, 1.0]
U_adj = TrialFESpace(V_adj, [VectorValue(d_boundary...), VectorValue(0, 0),VectorValue(0, 0)])
P_adj = TrialFESpace(Q_adj,0.0)


uh = interpolate(a-> VectorValue(rand(2)...),U)
ph = interpolate(a-> rand() ,P)

uhadj = interpolate(a-> VectorValue(rand(2)...),U_adj)
phadj = interpolate(a-> rand(),P_adj)





uh.free_values .= res0["uh"].free_values
ph.free_values .= res0["ph"].free_values

uhadj.free_values .= res0["uhadj"].free_values
phadj.free_values .= res0["phadj"].free_values

writevtk(Ω0, "Graduh", nsubcells=order, cellfields=["uh"=>uh, "∇(uh)"=>∇(uh) ])


Γ = BoundaryTriangulation(model0; tags="airfoil")
dΓ = Measure(Γ, order*2)
nΓ = get_normal_vector(Γ)

I1 = sum(∫(uhadj ⋅(∇(uh) ⋅ uh ⋅ (∇⋅(v_field))) )dΩ)
I2 = sum(∫(uhadj ⋅(∇(uh) ⋅ ( ∇(v_field) ⋅uh  ) ) )dΩ)
I3 = sum(∫(uhadj ⋅(∇(ph) ⋅ (∇⋅(v_field))) )dΩ)
I4 = sum(∫(2 * ν * ε(uh)⊙ε(uhadj)⋅(∇⋅(v_field)))dΩ)
I5 = sum(∫(phadj ⋅(  (∇⋅uh) ⋅ (∇⋅v_field)))dΩ)
I6 = sum(∫(phadj ⋅(  (∇(uh)) ⊙ (∇(v_field))))dΩ)

Ia1 = sum(∫(uhadj ⋅ (ν * ∇(uh) ⋅ nΓ) * (v_field ⋅ nΓ))dΓ)

# Pressure contribution  
Ia2 = sum(∫(uhadj ⋅ (ph * nΓ) * (v_field ⋅ nΓ))dΓ)

# Convective flux contribution
Ia3 = sum(∫(uhadj ⋅ ((uh ⋅ nΓ) * uh) * (v_field ⋅ nΓ))dΓ)

# Adjoint viscous stress
Ia4 = sum(∫((ν * ∇(uhadj) ⋅ nΓ) ⋅ uh * (v_field ⋅ nΓ))dΓ)

# Adjoint pressure
Ia5 = sum(∫((phadj * nΓ) ⋅ uh * (v_field ⋅ nΓ))dΓ)

boundary_sensitivity = sum(∫(uhadj ⋅ (ν * ∇(uh) ⋅ nΓ + ph * nΓ) * (v_field ⋅ nΓ))dΓ)


sum([I1,I2,I3,I4,I5,I6])
sum([Ia1, Ia2, Ia3,Ia4, Ia5])

G = sum(∫((-ph * uhadj + ν*(∇(uh) + transpose(∇(uh))) ⋅ uhadj) ⋅v_field)dΓ)


Gb = sum(∫( ((-ph  + 2ν * ε(uh)) ⋅ nΓ) ⋅ uhadj * (v_field ⋅ nΓ) )dΓ)


#writevtk(Ω0, "Γi", nsubcells=order, cellfields=["Gp"=>Gp, "Gu"=>Gu,  "∇(uh)"=>∇(uh) ])






# Convective terms
conv1 = transpose(∇(uh)) ⋅ uhadj ⊗ uh
conv2 = uhadj ⊗ (∇(uh) ⋅ uh)

# Pressure and continuity terms
pres = ∇(ph) ⊙ uhadj
cont = (∇⋅(uh)) ⊗ phadj

# Full T tensor
T = 2ν * (ε(uh) ⊙ ε(uhadj))+  + conv1 + conv2 + pres + cont


δJ = sum(∫( T ⊙ ∇(v_field) )dΩ) 


ADJ_grad[IIDX] = δJ


end



plot(ADJ_grad, label="adjoint gradient")


plot!(fd1, label = "finite difference gradient")

plot( - vcat(-reverse(ADJ_grad[21:end]), ADJ_grad[1:20]))



# Compute ∇v ⋅ n (normal derivative of velocity)
∇v_dot_n = ∇(uhadj) ⋅ nΓ

# First term: ν (n ⋅ ∇)v ⋅ ∇u
viscous_term = ν * (∇v_dot_n ⋅ ∇(uh))

# Second term: q * n ⋅ ∇u
q = phadj  # Or use adjoint pressure if appropriate
pressure_term = q * (nΓ ⋅  ∇(uh))

# Full integrand
Gb = (viscous_term - pressure_term) ⋅ v_field

δJ = -sum(∫(Gb)dΓ)
FDg[IIDX]



ADJ_grad
plot(FDg)

for (i,a) in enumerate(ADJ_grad)
  if i>20
  println(-a)
  else 
    println(a4)
  end
end
