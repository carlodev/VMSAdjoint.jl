using Revise
using AirfoilTools
using VMSAdjoint
using Gridap, GridapGmsh
using SegregatedVMSSolver.ParametersDef
using Plots
using JLD2

"""
Here we test different Sensitivity formulation coming from the priamal (time-averaged) and adjoint solution
"""

fname = "n0012.csv" #airfoil coordinates to load
AoA = 2.5
order, D = 2, 2
ν = 0.001

ap0 = get_airfoil_coordinates(joinpath(@__DIR__, fname))


control_px = collect(LinRange(0.05, 0.95, 20))

control_points = ControlPoints(control_px, control_px)


R = 0.15 #support radius. If is bigger, the deformation radius increases, try increasing it
RBFfun = RBFFunctionLocalSupport(RBF_CP4, R)


sprob = StabilizedProblem(VMS(order))

physicalp = PhysicalParameters(Re=1000, u_in=[1.0, 0.0])
timep = TimeParameters(dt=0.05, tF=10.0, time_window=(7.5, 10.0))

meshinfo = AirfoilMesh(AoA=0.0, meshref=2)
meshp = MeshParameters((1, 1), 2, meshinfo)
exportp = ExportParameters(printinitial=true, printmodel=true, name_tags=["airfoil"], fieldexport=[["uh", "ph", "friction"]])

solverp = SolverParameters(θ=1.0)

simparams = SimulationParameters(timep, physicalp, solverp, exportp)

airfoil_case = Airfoil(meshp, simparams, sprob)

rbfg = RBFGeometry(control_points, RBFfun)
rbfd0 = RBFDesign(rbfg, ap0)



#Select the solution to load
adj_sol0 = load("M/ADJ_SOL03.jld2")["adj_sol"]

#Loading the Solutions
am = adj_sol0.airfoil_model
uh,ph=adj_sol0.uh,adj_sol0.ph
uhadj,phadj=adj_sol0.uhadj,adj_sol0.phadj

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


reffe = ReferenceFE(lagrangian, VectorValue{D,Float64}, order)
VV0 = FESpace(am.model, reffe; conformity=:H1) 
function fx(x)
    return x
end

m0 = get_free_dof_values(interpolate_everywhere(fx, VV0))#space where to interpolate the variation


Jvol = zeros(length(δv))
Jbound = zeros(40)

for (IIDX, δ) in enumerate(δv)
    
    
    rbfd1 = perturb_DesignParameter(rbfd0, [IIDX], [δ])
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

    # Jbound[IIDX] = sum(∫( (T ⋅ nΓ) ⋅ v_field )dΓ) #it should be with -1, but the nΓ are pointing inward, so it should be ok
    Jbound[IIDX] =     sum(-ν* ∫( ((∇(uh) ⋅ tΓ) ⋅ nΓ) * ((∇(uhadj) ⋅ tΓ) ⋅ nΓ) ⋅ (v_field ⋅nΓ) )dΓ)


end

using Plots, LaTeXStrings
Plots.default(linewidth=2)

fd_CL = ([-0.13144596507117487,  -0.13144049689328133,-0.13142133204113804, -0.1314023025427092, -0.13138768434632345, -0.13137789568738115, -0.13137206566430665, -0.13136892753012633, -0.13136727395150294, -0.13136618863278962, -0.13136515002082605, -0.13136399015044142,-0.1313628193964369, -0.13136187795817136, -0.13136135770187202,  -0.13136154083666163, -0.13136252142175023, -0.1313639333878327, -0.13136581791491247, -0.1313665747711074,
-0.1313251131875929, -0.13132530677934243, -0.13132961301053342, -0.13133267190477554, -0.13133601782084586, -0.13133996905855333, -0.13134434732175732, -0.13134886622886952, -0.13135334579317356, -0.13135776040897643, -0.13136222610612688, -0.13136700121656528, -0.1313724027129784, -0.13137877509047755, -0.13138650213631226, -0.13139588932695623, -0.13140732515465886, -0.13142199456696485, -0.13144125040514898, -0.13146976870328622 ] .+ 0.13138352156761407)./δv

fd_CD = ([ 0.1216404544975305, 0.12164022538139854, 0.121639265678663, 0.12163864470430423, 0.12163829871723085, 0.12163808109652317, 0.12163789670091155, 0.12163769948157467, 0.12163746739003616, 0.1216371925430035, 0.12163687744654388, 0.12163653066004938, 0.12163616545553438, 0.12163579709642698, 0.12163544348892595, 0.12163511230041991, 0.1216348114898485, 0.1216345602199478, 0.12163435681023815, 0.12163427508648679, 0.12163706730354933, 0.12163706106603764, 0.12163663022813836, 0.1216362323852443, 0.12163596143757201, 0.12163582081275652, 0.12163578439796892, 0.12163582287820633, 0.12163591042708378, 0.12163602986232984, 0.12163617245015604, 0.12163633750628997, 0.12163652509682873, 0.12163673532585907, 0.12163697580637903, 0.12163722805357227, 0.12163748254041373, 0.12163779866991345, 0.12163817398616644, 0.12163895833965868] .- 0.12163386560377128) ./ δv


grad_c = -2 .* [Jbound[1:20]; - Jbound[21:end]]
fd_c = [fd_CD[1:20]; - fd_CD[21:end]]


plot(grad_c, label = "Adjoint", linecolor=:black)
scatter!(fd_c, label = "Finite Differences", markercolor=:black, markershape=:xcross, markersize=5.0, markerstrokewidth=2.0)
plot!(xlabel="Design variable β", ylabel = L"C_D"* " gradient (outward deformation >0)")

savefig("CD_grad.svg")
savefig("CD_grad.pdf")


#Plot airfoil
modelname00 = create_msh(meshinfo, rbfd0, physicalp, "MeshPerturb"; iter=IIDX + 100)
model00 = GmshDiscreteModel(modelname00)
am00 = AirfoilModel(model00, airfoil_case; am=nothing)

plot(aspect_ratio = 1)
plot!(am00.ap.xu, am00.ap.yu, linecolor=:black, label = false)
plot!(am00.ap.xl, am00.ap.yl, linecolor=:black, label=false)
scatter!(control_points.cu, rbfd0.cy.cu, markercolor=:black, markersize=3.0, label = false)
scatter!(control_points.cl, rbfd0.cy.cl, markercolor=:black, markersize=3.0, label=false)

dy = 0.05
annotate!(control_points.cu,rbfd0.cy.cu .+ dy , ["$t" for t in collect(1:20)], (8))

annotate!(control_points.cl,rbfd0.cy.cl .- dy , ["$t" for t in collect(21:40)], (8))

plot!(axis=false,       # removes entire axis (ticks, spines, labels)
grid=false,       # removes grid
framestyle=:none  # removes plot frame
)


savefig("Airfoil_DesignParameters.svg")
savefig("Airfoil_DesignParameters.pdf")