using Revise
using AirfoilTools
using VMSAdjoint
using Gridap, GridapGmsh
using SegregatedVMSSolver.ParametersDef
using Plots
using JLD2

fname = "n0012.csv" #airfoil coordinates to load
AoA = 4.0

ap0 = get_airfoil_coordinates(joinpath(@__DIR__, fname))


control_px = collect(LinRange(0.05,0.95,20))

control_points = ControlPoints(control_px,control_px)


R = 0.12 #support radius. If is bigger, the deformation radius increases, try increasing it
RBFfun = RBFFunctionLocalSupport(RBF_CP4, R)

rbfg = RBFGeometry(control_points,RBFfun)
rbfd0 = RBFDesign(rbfg, ap0)
rbfd1 = create_AirfoilDesign(rbfd0,w1)

rbfd2 = perturb_DesignParameter(rbfd1, [1,2], [0.02,-0.02])



Plots.default(linewidth=2)
plot(rbfd1.ap.xu,rbfd1.ap.yu, linecolor=:black, label="w1",aspect_ratio=:equal)
plot!(rbfd1.ap.xl,rbfd1.ap.yl, linecolor=:black, label=false)
plot!(rbfd2.ap.xu,rbfd2.ap.yu, linecolor=:red, label="w2")
plot!(rbfd2.ap.xl,rbfd2.ap.yl, linecolor=:red, label=false)


plot(w1, linecolor=:black, label="w1")
plot!(w2t, linecolor=:red, label="w2")

plot!(lb, linecolor=:black, linestyle=:dash, label="lower limit")
plot!(ub, linecolor=:black, linestyle=:dash, label="upper limit")
w2t = w1 - 0.01.*Jtot1



lb,ub = bounds_w(rbfd0,40)

(w1 - w0) ./ Jtot0

(w2 - w1) ./ Jtot1



adj_sol0 = load("results/ADJ_SOL0.jld2")["adj_sol"]
w0 = adj_sol0.βv
J0 = load("results/GJ0.jld2")["JtotD"]

Jtot0 = sum([J0[:J1s],J0[:J2s1],J0[:J2s2],J0[:J2s3],J0[:J2s4],J0[:J2s5]])

δ = 0.001 .* [ones(20); -ones(20)]
fd1 = [0.15061199763266886, 0.15033356437931247, 0.15054547342277966, 0.15068591045978938, 0.15076471239644818, 0.15080041866664554, 0.1508096963679196, 0.15080510797389515, 0.1507946732613078, 0.15078293221160388,
0.15077189029754565, 0.15076193722181774, 0.15075267805036696, 0.15074348134977245, 0.15073389029611, 0.15072321747440623, 0.15071101899152073, 0.15069786805225782, 0.1506835761120356, 0.1506788562428705,
0.15108644020718007,0.15103615756092043,0.15095179569982337, 0.15094136073416822,0.150937255072023, 0.15093100157630426,0.15092152790985264,0.15090871506491885,0.15089257261633635,0.1508732185434178,
0.15085056704768893, 0.1508240877338871, 0.15079282217996806,0.1507552649922497, 0.15070910031304607, 0.15065114021610662, 0.15057689818849648, 0.1504763556947541, 0.15033280702133686,0.1500987153528417  ] .- 0.1506570405150251
plot(J0[:J2s])


δ[1]=  δ[1]
plot(J0[:J2s])
scatter!(fd1 ./δ)

adj_sol1 = load("results/ADJ_SOL1.jld2")["adj_sol"]
w1 = adj_sol1.βv
J1 = load("results/J1.jld2")["JtotD"]
Jtot1 = sum([J1[:J1s],J1[:J2s1],J1[:J2s2],J1[:J2s3],J1[:J2s4],J1[:J2s5]])


adj_sol2 = load("results/ADJ_SOL2.jld2")["adj_sol"]
w2 = adj_sol2.βv
J2 = load("results/J2.jld2")["JtotD"]
Jtot2 = sum([J2[:J1s],J2[:J2s1],J2[:J2s2],J2[:J2s3],J2[:J2s4],J2[:J2s5]])


Plots.default(linewidth=2)
plot(Jtot2)



plot(w1)
plot!(Jtot1)
plot!(w2)


plot(Jtot0)
plot!(Jtot1)

plot([Jtot1[1:20];-Jtot1[21:end]])

J = J2

# plot(Jtot)
plot(J[:J1s] , label="geometric Gradient")
plot!(J[:J2s1], label="J2s1")
plot!(J[:J2s2], label="J2s2")
plot!(J[:J2s3],label="J2s3")
plot!(J[:J2s4],label="J2s4")
plot!(J[:J2s5],label="J2s5")



plot(adj_sol0.fgrad)
plot(gg, label="geometric gradient")
plot!(ag, label="adjoint gradient")
plot!(adj_sol1.fgrad)


w0 = [0.0344792, 0.0460489, 0.0526251, 0.057164, 0.0590081, 0.0599102, 0.0595747, 0.0583145, 0.05632, 0.0536866, 0.052162, 0.0487619, 0.0449802, 0.0409174, 0.0345058, 0.0301515, 0.0258337, 0.0196051, 0.0139143, 0.0076108, -0.0344792, -0.0460489, -0.0526251, -0.057164, -0.0590081, -0.0599102, -0.0595747, -0.0583145, -0.05632, -0.0536866, -0.052162, -0.0487619, -0.0449802, -0.0409174, -0.0345058, -0.0301515, -0.0258337, -0.0196051, -0.0139143, -0.0076108]

w1 = [0.03457127530617489, 0.046136270586943105, 0.05269126031165398, 0.05721274170936788, 0.059040787472723004, 0.05992826626849519, 0.059579785838018984, 0.0583083775015143, 0.05630453811675129, 0.05366361682147368, 0.0521332078897122, 0.0487288178600049, 0.04494418236096677, 0.04087961574405527, 0.034467170727691954, 0.03011295356461803, 0.025796157134569748, 0.019568885116587394, 0.013880264549259148, 0.007579380259992863, -0.03446295264610677, -0.04603585366282078, -0.052609693333660876, -0.057145268625741596, -0.05898618487236355, -0.05988539100135636, -0.05954732749296684, -0.058284916220587936, -0.05628861106161116, -0.053653867385525314, -0.05212846542976592, -0.048728199760316476, -0.044947138546486264, -0.04088602167181315, -0.034477484643123245, -0.030128320116068922, -0.025818746321743337, -0.019603071608441137, -0.013933338609764315, -0.00766253609842717]
