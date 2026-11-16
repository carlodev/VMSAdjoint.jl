using Plots

using AirfoilTools
using AirfoilTools.AirfoilCST
using VMSAdjoint
using Gridap, GridapGmsh
using SegregatedVMSSolver.ParametersDef
"""
Example taken from:
Sorgiovanni, G., Quadrio, M., Ponzini, R., 2016. A robust open-source adjoint optimization method for external aerodynamics. Politecnico di Milano, Milan.

Starting from a NACA0012 finding the derivaties with the adjoint method
"""


alpha= 0.0
fname = "../n0012.csv" #airfoil coordinates to load
AoA = 2.5
ap0 = get_airfoil_coordinates(joinpath(@__DIR__, fname))

control_px = collect(LinRange(0.05,0.95,20))
control_points = ControlPoints(control_px,control_px)
R = 0.15 
RBFfun = RBFFunctionLocalSupport(RBF_CP4, R)

rbfg = RBFGeometry(control_points,RBFfun)
rbfd = RBFDesign(rbfg, ap0)




cst1 = perturb_DesignParameter(rbfd, 2, 0.023)
Plots.default(lw=2)

plot(ap0.xu, ap0.yu,lw=4, lc=:black,label="original")
plot!(ap0.xl, ap0.yl,lw=4, lc=:black, label=false )
plot!(cst1.ap.xu, cst1.ap.yu, lc=:red,label="perturbed")
plot!(cst1.ap.xl, cst1.ap.yl, lc=:red, label=false )
plot!(ylabel="y/c", xlabel="x/c")
plot!(axis = true,
grid = true,
aspect_ratio = :equal,
legendfontsize = 12,tickfontsize = 12,guidefontsize = 12,
)

savefig("NACA0012_Perturbed_d2.svg")
savefig("NACA0012_Perturbed_d2.pdf")



plot(ap0.xu, ap0.yu,lw=4, lc=:black,label=false)
scatter!(control_px, rbfd.cy.cu, markershape=:circle, markercolor=:black, label="control points")


plot!(cst1.ap.xu, cst1.ap.yu, lc=:red,label=false)
plot!(axis = false,
grid = false,
aspect_ratio = :equal,
legendfontsize = 16,tickfontsize = 12,guidefontsize = 12,
)
plot!(xlims=[0.0, 0.25])

savefig("NACA0012_Perturbed_detail.svg")
savefig("NACA0012_Perturbed_detail.pdf")

