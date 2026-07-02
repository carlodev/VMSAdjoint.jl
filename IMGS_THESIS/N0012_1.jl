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



using Plots

cst1 = perturb_DesignParameter(rbfd, 25, 0.3)

plot(ap0.xu, ap0.yu, lc=:white)
plot!(ap0.xl, ap0.yl, lc=:white )

scatter!(control_px, rbfd.cy.cu, markershape=:circle, markercolor=:white)
scatter!(control_px, rbfd.cy.cl, markershape=:circle, markercolor=:white)

for i in eachindex(control_px)
    annotate!(control_px[i], rbfd.cy.cu[i]+ 0.01, text("$i", 8, :white, :bottom))
    annotate!(control_px[i], rbfd.cy.cl[i]-0.01, text("$(i + 20)", 8, :white, :top))
end


plot(cst1.ap.xu, cst1.ap.yu, lc=:white)
plot!(cst1.ap.xl, cst1.ap.yl, lc=:white )
plot!(axis = false,
grid = false,
legend = false,
background_color = :transparent,
foreground_color = :transparent,
ticks = nothing,
border = :none,
aspect_ratio = :equal
)

savefig("Non-Physical.svg")


x = collect(1:0.01:5)
y(x) = sin(2.5*x) + log(x)


plot(x,y.(x), lw=3, lc=:green1)

plot!(
legend = false,
background_color = :transparent,
foreground_color = :transparent,
aspect_ratio = :equal
)