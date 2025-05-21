using Revise
using AirfoilTools
using VMSAdjoint
using Gridap, GridapGmsh
using SegregatedVMSSolver.ParametersDef
using Plots
using JLD2


adj_sol = load("results/ADJ_SOL10.jld2")["adj_sol"]

ap = adj_sol.airfoil_model.ap
cp = adj_sol.pressure_distribution

Plots.default(linewidth=2, linecolor=:black)

plot(ap.xu[1:end-1], cp.uval[1:end-1], linecolor=:blue)
plot!(ap.xl[1:end-1], cp.lval[1:end-1], linecolor=:red)
yflip!()

using Trapz, Interpolations
xx = 0.0:0.001:1.0

cpu = linear_interpolation(ap.xu[1:end-1], cp.uval[1:end-1],extrapolation_bc=Line()).(xx)
cpl = linear_interpolation(ap.xl[3:end-1], cp.lval[3:end-1],extrapolation_bc=Line()).(xx)

trapz(xx, cpu-cpl)

