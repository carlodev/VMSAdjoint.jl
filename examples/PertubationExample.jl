
using Revise
using AirfoilTools
using VMSAdjoint
using Gridap, GridapGmsh
using SegregatedVMSSolver.ParametersDef

AoA = 15.0 #2.5
Re = 15_000
NW0 = 15 #number of w0 to parametrize top and bottom

cst0 = CST_NACA0012(;N=NW0, t=0.0)
cstd0 = AirfoilCSTDesign(cst0,1000)

cstd0.cstg.cstw.wu
cstd0.cstg.cstw.wl

cstd1 = perturb_DesignParameter(cstd0, 1, -0.08)
cstd1 = perturb_DesignParameter(cstd1, NW0+1, 0.08)

cstd1 = perturb_DesignParameter(cstd1, 2, 1.5)
cstd1 = perturb_DesignParameter(cstd1, NW0+2, 0.07)

cstd1 = perturb_DesignParameter(cstd1, 3, -0.07)
cstd1 = perturb_DesignParameter(cstd1, NW0+3, 0.07)

cstd1.cstg.cstw.wu[1]
cstd1.cstg.cstw.wl[1]

cstd1.cstg.cstw.wu

using Plots
Plots.default(aspect_ratio=:equal, linewidth=2)
plot(cstd0.ap.xu,cstd0.ap.yu, linecolor=:black, label="original")
plot!(cstd0.ap.xl,cstd0.ap.yl, linecolor=:black, label=false)

plot!(cstd1.ap.xu,cstd1.ap.yu, linecolor=:red, label="perturbed")
plot!(cstd1.ap.xl,cstd1.ap.yl, linecolor=:red, label=false)