using Revise
using AirfoilTools
using VMSAdjoint
using Gridap, GridapGmsh
using SegregatedVMSSolver.ParametersDef
using Plots
using JLD2





adj_sol = load("results/ADJ_SOL0.jld2")["adj_sol"]
am = adj_sol.airfoil_model

AoA = 15.0#2.5
Re = 15_000
NW0 = 15 #number of w0 to parametrize top and bottom

cst0 = CST_NACA0012(;N=NW0, t=0.0)

cstd0 = AirfoilCSTDesign(cst0,200)
modelname =create_msh(meshinfo,cstd0, physicalp ; iter = 200)
model = GmshDiscreteModel(modelname)
writevtk(model, "model_t")


am =  AirfoilModel(model, airfoil_case; am=nothing)


xx0 = maximum([minimum(am.ap.xu);minimum(am.ap.xl)] )+ leading_edge_cutoff
xx1 =minimum([maximum(am.ap.xu);maximum(am.ap.xl)] )- trailing_edge_cutoff
leading_edge_cutoff = 0.01
trailing_edge_cutoff = 0.01
tmin=0.005

@assert xx1 > xx0 "Airfoil x-coordinates don't span a valid range"

# Filter points to exclude the leading edge region

# For upper surface
valid_upper_idx = findall(x -> xx1  >= x >= xx0 , am.ap.xu)
xu_filtered = am.ap.xu[valid_upper_idx]
yu_filtered = am.ap.yu[valid_upper_idx]

# For lower surface
valid_lower_idx = findall(x ->xx1  >= x >= xx0 , am.ap.xl)
xl_filtered = am.ap.xl[valid_lower_idx]
yl_filtered = am.ap.yl[valid_lower_idx]




# Verify the filtered coordinates are sorted
@assert issorted(xu_filtered) "Upper surface x-coordinates not sorted after filtering"
@assert issorted(xl_filtered) "Lower surface x-coordinates not sorted after filtering"

# Create evaluation points (excluding leading edge)
xx = collect(LinRange(xx0 , xx1, 201))

# Interpolate surfaces
yu = linear_interpolation(xu_filtered, yu_filtered, extrapolation_bc=Line()).(xx)
yl = linear_interpolation(xl_filtered, yl_filtered, extrapolation_bc=Line()).(xx)

# Calculate thickness
Δy = yu - yl

# Check for minimum thickness violations, excluding endpoints for stability
thickness_violations = tmin .- Δy


penalty = sum((violation > 0) ? violation^2 : 0.0 for violation in thickness_violations)





 