###################################################################################
#OBJECTIVE FUNCTIONS
##################################################################################





###################################################################################
#Iterators
##################################################################################
"""
   compute_airfoil_forces(uh::SingleFieldFEFunction,ph::SingleFieldFEFunction,nΓ::OperationCellField,dΓ::GenericMeasure,ν::Float64)

It computes Drag and Lift over airfoil boundary. It takes into account pressure and velocity gradient.
It needs the normals pointings outward respect to the body.
"""

function compute_airfoil_forces(uh::SingleFieldFEFunction, ph::SingleFieldFEFunction, nΓ::OperationCellField, dΓ::GenericMeasure, ν::Float64)
    IForce = ∫(-ph ⋅ nΓ + ν * (∇(uh) + transpose(∇(uh))) ⋅ nΓ)dΓ #+ 
    D, L = sum(IForce)
    return D, L
end

"""
    compute_airfoil_coefficients(uh::SingleFieldFEFunction,ph::SingleFieldFEFunction,nΓ::OperationCellField,dΓ::GenericMeasure, physicalp::PhysicalParameters)

Compute the normaization of airfoil forces, obtaining CD and CL
"""
function compute_airfoil_coefficients(uh::SingleFieldFEFunction, ph::SingleFieldFEFunction, nΓ::OperationCellField, dΓ::GenericMeasure, physicalp::PhysicalParameters)
    @unpack c, u_in_mag, ν = physicalp

    q = 0.5 * c * u_in_mag

    D, L = compute_airfoil_forces(uh, ph, nΓ, dΓ, ν)
    CD = D / q
    CL = L / q
    println("----------")
    println("CL = $CL ; CD = $CD")
    println("----------")

    return CD, CL
end


"""
    obj_fun(am::AirfoilModel, mcase::AdjointProblem, uh_vec::AbstractVector,ph_vec::AbstractVector, fun)

Wrapper that evaluates `fun` the objective function. `fun` is a user-defined function.
See eg. `compute_drag`, `compute_lift` functions. It return fitnessval,E. 
`fitnessval` is the function value that need to be minimized. `E` is a value that we want to monitor.
"""
function obj_fun(am::AirfoilModel, vbcase::Airfoil, uh, ph, fun::Function)
    @sunpack order = vbcase
    @unpack model, params = am

    Γ = BoundaryTriangulation(model; tags="airfoil")
    dΓ = Measure(Γ, 2 * order)
    nΓ = -1 .* get_normal_vector(Γ)

    physicalp = vbcase.simulationp.physicalp
    CD, CL = compute_airfoil_coefficients(uh, ph, nΓ, dΓ, physicalp)

    thick_pen = thickness_penalty(am)
    fitnessval = fun([CD, CL])

    println("----------")
    println("fitnessval = $(fitnessval) ; thickness_penalty = $(thick_pen)")
    println("----------")

    return fitnessval + thick_pen, [CD, CL]
end

function thickness_penalty(am::AirfoilModel; tmin=0.005, α=1_000.0)
    
    leading_edge_cutoff = 0.01
    trailing_edge_cutoff = 0.01

    # Find the overall x-range of the airfoil
    xx0 = maximum([minimum(am.ap.xu);minimum(am.ap.xl)] )+ leading_edge_cutoff
    xx1 =minimum([maximum(am.ap.xu);maximum(am.ap.xl)] )- trailing_edge_cutoff
 

    @assert xx1 > xx0 "Airfoil x-coordinates don't span a valid range"

    # Filter points to exclude the leading and trailing edge regions

    # For upper surface
    valid_upper_idx = findall(x -> xx1  >= x >= xx0 , am.ap.xu)
    xu_filtered = am.ap.xu[valid_upper_idx]
    yu_filtered = am.ap.yu[valid_upper_idx]

    # For lower surface
    valid_lower_idx = findall(x ->xx1  >= x >= xx0 , am.ap.xl)
    xl_filtered = am.ap.xl[valid_lower_idx]
    yl_filtered = am.ap.yl[valid_lower_idx]

    # Verify the filtered coordinates are sorted
    @assert issorted(xu_filtered) "Upper surface x-coordinates not sorted"
    @assert issorted(xl_filtered) "Lower surface x-coordinates not sorted"

    # Create evaluation points
    xx01 = maximum([minimum(xu_filtered);minimum(xl_filtered)] )
    xx11 =minimum([maximum(xu_filtered);maximum(xl_filtered)] )
    xx = collect(LinRange(xx01, xx11, 201))

    
    # Interpolate surfaces
    yu = linear_interpolation(xu_filtered, yu_filtered).(xx)
    yl = linear_interpolation(xl_filtered, yl_filtered).(xx)

    # Calculate thickness
    Δy = yu - yl

    # Check for minimum thickness violations
    thickness_violations = tmin .- Δy
    violations_locations = findall(thickness_violations.>0)
    !isempty(violations_locations) && println("Violations Locations x= $(xx[violations_locations])")
    
    #Σviolations^2
    penalty = sum((violation > 0) ? violation^2 : 0.0 for violation in thickness_violations)


    if isnan(penalty) || isinf(penalty)
        @error "Invalid penalty value detected. Check airfoil coordinates."
        return α   # Return a large penalty for invalid configurations
    end

    return penalty * α
end

function dJobj_fun(fun::Function, CDCL::Vector{Float64})
    @assert length(CDCL) == 2
    ForwardDiff.gradient(x -> fun(x), CDCL)
end
