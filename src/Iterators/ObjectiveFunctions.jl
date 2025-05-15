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

function compute_airfoil_forces(uh::SingleFieldFEFunction,ph::SingleFieldFEFunction,nΓ::OperationCellField,dΓ::GenericMeasure,ν::Float64)
    IForce = ∫( -ph ⋅ nΓ+ ν* transpose(∇(uh)) ⋅ nΓ)dΓ #+ 
    D,L = sum(IForce)
    return D,L
end

"""
    compute_airfoil_coefficients(uh::SingleFieldFEFunction,ph::SingleFieldFEFunction,nΓ::OperationCellField,dΓ::GenericMeasure, physicalp::PhysicalParameters)

Compute the normaization of airfoil forces, obtaining CD and CL
"""
function compute_airfoil_coefficients(uh::SingleFieldFEFunction,ph::SingleFieldFEFunction,nΓ::OperationCellField,dΓ::GenericMeasure, physicalp::PhysicalParameters)
    @unpack c, u_in_mag, ν = physicalp
    
    q = 0.5 * c*u_in_mag

    D,L = compute_airfoil_forces(uh,ph,nΓ,dΓ,ν )
    CD = D/q
    CL = L/q
    println("----------")
    println("CL = $CL ; CD = $CD")
    println("----------")

    return CD,CL
end


"""
    obj_fun(am::AirfoilModel, mcase::AdjointProblem, uh_vec::AbstractVector,ph_vec::AbstractVector, fun)

Wrapper that evaluates `fun` the objective function. `fun` is a user-defined function.
See eg. `compute_drag`, `compute_lift` functions. It return fitnessval,E. 
`fitnessval` is the function value that need to be minimized. `E` is a value that we want to monitor.
"""
function obj_fun(am::AirfoilModel, mcase::AdjointProblem, uh,ph, fun)
    @unpack model, params = am
    # @unpack U,P = params
    Γ = BoundaryTriangulation(model; tags="airfoil")
    dΓ = Measure(Γ, 4)
    nΓ = -1 .* get_normal_vector(Γ)

    # uh=FEFunction(U,uh_vec)
    # ph=FEFunction(P,ph_vec)
    
    physicalp =  mcase.vbcase.simulationp.physicalp
    CD,CL = compute_airfoil_coefficients(uh,ph,nΓ,dΓ, physicalp)
    fitnessval = fun([CD,CL])

    return fitnessval, [CD,CL]
end


function dJobj_fun(fun::Function, CDCL::Vector{Float64})
    @assert length(CDCL) == 2
    ForwardDiff.gradient(x -> fun(x), CDCL)

end
