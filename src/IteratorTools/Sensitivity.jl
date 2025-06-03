
"""
    compute_sensitivity(model::DiscreteModel, params::Dict{Symbol,Any}, uh0,ph0,ϕu0, ϕp0; objective_function=compute_drag)

From the solution of the primal flow `uh0` `ph0`, and the adjoint flow `ϕu0` `ϕp0` it computes the senstivities according to the objective function.
J2 contributions are splitted and the gradients computed individually to avoid numerical cancellation
"""

function compute_sensitivity(am0,am1,  δ::Float64, simcase::Airfoil,  uh,ph,uhadj,phadj, J)

    @sunpack D,order,u_in, ν = simcase
    @unpack VV0, dΩ, reffe = am0.params

    function fx(x)
        return x
    end
    
    m0 = get_free_dof_values(interpolate_everywhere(fx,VV0))
    VV1 = FESpace(am1.model, reffe; conformity=:H1)
    m1 = get_free_dof_values(interpolate_everywhere(fx,VV1))
    vi = (m1-m0)./ δ

    v_field = FEFunction(VV0,vi)
   



    # Convective terms
    conv1 = transpose(∇(uh)) ⋅ uhadj ⊗ uh
    conv2 = uhadj ⊗ (∇(uh) ⋅ uh)

    # Pressure and continuity terms
    pres = ∇(ph) ⊙ uhadj
    cont = (∇⋅(uh)) ⊗ phadj

    # Full T tensor
    T = 2ν * (ε(uh) ⊙ ε(uhadj))+  + conv1 + conv2 + pres + cont
    J2 = sum(∫( T ⊙ ∇(v_field) )dΩ) 
    println("J2")

    return J2
    
end



