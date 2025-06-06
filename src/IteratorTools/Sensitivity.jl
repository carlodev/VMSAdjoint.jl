
"""
    compute_sensitivity(am0,am1,  δ::Float64, simcase::Airfoil,  uh,ph,uhadj,phadj)

From the solution of the primal flow uh,ph and the adjoint solution uhadj and phadj it computes the senstivities.
am0:: Airfoil model of un-pertubed geometry
am0:: Airfoil model of pertubed geometry
δ:: is the signed perturbation
"""

function compute_sensitivity(am0::AirfoilModel,am1::AirfoilModel,  δ::Float64, simcase::Airfoil,  uh,ph,uhadj,phadj)

    @sunpack D,order,u_in, ν = simcase
    @unpack VV0, dΩ, reffe = am0.params

    function fx(x)
        return x
    end
    
    m0 = get_free_dof_values(interpolate_everywhere(fx,VV0)) #it gives the nodes position in the un-perturbed mesh
    VV1 = FESpace(am1.model, reffe; conformity=:H1)
    m1 = get_free_dof_values(interpolate_everywhere(fx,VV1)) #it gives the nodes position in the perturbed mesh
    vi = (m1-m0)./ δ

    v_field = FEFunction(VV0,vi)
   

    # Convective terms
    conv1 = transpose(∇(uh)) ⋅ uhadj ⊗ uh
    conv2 = uhadj ⊗ (∇(uh) ⋅ uh)

    # Pressure and continuity terms
    pres = ∇(ph) ⊙ uhadj
    cont = (∇⋅(uh)) ⊗ phadj

    # Full T tensor
    T = 2ν * (ε(uh) ⊙ ε(uhadj))+  conv1 + conv2 + pres + cont
    J2 = sum(∫( T ⊙ ∇(v_field) )dΩ) 

    return J2
    
end