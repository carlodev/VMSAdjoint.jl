
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
    @unpack tΓ, nΓ, dΓ = am0.params

    function fx(x)
        return x
    end
    
    m0 = get_free_dof_values(interpolate_everywhere(fx,VV0)) #it gives the nodes position in the un-perturbed mesh
    VV1 = FESpace(am1.model, reffe; conformity=:H1)
    m1 = get_free_dof_values(interpolate_everywhere(fx,VV1)) #it gives the nodes position in the perturbed mesh
    vi = (m1-m0)./ δ

    v_field = FEFunction(VV0,vi)

    #-2 coming from the derivative of dCL/dβ or dCD/dβ; to be adjusted if more complex derivations
    J_sens = -2 .*    sum(-ν* ∫( ((∇(uh) ⋅ tΓ) ⋅ nΓ) * ((∇(uhadj) ⋅ tΓ) ⋅ nΓ) ⋅ (v_field ⋅nΓ) )dΓ)

    return J_sens
    
end