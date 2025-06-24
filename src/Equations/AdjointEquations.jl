####################################################################
#STEADY ADJOINT
####################################################################

function equations_adjoint( simcase::Airfoil,params::Dict{Symbol,Any},time_dep::Symbol)
    equations_adjoint( simcase,params,Val(time_dep))
end

"""
    adjoint_conservation(uh)

Steady conservation equations of the adjoint problem 
"""
function adjoint_conservation(uh)
    Rmadj(u, p) =  - ∇(p) - transpose(∇(u)) ⋅ uh - ∇(u)⋅uh
    Rcadj(u) = -1 .* (∇ ⋅ (u))
    return Rmadj, Rcadj
end

"""
    adjoint_steady_VMS(params::Dict{Symbol,Any})

It provides the set of Adjoint Equations
Symbol convention: O. Soto, & R. Lohner. (2004). On the Boundary Computation of Flow Sensitivities. https://doi.org/10.2514/6.2004-112
"""
function equations_adjoint( simcase::Airfoil,params::Dict{Symbol,Any},::Val{:steady})
    @sunpack ν, dt, D, θ = simcase
    @unpack uh,dΩ,Ω = params

    stab_coeff = SegregatedVMSSolver.Equations.compute_stab_coeff(simcase,params)
    Tm = SegregatedVMSSolver.Equations.momentum_stabilization(uh, stab_coeff, simcase)
    Tc = SegregatedVMSSolver.Equations.continuity_stabilization(uh, stab_coeff, simcase)

    Rmadj, Rcadj = adjoint_conservation(uh)

    TRm(u, p) = Tm * Rmadj(u, p)
    
    ADJBᴳ((u, p), (v, q)) = ∫(ν * ∇(v) ⊙ ∇(u) - q * (∇ ⋅ u))dΩ + ∫(v ⊙ Rmadj(u, p))dΩ
    
    ADJB_SUPG((u, p), (v, q)) = ∫( (- transpose(∇(v)) ⋅ uh - ∇(v)⋅uh- ∇(q)) ⊙ TRm(u, p))dΩ +
            -1 * ∫( Tc ⋅ (∇ ⋅ v) ⊙ Rcadj(u))dΩ

    ADJB_VMS1((u, p), (v, q)) = ∫((-uh ⋅ (∇(v))') ⊙ TRm(u, p))dΩ


    res_adj( (u, p), (v, q)) = ADJBᴳ((u, p), (v, q)) + ADJB_SUPG((u, p), (v, q)) + ADJB_VMS1((u, p), (v, q))
    
    rhs_adj((v,q)) = ∫(VectorValue(zeros(D)...) ⋅ v)dΩ

    return res_adj, rhs_adj
end



"""
    equations_adjoint( simcase::Airfoil,params::Dict{Symbol,Any},::Val{:unsteady})

Unsteady Adjoint Equations
"""
function equations_adjoint( simcase::Airfoil,params::Dict{Symbol,Any},::Val{:unsteady})
    @sunpack ν, dt, D, θ = simcase
    @unpack uh,dΩ,Ω = params

    res_adj_s, rhs_adj_s = equations_adjoint( simcase,params,Val(:steady))

    stab_coeff = SegregatedVMSSolver.Equations.compute_stab_coeff(simcase,params)
    Tm = SegregatedVMSSolver.Equations.momentum_stabilization(uh, stab_coeff, simcase)


    time_sign = 1.0
    m_adj(t, (u, p), (v, q)) =   time_sign *∫(u ⋅ v)dΩ -  ∫( Tm ⋅ (-uh ⋅ ∇(v) - ∇(v) ⋅ uh - ∇(q)) ⋅ u)dΩ


    res_adj(t, (u, p), (v, q)) = res_adj_s( (u, p), (v, q))
    
    rhs_adj(t,(v,q)) = rhs_adj_s((v,q))


    return m_adj, res_adj, rhs_adj
end
