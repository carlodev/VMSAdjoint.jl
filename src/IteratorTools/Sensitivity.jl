
"""
    compute_sensitivity(model::DiscreteModel, params::Dict{Symbol,Any}, uh0,ph0,ϕu0, ϕp0; objective_function=compute_drag)

From the solution of the primal flow `uh0` `ph0`, and the adjoint flow `ϕu0` `ϕp0` it computes the senstivities according to the objective function.
J2 contributions are splitted and the gradients computed individually to avoid numerical cancellation
"""
function  compute_sensitivity(simcase::Airfoil, am::AirfoilModel, d_boundary::Vector{Float64}, uh0,ph0,ϕu0, ϕp0, J::Function; i::Int64=1)
    @sunpack D,order,u_in, ν = simcase
    model = am.model
    V,Q = create_primal_spaces(model, simcase)

    u0 = VectorValue(u_in...)
    u_walls = VectorValue(zeros(D)...)
    p0 = 0.0

    U = TrialFESpace(V, [u0, u0, u_walls])
    P = TrialFESpace(Q, p0)


    V_adj,Q_adj = create_adjoint_spaces(model, simcase)
    U_adj = TrialFESpace(V_adj, [VectorValue(d_boundary...), VectorValue(zeros(D)...),VectorValue(zeros(D)...)])
    P_adj = TrialFESpace(Q_adj,0.0)
    
    
    uh = FEFunction(U,uh0.free_values)
    ph = FEFunction(P,ph0.free_values)
    ϕu = FEFunction(U_adj,ϕu0.free_values)
    ϕp = FEFunction(P_adj,ϕp0.free_values)

    Ω = Triangulation(model)
    dΩ = Measure(Ω, order*2)

    J1, _ = obj_fun(am, simcase, uh,ph, J)

    Γ = BoundaryTriangulation(model; tags="airfoil")
    dΓ = Measure(Γ, order*2)
    nΓ =  get_normal_vector(Γ) #-1, pointing outward

    J2_1 = sum(∫(ϕu⋅(transpose(∇(uh))⋅uh ))dΩ) 
    J2_2 = sum(∫(ϕu⋅(∇(ph)))dΩ)
    J2_3 = sum(∫(ϕp⋅(∇⋅(uh)))dΩ)
    J2_4 = ν*sum(∫((∇(ϕu)⊙∇(uh)))dΩ)
    J2_5 = -ν* sum(∫(ϕu⋅transpose(∇(uh))⋅ nΓ)dΓ)
    
    #export i
    dir = "MeshPerturb"
    mkpath(dir)
    fname = joinpath(dir,"model-$(i)")
    writevtk(Ω, fname,nsubcells=order,  cellfields=["uh"=>uh,"ph"=>ph, "uadj"=>ϕu,"padj"=>ϕp])

    
    return J1,(J2_1,J2_2,J2_3,J2_4,J2_5)
end