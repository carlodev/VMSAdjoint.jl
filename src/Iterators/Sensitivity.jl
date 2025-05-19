
"""
    compute_sensitivity(model::DiscreteModel, params::Dict{Symbol,Any}, uh0,ph0,ϕu0, ϕp0; objective_function=compute_drag)

From the solution of the primal flow `uh0` `ph0`, and the adjoint flow `ϕu0` `ϕp0` it computes the senstivities according to the objective function.
J2 contributions are splitted and the gradients computed individually to avoid numerical cancellation
"""
function  compute_sensitivity(simcase::Airfoil, am::AirfoilModel, d_boundary::Vector{Float64}, uh0,ph0,ϕu0, ϕp0, J::Function; i::Int64=1)
    @sunpack D,order,u_in, ν ,u_in = simcase
    model = am.model
    V,Q = create_primal_spaces(model, simcase)

    u0 = VectorValue(u_in...)
    u_walls = VectorValue(zeros(D)...)
    p0 = 0.0

    U = TrialFESpace(V, [u0, u0, u_walls])
    P = TrialFESpace(Q, p0)


    V_adj,Q_adj = create_adjoint_spaces(model, simcase)
    U_adj = TrialFESpace(V_adj, [VectorValue(d_boundary...), VectorValue(zeros(D)...),VectorValue(zeros(D)...)])
    P_adj = TrialFESpace(Q_adj)
    
    
    uh = FEFunction(U,uh0.free_values)
    ph = FEFunction(P,ph0.free_values)
    ϕu = FEFunction(U_adj,ϕu0.free_values)
    ϕp = FEFunction(P_adj,ϕp0.free_values)

    Ω = Triangulation(model)
    dΩ = Measure(Ω, order*2)

    J1, _ = obj_fun(am, simcase, uh,ph, J)

    Γ = BoundaryTriangulation(model; tags="airfoil")
    dΓ = Measure(Γ, order*2)
    nΓ =  -1 .*get_normal_vector(Γ) #-1, pointing outward

    J2_1 = sum(∫(ϕu⋅(transpose(∇(uh))⋅uh ))dΩ) 
    J2_2 = sum(∫(ϕu⋅(∇(ph)))dΩ)
    J2_3 = sum(∫(ϕp⋅(∇⋅(uh)))dΩ)
    J2_4 = ν*sum(∫((∇(ϕu)⊙∇(uh)))dΩ)
    J2_5 = -ν* sum(∫(ϕu⋅transpose(∇(uh))⋅ nΓ)dΓ)
    
    #export i
    writevtk(Ω, "MeshPerturb/model-$(100+i)", cellfields=["uh"=>uh,"ph"=>ph, "uadj"=>ϕu,"padj"=>ϕp])

    
    return J1,(J2_1,J2_2,J2_3,J2_4,J2_5)
end



# ###################################################################################
# #ITERATORS
# ##################################################################################
# function iterate_optimization(des_params, params; iter=0, objective_function=compute_drag)
    
#     modelname = create_msh(des_params; iter = iter,  mesh_ref=4)
#     model = GmshDiscreteModel(modelname)
#     uh,ph = solve_inc_primal(model, tagname; filename=joinpath("Results_primal", "res-$(iter)"))
#     _, S = obj_fun(model,params,uh, ph,objective_function)

# return model, S, (uh,ph)
# end



# function iterate_optimization(uh,ph,model, des_params, params; iter=0, detail=false, δ=0.02, α=0.01, objective_function=compute_drag)

# ϕu, ϕp = solve_inc_adj_s(model, uh,ph, params; filename="adj-$iter")#joinpath("Results","inc-adj-res-iter-$iter")

# J1ref,J2ref= compute_sensitivity(model,params, uh,ph,ϕu, ϕp; objective_function=objective_function)
# modelgrid0 = copy(model.grid.node_coordinates)
# Nc = get_designparameters_number(des_params)
# tag = get_designparameters_tags(des_params)

# vv = (tag .== "top") .- (tag .== "bottom")
# J1=zeros(Nc)
# J2=zeros(Nc)

# shift = vv.*δ
# for (i,ss) in enumerate(shift)
#         des_tmp = perturb_DesignParameters(des_params, i, ss)
#         modelname = create_msh(des_tmp; iter = i+100, mesh_ref=4,AoA=4.0)
#         model_tmp = GmshDiscreteModel(modelname)
#         model.grid.node_coordinates .= model_tmp.grid.node_coordinates 
#        J1tmp,J2tmp=compute_sensitivity(model,params, uh,ph,ϕu, ϕp; objective_function=objective_function)
#        J1[i] = J1tmp
#        J2[i] = J2tmp
#        model.grid.node_coordinates .= modelgrid0
# end



# J1s = (J1 .- J1ref)./ (δ) #Geometric gradient
# J2s =  (J2 .- J2ref)./(δ)
# Jtot = J1s +J2s #Total Gradient

# #Step, -1 because opposite direction of the gradient
# ΔD = ShiftUpdate(-vv.*α.*Jtot, tag)

# contr_new = perturb_DesignParameters(des_params, ΔD)

# modelname = create_msh(contr_new; iter = iter+1, mesh_ref=4, AoA=4.0)
# model = GmshDiscreteModel(modelname)
# model.grid.node_coordinates .= model.grid.node_coordinates 

# uh,ph = solve_inc_primal_s(model, params; filename=joinpath("Results_primal", "res-$(iter+1)"))

# fitnessval, S = obj_fun(model,params,uh, ph,objective_function)

# if detail
#     return model, (J1s,J2s,Jtot), (fitnessval,S), contr_new, (uh,ph)
# else
#     return model, Jtot, S, contr_new, (uh,ph)
# end


# end


