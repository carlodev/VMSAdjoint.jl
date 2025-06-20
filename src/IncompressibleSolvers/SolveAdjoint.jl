# Main entry point
function solve_inc_adj(am::AirfoilModel, simcase::Airfoil, d_bc::Vector{Float64}, timed::Symbol, uh, ph)
    solve_inc_adj(am, simcase,d_bc, Val(timed), uh, ph)
end

function solve_inc_adj(am::AirfoilModel, simcase::Airfoil, d_bc::Vector{Float64}, filename::String, timed::Symbol, uh, ph)
    solve_inc_adj(am, simcase,d_bc, filename, Val(timed), uh, ph)
end

function solve_inc_adj(am::AirfoilModel, simcase::Airfoil, d_bc::Vector{Float64}, filename::String,  ::Val{:steady}, uh, ph)
    solve_inc_adj_steady(am, simcase,d_bc, filename,  uh, ph)
end


function solve_inc_adj(am::AirfoilModel, simcase::Airfoil, d_bc::Vector{Float64}, ::Val{:steady}, uh, ph)
    filename = "inc-adj-steady"
    return solve_inc_adj_steady(am, simcase,d_bc, filename,uh,ph)
end




# function solve_inc_adj(am::AirfoilModel, simcase::Airfoil, ::Val{:unsteady}; uh0=nothing, ph0=nothing)
#     filename = "inc-primal-unsteady"
#     return solve_inc_primal_unsteady(am, simcase, filename, uh0, ph0)
# end




function create_adjoint_spaces(model, simcase::Airfoil)
    @sunpack  order, D = simcase
    reffe_u_adj = ReferenceFE(lagrangian, VectorValue{D,Float64}, order)
    V_adj = TestFESpace(model, reffe_u_adj, conformity=:H1, dirichlet_tags=["airfoil", "outlet","limits"], dirichlet_masks=[(true,true), (true,true),(false,true) ])
    reffe_p_adj = ReferenceFE(lagrangian, Float64, order)
    Q_adj = TestFESpace(model, reffe_p_adj, conformity=:H1, dirichlet_tags="inlet")
   
    return V_adj,Q_adj
end


# function solve_inc_adj_u(model, primal_sol_uh::Tuple, primal_sol_ph::Tuple, adjstart::Tuple, params::Dict{Symbol,Any}; filename="res-adj-unsteady")
#     @unpack D,order,t_endramp,t0,tf,θ,dt,Ω,d_boundary = params

#     uh0, UH = primal_sol_uh
#     ph0, PH = primal_sol_ph

#     ϕuh0, ϕph0 = adjstart
    
#     V_adj,Q_adj = create_adjoint_spaces(model, params)

  

#     uin(t) = (t < t_endramp) ? (1.0 - 1 .*(t_endramp-t)/t_endramp) : 1.0

#     d0(x,t) = d_boundary #D == 2 ? VectorValue(-1 .*uin(t), 0.0) :  VectorValue(-1 .*uin(t), 0.0, 0.0)
#     d0(t::Real) = x -> d0(x,t)
#     u_walls(x,t) = VectorValue(zeros(D)...)
#     u_walls(t::Real) = x -> u_walls(x,t)

#     U_adj = TransientTrialFESpace(V_adj, [d0, u_walls,u_walls])
#     P_adj = TrialFESpace(Q_adj, 0.0)

#     Y_adj = MultiFieldFESpace([V_adj, Q_adj])
#     X_adj = TransientMultiFieldFESpace([U_adj, P_adj])

#     ϕuh0 = interpolate(VectorValue(0, 0), U_adj(0.0),)
#     ϕph0 = interpolate(0.0, P_adj)

#     ϕxh0 = interpolate([ϕuh0, ϕph0], X_adj(0.0))
    
#     Nfields = length(UH)

#     ΦUH = [copy(UH[Nfields])]
#     ΦPH = [copy(PH[Nfields])]

#     copyto!(params[:uh].free_values, UH[Nfields])

#     uh = uh0
#     ph = ph0

#     updatekey(params, :uh,uh)
#     updatekey(params, :ph,ph)


#     m_adj, res_adj, rhs = eq_adjoint_unsteady(params)
#     op = TransientAffineFEOperator(m_adj, res_adj, rhs, X_adj, Y_adj)



#     ls = LUSolver()

#     θ_adj = 0.0
#     ode_solver = ThetaMethodBackw(ls, dt, θ_adj)
  

#     sol_adj = Gridap.solve(ode_solver, op, ϕxh0, t0, tf)

#     res_path = "Results_adj"
#     mkpath(res_path)

#     #Adjoint going backwards
#     createpvd(filename) do pvd
#         for (idx, (ϕxh, t)) in enumerate(sol_adj)
#             ϕuh = ϕxh[1]
#             ϕph = ϕxh[2]


#             pvd[t] = createvtk(Ω, joinpath(res_path, "$(filename)_$t" * ".vtu"), cellfields=["phi-uh" => ϕuh, "phi-ph" => ϕph,
#                 "uh" => uh, "ph" => ph])
#             push!(ΦUH, copy(ϕuh.free_values))
#             push!(ΦPH, copy(ϕph.free_values))

#             IDX = Nfields - idx + 1
#             println("Adjoint solved at time step $t")
#             copyto!(params[:uh].free_values, UH[IDX])
#             copyto!(params[:ph].free_values, PH[IDX])

#         end
#     end

#     return ΦUH, ΦPH
# end


function solve_inc_adj_steady(am::AirfoilModel, simcase::Airfoil,d_boundary::Vector{Float64}, filename, uh,ph)
    
    @unpack params,model = am
    @unpack Ω = params
    @sunpack order = simcase


    Γout = BoundaryTriangulation(model; tags="outlet")
    dΓout = Measure(Γout, order*2)
    nΓout =  get_normal_vector(Γout)

    Γlim = BoundaryTriangulation(model; tags="limits")
    dΓlim = Measure(Γlim, order*2)
    nΓlim =  get_normal_vector(Γlim)

    Γairfoil = BoundaryTriangulation(model; tags="airfoil")
    dΓairfoil = Measure(Γairfoil, order*2)
    nΓairfoil = get_normal_vector(Γairfoil)


    updatekey(params, :dΓout,dΓout)
    updatekey(params, :nΓout,nΓout)

    updatekey(params, :dΓlim,dΓlim)
    updatekey(params, :nΓlim,nΓlim)
    
    updatekey(params, :dΓairfoil,dΓairfoil)
    updatekey(params, :nΓairfoil,nΓairfoil)

    V_adj,Q_adj = create_adjoint_spaces(model, simcase)
    println("Adjoint Boundary condition value: $d_boundary")

    U_adj = TrialFESpace(V_adj, [VectorValue(d_boundary...), VectorValue(0, 0),VectorValue(0, 0)])
    P_adj = TrialFESpace(Q_adj,0.0)

    Y_adj = MultiFieldFESpace([V_adj, Q_adj])
    X_adj = MultiFieldFESpace([U_adj, P_adj])

    updatekey(params, :uh,uh)
    updatekey(params, :ph,ph)

    res_adj, rhs_adj = equations_adjoint(simcase,params,:steady)


    op_adj = AffineFEOperator(res_adj, rhs_adj, X_adj, Y_adj)

    ls = LUSolver()
    solver = LinearFESolver(ls)

    ϕu, ϕp = Gridap.solve(solver, op_adj)
    
    @info "Solving Steady Adjoint ..."
    res_path = "Results_adj"
    mkpath(res_path)

    if !isnothing(filename)
        writevtk(Ω, joinpath(res_path, "$(filename)" * ".vtu"), nsubcells=order, cellfields=["phi-u" => ϕu, "phi-p" => ϕp, "uh" => uh, "ph" => ph])
    end

    return ϕu, ϕp
end
