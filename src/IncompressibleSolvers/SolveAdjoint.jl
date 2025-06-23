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

function solve_inc_adj(am::AirfoilModel, simcase::Airfoil, d_bc::Vector{Float64}, filename::String,  ::Val{:unsteady}, uh, ph)
    solve_inc_adj_unsteady(am, simcase,d_bc, filename,  uh, ph)
end


function solve_inc_adj(am::AirfoilModel, simcase::Airfoil, d_bc::Vector{Float64}, ::Val{:unsteady}, uh, ph)
    filename = "inc-adj-unsteady"
    return solve_inc_adj_unsteady(am, simcase,d_bc, filename,uh,ph)
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



function solve_inc_adj_unsteady(am::AirfoilModel, simcase::Airfoil,d_boundary::Vector{Float64}, filename, uh,ph)
    
    @unpack params,model = am
    @sunpack order = simcase

    @sunpack D,order,t_endramp,t0,tF,θ,dt,u_in,time_window = simcase
    @sunpack M = simcase #here now M is the step to save .vtu files
    @unpack Ω,UH = params


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
    
    
    u0(x,t) = VectorValue(d_boundary...)
    u0(t::Real) = x -> u0(x,t)

    u_walls(x,t) = VectorValue(zeros(D)...) 
    u_walls(t::Real) = x -> u_walls(x,t)

    p0(x,t) = 0.0
    p0(t::Real) = x -> p0(x,t)


    U_adj = TransientTrialFESpace(V_adj, [u0, u_walls,u_walls])
    P_adj = TransientTrialFESpace(Q_adj,p0)

    Y_adj = MultiFieldFESpace([V_adj, Q_adj])
    X_adj = MultiFieldFESpace([U_adj, P_adj])

    uh0_adj = interpolate(u0(0.0), U_adj(0.0))
    ph0_adj = interpolate(p0(0.0), P_adj(0.0))
    copyto!(am.params[:uh].free_values, UH[end])


    xh0_adj = interpolate([uh0_adj, ph0_adj], X_adj(0.0))


    m_adj, res_adj, rhs_adj =  equations_adjoint( simcase, am.params,:unsteady)

    op_adj = TransientLinearFEOperator((res_adj, m_adj), rhs_adj, X_adj, Y_adj)

    ls = LUSolver()

    println("Solve Unsteady Adjoint")

    ode_solver = ThetaMethod(ls,dt,θ)

    sol = Gridap.solve(ode_solver, op_adj, t0, tF, xh0_adj)

    UH_ADJ = [copy(uh0_adj.free_values)]
    PH_ADJ = [copy(ph0_adj.free_values)]
    

    res_path = "Results_unsteady_primal"
    mkpath(res_path)

    time_vec = collect(t0:dt:tF)
    t_length = length(time_vec)

    createpvd(filename) do pvd
        pvd[t_length] = createvtk(Ω, nsubcells = order, joinpath(res_path, "$(filename)_$tF" * ".vtu"), cellfields=["uh-adj" => uh0_adj, "ph-adj" => ph0_adj])
        for (idx,(t, xhtn)) in enumerate(sol)
            ϕu = xhtn[1]
            ϕp = xhtn[2]
            
            push!(UH_ADJ, copy(ϕu.free_values))
            push!(PH_ADJ, copy(ϕp.free_values))

            idx_adj = Int(t_length - idx)
            t_adj = time_vec[idx_adj]

            idx_adj > 0 ? copyto!(am.params[:uh].free_values,UH[idx_adj]) : nothing

            println("Adjoint solved at time step $t_adj")
                   

            if true #mod(idx,M)==0
                pvd[t_adj] = createvtk(Ω, nsubcells = order, joinpath(res_path, "$(filename)_$t_adj" * ".vtu"), cellfields=["uh-adj" => ϕu, "ph-adj" => ϕp])
            end

        end
    end

    jldsave("UnsteadyAdjointFields.jld2"; UH_ADJ,PH_ADJ)

    time_window_adj = (time_window[1], tF - time_window[1])

    @assert time_window_adj[2] > time_window_adj[1] "Adjoint Time Window Averaging not consistent"

    avg_UH_ADJ,avg_PH_ADJ = time_average_fields(UH_ADJ,PH_ADJ, time_window_adj,dt, t0)

    uh0_adj.free_values .=  avg_UH_ADJ
    ph0_adj.free_values .=  avg_PH_ADJ
    

    return uh0_adj, ph0_adj
end
