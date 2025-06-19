# Main entry point
function solve_inc_primal(am::AirfoilModel, simcase::Airfoil, timed::Symbol; uh0=nothing, ph0=nothing)
    solve_inc_primal(am, simcase, Val(timed); uh0=uh0, ph0=ph0)
end

function solve_inc_primal(am::AirfoilModel, simcase::Airfoil, filename::String, timed::Symbol; uh0=nothing, ph0=nothing)
    solve_inc_primal(am, simcase, filename, Val(timed); uh0=uh0, ph0=ph0)
end

function solve_inc_primal(am::AirfoilModel, simcase::Airfoil,  filename::String, ::Val{:steady}; uh0=nothing, ph0=nothing)
    return solve_inc_primal_steady(am, simcase, filename,uh0,ph0)
end

function solve_inc_primal(am::AirfoilModel, simcase::Airfoil,  filename::String, ::Val{:unsteady}; uh0=nothing, ph0=nothing)
    return solve_inc_primal_unsteady(am, simcase, filename, uh0, ph0)
end



function solve_inc_primal(am::AirfoilModel, simcase::Airfoil, ::Val{:steady}; uh0=nothing, ph0=nothing)
    filename = "inc-primal-steady"
    return solve_inc_primal_steady(am, simcase, filename,uh0,ph0)
end

function solve_inc_primal(am::AirfoilModel, simcase::Airfoil, ::Val{:unsteady}; uh0=nothing, ph0=nothing)
    filename = "inc-primal-unsteady"
    return solve_inc_primal_unsteady(am, simcase, filename, uh0, ph0)
end


function rotation(n::VectorValue{2,Float64})
    n1, n2 = [n...] ./ norm(n)
    VectorValue(-n2, n1)
end

function create_primal_spaces(model, simcase::Airfoil)
    @sunpack tagname, order, D = simcase
    reffeᵤ = ReferenceFE(lagrangian, VectorValue{D,Float64}, order )
    V = TestFESpace(model, reffeᵤ, conformity=:H1, dirichlet_tags=["inlet", "limits", "airfoil"],  dirichlet_masks=[(true,true), (false,true),(true,true) ])
    reffeₚ = ReferenceFE(lagrangian, Float64, order)
    Q = TestFESpace(model, reffeₚ, conformity=:H1, dirichlet_tags=["outlet"])

    return V,Q
end

function solve_inc_primal_unsteady(am::AirfoilModel, simcase::Airfoil, filename, uh00, ph00)
    
    @sunpack D,order,t_endramp,t0,tF,θ,dt,u_in,time_window = simcase
    @sunpack M = simcase #here now M is the step to save .vtu files
    @unpack model = am

    V,Q = create_primal_spaces(model,simcase)

    uin(t) = (t < t_endramp) ? (1.0 - 1 .*(t_endramp-t)/t_endramp) : 1.0

    u0(x,t) = VectorValue(u_in...)
    u0(t::Real) = x -> u0(x,t)

    u_walls(x,t) = VectorValue(zeros(D)...) 
    u_walls(t::Real) = x -> u_walls(x,t)


    p0(x,t) = 0.0
    p0(t::Real) = x -> p0(x,t)

    U = TransientTrialFESpace(V, [u0, u0, u_walls])
    P = TransientTrialFESpace(Q, p0)

    Y = TransientMultiFieldFESpace([V, Q])
    X = TransientMultiFieldFESpace([U, P])

    setup_shape_opt_spaces!(D, order, model, am)

    
    degree = order * 2
    Ω = Triangulation(model)
    dΩ = Measure(Ω, degree)
    updatekey(am.params,:Ω,Ω)
    updatekey(am.params,:dΩ,dΩ)
    
    uh0 = interpolate(u0(0.0), U(0.0))
    ph0 = interpolate(p0(0.0), P(0.0))

    if isnothing(uh00) && t_endramp ==  t0
        #initialize with steady solution
        uh00,ph00 = solve_inc_primal_steady(am, simcase, nothing, uh00,uh00)
        uh0.free_values .=  uh00.free_values
        ph0.free_values .=  ph00.free_values
    end




    xh0 = interpolate([uh0, ph0], X(0.0))

    updatekey(am.params, :uh,uh0)
    m, res, rhs =  equations_primal( simcase, am.params,:unsteady)

    op = TransientLinearFEOperator((res, m), rhs, X, Y)

    ls = LUSolver()

    ode_solver = ThetaMethod(ls,dt,θ)

    sol = Gridap.solve(ode_solver, op, t0, tF, xh0)

    UH = [copy(uh0.free_values)]
    PH = [copy(ph0.free_values)]
    

    res_path = "Results_unsteady_primal"
    mkpath(res_path)

    createpvd(filename) do pvd
        pvd[0] = createvtk(Ω, nsubcells = order, joinpath(res_path, "$(filename)_0" * ".vtu"), cellfields=["uh" => uh0, "ph" => ph0])
        for (idx,(t, xhtn)) in enumerate(sol)
            uh = xhtn[1]
            ph = xhtn[2]
            push!(UH, copy(uh.free_values))
            push!(PH, copy(ph.free_values))
            println("Primal solved at time step $t")
                 
            copyto!(am.params[:uh].free_values,uh.free_values)
           
            if mod(idx,M)==0
                pvd[t] = createvtk(Ω, nsubcells = order, joinpath(res_path, "$(filename)_$t" * ".vtu"), cellfields=["uh" => uh, "ph" => ph])
            end

        end
    end

    avg_UH,avg_PH = time_average_fields(UH,PH, time_window,dt, t0)

    uh0.free_values .=  avg_UH
    ph0.free_values .=  avg_PH
    
    return uh0,ph0

end


function time_average_fields(UH,PH, time_window,dt::Float64, t0::Float64)
    # Number of stored snapshots (including t0)
    nt = length(UH)
    
    #Time to start the averaging
    t_0 = time_window[1]

    # Reconstruct time values (t0, t0+dt, ..., tF)
    times = t0:dt:(t0 + dt*(nt-1))

    # Find indices where t ≥ t01
    indices = findall(t -> t ≥ t_0, times)

    # Average the snapshots over the selected indices
    avg_UH = Statistics.mean(UH[indices])
    avg_PH =  Statistics.mean(PH[indices])
    return avg_UH,avg_PH

end

function solve_inc_primal_steady(am::AirfoilModel, simcase::Airfoil, filename, uh00,ph00)
    @sunpack D,order,t_endramp,t0,tf,θ,dt,u_in,petsc_options = simcase
    @sunpack matrix_freq_update = simcase #number of iterations in steady
    @unpack model = am
    V,Q = create_primal_spaces(model,simcase)

    n_inner_iter = matrix_freq_update

    u0 = VectorValue(u_in...)
    u_walls = VectorValue(zeros(D)...) 
    p0 = 0.0

    U = TrialFESpace(V, [u0, u0, u_walls])
    P = TrialFESpace(Q, [p0])
    
    updatekey(am.params,:U,U)
    updatekey(am.params,:P,P)

    Y = MultiFieldFESpace([V, Q])
    X = MultiFieldFESpace([U, P])
    
    setup_shape_opt_spaces!(D, order, model, am)

    degree = order*2
    Ω = Triangulation(model)
    dΩ = Measure(Ω, degree)
    updatekey(am.params,:Ω,Ω)
    updatekey(am.params,:dΩ,dΩ)



    ls = LUSolver()

    solver = LinearFESolver(ls)

    uh = interpolate(u0, U)
    ph = interpolate(0.0, P)
 
    if !isnothing(uh00)
        uh.free_values .=  uh00.free_values
     end
     
     if !isnothing(ph00)
         ph.free_values .=  ph00.free_values
     end

    for i = 1:n_inner_iter
        println(i)
        uh,ph = solve_steady_primal(uh,ph,X,Y,simcase, am.params,solver )
    end



    if !isnothing(filename)
        writevtk(Ω, nsubcells=order, filename,  cellfields=["uh" => uh, "ph" => ph])
    end
    

    return uh,ph
end

function solve_steady_primal(uh,ph,X,Y,simcase, params,solver )
    xh = interpolate([uh, ph], X)

    updatekey(params, :uh,uh)
    res, rhs = equations_primal( simcase,params,:steady)

    op = AffineFEOperator(res, rhs, X, Y)

  
    @info "solving steady problem..."

    # Gridap.solve!(xh, solver, op)
    Gridap.solve!(xh, solver, op)

    uh, ph = xh 
    return  uh, ph
end


function setup_shape_opt_spaces!(D::Int64, order::Int64, model, am::AirfoilModel)
    reffe = ReferenceFE(lagrangian, VectorValue{D,Float64}, order )
    VV0 = FESpace(model, reffe; conformity=:H1)

    Γ = BoundaryTriangulation(am.model; tags="airfoil")
    dΓ = Measure(Γ, order*2)
    nΓ = -get_normal_vector(Γ) #beacuse they point inward 
    tΓ = rotation ∘ nΓ

    updatekey(am.params,:reffe,reffe)
    updatekey(am.params,:VV0,VV0)
    updatekey(am.params,:Γ,Γ)
    updatekey(am.params,:dΓ,dΓ)
    updatekey(am.params,:nΓ,nΓ)
    updatekey(am.params,:tΓ,tΓ)
    
end