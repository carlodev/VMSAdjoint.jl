

function solve_adjoint_optimization(adjp::AdjointProblem)
    @unpack solver = adjp

    #create initial
    w_init = get_DesignParameters(adjp.adesign)

    #create bounds
    Ndes = length(w_init)
    lb,ub = bounds_w(adjp.adesign,  Ndes)


    f, ∇f! = make_f_and_∇f(adjp, Ndes)


    ls = LineSearches.BackTracking(
        c_1 = 1e-3,
        ρ_hi = 0.9,
        ρ_lo = 0.25,
        iterations = 1,
        order = 3,
        maxstep = 4e-2,  # ⬅️ controls max step size
        cache = nothing
    )    
    # L-BFGS optimizer with line search control
    result = optimize(f, ∇f!,lb,ub, w_init, Fminbox(LBFGS(linesearch = ls);))
    
    return true
end








function bounds_w(adesign::RBFDesign, Ndes::Int64)
    Nhalf = Int(Ndes ÷ 2)
    Δy = 0.005
    lb = [-Δy .* ones(Nhalf);-0.25.* ones(Nhalf)]
    ub =[0.25 .* ones(Nhalf); Δy.* ones(Nhalf)]
    lb[1] = Δy
    ub[Nhalf+1] = -Δy

    # lb = -0.25 .* ones(Ndes)
    # ub = 0.25 .* ones(Ndes)

    
    
    [@assert ub1>lb1 for (ub1,lb1) in zip(ub,lb)]

    return lb,ub
end


function bounds_w(adesign::AirfoilCSTDesign,  Ndes::Int64)
    Nhalf = Int(Ndes ÷ 2)
    nose = 0.08
    lb = [nose; nose; nose; -0.15.* ones(Nhalf-3); -1.5 .* ones(Nhalf)]
    ub = [1.5 .* ones(Nhalf); nose; -nose;-nose; 0.15 .* ones(Nhalf-3)]
    return lb,ub
end


mutable struct SharedCache
    fval::Float64
    grad::Vector{Float64}
    valid::Bool
    iter::Int64
    uh
    ph
    adj_bc::Vector{Float64}
    am
    Cp
    adjp::AdjointProblem
    CDCL::Vector{Float64}
end

function make_f_and_∇f(adjp::AdjointProblem, N::Int64)
    cache = SharedCache(NaN, zeros(N), true, -1, nothing, nothing, [0.0,0.0], nothing,nothing,adjp, [0.0,0.0])

    x_last = similar(zeros(N))

    f(x) = begin
        if !isequal(x, x_last)
            copy!(x_last, x)
            cache.fval, cache.iter, cache.uh, cache.ph, cache.adj_bc, cache.am, cache.Cp, cache.CDCL = eval_f(x, cache)  # compute both
        end
        return cache.fval
    end

    ∇f!(g, x) = begin
        println("g start = $g")

        if !isequal(x, x_last)
            copy!(x_last, x)
            cache.fval, cache.iter, cache.uh, cache.ph, cache.adj_bc, cache.am, cache.Cp = eval_f(x, cache)
        end
        eval_∇f!(g, x, cache)  # must recompute if f not called

        copyto!(cache.grad, g)
        println("g end = $g")

    end

    return f, ∇f!
end



function eval_f(w::Vector, cache::SharedCache)

    @unpack  iter, uh, ph, adjp, am= cache
    @unpack J,adesign, vbcase, timesol,solver = adjp

    println("w = $w")

    meshinfo = vbcase.meshp.meshinfo
    physicalp = vbcase.simulationp.physicalp


    iter = iter+1

    @info "Iteration $iter started"

    #create the new airfoil model from the weights w
    adesign = create_AirfoilDesign(adesign,w)
    modelname =create_msh(meshinfo,adesign, physicalp ; iter = iter)
    model = GmshDiscreteModel(modelname)
    writevtk(model, "model_$iter")
    am =  AirfoilModel(model, vbcase; am=am)


    #Solve Primal 
    if timesol==:steady
        filename = joinpath("Results_primal", "SOL_$(iter).vtu")
    else
           filename = "sol_$(iter)"
    end

    uh,ph = solve_inc_primal(am, vbcase, filename, timesol; uh0=uh,ph0=ph)    
 
    #### extract results from primal solution: Cp (and Cf)
    PressureCoefficient = get_aerodynamic_features(am,uh,ph)


    fval, CDCL = obj_fun(am, vbcase, uh,ph, J)
 
    #### Adjoint Boundary Conditions
    #use the CLCD value to set the Boundary Condition
    adj_bc = -dJobj_fun(J, CDCL)
    

    return fval, iter, uh, ph, adj_bc, am, PressureCoefficient, CDCL
end


function eval_∇f!(grad::Vector, w::Vector,  cache::SharedCache)
    @unpack  iter, uh, ph, adjp, am, adj_bc= cache
    @unpack J,adesign, vbcase, timesol,solver = adjp

    airfoil_case= vbcase
    uhadj,phadj = solve_inc_adj(am, airfoil_case, adj_bc, "inc-adj-steady-$iter", :steady, uh, ph)


    J1ref,(J2_1ref,J2_2ref,J2_3ref,J2_4ref,J2_5ref)= compute_sensitivity(airfoil_case, am,adj_bc, uh,ph,uhadj,phadj, J; i=iter)


    
    Ndes = length(w) #number of design parameters
    δ = solver.δ #0.01
    shift = CSTweights(Int(Ndes/2), δ)
    shiftv =   vcat(shift)

    
    J1, (J2_1,J2_2,J2_3,J2_4,J2_5) = iterate_perturbation(shiftv,adesign,am, airfoil_case,adj_bc,uh,ph,uhadj,phadj, J  )

    # shiftv = δ
    J1s = (J1 .- J1ref)./shiftv
    J2s1 = (J2_1 .- J2_1ref)./shiftv
    J2s2 = (J2_2 .- J2_2ref)./shiftv
    J2s3 = (J2_3 .- J2_3ref)./shiftv
    J2s4 = (J2_4 .- J2_4ref)./shiftv
    J2s5 = (J2_5 .- J2_5ref)./shiftv
    
    
    J2s = J2s1 + J2s2 + J2s3 + J2s4 + J2s5
    Jtot = J1s .+ J2s
    
    println("Geometric Gradient:")
    println(J1s)

    println("Adjoint Gradient:")
    println(J2s1)
    println(J2s2)
    println(J2s3)
    println(J2s4)
    println(J2s5)
    JtotD = Dict(:J1s=>J1s, :J2s1=>J2s1 , :J2s2=>J2s2 , :J2s3=>J2s3, :J2s4=>J2s4, :J2s5=>J2s5 )
    
    jldsave("results/J$(iter).jld2"; JtotD)

    grad[:] = Jtot 


    #update values iteration


    adj_sol = AdjSolution(iter, cache.fval, cache.CDCL, w, grad, am, adj_bc, cache.Cp, uh, ph )


    jldsave("results/ADJ_SOL$(iter).jld2"; adj_sol)
    @info "Iteration $iter completed"


end




function iterate_perturbation(shift::Vector{Float64},adesign::AirfoilDesign,am::AirfoilModel,airfoil_case::Airfoil,adj_bc::Vector{Float64}, uh,ph,uhadj,phadj, J::Function  )
    Ndes = length(shift)

    meshinfo = airfoil_case.meshp.meshinfo
    physicalp =airfoil_case.simulationp.physicalp
    
    J1=zeros(Ndes)    #Geometric term
    J2_1=zeros(Ndes)
    J2_2=zeros(Ndes)
    J2_3=zeros(Ndes)
    J2_4=zeros(Ndes)
    J2_5=zeros(Ndes)

    for (i,ss) in enumerate(shift)
        println("Perturbation Domain $i")

        adesign_tmp = perturb_DesignParameter(adesign, i, ss)

        modelname_tmp =create_msh(meshinfo,adesign_tmp, physicalp,"MeshPerturb"; iter = i+100)
        model_tmp = GmshDiscreteModel(modelname_tmp)
        am_tmp =  AirfoilModel(model_tmp, airfoil_case, am=am)

        J1tmp,(J2_1tmp,J2_2tmp,J2_3tmp,J2_4tmp,J2_5tmp)= compute_sensitivity(airfoil_case, am_tmp,adj_bc, uh,ph,uhadj,phadj, J; i=i+100)

        J1[i] = J1tmp
        J2_1[i] = J2_1tmp
        J2_2[i] = J2_2tmp
        J2_3[i] = J2_3tmp
        J2_4[i] = J2_4tmp
        J2_5[i] = J2_5tmp

    end
    
    return J1, (J2_1,J2_2,J2_3,J2_4,J2_5)
end
