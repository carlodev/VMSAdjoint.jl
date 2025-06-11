
struct AirfoilNormals
    nu::Vector{Vector{Float64}}
    nl::Vector{Vector{Float64}}
end

struct AirfoilModel
    ap::AirfoilPoints
    model::Gridap.Geometry.UnstructuredDiscreteModel
    an::AirfoilNormals
    params::Dict{Symbol,Any}
end

struct AirfoilScalar
    uval::Vector{Float64}
    lval::Vector{Float64}
end


function AirfoilScalar(am::AirfoilModel)
    nu = length(am.ap.xu)
    nl = length(am.ap.xl)
    return AirfoilScalar(zeros(nu), zeros(nl))
end

@with_kw struct DesignBounds
    upper::Float64=0.5
    lower::Float64=-0.5
    Δy::Float64=0.005
end


@with_kw struct AdjSolver
    max_iter::Int64 = 10 #maximum number of adjoint iterations
    tol::Float64 = 2.5e-2 #tolerance convergence
    αg::Float64 = 2.0 #alphaguess, reduce if the geometry is chaning too fast
    δ::Float64=0.0001 #Perturbation of the design parameters for finite differences
    bounds::DesignBounds = DesignBounds()
end


@with_kw struct AirfoilMesh <:MeshInfo
    AoA::Real #Angle of Attack - degrees
    meshref::Int64=1
    folder::String="MeshFiles"
end

struct AdjointProblem
    adesign::AirfoilDesign
    vbcase::VelocityBoundaryCase
    solver::AdjSolver
    timesol::Symbol
    JJfact::Tuple #objective function
    function AdjointProblem(   adesign::AirfoilDesign,vbcase::VelocityBoundaryCase,solver::AdjSolver,timesol::Symbol,JJfact::Tuple)
        if timesol ∉ (:steady, :unsteady)
            throw(ArgumentError("Invalid timesol: $timesol. Must be :steady, :unsteady"))
        end
        new(adesign, vbcase,solver,timesol,JJfact)
    end
end

function AdjointProblem(    adesign::AirfoilDesign,vbcase::VelocityBoundaryCase,timesol::Symbol, JJfact::Tuple)
    solver=AdjSolver()
return     AdjointProblem(    adesign,vbcase,solver,timesol,JJfact)
end

"""
    AdjSolution

It stores the results of the adjoint iteration.
fitnessval : value of the fitness function, it should go towards zero as the iterations are increasing
βv : value of the design parameters (e.g. CST weights) - it has been kept general
fgrad : ∂J/∂βv derivative of the objective function with respect to the design variable βv
airfoil_model : we store the arifoil model at this
"""
struct AdjSolution
    iter::Int64
    fitnessval::Float64 
    CDCL::Vector{Float64}
    βv::Vector{Float64} 
    fgrad::Vector{Float64} 
    airfoil_model::AirfoilModel
    bc_adj::Vector{Float64} 
    pressure_distribution::AirfoilScalar
    uh
    ph
    uhadj
    phadj
end



# AdjSolution(; iter=-1, fitnessval=Inf, CDCL=[0.0,0.0],βv=Float64[],fgrad=Float64[], bc_adj=[0.0,0.0]; pressure_distribution=AirfoilScalar() )





