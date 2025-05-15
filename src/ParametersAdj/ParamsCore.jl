
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




@with_kw mutable struct AdjBC
    tagname::String="airfoil" #for multiple tags, we can put this as Vector{String}
    bc::Vector #Adjoint BC on the tagname; it depends on what we want to optimize
end

@with_kw struct AdjSolver
    max_iter::Int64 = 10 #maximum number of adjoint iterations
    tol::Float64 = 2.5e-2 #tolerance convergence
    opt_alg::Symbol = :LD_LBFGS #OPTIMIZATION algorithm
end

@with_kw struct AirfoilMesh <:MeshInfo
    AoA::Real #Angle of Attack - degrees
    cstdesign::AirfoilCSTDesign
    meshref::Int64=1
    folder::String="MeshFiles"
end

struct AdjointProblem
    bc::AdjBC
    vbcase::VelocityBoundaryCase
    solver::AdjSolver
end

function AdjointProblem(    bc::AdjBC,vbcase::VelocityBoundaryCase)
    solver=AdjSolver()
return     AdjointProblem(    bc,vbcase,solver)
end

"""
    AdjSolution

It stores the results of the adjoint iteration.
fitnessval : value of the fitness function, it should go towards zero as the iterations are increasing
βv : value of the design parameters (e.g. CST weights) - it has been kept general
fgrad : ∂J/∂βv derivative of the objective function with respect to the design variable βv
airfoil_model : we store the arifoil model at this
"""
@with_kw mutable struct AdjSolution
    iter::Int64 = 0
    i::Int64=0
    fitnessval::Vector{Float64}
    βv::Vector{Float64}
    fgrad::Vector{Float64}
    airfoil_model::AirfoilModel
    pressure_distribution::AirfoilScalar
end

function AdjSolution(adjp::AdjointProblem, am::AirfoilModel)
    max_iter = adjp.solver.max_iter
    nβv = length(length(adjp.vbcase.meshp.meshinfo.cstdesign.cstg.cstw))
    fitnessval = zeros(max_iter)
    βv = zeros(nβv)
    fgrad = zeros(nβv)
    pressure_distribution = AirfoilScalar(am)
    AdjSolution(fitnessval=fitnessval,βv=βv,fgrad=fgrad,airfoil_model=am, pressure_distribution=pressure_distribution)

end

struct FieldsSol
    U
    V
    P
    Q
    uh_vec::Vector{VectorValue}
    ph_vec::Vector{Float64}
end




