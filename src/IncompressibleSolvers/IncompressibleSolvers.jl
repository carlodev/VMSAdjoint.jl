module IncompressibleSolvers

using SegregatedVMSSolver.ParametersDef

using VMSAdjoint
using VMSAdjoint.Equations
using VMSAdjoint.Interfaces
using VMSAdjoint.ParametersAdj

using Gridap
using Gridap.Algebra
using Gridap.FESpaces

using Parameters

using JLD2


export solve_inc_primal
include("SolvePrimal.jl")

# export solve_inc_adj_u
# export solve_inc_adj_s
# include("SolveAdjoint.jl")

# export solve_inc_direct_differentiation_u
# export solve_inc_direct_differentiation_s
# include("SolveDirectDifferentiation.jl")

end