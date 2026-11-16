module IncompressibleSolvers

using SegregatedVMSSolver.ParametersDef

using VMSAdjoint
using VMSAdjoint.Equations
using VMSAdjoint.Interfaces
using VMSAdjoint.ParametersAdj

using Gridap
using Gridap.Algebra
using Gridap.FESpaces
using Gridap.MultiField: BlockMultiFieldStyle, ConsecutiveMultiFieldStyle

using GridapPETSc

using GridapSolvers
using GridapSolvers.BlockSolvers: BlockTriangularSolver, BiformBlock, LinearSystemBlock
using GridapSolvers.LinearSolvers: FGMRESSolver

using Statistics
using Parameters
using LinearAlgebra: norm, lu

using JLD2


export GPULUSolver
export SolverContext
export linear_solver
export with_petsc
include("LinearSolverBackends.jl")

export march_frozen_theta
include("TransientSolve.jl")

export solve_inc_primal
export tangential_normal_derivative
include("SolvePrimal.jl")

export solve_inc_adj
include("SolveAdjoint.jl")

# export solve_inc_direct_differentiation_u
# export solve_inc_direct_differentiation_s
# include("SolveDirectDifferentiation.jl")

end