module IncompressibleSolvers


using LinearAlgebra
using Parameters
using PartitionedArrays
using MPI
using SegregatedVMSSolver

using SegregatedVMSSolver.ParametersDef

using VMSAdjoint.Equations
using VMSAdjoint.Interfaces
using VMSAdjoint.ParametersAdj

using Gridap
using Gridap.Algebra
using Gridap.FESpaces
using Gridap.Arrays
using Gridap.CellData

using GridapDistributed
using GridapDistributed.MultiField
using GridapDistributed.CellData


using GridapPETSc

using Statistics
using JLD2


export solve_inc_primal
include("SolvePrimal.jl")

export solve_inc_adj
include("SolveAdjoint.jl")

# export solve_inc_direct_differentiation_u
# export solve_inc_direct_differentiation_s
# include("SolveDirectDifferentiation.jl")

end