"""
    VMSAdjointCUDSSExt

Package extension implementing [`VMSAdjoint.GPULUSolver`](@ref) on top of
NVIDIA cuDSS (via CUDSS.jl). Loaded automatically when `CUDA` and `CUDSS`
are present alongside `VMSAdjoint`.

The Gridap-assembled `SparseMatrixCSC` is converted to a device CSR matrix;
cuDSS performs the sparsity analysis once and only refactorizes numerically
when `numerical_setup!` is called with new matrix values (same sparsity
pattern — guaranteed here because Gridap reassembles on a fixed mesh/space).
"""
module VMSAdjointCUDSSExt

using VMSAdjoint
using VMSAdjoint.IncompressibleSolvers: GPULUSolver
using Gridap
using Gridap.Algebra
using SparseArrays
using LinearAlgebra
using CUDA
using CUDA.CUSPARSE
using CUDSS

struct GPULUSymbolicSetup <: Gridap.Algebra.SymbolicSetup end

Gridap.Algebra.symbolic_setup(::GPULUSolver, ::AbstractMatrix) = GPULUSymbolicSetup()

mutable struct GPULUNumericalSetup{S<:CudssSolver,V<:CuVector} <: Gridap.Algebra.NumericalSetup
    solver::S
    x_gpu::V
    b_gpu::V
end

function Gridap.Algebra.numerical_setup(::GPULUSymbolicSetup, mat::AbstractMatrix)
    A_gpu = CuSparseMatrixCSR(SparseMatrixCSC{Float64,Int}(mat))
    n = size(mat, 1)
    x_gpu = CUDA.zeros(Float64, n)
    b_gpu = CUDA.zeros(Float64, n)

    solver = CudssSolver(A_gpu, "G", 'F') # general (non-symmetric) matrix, full storage
    cudss("analysis", solver, x_gpu, b_gpu)
    cudss("factorization", solver, x_gpu, b_gpu)

    return GPULUNumericalSetup(solver, x_gpu, b_gpu)
end

function Gridap.Algebra.numerical_setup!(ns::GPULUNumericalSetup, mat::AbstractMatrix)
    A_gpu = CuSparseMatrixCSR(SparseMatrixCSC{Float64,Int}(mat))
    cudss_set(ns.solver, A_gpu)
    cudss("refactorization", ns.solver, ns.x_gpu, ns.b_gpu)
    return ns
end

function Gridap.Algebra.solve!(x::AbstractVector, ns::GPULUNumericalSetup, b::AbstractVector)
    copyto!(ns.b_gpu, b)
    cudss("solve", ns.solver, ns.x_gpu, ns.b_gpu)
    copyto!(x, ns.x_gpu)
    return x
end

end
