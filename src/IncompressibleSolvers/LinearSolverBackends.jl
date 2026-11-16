#############################################################################
# Linear-solver backend selection
#############################################################################
#
# Four backends, selected by `AdjSolver(linear_backend = ...)`:
#
#   :cpu    UMFPACK LU (Gridap `LUSolver`). Direct, serial, robust. Default.
#   :gpu    cuDSS LU on an NVIDIA GPU (package extension, needs `using CUDA, CUDSS`).
#   :petsc  PETSc KSP driven by `SolverParameters.petsc_options`. Requires the
#           whole solve to run inside `with_petsc(...)`.
#   :block  GridapSolvers FGMRES preconditioned by a block-triangular
#           velocity/pressure splitting with a Cahouet-Chabard pressure Schur
#           approximation. Serial today, distributed-ready.
#
# Backends that need to know about the FE spaces (`:block`) receive a
# `SolverContext`; the others ignore it.
#############################################################################

"""
    SolverContext(; P, Q, dΩ, ν, τ)

Everything a linear-solver backend may need beyond the assembled matrix:

- `P`, `Q`  — pressure trial/test space (to assemble the Schur approximation)
- `dΩ`      — volume measure
- `ν`       — kinematic viscosity
- `τ`       — effective time scale of the stage operator: `θ*dt` for a transient
              solve, `Inf` for a steady one (see [`schur_biform`](@ref)).

The default (all-`nothing`) context is what the backends that ignore it get.
"""
@with_kw struct SolverContext
    P::Any = nothing
    Q::Any = nothing
    dΩ::Any = nothing
    ν::Float64 = 0.0
    τ::Float64 = Inf
end

"Backends whose matrices must be assembled block-wise (velocity/pressure blocks)."
needs_block_assembly(backend::Symbol) = backend === :block

"""
    multifield_style(backend)

`MultiFieldStyle` for the velocity/pressure product space: block-wise for the
`:block` backend (GridapSolvers block preconditioners index the blocks
directly), Gridap's default interleaved layout otherwise.
"""
multifield_style(backend::Symbol) =
    needs_block_assembly(backend) ? BlockMultiFieldStyle() : ConsecutiveMultiFieldStyle()

#############################################################################
# GPU (cuDSS)
#############################################################################

"""
    GPULUSolver()

Gridap `LinearSolver` that factorizes the system on an NVIDIA GPU with the
cuDSS sparse direct solver (LU for general matrices). The assembled CPU matrix
is copied to the device at each (re)factorization; the sparsity analysis is
reused across refactorizations, so transient solves only pay the numerical
refactorization.

The implementation lives in the `VMSAdjointCUDSSExt` package extension:
it activates when `CUDA` and `CUDSS` are loaded, i.e.

    using CUDA, CUDSS
    using VMSAdjoint

Without them, any attempt to solve with this backend raises an informative error.
"""
struct GPULUSolver <: Gridap.Algebra.LinearSolver end

function Gridap.Algebra.symbolic_setup(::GPULUSolver, ::AbstractMatrix)
    error("GPULUSolver requires the CUDSS extension: run `using CUDA, CUDSS` " *
          "(with a CUDA-capable GPU) before solving with linear_backend = :gpu.")
end

#############################################################################
# PETSc
#############################################################################

"""
    with_petsc(f, simcase::Airfoil)

Run `f()` inside a PETSc session configured with
`simcase.simulationp.solverp.petsc_options`. Mandatory around any solve that
uses `linear_backend = :petsc`, because `PETScLinearSolver` may only be built
and used while PETSc is initialized:

    with_petsc(simcase) do
        solve_adjoint_optimization(adjp)
    end

The option string is the same one `SegregatedVMSSolver` uses, so options tuned
there transfer unchanged.
"""
function with_petsc(f::Function, simcase::Airfoil)
    @sunpack petsc_options = simcase
    GridapPETSc.with(args=split(petsc_options)) do
        f()
    end
end

#############################################################################
# Block-preconditioned Krylov (GridapSolvers)
#############################################################################

"""
    schur_biform(ctx::SolverContext)

Cahouet-Chabard-style approximation of the pressure Schur complement
`S = B A⁻¹ Bᵀ` of the velocity/pressure system, as a single bilinear form:

    ŝ(p,q) = τ ∫ ∇p⋅∇q  +  (1/ν) ∫ p q

The two limits it interpolates between:

- `τ = θ·dt → 0` (small time step): the velocity block is mass-dominated,
  `A ≈ M/(θ dt)`, hence `S ≈ θ dt · Lp` — the pressure Laplacian term.
- `τ → ∞` (steady problem): the velocity block is diffusion-dominated,
  `A ≈ νK`, hence `S ≈ Mp/ν` — the pressure mass term.

Summing the two operators (rather than their inverses, as in the original
Cahouet-Chabard preconditioner) costs one solve per application instead of two
and is the variant normally used in practice; it is spectrally equivalent to
`S` uniformly in `dt` for Stokes and degrades gracefully with Reynolds number.
"""
function schur_biform(ctx::SolverContext)
    @unpack dΩ, ν, τ = ctx
    # a steady solve has no mass term: drop the Laplacian contribution rather
    # than multiplying it by Inf
    τ_eff = isfinite(τ) ? τ : 0.0
    return (p, q) -> ∫(τ_eff * (∇(p) ⋅ ∇(q)) + (1 / ν) * (p * q))dΩ
end

"""
    block_schur_solver(ctx::SolverContext; kwargs...)

FGMRES preconditioned by an upper block-triangular velocity/pressure solver:

    │ A  B │        │ Ã  B │⁻¹
    │ Bᵀ C │   ≈    │ 0  S̃ │

with `Ã` the assembled velocity block (solved by `solver_u`) and `S̃` the
[`schur_biform`](@ref) approximation of the pressure Schur complement (solved
by `solver_p`). The off-diagonal `B` block is taken straight from the assembled
system (`LinearSystemBlock`), so the preconditioner tracks the current
linearization for free.

Both block solves default to `LUSolver()`. That already beats the monolithic LU
— each block is much smaller and far better conditioned than the saddle-point
system — while keeping the method as robust as a direct solver. For the
distributed case, substituting `PETScLinearSolver()` or `GMGLinearSolver` per
block is the only change needed; the block structure is unchanged.
"""
function block_schur_solver(ctx::SolverContext;
    solver_u=LUSolver(), solver_p=LUSolver(),
    m::Int=30, rtol::Float64=1e-8, maxiter::Int=200, verbose::Bool=false)

    isnothing(ctx.Q) && error(
        "linear_backend = :block needs a SolverContext carrying the pressure spaces; " *
        "it cannot be built from `linear_solver(:block)` alone.")

    blocks = [LinearSystemBlock() LinearSystemBlock();
        LinearSystemBlock() BiformBlock(schur_biform(ctx), ctx.P, ctx.Q)]
    precond = BlockTriangularSolver(blocks, [solver_u, solver_p])

    return FGMRESSolver(m, precond; rtol=rtol, maxiter=maxiter, verbose=verbose)
end

#############################################################################
# Dispatch
#############################################################################

const LINEAR_BACKENDS = (:cpu, :gpu, :petsc, :block)

"""
    linear_solver(backend::Symbol, [ctx::SolverContext])

Build the Gridap `LinearSolver` of the requested backend:
- `:cpu`   -> `LUSolver()` (UMFPACK, default)
- `:gpu`   -> [`GPULUSolver`](@ref) (NVIDIA cuDSS; requires `using CUDA, CUDSS`)
- `:petsc` -> `PETScLinearSolver()` (requires an enclosing [`with_petsc`](@ref))
- `:block` -> [`block_schur_solver`](@ref) (needs a populated `ctx`)
"""
linear_solver(backend::Symbol) = linear_solver(backend, SolverContext())

function linear_solver(backend::Symbol, ctx::SolverContext)
    backend === :cpu && return LUSolver()
    backend === :gpu && return GPULUSolver()
    backend === :petsc && return GridapPETSc.PETScLinearSolver()
    backend === :block && return block_schur_solver(ctx)
    throw(ArgumentError("linear_backend must be one of $LINEAR_BACKENDS, got :$backend"))
end
