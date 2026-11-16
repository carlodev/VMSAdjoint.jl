#############################################################################
# Frozen-operator θ-method driver (assembly + factorization reuse)
#############################################################################
#
# Both the primal and the adjoint VMS systems are *linear* once the
# linearization field `params[:uh]` is fixed, so the θ-method stage matrix
#
#     A = M + θ·dt·K(uh)
#
# only changes when `uh` changes. Re-assembling and re-factorizing `A` at every
# time step — what a plain `TransientLinearFEOperator` does — is the dominant
# cost of an unsteady run and is almost entirely avoidable.
#
# Gridap already knows how to reuse a factorization: `ThetaMethod` on a
# `TransientLinearFEOperator` declared with `constant_forms = (true, true)`
# sets `reuse = true` on the stage operator, which skips both the Jacobian
# assembly and `numerical_setup!` on every step after the first (the residual
# is still re-assembled every step, so time-dependent right-hand sides and
# Dirichlet data remain correct).
#
# Declaring the forms constant for the *whole* run would freeze `uh` at `t0`,
# which is wrong. So instead we march in chunks of `matrix_freq_update` steps:
# the operator is rebuilt (and the factorization recomputed) once per chunk
# around the latest state, and reused within it. This is the same
# `matrix_freq_update` trade-off `SegregatedVMSSolver` makes, expressed
# entirely through the public Gridap ODE API — no reimplementation of the
# θ-method, of the Dirichlet lifting, or of the transient FE spaces.
#
# Cost per chunk of N steps: 1 assembly + 1 factorization + N residual
# assemblies + N back-substitutions, instead of N of each.
#############################################################################

"""
    march_frozen_theta(make_operator, lsolver, X, xh0, t0, tF, dt, θ, freq;
                       on_step!, refresh!, label)

March a linear transient problem from `t0` to `tF` with the θ-method, rebuilding
the operator every `freq` steps and reusing its factorization in between.

- `make_operator()` -> `TransientLinearFEOperator`, built with
  `constant_forms = (true, true)` around the *current* linearization field.
  Called once per chunk.
- `lsolver`         -> the Gridap `LinearSolver` (see [`linear_solver`](@ref)).
- `X`               -> transient trial space; `X(t)` must be valid at every
  chunk boundary.
- `xh0`             -> initial `FEFunction` on `X(t0)`.
- `on_step!(step, t, xhtn)` -> called after every step, with the 1-based global
  step index, the physical time and the solution `FEFunction`. Note `xhtn` is
  reused across steps: copy anything you keep.
- `refresh!(step, state)` -> called at every chunk boundary with the free-value
  vector of the latest state. This is where the caller updates the linearization
  field that `make_operator` closes over.
- `label`           -> string used in the per-chunk timing log.

Returns the final state as an `FEFunction` on `X(tF)`.

`freq = 1` reproduces the fully-updated (previous) behaviour exactly, at the
previous cost; it is the reference to validate any larger value against.
"""
function march_frozen_theta(make_operator::Function, lsolver, X, xh0,
    t0::Float64, tF::Float64, dt::Float64, θ::Float64, freq::Int;
    on_step!::Function, refresh!::Function, label::String="transient")

    nsteps = max(0, round(Int, (tF - t0) / dt))
    freq = max(1, freq)
    ode_solver = ThetaMethod(lsolver, dt, θ)

    state = copy(get_free_dof_values(xh0))
    step = 0
    nchunks = cld(nsteps, freq)
    chunk = 0

    while step < nsteps
        chunk += 1
        nchunk = min(freq, nsteps - step)
        # accumulate times from the step index, never by adding dt repeatedly,
        # so chunk boundaries stay exactly on the global time grid
        t_start = t0 + step * dt
        t_end = t0 + (step + nchunk) * dt
        first_step = step + 1

        chunk_t0 = time()
        op = make_operator()
        xh_start = FEFunction(X(t_start), state)

        # Gridap assembles and factorizes lazily, on the first step of the
        # `Gridap.solve` iteration, so the cost split that matters is
        # first-step (assembly + factorization + solve) versus the reused
        # steps that follow. Timing `make_operator` alone would report ~0.
        first_step_time = 0.0
        chunk_time = @elapsed begin
            for (t, xhtn) in Gridap.solve(ode_solver, op, t_start, t_end, xh_start)
                step += 1
                on_step!(step, t, xhtn)
                copyto!(state, get_free_dof_values(xhtn))
                step == first_step && (first_step_time = time() - chunk_t0)
            end
        end

        reused = nchunk - 1
        reused_msg = reused > 0 ?
                     ", $(reused) reused steps $(round(chunk_time - first_step_time, digits=2))s " *
                     "($(round((chunk_time - first_step_time) / reused, digits=3))s/step)" : ""
        @info "$label chunk $chunk/$nchunks: steps $first_step..$step " *
              "(t $(round(t_start, digits=6)) -> $(round(t_end, digits=6))), " *
              "first step (assembly+factorization+solve) $(round(first_step_time, digits=2))s" *
              reused_msg

        refresh!(step, state)
    end

    return FEFunction(X(t0 + nsteps * dt), state)
end

"""
    solver_context(simcase, params, P, Q; transient)

Build the [`SolverContext`](@ref) the `:block` backend needs. `τ = θ·dt` for a
transient stage operator, `Inf` for a steady one.
"""
function solver_context(simcase::Airfoil, params::Dict{Symbol,Any}, P, Q; transient::Bool)
    @sunpack ν, dt, θ = simcase
    @unpack dΩ = params
    return SolverContext(P=P, Q=Q, dΩ=dΩ, ν=ν, τ=transient ? θ * dt : Inf)
end
