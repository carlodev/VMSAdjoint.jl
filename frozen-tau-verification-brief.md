# Handoff brief — quantifying the frozen-τ approximation in the VMS adjoint

**Purpose of this document.** It is a self-contained specification for one task, written to be
handed to a Claude session that has the VMS *adjoint* solver package available. It carries all
the context needed: you should not have to read the paper to act on it.

**Status of what follows.** The derivation, the discrete framing and the implementation plan
were worked out against the paper source and against the *primal* solver
(`SegregatedVMSSolver.jl`, `~/.julia/dev/SegregatedVMSSolver.jl`, branch `limiter`). The
**adjoint package was not available** in the session that produced this brief — it is not in
`~/.julia/dev`, not in `~/Documents`, and no local file mentioned `adjoint`. So every code
snippet below is written in `SegregatedVMSSolver.jl` idioms and must be adapted to the adjoint
package's actual structure. Treat names, spaces and loop placement as suggestions; treat the
mathematics as verified.

---

## 1. The task

Paper: *"Unsteady adjoint variational multiscale method for airfoil applications"*, source in
`c:\Users\mail\Documents\Articles\VMS-Adjoint\v1\`. The relevant section is
`Adjoint_Formulation.tex`, `\subsection{The frozen stabilization-parameter approximation}`
(`\label{sec:adjoint-frozen-tau}`). It currently ends on this TODO:

> report, for at least one design variable and one mesh, the ratio ‖δB^τ‖/‖B^†‖ and show that
> it is small compared with the gradient error measured in Section~\ref{sec:results-gradient}.
> Without this the Taylor-remainder test alone can be read as evidence that the frozen-τ adjoint
> is simply not the derivative of the discrete objective.

Deliver: the numbers, and the paragraph that replaces the TODO.

A closely related second TODO, currently commented out in `Sensitivity.tex` (around line 58,
`sec:sensitivity-mesh`), asks the same question about the *mesh-metric* dependence of τ. The
same implementation answers both — see §8.

---

## 2. Notation and the primal formulation

Paper macros: `\dg` = †, `\ua` = φ_u (adjoint velocity), `\pa` = φ_p (adjoint pressure),
`\va` = v† (adjoint weighting function), `\Rm` = R_m, `\Rc` = R_c, `\bbeta` = β (design vars).

Primal weak form (eq. `primal-weakform`), incompressible NS, equal-order T_k/T_k triangles:

```
B(v_h, q_h; u_h, p_h) = B_G + B_STAB + B_VMS1 + B_VMS2 = 0
```

with resolved-scale residuals (eq. `resolved-residuals`)

```
R_m = ∂u_h/∂t + u_h·∇u_h + ∇p_h − ν Δu_h
R_c = ∇·u_h
```

and the τ-weighted groups

```
B_STAB = ∫ τ_M (u_h·∇v_h)·R_m      (SUPG)
       + ∫ τ_M (∇q_h)·R_m          (PSPG)
       + ∫ τ_C (∇·v_h) R_c         (LSIC)

B_VMS1 = ∫ (u_h·∇v_hᵀ)·(τ_M R_m)                  (cross stress)
B_VMS2 = −∫ ∇v_h : (τ_M R_m ⊗ τ_M R_m)            (subgrid Reynolds stress)
```

Stabilization parameters (eq. `tau`):

```
τ_M = (4/Δt² + u_h·G u_h + C_I ν² G:G)^(−1/2)
τ_C = (τ_M g·g)^(−1)
```

G is the element metric tensor, g the associated vector, both from the parent-element mapping.
**τ_M and τ_C depend on the state through u_h** — that is the entire subject of this brief.

Adjoint weak form (eq. `adjoint-weakform`): find (φ_u, φ_p) such that

```
B†(t; v†, q†; φ_u, φ_p) = B_G† + B_STAB† + B_VMS1† + B_VMS2† = −∂I_C/∂U [v†, q†]   ∀ (v†, q†)
```

The transposition treats τ_M, τ_C as constants. That is the approximation to be quantified.

Objectives: time-averaged C_L, C_D, C_D/C_L, and an inverse-design ½(C̄_L − C̄_L*)².
Adjoint is linear, integrated **backward** in time, reading the stored primal trajectory.

---

## 3. What the ratio actually is (the discrete framing)

This is the framing that makes the result publishable — use it in the write-up.

Let `A = ∂R/∂U` be the exact tangent of the discrete primal residual and `Ã = A − E` the
frozen-τ tangent that is actually transposed, `E` being the block coming from differentiating
τ_M and τ_C. The solver produces `Ψ̃` from

```
Ãᵀ Ψ̃ = −∂I_C/∂U
```

Substituting into the *exact* adjoint equation:

```
Aᵀ Ψ̃ + ∂I_C/∂U = (Ã + E)ᵀ Ψ̃ + ∂I_C/∂U = Eᵀ Ψ̃
```

So **δB^τ is exactly the residual of the true adjoint equation evaluated at the computed adjoint
field**, and ‖δB^τ‖/‖B†‖ is a relative residual — a standard, well-defined object, not a loose
term-by-term size comparison. Everything entering it is already in memory during the backward
sweep, so the ratio costs no extra solve.

### Three conceptual points to preserve in the paper text

1. **It is not the size of ∂τ/∂u that is being measured.** `∂τ_M/∂u_h = −τ_M³ G u_h` is not
   small on its own — on stretched near-wall elements it is O(1) in the relevant units. What
   makes the dropped contribution small is that it always appears multiplied by a factor
   containing `R_m` (or `R_c`). The smallness comes from *consistency of the discretization*,
   not from τ being insensitive to the state. This is why the approximation is qualitatively
   weaker than frozen turbulence in a RANS adjoint, where the dropped term carries ν_t and does
   not vanish under refinement.

2. **The denominator is the whole adjoint operator applied to Ψ**, which equals the adjoint
   right-hand side — not "the other terms" in a loose sense.

3. **A small operator residual does not by itself imply a small gradient error.** The map from
   adjoint residual to gradient error is `A^{−T}` followed by the sensitivity kernel; an
   ill-conditioned adjoint could amplify a 1e−4 residual into a 1e−2 gradient error. §7 measures
   that amplification instead of assuming it away. Without §7 you have shown the term is small
   *in the equation*; with §7 you have shown it is small *in what the optimizer consumes*, which
   is what the TODO demands.

Also worth stating explicitly in the paper: the ratio itself does **not** depend on the design
variable — it is a property of the adjoint field, one number per objective per mesh per time
level. The design-variable dependence the TODO asks for appears only once the defect is
propagated to the gradient (§7).

---

## 4. Where the approximation lives in the code

In the primal package, `src/Commons/Equations/StabilizationOperations.jl`, lines 1–2:

```julia
val_u(x) = x
val_u(x::Gridap.Fields.ForwardDiff.Dual) = x.value
```

used inside `momentum_stabilization` (~line 81) as

```julia
uu_new = VectorValue(val_u.(uun)...)
```

That line **is** the frozen-τ approximation: it strips the dual part of `u` before τ_M is
formed, so ForwardDiff never differentiates τ. `compute_dt_eff` calls `val_u` too, for the
limiters. Expect the same construct in the adjoint package; find it first, because §9 exploits
it for a near-free correctness proof.

---

## 5. Closed form of δB^τ

The paper writes only the SUPG/PSPG pair (eq. `tau-variation-terms`) plus "analogous τ_C and VMS
terms". The complete expression follows.

Differentiating eq. `tau` at fixed mesh, with G symmetric:

```
δτ_M = −τ_M³ (G u_h · δu)
δτ_C = +(τ_M / (g·g)) (G u_h · δu)          [ = τ_C² (g·g) τ_M³ (G u_h · δu), as in the paper ]
```

**Mechanical rule.** Take each primal τ-weighted integrand, put `(v_h, q_h) → (φ_u, φ_p)`, and
substitute, with the scalar `a := (G u_h)·v†`:

```
τ_M  →  −τ_M³ · a          τ_C  →  +(τ_M/(g·g)) · a          τ_M²  →  −2 τ_M⁴ · a
```

Term by term:

| primal group | τ-variation contribution |
|---|---|
| SUPG  | `−τ_M³ a [(u_h·∇φ_u)·R_m]` |
| PSPG  | `−τ_M³ a [(∇φ_p)·R_m]` |
| LSIC  | `+(τ_M/(g·g)) a (∇·φ_u) R_c` |
| VMS1  | `−τ_M³ a [(u_h·∇φ_uᵀ)·R_m]` |
| VMS2  | `+2 τ_M⁴ a [∇φ_u : (R_m ⊗ R_m)]` |

Collecting:

```
δB^τ(v†) = ∫_Ω v† · s^τ dΩ ,          s^τ := Θ · (G u_h)

Θ = −τ_M³ [ (u_h·∇φ_u + u_h·∇φ_uᵀ + ∇φ_p) · R_m ]     ← SUPG + VMS1 + PSPG
    + 2 τ_M⁴ [ ∇φ_u : (R_m ⊗ R_m) ]                    ← VMS2
    + (τ_M/(g·g)) (∇·φ_u) R_c                          ← LSIC
```

δB^τ has **no q† component** — δτ depends on δu only — so it lives entirely in the velocity
block. This confirms the paper's statement that the neglected terms "act as zeroth-order
reaction sources in the adjoint momentum equation": `s^τ` is literally that source.

### Implementation traps

1. **R_m must include ∂u_h/∂t**, evaluated with the same θ-scheme difference the solver uses
   (θ = 0.55 in this work). Dropping it understates the defect by the dominant term in an
   unsteady flow.
2. **dt_limiter.** With `AdvectiveLimiter` / `CombinedLimiter`,
   `dt_eff = C_adv/√(u·G·u)` is itself state-dependent, so `∂τ_M/∂u` picks up a second
   contribution and the closed form above is incomplete. The default in `TensorFormulation` is
   `NoLimiter()`; confirm which the paper runs used and state it.
3. **VMS2 may not exist in the primal.** `git grep` for `⊗` / `VMS2` across *every local branch*
   of `SegregatedVMSSolver.jl` returns nothing: `segregated_equations` assembles `Auu_vms1`,
   `Tuu_vms1`, `Aup_vms1` but no subgrid Reynolds-stress term. Eqs. `BVMS2` / `adj-BVMS2` may be
   formulation-only. **Verify this against the adjoint package.** If VMS2 is not assembled, drop
   its line from Θ — and the paper needs a sentence saying the term is formulated but not active
   in the reported runs.
4. **Skew terms.** `segregated_equations` carries `skewcoeff`-weighted τ_M terms. `skew` is
   `false` for `VMS()` by default but `true` for `SUPG()`. If a run uses skew, those integrands
   contribute to δB^τ as well, by the same substitution rule.

---

## 6. Steps 1–3: assemble the ratio

### Defect source

```julia
"""
    tau_defect_source(uh, ∂tuh, ph, φu, φp, stab_coeff, simcase; vms2=false)

Vector field s^τ with δB^τ(v†) = ∫ v†⋅s^τ dΩ: the source that the exact linearization of
τ_M, τ_C adds to the adjoint momentum equation and that the frozen-τ formulation drops.
"""
function tau_defect_source(uh, ∂tuh, ph, φu, φp, stab_coeff, simcase; vms2=false)
    @unpack G, gg = stab_coeff
    @sunpack ν, D  = simcase

    τm = momentum_stabilization(uh, stab_coeff, simcase)

    Rm = ∂tuh + (cconv ∘ (uh, ∇(uh))) + ∇(ph)      # + (-ν*Δ(uh)) for order ≥ 2
    Rc = ∇ ⋅ uh

    # mirror EXACTLY the contractions used in the primal SUPG / VMS1 / PSPG integrands
    wM = (cconv ∘ (uh, ∇(φu))) + (∇(φu) ⋅ uh) + ∇(φp)

    Θ = -(τm^3) * (wM ⋅ Rm) + (τm / gg) * ((∇ ⋅ φu) * Rc)
    vms2 && (Θ = Θ + 2*(τm^4) * (∇(φu) ⊙ (Rm ⊗ Rm)))

    return Θ * (G ⋅ uh)
end
```

`cconv(u, ∇v) = (∇v)'⋅u` is the primal package's convective wrapper. The two convective
contractions in `wM` correspond to SUPG (`u_h·∇φ_u`) and VMS1 (`u_h·∇φ_uᵀ`) respectively —
match them to whatever index convention the adjoint package's own VMS1 term uses, do not trust
the transcription blindly.

### Inside the backward loop

```julia
sτ  = tau_defect_source(uh_n, ∂tuh_n, ph_n, φu_n, φp_n, stab_coeff, simcase)
δbτ = assemble_vector(v -> ∫(v ⋅ sτ)dΩ, V_adj)   # velocity test space, adjoint BCs
```

### Denominator

Take `b†` as the **velocity block** of `Aᵀ Ψ̃`, i.e. the assembled adjoint right-hand side
`−∂I_C/∂U`. Restricting to the velocity block is apples-to-apples (δb^τ touches nothing else)
and is the conservative choice — using the full vector only shrinks the ratio.

Assemble `b†` two ways once, as a free implementation check: (a) directly as the adjoint RHS,
(b) term-by-term from eq. `adjoint-weakform` with the computed Ψ. They must agree to solver
tolerance. If they do not, the adjoint assembly has a bug and the ratio is meaningless.

### Norm

Plain ℓ₂ of coefficient vectors is not a norm on the dual space and a referee can say so. Use
the mass-matrix dual norm — one CG solve per level, negligible cost:

```julia
M = assemble_matrix((u,v) -> ∫(u ⋅ v)dΩ, U_adj, V_adj)   # once, outside the loop
dualnorm(b) = sqrt(dot(b, M \ b))
ετ_n = dualnorm(δbτ) / dualnorm(b†)
```

### What to log

CSV per time level: `t, ‖δbτ‖, ‖b†‖, ετ`. Report `max_n ετ`, the window RMS, and plot ετ
against `t/ftt`.

**Run it on both meshes** — the baseline and the coarse-leading-edge mesh already used in
`sec:results-mesh`. A decreasing ετ under refinement is the observable consequence of the
paper's "vanishes at the rate of the local truncation error" argument, and turns an assertion
into a measurement. A third, systematically refined mesh would make it a convergence statement.

---

## 7. Step 4: propagate to the gradient (the number the TODO really wants)

The ratio bounds the operator defect. The referee cares about the gradient. One extra **linear**
backward sweep gives it, reusing the stored primal trajectory:

1. Sweep the same adjoint operator backward with terminal condition `δφ_u(t_F) = 0`, homogeneous
   datum on Γ_a, and volume source `−δb^τ(tⁿ)` injected at each level in place of `−∂I_C/∂U`.
   Result: `δΨ`.
2. Feed `δΨ` through the **existing** sensitivity routine (eq. `pIpb` / `gradient-boundary`),
   same primal, zeroing the explicit `∂I_C/∂β` part (it does not involve Ψ). Result: `δg`.
3. Report `‖δg‖₂/‖g‖₂`, `max_k |δg_k|/|g_k|` and the control point where it occurs, and
   `1 − cos∠(g, g+δg)`.

That answers "for at least one design variable" for all forty at once, and `‖δg‖/‖g‖` is
directly comparable with the finite-difference scatter in `sec:results-gradient`.

**Caveat to state in one sentence in the paper:** the correction sweep uses `Ã`, not `A`, so
`δg` is the leading term of a Neumann series and the residual error is `O(ε_τ²)` — which is
precisely why the estimate is sharp once ε_τ ≪ 1.

**Time-level bookkeeping.** In a segregated solver τ is built on `u_adv`, an extrapolation of
previous time steps, not on the current unknown. Within a step τ is genuinely constant; across
steps it depends on the state history. So the defect computed at step *n* must be injected into
the adjoint equation at the step(s) `u_adv^n` was extrapolated from, split by the extrapolation
weights. Getting this wrong shifts δg by a step. §9 catches it.

---

## 8. Bonus: the same code discharges the mesh-metric TODO

The commented-out eq. `gradient-metric` in `Sensitivity.tex` (~line 49–58) has identical
structure; only the chain rule changes:

```
∂τ_M/∂G : δG = −½ τ_M³ [ (u_h ⊗ u_h) : δG + 2 C_I ν² G : δG ]
```

and `∂G/∂β_k` needs no flow solve: morph the mesh by `ε e_k` with the same RBF operator,
recompute G via `compute_G` (`StabilizationParameters.jl`), and difference. Reuse Θ verbatim
with `(G u_h)·v†` replaced by that contraction. One implementation, both neglected
contributions.

---

## 9. Step 5: prove the δB^τ implementation is right (nearly free)

Make the freeze switchable:

```julia
const FREEZE_TAU = Ref(true)
val_u(x) = x
val_u(x::Gridap.Fields.ForwardDiff.Dual) = FREEZE_TAU[] ? x.value : x
```

Then, at one time level with a random state perturbation `δU`, form the AD directional
derivative of the primal weak form contracted with Ψ, once with `FREEZE_TAU[] = true` and once
`false`. The difference **must equal** the closed-form `δB^τ(δu)` to machine precision.

This validates the hand derivation, the contraction conventions, and — if performed on the fully
discrete operator rather than one spatial form — the time-level bookkeeping of §7. It belongs in
`test/CommonsTests/EquationsTests/` as a permanent regression test.

---

## 10. What to report, and suggested replacement text

Suggested table (one row per mesh, one block per objective):

| mesh | max ε_τ | rms ε_τ | ‖δg‖/‖g‖ | max_k \|δg_k\|/\|g_k\| | 1 − cos∠ |
|---|---|---|---|---|---|

Replacement for the TODO paragraph in `sec:adjoint-frozen-tau`:

```latex
Second, the magnitude of the neglected contribution can be bounded a~posteriori,
because all of its factors are available once the primal and adjoint solutions
have been computed. Assembling~\eqref{eq:tau-variation-terms} against the
computed adjoint field at every time level of the backward sweep gives the
residual of the exact adjoint equation, whose size relative to the adjoint
right-hand side is $\|\delta B^{\tau}\|/\|B^{\dg}\| = \num{}$ at most and
$\num{}$ in root-mean-square over the averaging window on the baseline mesh,
falling to $\num{}$ under refinement, as the residual argument above predicts.
Propagating this defect through one additional linear adjoint sweep and through
the sensitivity kernel~\eqref{eq:gradient-boundary} bounds its effect on the
gradient itself at $\|\delta\bm{g}\|/\|\bm{g}\| = \num{}$, with a largest
component-wise contribution of \num{}\% at control point \num{}; both are well
below the finite-difference scatter reported in
Section~\ref{sec:results-gradient}.
```

The full δB^τ expression of §5 should go into `Appendix_Adjoint_Operators.tex`, which already
collects the term-by-term transposition and explicitly defers this correction to
`sec:adjoint-frozen-tau`.

---

## 11. Unrelated discrepancy noticed in passing

`StabilizationOperations.jl` (~line 84) computes `τ₁ = Ci[1]*(2/dt_eff)^2` with default
`Ci = [4, 36]`, i.e. `16/Δt²`, while eq. `tau` in the paper states `4/Δt²`. It does not affect
this task — τ₁ is state-independent, so `∂τ_M/∂u` is unchanged — but the paper and the code
disagree on the formula as written. Worth reconciling before submission.
