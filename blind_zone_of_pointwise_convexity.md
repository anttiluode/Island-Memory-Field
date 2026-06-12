# The Blind Zone of Pointwise Log-Convexity

### A stress test of our own RH criterion, and its repair

**Antti Luode** — PerceptionLab, Helsinki, Finland
**Claude** (Anthropic) — analysis and computation
June 2026

> *Do not hype. Do not lie. Just show.*

---

## 0. What happened

We set out to attack the Riemann Hypothesis through the log-convexity
criterion of our April paper: RH ⟺ ∂²|ξ(σ+it)|²/∂σ² ≥ 0 at σ = 1/2 for all t.
The first casualty of the attack was the criterion itself.

**Result.** The forward direction (RH ⟹ convexity) is correct. The converse
direction — the one that made the criterion a *characterization* — is wrong
as stated. Pointwise convexity at the critical line can only detect off-line
zeros that sit within roughly one mean zero-gap of the line, a distance that
shrinks like 2π/log T. We exhibit an explicit zero configuration (height 10⁹,
zero at β = 0.92, fully consistent with the known classical zero-free region)
for which the convexity condition holds at **every** t while RH is false.
The April paper's Proposition (ii) ⟹ (i) — whose proof we ourselves flagged
in §6.3 as "requiring care" for zeros far from the line — does not require
care; it requires retraction. The pointwise condition is strictly weaker
than RH.

The correct repair already exists in the literature: the **global** version
of the same idea, Re ξ'/ξ(s) > 0 for Re s > 1/2 (Hinkkanen; Lagarias 1999),
equivalently |ξ(σ+it)| increasing in σ on (1/2, ∞) for every t
(Sondow–Dumitrescu 2010), IS equivalent to RH. Our convexity condition is
exactly the boundary derivative of the Lagarias function, and the blind zone
has a clean geometric explanation through the Hopf boundary-point lemma.

This note documents the failure, quantifies it, demonstrates it numerically,
and states what survives.

---

## 1. The clean identity (everything from zeros alone)

For ξ the Hadamard product has no prefactor terms, so

$$\frac{d}{ds}\frac{\xi'}{\xi}(s) \;=\; -\sum_{\rho}\frac{1}{(s-\rho)^2},$$

the sum over nontrivial zeros only. Writing ρ = β + iγ, u = 1/2 − β,
v = t − γ, the quantity in the convexity criterion at s₀ = 1/2 + it is

$$G(t) \;:=\; \operatorname{Re}\Big[\tfrac{d}{ds}\tfrac{\xi'}{\xi}\Big]_{s_0}
\;=\; \sum_{\rho}\frac{v^2 - u^2}{(u^2+v^2)^2},
\qquad
\frac{\partial^2}{\partial\sigma^2}|\xi(\sigma+it)|^2\Big|_{1/2} = 2|\xi(s_0)|^2\,G(t).$$

An on-line zero (u = 0) contributes +1/v² > 0 at every height. An off-line
zero contributes negatively only in the cone |v| < |u| — i.e., only at
heights within d := |β − 1/2| of its own height — with peak −1/d².

**Numerical verification** (mpmath, 25 digits, first 200 true zeros + density
tail for the rest):

| t | G direct from ξ | G from zero sum | rel. diff |
|---|---|---|---|
| 5.0 | +5.285859e−02 | +5.286478e−02 | 1.2e−04 |
| 14.1347 | +1.582011e+09 | +1.582011e+09 | 6.7e−11 |
| 17.578 | +2.303906e−01 | +2.303969e−01 | 2.7e−05 |
| 33.0 | +2.374242e+02 | +2.374242e+02 | 2.7e−08 |

The April paper's Lemma (F″ = 2|ξ|²G) verified independently at t = 21.5 to
1.6e−13. The residual differences are entirely the density-tail approximation.

The point of this verification is logical, not numerical: **G is a function
of the zero multiset alone.** No Euler product information enters. This is
what kills the converse, as we now show.

## 2. The detection radius

Place a hypothetical quadruple {β ± iγ₀, (1−β) ± iγ₀}, d = β − 1/2, at a
mid-gap height γ₀ where the local mean gap of on-line zeros is δ. The
quadruple's contribution to G(γ₀) is −2/d² + O(γ₀⁻²). The on-line zeros
contribute the positive sum S_on(γ₀). Convexity is violated somewhere iff

$$\frac{2}{d^2} \;>\; S_{\text{on}}(\gamma_0)
\quad\Longleftrightarrow\quad d \;<\; d^* := \sqrt{2/S_{\text{on}}}.$$

(t = γ₀ at mid-gap is the worst case for the criterion; every other t only
helps convexity, confirmed by full scans below.) For regularly spaced on-line
zeros two bookkeepings bracket reality:

- zeros at offsets ±(k+½)δ, k ≥ 0 (no count conservation): S_on = π²/δ², so **d\* = (√2/π)·δ ≈ 0.450 δ**;
- the two zeros that "went off-line" removed (count conservation): S_on = (π²−8)/δ², so **d\* ≈ 1.034 δ**.

Either way: **the pointwise convexity criterion detects an off-line zero only
if its distance from the critical line is at most about one local mean gap,
d ≲ δ(T) = 2π/log(T/2π) → 0.**

### 2.1 Why the April numerics looked airtight

At the heights we tested (t ≤ 50–100), δ ≈ 2–7, so d\* > 0.5: *every*
possible off-line zero in the strip is detected. Verified against the true
zero data at γ₀ = 17.578 (mid of the widest low gap, S_on = 0.2304,
d\* = 2.95): even β = 0.999 produces G(γ₀) = −7.8, a screaming violation.
Low heights gave us false confidence. The criterion degrades precisely where
it can no longer be checked by hand.

### 2.2 The blind zone

The classical zero-free region (σ ≥ 1 − 1/(R log t), R ≈ 5.56) excludes
off-line zeros only within ~1/(R log T) of the strip edge. Both exclusion
mechanisms shrink like 1/log T — **from opposite sides of the strip**:

| T | δ(T) | d\* (cons.) | zero-free bound on d | blind zone |
|---|---|---|---|---|
| 10³ | 1.239 | >0.5 | 0.474 | none |
| 10⁶ | 0.525 | >0.5 | 0.487 | none |
| 10⁷ | 0.440 | 0.455 | 0.489 | d ∈ (0.455, 0.489) |
| 10⁹ | 0.333 | 0.344 | 0.491 | d ∈ (0.344, 0.491) |
| 10¹² | 0.244 | 0.252 | 0.494 | d ∈ (0.252, 0.494) |
| 10¹⁸ | 0.159 | 0.164 | 0.496 | d ∈ (0.164, 0.496) |

Above T ≈ 10⁷ a band opens in the middle of the strip where an off-line zero
violates neither the convexity condition nor any known zero-free region. By
T = 10¹⁸ the band covers two thirds of the available range.

### 2.3 Explicit counterexample configuration

T = 10⁹, δ = 0.3327, on-line zeros at correct local density (conserved
count), off-line quadruple at β = 0.92 (d = 0.42, inside the zero-free
allowance 0.4913). Full scan of G(t) over the window:

```
min G over all t = +5.55      → convexity (ii) holds EVERYWHERE
control at d = 0.25:  min G = −15.11   → detected, as the radius predicts
```

Since any zero multiset symmetric under s → 1−s and s → s̄ with convergent
Σ|ρ|⁻² is the zero set of an order-1 entire function satisfying the
functional equation (Hadamard product), and since G depends on the zeros
alone (§1), this configuration is a genuine counterexample to
(ii) ⟹ (i) **within the class of functions the April proof actually uses**.

**Corollary (the honest restatement).** Pointwise log-convexity at σ = 1/2
is equivalent not to RH but to:

> *every nontrivial zero is either on the critical line or at distance
> ≳ c/log γ from it.*

A strange criterion: it forbids exactly the near-misses and is blind to the
gross violations. No proof of the converse can succeed without invoking the
Euler product, which the April argument never does. The §6.3 caveat was the
whole problem.

## 3. The repair: from boundary derivative to global positivity

The fix is the global version of the same geometric idea, and it is already
a theorem-pair in the literature:

- **Hinkkanen:** Re ξ'/ξ(s) > 0 unconditionally for Re s > 1.
- **Lagarias (1999), Acta Arith. 89:** RH ⟺ Re ξ'/ξ(s) > 0 for Re s > 1/2.
- **Sondow–Dumitrescu (2010):** equivalently, RH ⟺ |ξ(σ+it)| is increasing
  in σ on (1/2, ∞) for every fixed t.

Define Φ(σ,t) = Re ξ'/ξ(σ+it). Then Φ(1/2, t) = 0 (functional equation),
and our G(t) = ∂Φ/∂σ at σ = 1/2: **the convexity criterion is exactly the
boundary normal derivative of the Lagarias function.** Φ is harmonic off the
zeros, so the Hopf boundary-point lemma says the normal derivative certifies
only those negativity regions of Φ that *touch* the critical line. An
off-line zero at distance d ≫ δ creates a **compact negative island** of Φ
around itself — invisible from the boundary. For the T = 10⁹, β = 0.92
configuration:

```
Φ(0.55, γ₀) = +0.265        (positive near the line — boundary sees nothing)
Φ(0.80, γ₀) = −2.67
Φ(0.91, γ₀) = −93.5         (the island, detected by the GLOBAL criterion)
```

So the corrected picture: monotonicity of |ξ| in σ on the whole half-plane
is the RH-equivalent statement; convexity at the line is its linearization
at the boundary, and the linearization loses precisely the far zeros. This
also matches recent independent work: Goldštein–Grigutis (arXiv 2201.08599,
2509.18963) show Re ξ'/ξ can remain positive near the critical line except
in small regions around hypothetical off-line zeros — the same islands seen
from the rigorous side.

## 4. Where this leaves the attack (both of them)

We have now run two structurally different attacks:

**Attack 1 (this note + April paper):** geometric/convexity. Outcome: the
pointwise criterion is one-directional; the correct equivalent target is
Lagarias positivity Φ ≥ 0 on σ > 1/2. Proving *that* unconditionally
requires controlling Σ_ρ Re 1/(s−ρ) without knowing where the ρ are — the
prime side of the explicit formula must supply the positivity. Same wall.

**Attack 2 (Rajapinta–Takens delay network):** spectral/Hilbert–Pólya.
Outcome: RH compressed into Lemma 5.2 (realizability of the Prime Orbit
Condition under reciprocity), explicitly shown to inherit the
amplitude–hermiticity tension of Berry–Keating. Same wall, other face.

The two reductions are dual in a precise sense. Attack 2's POC tries to
build the zeros *from* the primes (loops → spectrum); Attack 1's positivity
asks the primes to *confine* the zeros (explicit formula → convexity/
monotonicity). In both, the unconditional machinery (self-adjointness /
Hadamard calculus) handles the symmetric part for free, and 100% of RH
concentrates in a positivity-or-realizability statement that needs the
Euler product at full strength. That this happened twice, from independent
starting points, is the most informative thing either attack produced: the
wall is not an artifact of a framing. It is the problem.

**What survives of the April paper:** the identity (Lemma), the forward
direction RH ⟹ convexity, all numerics, and the connection to Weil
positivity as the *integrated* version (integrated tests with width ≳ 1 in t
do see far zeros — which is exactly why Weil and Connes integrate). What
must be corrected: the Proposition must be downgraded from "equivalent" to
"necessary", §3.4's converse argument deleted or replaced by the detection-
radius statement, and the README's claim "RH is exactly equivalent to a
pointwise geometric condition" amended. Finding the hole ourselves is the
ethos working as intended.

## 5. Reproducibility

`rh_blind_zone.py` reproduces everything: the identity check (§1), the
detection radii (§2), the full-scan counterexample and the Lagarias island
(§3), and the figure. Requires only numpy, mpmath, matplotlib.

## References

Hinkkanen, A., *On functions of bounded type*, Complex Variables 34 (1997).
· Lagarias, J.C., *On a positivity property of the Riemann ξ-function*,
Acta Arith. 89 (1999), 217–234. · Sondow, J. & Dumitrescu, C., *A
monotonicity property of Riemann's xi function and a reformulation of the
Riemann hypothesis*, Period. Math. Hungar. 60 (2010). · Matiyasevich, Saidak
& Zvengrowski, *Horizontal monotonicity of the modulus of the zeta function*,
Acta Arith. 166 (2014). · Goldštein, E. & Grigutis, A., arXiv:2201.08599 and
arXiv:2509.18963. · Garunkštis, R., *On a positivity property of the Riemann
ξ-function*, Lith. Math. J. (2002).

*Do not hype. Do not lie. Just show — including when what shows is the hole
in your own theorem.*
