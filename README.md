# Island Memory & The RH Blind Zone

> *Do not hype. Do not lie. Just show.*

A mathematical exploration of the **blind zone of pointwise log-convexity** on the Riemann critical line, alongside a continual-learning neural network architecture (**IslandNet**) constructed directly from the complex-analytic machinery that survives the retraction.

---

## 1. Overview & Lineage

This repository contains two interconnected tracks of research spanning analytic number theory and alternative machine learning primitives:

1.  **The Blind Zone Analysis (`rh_blind_zone.py` / `blind_zone_of_pointwise_convexity.md`):** A rigorous stress test and subsequent retraction of a local, pointwise log-convexity criterion for the Riemann Hypothesis (RH). We show that while the forward direction ($RH \implies \text{convexity}$) holds, the converse fails because local convexity cannot detect off-line zeros outside a narrow detection radius $\delta(T) \approx 2\pi/\log T$. 
2.  **The Island Memory Network (`island_net.py`):** A functional neural architecture built directly on the piece of the blind-zone mathematics that transfers cleanly to engineering: the **influence kernel of a complex pole**. It implements a continual-learning field where memory slots suppress interference entirely through geometry rather than parameter freezing or explicit gating.

---

## 2. Core Mathematical Architecture

### The Blind Zone Check
The pointwise convexity metric evaluated along the critical line ($\sigma = 1/2$) is given by:

$$G(t) = \frac{\partial^2}{\partial \sigma^2} \left| \xi\left(\frac{1}{2} + it\right) \right|^2$$

Through direct high-precision verification against the Hadamard product over the non-trivial zeros $\rho = \beta + i\gamma$, we confirm that local convexity is structurally blind to zeros sitting far from the critical line at high $T$. A global criteria—such as the Lagarias version $\text{Re}(\xi'/\xi(s)) > 0$ for $\text{Re}(s) > 1/2$—is required to catch these anomalous states.

### The Island Memory Field
While local characterization of the zeta function has a blind zone, the underlying complex-variable formulation offers an elegant primitive for handling catastrophic forgetting. In `island_net.py`, we define a continuous field $F(s)$ driven by complex poles (memory slots) $\rho_k$:

$$F(s) = \sum_{k} \frac{C_k}{s - \rho_k}$$

* **Working Memory:** Read out on the boundary line $\text{Re}(s) = 0$. An input vector $x$ maps to positions $t_j(x)$ on the line, evaluating the real and imaginary components of $F(it_j(x))$.
* **Archival (Deep Storage):** Moving a task's memory slots to a deeper coordinate ($\text{Re}(\rho_k) = d_{\text{cold}}$) naturally scales down its boundary contribution by a factor of $\sim 1/d$. 
* **Interference Suppression:** Cross-task interference is mitigated geometrically via the pole's adjoint kernel gradient, preserving older task spaces without hard gating masks.

---

## 3. Repository Structure

├── blind_zone_of_pointwise_convexity.md  # Detailed mathematical analysis and retraction statement├── rh_blind_zone.py                       # High-precision verification engine (mpmath) & counterexample verification├── island_net.py                         # PyTorch/NumPy core implementation of the Island Memory Network├── blind_zone.png                        # Plot displaying detection radii thresholds and the blind zone space└── island_memory.png                     # Diagram verifying task preservation across sequential optimization
---

## 4. Getting Started

### Prerequisites
The math validation requires high-precision float support (`mpmath`), while the baseline modeling requires standard visualization and science stacks.

```bash
pip install numpy scipy torch matplotlib mpmath
```

# Running the Verification Engines

To verify the identity match for $G(t)$, the shrinking detection radius, and visualize the blind-zone gap:

```bash
python rh_blind_zone.py
```


To run the Island Memory network profile and observe geometric task preservation without parameter constraints:

```bash
python island_net.py
```


# 5. Verified Scorecard

| Experiment / Metric | What is Measured / Exhibited | Status |
| :---- | :---- | :---- |
| **Identity Verification** | Direct $\\xi(s)$ derivatives vs. Hadamard zero-sum calculation | **Match verified** at high precision ($10^{-25}$) |
| **Counterexample Scan** | Synthetic zero at $T=10^9, \\beta=0.92$ leaves $G(t) \> 0$ everywhere | **Confirmed Blind** (Pointwise condition fails to witness) |
| **Global Catch** | Lagarias condition $\\text{Re}(\\xi'/\\xi) \> 0$ on the same configuration | **Success** (Correctly flags the off-line zero) |
| **Continual Field Retention** | Sequential task performance decay rate as slots shift to $d\_{\\text{cold}}$ | **Geometric suppression verified** |

# 6. Project Stance

This repository represents an exercise in transparent scientific computing. When the local equivalence criterion was shown numerically to be incomplete, 
the property was downgraded from a characterization to a necessary condition, the analytical limit was computed, and the underlying structural representation 
was translated into an unconstrained machine learning memory primitive.

