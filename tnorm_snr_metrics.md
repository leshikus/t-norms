# T-Norm Quality Metrics: Signal-to-Noise Perspective

Based on: Ryzhov A.P., Fedotov F.A. — *Fuzzy Risk Aggregation on Trees Using T-norms*,
Lomonosov Readings 2026, Moscow State University, April 1, 2026.

---

## Setup

A directed tree with leaves $x_1,\ldots,x_n \in [0,1]$ (higher = better).
All internal nodes use the same strict Archimedean t-norm:
$$T(x,y) = g^{-1}(g(x)+g(y))$$

The root value collapses to:
$$R(x_1,\ldots,x_n) = g^{-1}\!\left(\sum_{i=1}^n g(x_i)\right)$$

**Profile:** $h(x) := -g'(x) > 0$, so $g(x) = \int_x^1 h(u)\,du$.

**Heterogeneous configuration:** $(t, a, \ldots, a)$ with $t < a$,
where $a$ is the background level and $t$ is a single incident leaf.

---

## SNR Metrics

### Incident sensitivity — "signal"
$$I_n(t,a) = \frac{\partial R_n}{\partial t} = \frac{h(t)}{h(R_n(t,a))}$$

### Background sensitivity — "noise"
$$B_n(t,a) = \frac{\partial R_n}{\partial a} = (n-1)\frac{h(a)}{h(R_n(t,a))}$$

### Selectivity — **SNR analogue**
$$\Gamma_n(t,a) = \frac{I_n(t,a)}{B_n(t,a)} = \frac{h(t)}{(n-1)\,h(a)}$$

> Measures local contrast between incident and background.
> Structurally analogous to **d′** in signal detection theory.

### Detection — finite root drop under incident
$$D_n(t;a) = R_n(a) - R_n(t,a)$$

---

## Optimization Problem

Minimize worst-case noise sensitivity subject to signal constraints:

$$\min_h \sup_{(t,a)\in\Omega^+} B_n(t,a)$$

subject to:
- $D_n(t;a) \geq D^*$ — minimum detection (signal floor)
- $\Gamma_n(t,a) \geq \Gamma^*$ — minimum selectivity (SNR floor)
- $J_\text{reg}(h) \leq M$ — profile regularity
- $h$ non-increasing — global robustness (Lipschitz $L=1$ in $\ell^1$)

Working domain: $\Omega^+ = \{(t,a): a\in[a_-,a_+],\, t\in[t_{\min},\eta a]\}$.

---

## Which T-Norm Maximizes $\Gamma_n$?

### Among all t-norms: $T_{\min}(x,y) = \min(x,y)$

**Proof.** For $(t,a,\ldots,a)$ with $t < a$: $R_n = t$, so $\partial R/\partial t = 1$,
$\partial R/\partial a = 0$, hence $\Gamma_n = 1/0 = +\infty$. $\square$

But $T_{\min}$ is not strict Archimedean — no generator exists.

### Among strict Archimedean t-norms: power profile $h(x) = Cx^{-\alpha}$

$$\Gamma_n(t,a) = \frac{1}{n-1}\left(\frac{a}{t}\right)^\alpha$$

**Monotonicity in $\alpha$:** Let $r = a/t > 1$. Then:
$$\frac{\partial \Gamma_n}{\partial \alpha} = \frac{r^\alpha \ln r}{n-1} > 0$$
since $r > 1 \Rightarrow \ln r > 0$. So $\Gamma_n$ is strictly increasing in $\alpha$. $\square$

| $\alpha$ | T-norm | $\Gamma_n(t,a)$ |
|---|---|---|
| 1 | Product $T_P = xy$ | $a\,/\,(n-1)t$ |
| 2 | Power | $(1/(n-1))(a/t)^2$ |
| $\alpha\to\infty$ | $\to T_{\min}$ | $\to+\infty$ |

### Why $\alpha\to\infty$ is infeasible

The regularizer at $\beta=0$:
$$J_\text{reg}(\alpha,0) = \alpha\ln\frac{1}{t_{\min}} \leq M \implies \alpha \leq \frac{M}{\ln(1/t_{\min})} < +\infty$$

This bounds $\alpha$ from above, making the 1D optimization over $\alpha$ non-trivial.

---

## Tradeoff Summary

| | $T_{\min}$ | Power, finite $\alpha$ | $\alpha\to\infty$ |
|---|---|---|---|
| $\Gamma_n$ | $+\infty$ | $(a/t)^\alpha/(n-1)$ | $+\infty$ |
| $J_\text{reg}$ | undefined | $\alpha\ln(1/t_{\min})$ | $+\infty$ — **violates** $\leq M$ |
| Class | not strict Archimedean | strict Archimedean | boundary |

---

## Foundations in Literature

| Result | Source |
|---|---|
| Generator representation $T(x,y)=g^{-1}(g(x)+g(y))$ | Klement, Mesiar, Pap — *Triangular Norms*, Kluwer, 2000 |
| Properties of strict/Archimedean t-norms | Klement, Mesiar, Pap — Position Papers I–III, *FSS* 2004 |
| Generator reconstruction from partial derivatives | [IEEE Xplore, IPMU'08](https://ieeexplore.ieee.org/document/5489144/) |
| $\partial R/\partial x_k = h(x_k)/h(R)$ | Implicit differentiation of $g(R)=\sum g(x_i)$ — corollary of above |

### Original contributions (Ryzhov & Fedotov 2026)

- Metrics $\Gamma_n$, $D_n$, $B_n$ as t-norm quality characteristics
- Power profile as canonical selectivity invariant
- Inverse aggregator selection problem
- Reduction to 1D optimization via estimation bounds

---

## Standard SNR Metrics for Reference

| Metric | Formula | Domain |
|---|---|---|
| SNR | $\mu_s/\sigma_n$ or $P_s/P_n$ | general |
| SNR (dB) | $10\log_{10}(P_s/P_n)$ | engineering |
| d′ (d-prime) | $(\mu_s-\mu_n)/\sigma$ | signal detection theory |
| Mahalanobis distance | $\sqrt{(\mu_s-\mu_n)^T\Sigma^{-1}(\mu_s-\mu_n)}$ | multivariate |
| AUC ROC | area under ROC curve | classification |

$\Gamma_n$ is closest in spirit to **d′**: both measure local contrast between signal and background,
normalized by the background level.
