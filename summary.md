# Summary: Fuzzy Risk Aggregation on Trees Using t-norms

**Source documents (Russian):**
- `orig/Нечёткая_агрегация_рисков_на_деревьях.pdf` — Main thesis draft (35 pp.)
- `orig/Презентация_последняя.pdf` — Conference presentation, Lomonosov Readings 2026, MSU, April 1 2026. Authors: Ryzhov A.P., Fedotov F.A.
- `orig/Черновой вариант.pdf` — Extended draft of the thesis paper (58 pp.)

---

## Overview

The work develops a mathematical framework for aggregating fuzzy risk scores in hierarchical (tree-structured) models. A rooted tree represents an organization or system: leaves hold atomic risk scores $x_i \in [0,1]$ (higher = safer), and each internal node aggregates its children using the same binary operator. The goal is to understand how leaf-level changes propagate to the root score $R$, and — crucially — how to *design* the aggregation operator to achieve desired sensitivity/robustness trade-offs.

---

## 1. The Model

A **t-norm** $T:[0,1]^2\to[0,1]$ is used at every internal node. A t-norm satisfies:
- **Commutativity**: $T(x,y)=T(y,x)$
- **Associativity**: $T(T(x,y),z)=T(x,T(y,z))$
- **Monotonicity**: $x_1\le x_2 \Rightarrow T(x_1,y)\le T(x_2,y)$
- **Neutral element**: $T(x,1)=x$

The work focuses on **strict Archimedean t-norms**, satisfying additionally $T(x,x)<x$ for $x\in(0,1)$ and strict monotonicity. This class admits the **generator representation**:
$$T(x,y)=g^{-1}(g(x)+g(y))$$
where $g:(0,1]\to[0,+\infty)$ is continuous, strictly decreasing, $g(1)=0$, $g(0^+)=+\infty$.

**Classical examples:**
| t-norm | Generator $g(x)$ | Notes |
|---|---|---|
| Product $T_\Pi(x,y)=xy$ | $-\ln x$ | Strict Archimedean |
| Łukasiewicz $T_L(x,y)=\max(0,x+y-1)$ | $1-x$ | Archimedean, NOT strict (nilpotent) |
| Gödel/Minimum $T_{\min}(x,y)=\min(x,y)$ | — | NOT Archimedean |
| Dombi $T_\lambda$, $\lambda>0$ | $\left(\frac{1-x}{x}\right)^\lambda$ | Strict Archimedean, 1-parameter family |

**Scale ambiguity**: $g$ and $cg$ ($c>0$) generate the same t-norm, so generators are unique only up to a positive scalar.

---

## 2. Tree Collapse Formula

**Key theorem**: If the same strict Archimedean t-norm with generator $g$ is used at every internal node, the root value depends only on the *multiset* of leaf values — not on tree topology:

$$R(x_1,\ldots,x_n) = g^{-1}\!\left(\sum_{i=1}^n g(x_i)\right)$$

This "collapse" to a sum in generator scale is the central technical tool.

---

## 3. Sensitivity Analysis

The generator framework yields a clean partial derivative formula. Let $h(x):=-g'(x)>0$ (the **profile** of the generator). Then:

$$\frac{\partial R}{\partial x_k} = \frac{h(x_k)}{h(R)}$$

### Homogeneous point metrics

At the uniform configuration $(a,\ldots,a)$, define $R_n(a)=g^{-1}(ng(a))$ and:
- **Single-leaf sensitivity**: $\sigma_n(a) = h(a)/h(R_n(a))$
- **Synchronous noise sensitivity**: $\rho_n(a) = n\cdot\sigma_n(a)$

These two metrics are *not independent* — they always satisfy $\rho_n = n\sigma_n$. This means the homogeneous setting cannot distinguish "incident vs. background noise."

---

## 4. Non-Homogeneous (Incident) Scenario

To distinguish an incident from background noise, the work introduces the configuration $(t,a,\ldots,a)$ with $0<t<a<1$: one "incident" leaf at level $t$, all others at background level $a$. The root becomes:
$$R(t,a)=g^{-1}(g(t)+(n-1)g(a))$$

Three complementary metrics are defined:

| Metric | Definition | Meaning |
|---|---|---|
| **Detection** $D_n(t;a)$ | $R_n(a)-R(t,a)>0$ | Finite drop of root score due to incident |
| **Incident sensitivity** $I_n(t,a)$ | $\|\partial R/\partial t\|=h(t)/h(R(t,a))$ | Local reaction to incident leaf |
| **Background sensitivity** $B_n(t,a)$ | $\|\partial R/\partial a\|=(n-1)h(a)/h(R(t,a))$ | Local reaction to background drift |
| **Selectivity** $\Gamma_n(t,a)$ | $I_n/B_n = h(t)/[(n-1)h(a)]$ | Signal-to-noise ratio |

**Crucial property of selectivity**: The denominator $h(R(t,a))$ cancels out, so $\Gamma_n(t,a)$ depends only on the *shape* of the profile $h$, not on the root value or generator scale.

---

## 5. Inverse Design Problem

The forward analysis motivates the **inverse problem**: given desired sensitivity/robustness requirements, *design* a generator $g$ (and hence a t-norm) that satisfies them.

**Working domain**: $\Omega=\{(t,a): a\in[a_-,1),\, t\in[t_{\min},a)\}$. A refined version $\Omega_\eta$ adds the constraint $t\le\eta a$ ($\eta\in(0,1)$) to ensure a non-trivial incident depth.

**Objective**: Minimize worst-case background sensitivity $\sup_\Omega B_n(t,a)$.

**Constraints**:
1. Detection floor: $\inf_\Omega D_n(t;a)\ge D_*$
2. Selectivity floor: $\inf_\Omega \Gamma_n(t,a)\ge\Gamma_*$
3. Regularity bound: $J_{\text{reg}}(h):=\int_{t_{\min}}^1|(\log h)'|\,dx\le M$
4. Global robustness: profile $h$ non-increasing ensures $0\le\partial R/\partial x_k\le 1$ and $|R(x)-R(y)|\le\sum|x_k-y_k|$ (1-Lipschitz in $\ell_1$)

This is an **infinite-dimensional variational problem** (the unknown is a function $h$), justified in spirit but not directly computable.

---

## 6. Parametric Families

### Classical families (benchmarks)

**Product ($q=1$ case)**:
$$R_\Pi(t,a)=ta^{n-1},\quad D_n=a^{n-1}(a-t),\quad \Gamma_n=\frac{a}{(n-1)t}$$
No free parameters — fixed trade-off.

**Dombi family** (parameter $\lambda>0$):
$$R_\lambda(t,a)=\frac{1}{1+S_\lambda(t,a)^{1/\lambda}},\quad \Gamma_n^{(\lambda)}=\frac{1}{n-1}\left(\frac{1-t}{1-a}\right)^{\lambda-1}\!\!\left(\frac{a}{t}\right)^{\lambda+1}$$
One parameter $\lambda$ controls everything simultaneously — selectivity strictly increases with $\lambda$, but detection $D_n$ does not improve monotonically. **Limitation**: a single parameter cannot independently control noise, detection, and selectivity.

### Novel parametric family (own construction)

To decouple the roles of parameters, the authors introduce a **piecewise power-law profile**:
$$h_{p,q,\tau}(x)=\begin{cases}x^{-p}, & 0<x\le\tau\\ \tau^{q-p}x^{-q}, & \tau<x\le 1\end{cases}$$
with $p>1$, $q>0$, $\tau\in(0,t_{\min}]$.

**Parameter roles**:
- $p$ — singularity strength near zero (ensures $g(0^+)=+\infty$, i.e., strict Archimedeanism)
- $q$ — profile slope in the working zone $[t_{\min},1]$ → controls selectivity and noise
- $\tau$ — switching point between the two regimes

**Key results** (when $\tau\le t_{\min}$, so the working zone uses only the right branch):

$$\Gamma_n^{(\theta)}(t,a) = \frac{1}{n-1}\left(\frac{a}{t}\right)^q \qquad \text{(selectivity depends only on } q\text{)}$$

$$J_{\text{reg}}(\theta) = q\ln\frac{1}{t_{\min}} \qquad \text{(regularizer depends only on } q\text{)}$$

$$J_{\text{noise}}(\theta) = (n-1)\eta^{-q} \quad (q\ne 1), \qquad (n-1)\eta \quad (q=1)$$

All three working-zone functionals reduce to functions of $q$ alone, giving a clean 1D optimization on the working domain.

---

## 7. Compatibility Theorems

### Selectivity vs. regularization

Conditions $J_\Gamma\ge\Gamma_*$ and $J_{\text{reg}}\le M$ are jointly satisfiable **if and only if**:
$$\frac{\ln((n-1)\Gamma_*)}{\ln(1/\eta)} \le \frac{M}{\ln(1/t_{\min})}$$

This gives a direct quantitative condition on problem parameters $(\Gamma_*, M, \eta, t_{\min})$ before any specific profile is chosen.

### Selectivity vs. noise

For $q\ne 1$, conditions $J_\Gamma\ge\Gamma_*$ and $J_{\text{noise}}\le B_*$ are jointly satisfiable **if and only if**:
$$(n-1)^2\Gamma_* \le B_*$$

This shows that increasing selectivity necessarily increases worst-case noise — the trade-off is *fundamental*, not just a limitation of the parametric family.

### Existence of optimum

Under the above compatibility conditions, the finite-dimensional optimization problem on the compact box $\Theta_0=[p_-,p_+]\times[q_-,q_+]\times[\tau_-,\tau_+]$ has a solution by Weierstrass's theorem (all functionals are continuous in $\theta$).

---

## 8. Design Axioms Approach (thesis draft)

The thesis draft explores an alternative approach: deriving specific generator classes from *design axioms* on sensitivity decay.

### Scale-invariance axiom → Schweizer–Sklar family
Require that marginal sensitivity decays as a power law when total risk is scaled. This uniquely forces:
$$f(s)=(1+s)^{-\alpha}, \quad g(x)=x^{-1/\alpha}-1, \quad T(x,y)=(x^{-\lambda}+y^{-\lambda}-1)^{-1/\lambda}$$
giving the Schweizer–Sklar t-norm. Metrics: $\sigma_n(a_0)=\left(\frac{2}{n+1}\right)^{\alpha+1}$ — polynomial decay in $n$.

### Shift-invariance axiom → Product t-norm
Require that adding a fixed increment to total risk reduces marginal sensitivity by a constant factor (memoryless property). This uniquely forces $f(s)=e^{-ks}$, hence $g(x)=-\frac{1}{k}\ln x$, i.e., the **product t-norm**. Metrics: $\sigma_n(a_0)=a_0^{n-1}$ — exponential decay in $n$.

This axiomatic approach shows that "natural" behavioral requirements uniquely determine the generator class.

---

## 9. Main Conclusions

1. **Tree collapse**: For any strict Archimedean t-norm applied uniformly, the root is $R=g^{-1}(\sum g(x_i))$ — tree topology is irrelevant.

2. **Selectivity formula**: $\Gamma_n(t,a)=h(t)/[(n-1)h(a)]$ depends only on the profile shape, not on absolute generator values. This is the central design handle.

3. **Homogeneous metrics are insufficient**: $\rho_n=n\sigma_n$ always — incident and noise cannot be separated in the homogeneous regime. The non-homogeneous scenario $(t,a,\ldots,a)$ is essential.

4. **Fundamental trade-off**: In the proposed parametric family (and more generally), improving selectivity increases worst-case noise sensitivity. The bound $(n-1)^2\Gamma_*\le B_*$ is a necessary condition for feasibility.

5. **Power-law profile as canonical etalon**: Selectivity depends only on $t/a$ if and only if $h(x)=Cx^{-\alpha}$. This gives the Schweizer–Sklar family as the canonical reference.

6. **Existence of optimal generator**: Under compatibility conditions, the finite-dimensional inverse problem has a solution. Selectivity and regularization are controlled by a single parameter $q$, enabling explicit compatibility analysis.
