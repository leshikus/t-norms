# T-norms with Generators in Elementary Functions

A t-norm $T:[0,1]^2\to[0,1]$ is **Archimedean** if $T(x,x)<x$ for all $x\in(0,1)$.
Every continuous Archimedean t-norm has an **additive generator** $g:(0,1]\to[0,\infty)$:
$$T(x,y)=g^{-1}\!\bigl(\min(g(x)+g(y),\,g(0))\bigr),$$
where $g$ is continuous, strictly decreasing, $g(1)=0$.

- **Strict** ($g(0^+)=+\infty$): $T(x,y)=g^{-1}(g(x)+g(y))$.
- **Nilpotent** ($g(0)<\infty$): $T(x,y)=\max(0,\ldots)$; there exist $x,y>0$ with $T(x,y)=0$.

Generators are unique up to a positive scalar ($g$ and $cg$ generate the same $T$).

**Elementary functions** here means: algebraic combinations of polynomials, $e^x$, $\ln x$, trigonometric and inverse trigonometric functions, and their compositions.

---

## Complete table of named families

| Family | Param. | Generator $g(x)$ | T-norm $T(x,y)$ | Type |
|--------|--------|-------------------|-----------------|------|
| **Łukasiewicz** | — | $1-x$ | $\max(0,x+y-1)$ | nilpotent |
| **Product** | — | $-\ln x$ | $xy$ | strict |
| **Schweizer–Sklar** | $\lambda\in\mathbb{R}$ | $1-x^\lambda$ ($\lambda>0$); $x^\lambda-1$ ($\lambda<0$) | $\bigl[\max(0,x^\lambda+y^\lambda-1)\bigr]^{1/\lambda}$ | nilp./strict |
| **Yager** | $\lambda\ge 1$ | $(1-x)^\lambda$ | $\max\!\bigl(0,1-[(1-x)^\lambda+(1-y)^\lambda]^{1/\lambda}\bigr)$ | nilpotent |
| **Aczel–Alsina** | $\lambda\ge 1$ | $(-\ln x)^\lambda$ | $\exp\!\bigl(-[(-\ln x)^\lambda+(-\ln y)^\lambda]^{1/\lambda}\bigr)$ | strict |
| **Dombi** | $\lambda>0$ | $\left(\dfrac{1-x}{x}\right)^{\!\lambda}$ | $\dfrac{1}{1+\left[\!\left(\dfrac{1-x}{x}\right)^{\!\lambda}+\left(\dfrac{1-y}{y}\right)^{\!\lambda}\right]^{1/\lambda}}$ | strict |
| **Frank** | $s>0,s\ne 1$ | $\ln\dfrac{s-1}{s^x-1}$ | $\log_s\!\left(1+\dfrac{(s^x-1)(s^y-1)}{s-1}\right)$ | strict |
| **Hamacher** | $\gamma>0$ | $\ln\dfrac{\gamma+(1-\gamma)x}{x}$ | $\dfrac{xy}{\gamma+(1-\gamma)(x+y-xy)}$ | strict |
| **Sugeno–Weber** | $\lambda>-1$ | $\ln\dfrac{1+\lambda}{1+\lambda x}$ | $\max\!\left(0,\dfrac{x+y-1+\lambda xy}{1+\lambda}\right)$ | nilpotent |

Additionally:
- **Min** $T_M(x,y)=\min(x,y)$: not Archimedean, no generator; arises as a limit of all strict families as $\lambda\to\infty$ (Aczel–Alsina, Dombi, Frank).

---

## Derivations

### Schweizer–Sklar

**For $\lambda<0$** (strict): $g_\lambda(x)=x^\lambda-1$.
$g_\lambda(1)=0$ ✓; $g_\lambda(0^+)=+\infty$ ✓.
$g_\lambda^{-1}(s)=(s+1)^{1/\lambda}$. Result: $T=(x^\lambda+y^\lambda-1)^{1/\lambda}$.

**For $\lambda>0$** (nilpotent): $g_\lambda(x)=1-x^\lambda$.
$g_\lambda(0)=1$ (finite) ✓.
$g_\lambda^{-1}(s)=(1-s)^{1/\lambda}$. Result: $T=\max(0,(x^\lambda+y^\lambda-1)^{1/\lambda})$.

Both cases give the same formula, unified as above.

### Yager

$g_\lambda(x)=(1-x)^\lambda$, $g_\lambda(0)=1$ (nilpotent), $g_\lambda^{-1}(s)=1-s^{1/\lambda}$.
$$T=\max\!\bigl(0,\,1-[(1-x)^\lambda+(1-y)^\lambda]^{1/\lambda}\bigr).$$

### Aczel–Alsina

$g_\lambda(x)=(-\ln x)^\lambda$, $g_\lambda(0^+)=+\infty$ (strict), $g_\lambda^{-1}(s)=e^{-s^{1/\lambda}}$.
$$T=\exp\!\bigl(-[(-\ln x)^\lambda+(-\ln y)^\lambda]^{1/\lambda}\bigr).$$

### Dombi

$g_\lambda(x)=\left(\frac{1-x}{x}\right)^\lambda$, $g_\lambda^{-1}(s)=\frac{1}{1+s^{1/\lambda}}$. See `dombi.md` for full analysis.

### Frank

$g_s(x)=\ln\frac{s-1}{s^x-1}$ for $s>1$ (and similarly for $s\in(0,1)$ by continuity). Verify $g_s(1)=\ln 1=0$ ✓; $g_s(0^+)=+\infty$ ✓ (since $s^0-1\to0^+$).

$g_s^{-1}(u)=\log_s(1+(s-1)e^{-u})$.

$$T_s(x,y)=g_s^{-1}(g_s(x)+g_s(y))=\log_s\!\left(1+\frac{(s^x-1)(s^y-1)}{s-1}\right).$$

### Hamacher

Try $g_\gamma(x)=\ln\frac{\gamma+(1-\gamma)x}{x}$.
$g_\gamma(1)=\ln 1=0$ ✓; $g_\gamma(0^+)=+\infty$ ✓ (for $\gamma>0$).
$g_\gamma^{-1}(u)=\frac{\gamma}{e^u-(1-\gamma)}$.

$$T=\frac{\gamma}{\frac{(\gamma+(1-\gamma)x)(\gamma+(1-\gamma)y)}{xy}-(1-\gamma)}$$
$$=\frac{\gamma xy}{(\gamma+(1-\gamma)x)(\gamma+(1-\gamma)y)-(1-\gamma)xy}.$$

Expanding the denominator: $\gamma^2+\gamma(1-\gamma)(x+y)+(1-\gamma)^2xy-(1-\gamma)xy = \gamma[\gamma+(1-\gamma)(x+y-xy)]$. Thus:
$$T_\gamma^H(x,y)=\frac{xy}{\gamma+(1-\gamma)(x+y-xy)}. \qquad \checkmark$$

### Sugeno–Weber

**Goal**: find $g$ with $g^{-1}(g(x)+g(y))=\frac{x+y-1+\lambda xy}{1+\lambda}$.

Set $f=g^{-1}$. The equation $f(u+v)=\frac{f(u)+f(v)-1+\lambda f(u)f(v)}{1+\lambda}$ must hold.

Let $F(u)=1+\lambda f(u)$. Then:
$$F(u+v)=1+\lambda f(u+v)=\frac{(1+\lambda f(u))(1+\lambda f(v))}{1+\lambda}=\frac{F(u)F(v)}{1+\lambda}.$$

This is a **multiplicative Cauchy equation**. The continuous solution with $f(0)=1$ (i.e., $F(0)=1+\lambda$) is $F(u)=(1+\lambda)r^u$ for some $r\in(0,1)$.

Choosing $r=e^{-1}$: $f(u)=\frac{(1+\lambda)e^{-u}-1}{\lambda}$, so $g(x)=f^{-1}(x)$ gives:
$$\boxed{g_\lambda^{SW}(x)=\ln\frac{1+\lambda}{1+\lambda x}.}$$
$g(1)=0$ ✓; $g(0)=\ln(1+\lambda)<\infty$ (nilpotent) ✓.

Verify: $g(x)+g(y)=\ln\frac{(1+\lambda)^2}{(1+\lambda x)(1+\lambda y)}$, so
$e^{-(g(x)+g(y))}=\frac{(1+\lambda x)(1+\lambda y)}{(1+\lambda)^2}$, and
$$T=\frac{(1+\lambda)\cdot\frac{(1+\lambda x)(1+\lambda y)}{(1+\lambda)^2}-1}{\lambda}=\frac{(1+\lambda x)(1+\lambda y)-(1+\lambda)}{\lambda(1+\lambda)}=\frac{x+y-1+\lambda xy}{1+\lambda}. \qquad \checkmark$$

---

## Limit cases and special values

| Condition | Result |
|-----------|--------|
| Schweizer–Sklar $\lambda=1$ | $T=(x+y-1)_+$ = Łukasiewicz |
| Schweizer–Sklar $\lambda=0$ | $T=xy$ = Product |
| Schweizer–Sklar $\lambda=-1$ | $T=xy/(x+y-xy)$ = Hamacher ($\gamma=0$) = Dombi ($\lambda=1$) |
| Schweizer–Sklar $\lambda\to-\infty$ | $T\to\min(x,y)$ |
| Yager $\lambda=1$ | Łukasiewicz |
| Yager $\lambda\to\infty$ | $T\to\min(x,y)$ |
| Aczel–Alsina $\lambda=1$ | $g(x)=-\ln x$, $T=xy$ = Product |
| Aczel–Alsina $\lambda\to\infty$ | $T\to\min(x,y)$ |
| Dombi $\lambda\to\infty$ | $T\to\min(x,y)$ |
| Dombi $\lambda\to 0^+$ | $T\to$ Drastic product |
| Frank $s=1$ (limit) | $T=xy$ = Product |
| Frank $s\to 0^+$ | $T\to T_L$ = Łukasiewicz |
| Frank $s\to\infty$ | $T\to\min(x,y)$ |
| Hamacher $\gamma=1$ | $g(x)=-\ln x$, $T=xy$ = Product |
| Hamacher $\gamma\to 0$ | $g(x)\to(1-x)/x$ (rescaled), $T\to xy/(x+y-xy)$ |
| Hamacher $\gamma=2$ | Einstein product $T=xy/(1+(1-x)(1-y))$ |
| Sugeno–Weber $\lambda=0$ | $g(x)\to(1-x)$ (rescaled), $T=\max(0,x+y-1)$ = Łukasiewicz |
| Sugeno–Weber $\lambda\to\infty$ | $g(x)\to -\ln x$ (rescaled), $T\to xy$ = Product |

---

## Relationships between families

The families are not disjoint — several coincide at specific parameter values:

```
Min ──────────────────────────────────── (limit of all strict families)
 │
 ├─ Frank (s→∞)    ←─ contains {Łukasiewicz, Product, Min} as limits
 │      └─ Product (s=1), Łukasiewicz (s→0)
 │
 ├─ Schweizer–Sklar (λ→−∞)
 │      └─ Product (λ=0), Łukasiewicz (λ=1), Hamacher γ=0 (λ=−1)
 │
 ├─ Aczel–Alsina (λ→∞)
 │      └─ Product (λ=1)
 │
 └─ Dombi (λ→∞)
        └─ Hamacher γ=0 (λ=1)
```

The **Frank family** is unique in being the only family satisfying $T(x,y)+S(x,y)=x+y$ (Frank 1979), where $S$ is the dual t-conorm.

The **Hamacher family** is identical to the **Ali–Mikhail–Haq (AMH)** family under the reparametrization $\theta=1-\gamma$:
$$T^{AMH}_\theta(x,y)=\frac{xy}{1-\theta(1-x)(1-y)}, \qquad \theta\in(-\infty,1).$$

The **Schweizer–Sklar** family at $\lambda<0$ corresponds to the **Clayton copula** (under the substitution $U=-\ln X$): $C_\theta(u,v)=(u^{-\theta}+v^{-\theta}-1)^{-1/\theta}$ with $\theta=-\lambda>0$.

The **Aczel–Alsina** family corresponds to the **Gumbel–Hougaard copula** under the same substitution.

---

## Profile $h=-g'$ for each family

The profile $h=-g'$ controls sensitivity via $\partial R/\partial x_k = h(x_k)/h(R)$.

| Family | $h(x) = -g'(x)$ | Monotone? |
|--------|-----------------|-----------|
| Łukasiewicz | $1$ (constant) | yes (flat) |
| Product | $1/x$ | decreasing ✓ |
| Schweizer–Sklar $\lambda<0$ | $-\lambda x^{\lambda-1}$ | decreasing iff $\lambda<0$ (always here) ✓ |
| Schweizer–Sklar $\lambda>0$ | $\lambda x^{\lambda-1}$ | increasing ✗ (robustness fails) |
| Yager | $\lambda(1-x)^{\lambda-1}$ | decreasing iff $\lambda\ge 1$ ✓ |
| Aczel–Alsina | $\lambda(-\ln x)^{\lambda-1}/x$ | decreasing iff $\lambda\ge 1$ ✓ |
| Dombi | $\lambda(1-x)^{\lambda-1}/x^{\lambda+1}$ | decreasing iff $\lambda\ge 1$ ✓ (see `dombi.md`) |
| Frank | $\frac{s^x\ln s}{(s^x-1)}$ (up to scale) | decreasing (for $s>1$) ✓ |
| Hamacher | $\frac{1}{x(\gamma+(1-\gamma)x)}$ | decreasing (for $\gamma\ge 0$) ✓ |
| Sugeno–Weber | $\frac{\lambda}{1+\lambda x}$ | decreasing ✓ |

A non-increasing profile guarantees $|\partial R/\partial x_k|\le 1$ and 1-Lipschitz stability.

---

## Selectivity formula for each family

In the incident scenario $(t,a,\ldots,a)$, $\Gamma_n = h(t)/[(n-1)h(a)]$:

| Family | $\Gamma_n(t,a)$ |
|--------|-----------------|
| Product | $\dfrac{a}{(n-1)t}$ |
| Schweizer–Sklar $\lambda<0$ | $\dfrac{1}{n-1}\left(\dfrac{a}{t}\right)^{1-\lambda}$ |
| Yager ($\tau=1-t$, $\alpha=1-a$) | $\dfrac{1}{n-1}\left(\dfrac{1-t}{1-a}\right)^{\lambda-1}$ |
| Dombi | $\dfrac{1}{n-1}\left(\dfrac{1-t}{1-a}\right)^{\!\lambda-1}\!\!\left(\dfrac{a}{t}\right)^{\!\lambda+1}$ |
| Aczel–Alsina | $\dfrac{1}{n-1}\left(\dfrac{\ln(1/t)}{\ln(1/a)}\right)^{\lambda-1}\cdot\dfrac{\ln(1/a)}{\ln(1/t)}\cdot\dfrac{a}{t} = \dfrac{1}{n-1}\left(\dfrac{\ln(1/t)}{\ln(1/a)}\right)^{\lambda-1}\dfrac{a}{t}$ |
| Hamacher | $\dfrac{a(\gamma+(1-\gamma)a)}{(n-1)t(\gamma+(1-\gamma)t)}$ |
| Sugeno–Weber | $\dfrac{1+\lambda a}{(n-1)(1+\lambda t)}$ |

For **power-law profile** $h(x)=Cx^{-q}$ (the case $\Gamma$ depends only on ratio $a/t$): $\Gamma_n=\frac{1}{n-1}(a/t)^q$. This uniquely characterizes **Schweizer–Sklar** ($q=1-\lambda$).

---

## Novel: piecewise power-law family (from the studied papers)

The papers under study introduce a 3-parameter family:

$$h_{p,q,\tau}(x)=\begin{cases}x^{-p}, & 0<x\le\tau\\ \tau^{q-p}\,x^{-q}, & \tau<x\le 1\end{cases}$$

with $p>1$, $q>0$, $\tau\in(0,t_{\min}]$. The generator is obtained by integration:

$$g(x)=\int_x^1 h_{p,q,\tau}(u)\,du \quad \text{(elementary, piecewise power-law)}.$$

When $\tau\le t_{\min}$, only the right branch $h(x)=\tau^{q-p}x^{-q}$ is active on the working zone, giving:
- $\Gamma_n(t,a) = \frac{1}{n-1}(a/t)^q$ — depends only on $q$
- $J_{\mathrm{reg}} = q\ln(1/t_{\min})$ — depends only on $q$
- Underlying t-norm is **Schweizer–Sklar** on the working zone (up to a scale change in $g$)

The $p$ parameter governs only the behavior near $x=0$ (ensuring $g(0^+)=+\infty$ for strict Archimedeanism), while $q$ controls all design-relevant behavior.

---

## On completeness

There is **no theorem** characterizing all t-norms with elementary generators — the set is dense (any sufficiently smooth $g$ that is positive, decreasing, with $g(1)=0$ and $g(0^+)=+\infty$ defines a valid t-norm, and one can always add $\sin$-perturbations to get elementary $g$ with non-elementary $T$).

The 9 families above represent **all standard named families** in the literature (cf. Klement–Mesiar–Pap, *Triangular Norms*, 2000) where both $g$ **and** $T$ admit compact elementary closed forms. The key structural reason each works: the functional equation $f(u+v) = \Phi(f(u), f(v))$ for a rational or exponential $\Phi$ admits power-law or exponential solutions $f$.

---

## Holomorphicity of generators

### At $x=1$

$g(1)=0$ for all families, so the relevant question is whether $g$ extends to a complex-analytic function in a neighborhood of $x=1$.

The three families that fail for non-integer $\lambda$ are **Yager**, **Aczel–Alsina**, and **Dombi** — all share the local form $(1-x)^\lambda$ near $x=1$, which has a branch point there:

$$\underbrace{(1-x)^\lambda}_{\text{Yager}}, \quad \underbrace{(-\ln x)^\lambda \approx (1-x)^\lambda}_{\text{Aczel–Alsina}}, \quad \underbrace{\left(\tfrac{1-x}{x}\right)^\lambda \approx (1-x)^\lambda}_{\text{Dombi}}$$

By contrast, Schweizer–Sklar has $x^\lambda = e^{\lambda\ln x}$, whose branch point is at $x=0$, not $x=1$, so it is holomorphic at $x=1$ for all $\lambda$.

All other families (Łukasiewicz, Product, Frank, Hamacher, Sugeno–Weber) are holomorphic at $x=1$ unconditionally.

### At $x=0$

Strict generators satisfy $g(0^+)=+\infty$, so they all have singularities at $x=0$ (logarithmic or pole/branch-point type). Among nilpotent generators:

- **Łukasiewicz**, **Yager**, **Sugeno–Weber**: holomorphic at $x=0$ ✓
- **Schweizer–Sklar** $\lambda>0$: $g(x)=1-x^\lambda$ has a branch point at $x=0$ for non-integer $\lambda$ ✗

### Logit factorization: $g = h \circ b$ with holomorphic $h$

All generators can be expressed as $h(b(x))$ with a **single fixed** $b$ and holomorphic $h$:

$$\boxed{b(x) = \ln\frac{x}{1-x} \quad\text{(logit)}, \qquad a = \mathrm{id}.}$$

**Why.** Every generator is analytic on the open interval $(0,1)$ — all $\ln$-arguments are positive, all bases of $(\cdot)^\lambda$ are in $(0,1)$. The singularities live only at the boundary. The logit pushes both endpoints to $\pm\infty$:
$$b(0^+)=-\infty, \qquad b(1^-)=+\infty,$$
so its inverse $b^{-1}=\sigma$ (sigmoid, $\sigma(u)=e^u/(1+e^u)$) maps all of $\mathbb{R}$ into $(0,1)$. Then $h=g\circ\sigma$ has no singularity at any finite point:

| Family | $h(u) = g(\sigma(u))$ | Holomorphic on $\mathbb{R}$? |
|--------|----------------------|------------------------------|
| Łukasiewicz | $(1+e^u)^{-1}$ | ✓ |
| Product | $\ln(1+e^{-u})$ (softplus) | ✓ |
| Schweizer–Sklar | $1-(1+e^{-u})^{-\lambda}$ or $(1+e^{-u})^{-\lambda}-1$ | ✓ |
| Yager | $(1+e^u)^{-\lambda}$ | ✓ ($1+e^u>0$) |
| Aczel–Alsina | $(\ln(1+e^{-u}))^\lambda$ | ✓ ($\ln(1+e^{-u})>0\ \forall u$) |
| Dombi | $e^{-\lambda u}$ | ✓ (entire) |
| Frank | $\ln\frac{s-1}{s^{\sigma(u)}-1}$ | ✓ ($s^{\sigma(u)}\ne 1$ for finite $u$) |
| Hamacher | $\ln\frac{\gamma+(1-\gamma)\sigma(u)}{\sigma(u)}$ | ✓ |
| Sugeno–Weber | $\ln\frac{1+\lambda}{1+\lambda\sigma(u)}$ | ✓ |

The critical case is Aczel–Alsina: $1+e^{-u}>1$ for all real $u$, so $\ln(1+e^{-u})>0$ everywhere and $(\cdot)^\lambda$ never encounters its branch point at $0$.

**General principle.** Any bijection $b:(0,1)\to\mathbb{R}$ with analytic inverse suffices (e.g.\ $b(x)=\tan(\pi(x-\tfrac{1}{2}))$). The logit is the canonical choice. The key property is that both boundary singularities ($x=0$ and $x=1$) are pushed to $\pm\infty$, leaving $h$ singularity-free on $\mathbb{R}$.

### Converse: does holomorphic $h$ produce a valid generator?

**No.** Holomorphicity is neither necessary nor sufficient. The implication runs only one way:
$$g \text{ valid generator} \implies h = g\circ\sigma \text{ holomorphic on } \mathbb{R}.$$
The converse fails: e.g.\ $h(u)=e^u$ is entire, but $g(x)=x/(1-x)$ is strictly *increasing*.

For $g = h\circ\mathrm{logit}$ to be a valid additive generator, $h:\mathbb{R}\to[0,\infty)$ must satisfy:

| Condition on $h$ | Corresponds to |
|---|---|
| $h$ strictly decreasing | $g$ strictly decreasing |
| $h(u)>0$ for all finite $u$ | $g>0$ on $(0,1)$ |
| $\lim_{u\to+\infty}h(u)=0$ | $g(1)=0$ |
| $\lim_{u\to-\infty}h(u)\in(0,+\infty]$ | $g(0^+)>0$ (nilpotent or strict) |

These four conditions are necessary and sufficient (holomorphicity is independent extra regularity; continuity alone suffices).

**Entireness forces strictness.** If $h$ is entire (holomorphic on all of $\mathbb{C}$) and satisfies the four conditions, Liouville's theorem forces $h$ to be unbounded, hence $h(-\infty)=+\infty$, giving a **strict** generator. Nilpotent generators require $h(-\infty)<\infty$, so their $h$ cannot be entire: Łukasiewicz $h=(1+e^u)^{-1}$ has poles at $u=i\pi(2k+1)$; Yager, Sugeno–Weber similarly.

### Power series $h(u)=\sum_{n\ge 0}a_n u^n$ around $u=0$

$u=0$ corresponds to $x=\sigma(0)=\tfrac{1}{2}$. Use the sigmoid expansion
$$\sigma(u)=\tfrac{1}{2}+\tfrac{1}{4}u-\tfrac{1}{48}u^3+\tfrac{1}{480}u^5-\cdots$$

**Łukasiewicz** $h(u)=(1+e^u)^{-1}$:
$$\tfrac{1}{2}-\tfrac{1}{4}u+\tfrac{1}{48}u^3-\tfrac{1}{480}u^5+\cdots$$
Even-degree terms $\ge 2$ vanish ($h(u)-\tfrac12$ is odd). General formula:
$$a_n = \frac{(-1)^n}{n!}\,\eta(-n), \qquad \eta(s)=(1-2^{1-s})\zeta(s)\text{ (Dirichlet eta)}.$$

**Product** $h(u)=\ln(1+e^{-u})$:
$$\ln 2-\tfrac{1}{2}u+\tfrac{1}{8}u^2-\tfrac{1}{192}u^4+\tfrac{1}{2880}u^6-\cdots$$
Odd-degree terms $\ge 3$ vanish (since $h'=-h_L$ and $h_L-\tfrac12$ is odd). Using $h_P'=-h_L$:
$$a_n = \frac{(-1)^{n+1}}{n!}\,\eta(n-1), \quad n\ge 1;\qquad a_0=\ln 2.$$

**Dombi** $h(u)=e^{-\lambda u}$:
$$1-\lambda u+\tfrac{\lambda^2}{2}u^2-\tfrac{\lambda^3}{6}u^3+\cdots$$
$$\boxed{a_n = \frac{(-\lambda)^n}{n!}.}$$

**Hamacher** $h(u)=\ln(1+\gamma e^{-u})$:
$$\ln(1+\gamma)-\frac{\gamma}{1+\gamma}u+\frac{\gamma}{2(1+\gamma)^2}u^2-\frac{\gamma(1-\gamma)}{6(1+\gamma)^3}u^3+\cdots$$
Exact formula via polylogarithm:
$$\boxed{a_n = \frac{(-1)^{n+1}}{n!}\operatorname{Li}_{1-n}(-\gamma),\quad n\ge 1;}\qquad a_0=\ln(1+\gamma).$$
where $\operatorname{Li}_s(z)=\sum_{k=1}^\infty z^k/k^s$. Equivalently, define $f_m(\gamma)=\sum_{k=1}^\infty(-1)^{k-1}\gamma^k k^m$ with recursion $f_m=\gamma f_{m-1}'$ and $f_0=\gamma/(1+\gamma)$; then $a_n=(-1)^n f_{n-1}(\gamma)/n!$. Product is the $\gamma=1$ special case.

**Sugeno–Weber** $h(u)=\ln\tfrac{1+\lambda}{1+\lambda\sigma(u)}$:
$$\ln\tfrac{2(1+\lambda)}{2+\lambda}-\frac{\lambda}{2(2+\lambda)}u+\frac{\lambda^2}{8(2+\lambda)^2}u^2-\cdots$$
Note $a_2/a_1 = -a_1/2$ (geometric-like start); the full series follows from $h(u) = \ln(1+\lambda)+\ln(1+e^u)-\ln(1+(1+\lambda)e^u)$, i.e., a difference of two Hamacher-type series with $\gamma=1$ and $\gamma=1+\lambda$.

**Schweizer–Sklar** $(\lambda>0)$, $h(u)=1-\sigma(u)^\lambda$; write $\sigma=\tfrac{1}{2}(1+\frac{u}{2}-\frac{u^3}{24}+\cdots)$:
$$\bigl(1-2^{-\lambda}\bigr)-\frac{\lambda}{2^{\lambda+1}}u-\frac{\lambda(\lambda-1)}{2^{\lambda+3}}u^2-\frac{\lambda^2(\lambda-3)}{3\cdot 2^{\lambda+4}}u^3-\cdots$$
In general $a_n=-2^{-\lambda}[\sigma^\lambda]_n$ where $[\sigma^\lambda]_n$ is the $u^n$ coefficient of $\exp(\lambda\ln\sigma(u))$ obtained by composing the $\ln\sigma$ series with $\exp$. No closed elementary form for general $n$.

**Yager** $h(u)=(1+e^u)^{-\lambda}=2^{-\lambda}(1+e^u)^{-\lambda}/2^{-\lambda}$:
$$2^{-\lambda}\Bigl(1-\tfrac{\lambda}{2}u+\tfrac{\lambda(\lambda-1)}{8}u^2+\tfrac{\lambda^2(3-\lambda)}{48}u^3+\cdots\Bigr)$$
Coefficients follow from $(1+w)^{-\lambda}$ with $w=(e^u-1)/2=\frac{u}{2}+\frac{u^2}{4}+\cdots$.

**Aczel–Alsina** $h(u)=(\ln(1+e^{-u}))^\lambda=h_P(u)^\lambda$; write $h_P=\ln 2+\Delta$ with $\Delta=-\frac{u}{2}+\frac{u^2}{8}-\cdots$:
$$(\ln 2)^\lambda-\frac{\lambda(\ln 2)^{\lambda-1}}{2}u+\frac{\lambda(\lambda-1+\ln 2)(\ln 2)^{\lambda-2}}{8}u^2-\cdots$$
General $a_n$ via $h_P^\lambda=(\ln 2)^\lambda(1+\Delta/\ln 2)^\lambda$, expanding by binomial series in $\Delta/\ln 2$.

---

## References

- **Schweizer, B., Sklar, A.** (1963). Associative functions and statistical triangle inequalities. *Publ. Math. Debrecen* 8, 169–186.
- **Frank, M.J.** (1979). On the simultaneous associativity of $F(x,y)$ and $x+y-F(x,y)$. *Aequationes Mathematicae* 19, 194–226.
- **Dombi, J.** (1982). A general class of fuzzy operators, the DeMorgan class of fuzzy operators and fuzziness measures induced by fuzzy operators. *Fuzzy Sets and Systems* 8(2), 149–163.
- **Aczél, J., Alsina, C.** (1982). Characterizations of some classes of quasilinear functions with applications to triangular norms and to synthesizing judgements. *Methods of Operations Research* 48, 3–22.
- **Yager, R.R.** (1980). On a general class of fuzzy connectives. *Fuzzy Sets and Systems* 4(3), 235–242.
- **Hamacher, H.** (1978). Über logische Verknüpfungen unscharfer Aussagen und deren zugehörige Bewertungsfunktionen. *Progress in Cybernetics and Systems Research* 3, 276–288.
- **Sugeno, M., Weber, S.** (1993). A note on an open problem and its solution. *Fuzzy Sets and Systems* 54(2), 195–196.
- **Klement, E.P., Mesiar, R., Pap, E.** (2000). *Triangular Norms*. Kluwer Academic Publishers, Dordrecht.
- **Ryzhov, A.P., Fedotov, F.A.** (2026). Fuzzy Risk Aggregation on Trees Using t-norms. MSU preprint (source documents in this repository).

*Note: citations from training knowledge (cutoff August 2025); verify details before citing.*
