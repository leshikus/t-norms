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
