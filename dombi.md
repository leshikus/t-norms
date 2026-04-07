# Convexity of the Dombi Generator

## Setup

The **Dombi t-norm** with parameter $\lambda > 0$ has additive generator

$$g_\lambda(x) = \left(\frac{1-x}{x}\right)^\lambda, \qquad x \in (0,1].$$

It satisfies $g_\lambda(1)=0$ and $g_\lambda(0^+)=+\infty$, confirming it is a valid strict Archimedean generator.

---

## First derivative (profile)

$$g_\lambda'(x) = -\lambda\,\frac{(1-x)^{\lambda-1}}{x^{\lambda+1}}$$

The profile $h_\lambda(x) := -g_\lambda'(x) = \lambda\,(1-x)^{\lambda-1} x^{-(\lambda+1)} > 0$ for all $x\in(0,1)$. $\checkmark$

---

## Second derivative

Differentiating $g'$ by the product rule:

$$g_\lambda''(x) = \lambda\,(1-x)^{\lambda-2}\,x^{-(\lambda+2)}\,\bigl(\lambda+1-2x\bigr).$$

**Derivation.** Write $g'(x) = -\lambda\,u(x)\,v(x)$ with $u=(1-x)^{\lambda-1}$, $v=x^{-(\lambda+1)}$.

$$u' = -(\lambda-1)(1-x)^{\lambda-2}, \qquad v' = -(\lambda+1)x^{-(\lambda+2)}.$$

$$g'' = -\lambda(u'v + uv') = -\lambda\!\left[-(\lambda-1)(1-x)^{\lambda-2}x^{-(\lambda+1)} - (\lambda+1)(1-x)^{\lambda-1}x^{-(\lambda+2)}\right]$$

$$= \lambda(1-x)^{\lambda-2}x^{-(\lambda+2)}\!\left[(\lambda-1)x + (\lambda+1)(1-x)\right]$$

$$= \lambda(1-x)^{\lambda-2}x^{-(\lambda+2)}\,(\lambda+1-2x). \qquad \square$$

---

## Sign analysis

For $x \in (0,1)$ the prefactor $\lambda(1-x)^{\lambda-2}x^{-(\lambda+2)}$ is **strictly positive** (regardless of $\lambda>0$), so

$$\operatorname{sgn} g_\lambda''(x) = \operatorname{sgn}(\lambda+1-2x).$$

The inflection point is at

$$x^* = \frac{\lambda+1}{2}.$$

| Region | $g_\lambda''$ | Shape of $g_\lambda$ |
|--------|--------------|----------------------|
| $x < x^*$ | $> 0$ | convex |
| $x = x^*$ | $= 0$ | inflection |
| $x > x^*$ | $< 0$ | concave |

---

## Convexity criterion

$g_\lambda$ is **convex on the entire interval $(0,1)$** if and only if $x^* \ge 1$:

$$\frac{\lambda+1}{2} \ge 1 \iff \boxed{\lambda \ge 1}.$$

For $0 < \lambda < 1$ the inflection point $x^* = (\lambda+1)/2 \in (\tfrac{1}{2}, 1)$ lies inside $(0,1)$, so $g_\lambda$ is neither convex nor concave globally.

---

## Connection to the profile and global robustness

Because $h_\lambda = -g_\lambda'$,

$$h_\lambda'(x) = -g_\lambda''(x) \propto -({\lambda+1-2x}).$$

| $\lambda$ | Profile $h_\lambda$ | $g_\lambda$ | Global robustness (1-Lipschitz in $\ell_1$) |
|-----------|--------------------|-----------|--------------------------------------------|
| $\lambda \ge 1$ | non-increasing on $(0,1)$ | convex | **yes** |
| $0 < \lambda < 1$ | has a local minimum at $x^*$; increases for $x > x^*$ | non-convex | **no** |

The 1-Lipschitz bound $|\partial R/\partial x_k| \le 1$ requires $h$ non-increasing (see summary §6), which is equivalent to $g$ being convex. So:

> **The Dombi t-norm provides a globally robust aggregator if and only if $\lambda \ge 1$.**

---

## Selectivity

The selectivity $\Gamma_n(t,a) = h_\lambda(t)/[(n-1)h_\lambda(a)]$ in the incident scenario $(t,a,\ldots,a)$:

$$\Gamma_n(t,a) = \frac{1}{n-1}\left(\frac{1-t}{1-a}\right)^{\!\lambda-1}\!\!\left(\frac{a}{t}\right)^{\!\lambda+1}.$$

Note $h(R)$ cancels, so $\Gamma_n$ is independent of the root value and of generator scale.

**Monotonicity in $\lambda$:** Taking the log,

$$\ln\Gamma_n = -\ln(n-1) + (\lambda-1)\ln\frac{1-t}{1-a} + (\lambda+1)\ln\frac{a}{t}$$

$$= -\ln(n-1) + \lambda\ln\frac{a(1-t)}{t(1-a)} + \ln\frac{a(1-a)}{t(1-t)}.$$

Since $t < a$ implies $\frac{a(1-t)}{t(1-a)} > 1$, the coefficient of $\lambda$ is positive, so **selectivity is strictly increasing in $\lambda$**. Higher $\lambda$ always improves signal-to-noise ratio, at fixed $(t,a)$.

---

## Detection

The inverse generator is $g_\lambda^{-1}(s) = \bigl(1+s^{1/\lambda}\bigr)^{-1}$. Define

$$A = \left(\frac{1-a}{a}\right)^\lambda, \qquad T = \left(\frac{1-t}{t}\right)^\lambda.$$

The homogeneous root and the incident root are

$$R_n(a) = g_\lambda^{-1}(nA) = \frac{a}{a + n^{1/\lambda}(1-a)},$$

$$R(t,a) = g_\lambda^{-1}\!\bigl(T+(n-1)A\bigr) = \frac{1}{1+\bigl(T+(n-1)A\bigr)^{1/\lambda}}.$$

Detection is the finite drop of the root score:

$$D_n(t;a) = R_n(a) - R(t,a) = \frac{a}{a+n^{1/\lambda}(1-a)} - \frac{1}{1+\bigl(T+(n-1)A\bigr)^{1/\lambda}}.$$

Since $t < a$ implies $T > A$, we have $T+(n-1)A > nA$, hence $R(t,a) < R_n(a)$ and $D_n > 0$ as expected.

**Closed form at uniform background $a$.** Write $\delta = n^{1/\lambda}(1-a)/a$ so $R_n(a)=1/(1+\delta)$. There is no further simplification unless $T$ and $A$ share a rational exponent, so the formula above is the natural closed form.

---

## Regularity

The log-derivative of the profile is

$$(\log h_\lambda)'(x) = \frac{d}{dx}\bigl[\log\lambda + (\lambda-1)\log(1-x) - (\lambda+1)\log x\bigr] = -\frac{\lambda+1}{x} - \frac{\lambda-1}{1-x} = -\frac{\lambda+1-2x}{x(1-x)}.$$

(Partial-fraction check: $\frac{\lambda+1-2x}{x(1-x)} = \frac{\lambda+1}{x} + \frac{\lambda-1}{1-x}$.)

**Antiderivative** (without absolute value): $F(x) = (\lambda+1)\ln x - (\lambda-1)\ln(1-x)$.

### Case $\lambda \ge 1$ (profile non-increasing, no sign change)

$|(\log h)'| = \frac{\lambda+1-2x}{x(1-x)}$ on all of $(0,1)$, so over a working zone $[t_{\min}, a_+]$ with $a_+<1$:

$$\boxed{J_{\mathrm{reg}} = (\lambda+1)\ln\frac{a_+}{t_{\min}} + (\lambda-1)\ln\frac{1-t_{\min}}{1-a_+}.}$$

- **$\lambda = 1$:** $J_{\mathrm{reg}} = 2\ln(a_+/t_{\min})$; taking $a_+\to 1^-$ gives $J_{\mathrm{reg}} = 2\ln(1/t_{\min}) < \infty$.
- **$\lambda > 1$:** the term $(\lambda-1)\ln\frac{1}{1-a_+}\to+\infty$ as $a_+\to 1^-$, so **$J_{\mathrm{reg}}=+\infty$ when the zone extends to 1**. For any fixed $a_+<1$ the integral is finite.

### Case $0 < \lambda < 1$ (inflection at $x^* = (\lambda+1)/2 < 1$)

$h_\lambda$ decreases on $(0, x^*)$ and increases on $(x^*, 1)$, so the absolute value must be split. For a zone $[t_{\min}, a_+]$ with $t_{\min} < x^* < a_+$:

$$J_{\mathrm{reg}} = 2F(x^*) - F(t_{\min}) - F(a_+), \qquad x^* = \frac{\lambda+1}{2},$$

$$F(x^*) = (\lambda+1)\ln\frac{\lambda+1}{2} - (\lambda-1)\ln\frac{1-\lambda}{2}.$$

If $a_+ \le x^*$ (zone stays in the convex part) the same formula as the $\lambda\ge 1$ case applies without splitting.

### Divergence summary

| $\lambda$ | $J_{\mathrm{reg}}$ as $a_+\to 1^-$ |
|-----------|--------------------------------------|
| $\lambda = 1$ | $2\ln(1/t_{\min})$ — **finite** |
| $\lambda > 1$ | $+\infty$ (logarithmic blow-up from the $\frac{\lambda-1}{1-x}$ term) |
| $0 < \lambda < 1$ | $+\infty$ (same cause near $x=1$) |

The Dombi generator has finite total regularity **only at $\lambda=1$**, coinciding exactly with the convexity boundary. For $\lambda\ne 1$, $J_{\mathrm{reg}}$ is finite on any compact subinterval $[t_{\min}, a_+]\subset(0,1)$ but diverges as the zone approaches 1.

---

## Summary

| Condition | Result |
|-----------|--------|
| $\lambda \ge 1$ | $g_\lambda$ convex, $h_\lambda$ non-increasing, robustness guaranteed |
| $\lambda = 1$ | reduces to the **product t-norm** ($g(x)=(1-x)/x$, which equals $1/x - 1$; agrees with $-\ln x$ up to reparametrization only at $x=1$) — actually a distinct t-norm, but shares the convexity boundary |
| $0 < \lambda < 1$ | $g_\lambda$ non-convex, inflection at $x^*=(\lambda+1)/2$, robustness fails |

The inflection point formula $x^* = (\lambda+1)/2$ gives a precise diagnostic: for any desired working interval $[t_{\min}, 1]$, convexity holds throughout that interval whenever $\lambda \ge 2t_{\min}-1$ (a weaker, domain-restricted condition that is easier to satisfy than $\lambda \ge 1$).
