# Aczel–Alsina T-norm

## Generator and T-norm

$$g_\lambda(x) = (-\ln x)^\lambda, \qquad \lambda \ge 1$$

- $g_\lambda(1) = 0$ ✓
- $g_\lambda(0^+) = +\infty$ ✓ (strict)
- $g_\lambda^{-1}(s) = e^{-s^{1/\lambda}}$

$$T_\lambda(x,y) = \exp\!\bigl(-[(-\ln x)^\lambda + (-\ln y)^\lambda]^{1/\lambda}\bigr)$$

## Profile $h = -g'$

$$g'(x) = \lambda(-\ln x)^{\lambda-1} \cdot \left(-\frac{1}{x}\right)$$

$$\boxed{h_a(x) = \frac{\lambda\,(-\ln x)^{\lambda-1}}{x}}$$

## Monotonicity of $h$

$$\frac{h'(x)}{h(x)} = \frac{1}{x}\left(\frac{\lambda-1}{\ln(1/x)} - 1\right)$$

Since $\ln(1/x) > 0$ for $x \in (0,1)$, $h$ is globally decreasing iff $\lambda \ge 1$.

## Special cases

| $\lambda$ | $h(x)$ | T-norm |
|-----------|---------|--------|
| $1$ | $1/x$ | Product |
| $\to\infty$ | concentrates near $x=0$ | $\min(x,y)$ |

## Selectivity

$$\Gamma_n(t,a) = \frac{h(t)}{(n-1)\,h(a)} = \frac{1}{n-1}\left(\frac{\ln(1/t)}{\ln(1/a)}\right)^{\lambda-1}\frac{a}{t}$$

---

## Variational characterization

**Problem.** Find a functional $\Delta[h]$ of the form

$$\Delta[h] = \inf_{0 \le t < a \le 1} \delta_1\!\left(\left|\frac{h(t)}{h(a)}\right|\right)$$

that is maximized over all profiles $h$ at $h = h_a$.

### Change of variables

Set $u = -\ln t$, $v = -\ln a$ ($u > v > 0$), and define the **rescaled profile**:

$$\tilde{h}(u) = e^{-u}\,h(e^{-u})$$

Then:
$$\frac{h(t)}{h(a)} = e^{u-v}\,\frac{\tilde{h}(u)}{\tilde{h}(v)}$$

For $h = h_\lambda$: $\tilde{h}(u) = \lambda u^{\lambda-1}$, so $\tilde{h}(u)/\tilde{h}(v) = (u/v)^{\lambda-1}$.

### Characterizing condition

The **log-log secant** of $\tilde{h}$ at a pair $u > v > 0$ is:

$$S(u,v) := \frac{\ln\tilde{h}(u) - \ln\tilde{h}(v)}{\ln(u/v)}$$

$h_\lambda$ is the unique decreasing profile (up to a positive scalar) for which

$$S(u,v) = \lambda - 1 \quad \text{for all } u > v > 0$$

i.e., $\ln\tilde{h}$ is affine in $\ln u$. In $(t,a)$ coordinates:

$$\ln\frac{h(t)}{h(a)} = \ln\frac{a}{t} + (\lambda-1)\ln\frac{\ln(1/t)}{\ln(1/a)}$$

From the ratio $h(t)/h(a)$ and the pair $(t,a)$ one recovers $S$ as:

$$S(u,v) = \frac{\ln(h(t)/h(a)) - \ln(a/t)}{\ln(\ln(1/t)/\ln(1/a))}$$

### Variational functional

$$\boxed{\Delta[h] = -\sup_{0 < t < a \le 1}\bigl(S(u,v) - (\lambda-1)\bigr)^2}$$

### Why $h_\lambda$ is the unique maximizer

$\Delta[h] \le 0$ for all admissible $h$, with $\Delta[h] = 0$ iff $S(u,v) = \lambda - 1$ for all pairs, i.e. iff

$$\tilde{h}(u) = e^C u^{\lambda-1} \implies h(x) = e^C\,\frac{(-\ln x)^{\lambda-1}}{x} \propto h_\lambda(x)$$

So $h_\lambda$ is the **unique maximizer** of $\Delta$, up to the positive scalar freedom in the generator.
