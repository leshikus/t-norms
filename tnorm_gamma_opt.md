# Optimal Profile for Worst-Case Selectivity

## Problem

Find $h: [0,1] \to (0,\infty)$ with $h(0) = 1$ that maximizes:
$$\min_{0 \leq t < a \leq 1} \Gamma_n(t,a) = \frac{h(t)}{(n-1)\,h(a)}$$

## Upper Bound

For any $h$, fix $a$ and let $t \to a^-$. By continuity, $h(t)/h(a) \to 1$, so:
$$\inf_{t<a} \frac{h(t)}{h(a)} \leq 1 \quad \Longrightarrow \quad \min\Gamma_n \leq \frac{1}{n-1}$$

## Achieving the Bound

If $h$ is **non-increasing**, then $t < a \Rightarrow h(t) \geq h(a)$, so $h(t)/h(a) \geq 1$ everywhere. The infimum equals 1 and the bound is tight.

If $h$ increases anywhere — $\exists\, t_0 < a_0$ with $h(t_0) < h(a_0)$ — then $\Gamma_n(t_0, a_0) < 1/(n-1)$, which is suboptimal.

## Canonical Solution

The unique $h$ that **achieves** the minimum (not just the infimum) at $\min \Gamma = 1/(n-1)$:

$$\boxed{h(x) = 1}$$

This gives:
$$g(x) = \int_x^1 1\,du = 1-x, \qquad T(x,y) = \max(0,\,x+y-1)$$

the **Łukasiewicz t-norm** — the only profile with uniform selectivity across all $(t,a)$ configurations.

## Summary

| $h$ | $\inf_{t<a} \Gamma_n$ | Minimum achieved? |
|---|---|---|
| Constant $= 1$ | $\tfrac{1}{n-1}$ | Yes |
| Strictly decreasing, $h(0)=1$ | $\tfrac{1}{n-1}$ | No (infimum only) |
| Has any increasing part | $< \tfrac{1}{n-1}$ | — |

Any strictly decreasing $h$ achieves the same infimum $1/(n-1)$ but never attains it; $h \equiv 1$ is the unique maximizer of the actual minimum.
