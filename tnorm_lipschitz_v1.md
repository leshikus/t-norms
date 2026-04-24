# T-Norm Quality Metrics: Signal-to-Noise Perspective

Based on: Ryzhov A.P., Fedotov F.A. — *Fuzzy Risk Aggregation on Trees Using T-norms*,
Lomonosov Readings 2026, Moscow State University, April 1, 2026.

---

## Setup

A directed tree with leaves $x_1,\ldots,x_n \in [0,1]$ (higher = better).
All internal nodes use the same strict Archimedean t-norm:
$$T(x,y) = g^{-1}(g(x)+g(y))$$

**Profile:** $h(x) := -g'(x) > 0$, so $g(x) = \int_x^1 h(u)\,du$.

---

## SNR Metrics

### Incident sensitivity — "signal"
$$I_n(t,a) = \frac{\partial R_n}{\partial t} = \frac{h(t)}{h(R_n(t,a))}$$

### Background sensitivity — "noise"
$$B_n(t,a) = \frac{\partial R_n}{\partial a} = (n-1)\frac{h(a)}{h(R_n(t,a))}$$

### Selectivity — **SNR analogue**
$$\Gamma_n(t,a) = \frac{I_n(t,a)}{B_n(t,a)} = \frac{h(t)}{(n-1)\,h(a)}$$

### Worst-case selectivity
The adversarial regime is $t < a$ (the incident node scores lower than the background), so
$$\Gamma_n^* = \inf_{0 \le t < a \le 1} \Gamma_n(t,a) = \frac{1}{n-1}\inf_{0 \le t < a \le 1}\frac{h(t)}{h(a)}.$$

---

## Lipschitz continuity of the t-norm

The partial derivative of $T$ with respect to its first argument:
$$\frac{\partial T}{\partial x}(x,y) = \frac{h(x)}{h(T(x,y))}$$

Since $T(x,y) \le x$ and $h > 0$, the condition $\frac{\partial T}{\partial x} \le 1$ reduces to
$$h(x) \le h(T(x,y)),$$
which holds whenever $h$ is **non-increasing** (equivalently, $g$ is convex, $g'' \ge 0$).

Define the **Lipschitz constant** of $T$:
$$L_T = \sup_{x,y \in (0,1)} \frac{h(x)}{h(T(x,y))}.$$

**Proposition.** $T$ is $L_T$-Lipschitz in each argument, hence globally:
$$|T(x,y) - T(x',y')| \le L_T\bigl(|x - x'| + |y - y'|\bigr).$$

$L_T \le 1$ iff $h$ is non-increasing (equivalently, $g$ is convex, $g'' \ge 0$).

By induction the aggregate $R_n$ inherits the same constant in each leaf variable.

**Connection to SNR.** The signal sensitivity $I_n(t,a) = h(t)/h(R_n)$ is exactly the local Lipschitz constant of $R_n$ with respect to $t$; $L_T$ is its global supremum.

---

## Reciprocity of $L_T$ and $\Gamma_n^*$

For a **strict** t-norm, $T(x,\cdot):(0,1)\to(0,x)$ is surjective (continuous, maps $y\to 0^+$ to $T\to 0^+$ and $y\to 1^-$ to $T\to x^-$). Therefore the sup over $y$ equals the sup over all $z\in(0,x)$:
$$L_T = \sup_{x,y\in(0,1)}\frac{h(x)}{h(T(x,y))} = \sup_{0<z<x<1}\frac{h(x)}{h(z)}.$$

The right-hand side is the reciprocal of $\inf_{0<t<a<1} h(t)/h(a) = (n-1)\,\Gamma_n^*$, so:

$$\boxed{L_T\cdot\Gamma_n^* = \frac{1}{n-1}.}$$

**Corollary (convex t-norms).** When $g$ is convex ($h$ non-increasing), $L_T = 1$ and the bound is saturated: $\Gamma_n^* = \frac{1}{n-1}$, the maximum attainable value. Any deviation from convexity (increasing region in $h$) drives $L_T > 1$ and $\Gamma_n^* < \frac{1}{n-1}$ proportionally.
