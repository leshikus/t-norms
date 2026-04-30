# T-Norm Quality Metrics: Signal-to-Noise Perspective

Based on: Ryzhov A.P., Fedotov F.A. — *Fuzzy Risk Aggregation on Trees Using T-norms*,
Lomonosov Readings 2026, Moscow State University, April 1, 2026.

---

## Setup

Consider a directed tree whose $n$ leaves take values $x_1,\ldots,x_n \in [0,1]$ (higher is better). Every internal node aggregates its children using the same strict Archimedean t-norm:
$$T(x,y) = g^{-1}(g(x)+g(y))$$

**Profile:** $h(x) := -g'(x) > 0$, so $g(x) = \int_x^1 h(u)\,du$.

Write $R_n(t, a)$ for the root value when one designated leaf (the *incident* node) takes value $t$ and all remaining $n-1$ leaves share the common *background* value $a$.

---

## SNR Metrics

### Incident sensitivity — signal
$$I_n(t,a) = \frac{\partial R_n}{\partial t} = \frac{h(t)}{h(R_n(t,a))}$$

### Background sensitivity — noise
$$B_n(t,a) = \frac{\partial R_n}{\partial a} = (n-1)\frac{h(a)}{h(R_n(t,a))}$$

### Selectivity — SNR analogue
$$\Gamma_n(t,a) = \frac{I_n(t,a)}{B_n(t,a)} = \frac{h(t)}{(n-1)\,h(a)}$$

### Worst-case selectivity
The adversarial regime is $t < a$ (the incident node scores lower than the background), giving
$$\Gamma_n^* = \inf_{0 \le t < a \le 1} \Gamma_n(t,a) = \frac{1}{n-1}\inf_{0 \le t < a \le 1}\frac{h(t)}{h(a)}.$$

---

## Lipschitz continuity of the t-norm

The partial derivative of $T$ with respect to its first argument is
$$\frac{\partial T}{\partial x}(x,y) = \frac{h(x)}{h(T(x,y))}.$$

Since $T(x,y) \le x$ and $h > 0$, the condition $\frac{\partial T}{\partial x} \le 1$ reduces to
$$h(x) \le h(T(x,y)),$$
which holds whenever $h$ is **non-increasing** (equivalently, $g$ is convex, $g'' \ge 0$).

Define the **Lipschitz constant** of $T$:
$$L_T = \sup_{x,y \in (0,1)} \frac{h(x)}{h(T(x,y))}.$$

**Proposition.** $T$ is $L_T$-Lipschitz in each argument, and therefore globally:
$$|T(x,y) - T(x',y')| \le L_T\bigl(|x - x'| + |y - y'|\bigr).$$

$L_T \le 1$ if and only if $h$ is non-increasing (equivalently, $g'' \ge 0$), and by induction the aggregate $R_n$ inherits the same Lipschitz constant in each leaf variable.

**Connection to SNR.** The signal sensitivity $I_n(t,a) = h(t)/h(R_n(t,a))$ is the local Lipschitz constant of $R_n$ with respect to $t$; $L_T$ is its global supremum.

---

## Reciprocity of $L_T$ and $\Gamma_n^*$

For a strict t-norm, $T(x,\cdot):(0,1)\to(0,x)$ is surjective: as $y$ ranges over $(0,1)$, $T(x,y)$ is continuous, tends to $0^+$ as $y\to 0^+$, and tends to $x^-$ as $y\to 1^-$. Hence the supremum over $y$ ranges over all $z\in(0,x)$:
$$L_T = \sup_{x,y\in(0,1)}\frac{h(x)}{h(T(x,y))} = \sup_{0<z<x<1}\frac{h(x)}{h(z)}.$$

Relabeling $(z,x)\mapsto(t,a)$ with $t<a$ gives $L_T = \sup_{0<t<a<1} h(a)/h(t)$, which is the reciprocal of $\inf_{0<t<a<1} h(t)/h(a) = (n-1)\,\Gamma_n^*$. Therefore:

$$\boxed{L_T\cdot\Gamma_n^* = \frac{1}{n-1}.}$$

**Corollary (convex generators).** When $g$ is convex ($h$ non-increasing), $L_T = 1$ and $\Gamma_n^* = \frac{1}{n-1}$, the maximum attainable value. Any deviation from convexity — an increasing region in $h$ — forces $L_T > 1$ and depresses $\Gamma_n^*$ proportionally.
