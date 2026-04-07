"""
Numerical algorithm for the q_{alpha,beta} generator family.

The shape function (log-derivative of profile) is:
    q_{alpha,beta}(x) = alpha + beta * exp(-( (-ln x - u0) / tau )^2)

The profile is h_{alpha,beta}(x) = x^{-alpha} * exp(beta * A(u0, tau, -ln x)),
where A(u0, tau, u) = integral_0^u exp(-((v - u0)/tau)^2) dv  (scaled erf integral).

Generator (no closed form): g(x) = integral_x^1 h(u) du.
Root: R_n(t, a) = g^{-1}(g(t) + (n-1)*g(a)).

Selectivity and regularizer have exact analytic formulas (from paper).
Detection and noise require numerical g and g^{-1}.
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import erf
import warnings

# ---------------------------------------------------------------------------
# Working domain defaults (can be overridden via DomainParams)
# ---------------------------------------------------------------------------

class DomainParams:
    """Working domain Omega^+ = {(t, a): a in [a_min, a_max], t in [t_min, eta*a]}."""
    def __init__(self, t_min=0.05, a_min=0.3, a_max=0.9, eta=0.5, n=5):
        self.t_min = t_min
        self.a_min = a_min
        self.a_max = a_max
        self.eta   = eta        # t <= eta * a
        self.n     = n          # number of leaves


class QabParams:
    """Parameters of the q_{alpha,beta} family and domain geometry."""
    def __init__(self, alpha, beta, domain: DomainParams):
        self.alpha  = alpha
        self.beta   = beta
        self.domain = domain
        # u0 and tau pinned to domain so that the Gaussian bump sits over the
        # working zone in log-space.  A natural choice: centre at the midpoint
        # of [-ln a_max, -ln t_min] with half-width equal to that half-range.
        u_lo = -np.log(domain.a_max)   # ~ 0  (a close to 1)
        u_hi = -np.log(domain.t_min)   # > 0  (t_min small)
        self.u0  = 0.5 * (u_lo + u_hi)
        self.tau = 0.5 * (u_hi - u_lo)


# ---------------------------------------------------------------------------
# Core analytic building blocks
# ---------------------------------------------------------------------------

def A_func(u0: float, tau: float, u: float) -> float:
    """
    A(u) = integral_0^u exp(-((v - u0)/tau)^2) dv
         = (tau * sqrt(pi) / 2) * (erf((u - u0)/tau) + erf(u0/tau))
    """
    return 0.5 * tau * np.sqrt(np.pi) * (erf((u - u0) / tau) + erf(u0 / tau))


def A_inf(u0: float, tau: float) -> float:
    """A(+inf) = (tau*sqrt(pi)/2)*(1 + erf(u0/tau))."""
    return 0.5 * tau * np.sqrt(np.pi) * (1.0 + erf(u0 / tau))


def profile(x: float, p: QabParams) -> float:
    """h_{alpha,beta}(x) = x^{-alpha} * exp(beta * A(-ln x))."""
    u = -np.log(x)
    return x ** (-p.alpha) * np.exp(p.beta * A_func(p.u0, p.tau, u))


# ---------------------------------------------------------------------------
# Generator  g(x) = integral_x^1 h(u) du   (numerical, log-space)
#
# Change of variables u = -ln(x):  h(e^{-v}) |dx/dv| = exp((alpha-1)*v + beta*A(v))
# So g(x) = integral_0^{-ln x}  phi(v) dv,   phi(v) = exp((alpha-1)*v + beta*A(v))
# This avoids the x^{-alpha} singularity at x=0 entirely.
# ---------------------------------------------------------------------------

def _phi(v: float, p: QabParams) -> float:
    """Integrand in log-space: phi(v) = exp((alpha-1)*v + beta*A(u0,tau,v))."""
    return np.exp((p.alpha - 1) * v + p.beta * A_func(p.u0, p.tau, v))


def generator(x: float, p: QabParams) -> float:
    """g(x) = integral_0^{-ln x} phi(v) dv.  Returns +inf for x -> 0."""
    if x <= 0.0:
        return np.inf
    if x >= 1.0:
        return 0.0
    U = -np.log(x)   # upper limit of integration
    val, _ = quad(_phi, 0.0, U, args=(p,), limit=200, epsabs=1e-10, epsrel=1e-10)
    return val


def generator_vec(xs: np.ndarray, p: QabParams) -> np.ndarray:
    return np.array([generator(x, p) for x in xs])


# ---------------------------------------------------------------------------
# Inverse generator  g^{-1}(s)  via Brent root-finding
# ---------------------------------------------------------------------------

def inv_generator(s: float, p: QabParams,
                  x_lo: float = 1e-9, x_hi: float = 1.0 - 1e-10) -> float:
    """Solve g(x) = s for x in (0, 1)."""
    if s <= 0.0:
        return 1.0
    if s >= generator(x_lo, p):
        return x_lo  # saturate at boundary
    f = lambda x: generator(x, p) - s
    return brentq(f, x_lo, x_hi, xtol=1e-12, rtol=1e-12)


# ---------------------------------------------------------------------------
# Root value R_n(t, a) and derived metrics
# ---------------------------------------------------------------------------

def root(t: float, a: float, p: QabParams) -> float:
    """R_n(t, a) = g^{-1}(g(t) + (n-1)*g(a))."""
    n = p.domain.n
    s = generator(t, p) + (n - 1) * generator(a, p)
    return inv_generator(s, p)


def root_homogeneous(a: float, p: QabParams) -> float:
    """R_n(a) = g^{-1}(n * g(a))."""
    n = p.domain.n
    s = n * generator(a, p)
    return inv_generator(s, p)


def detection(t: float, a: float, p: QabParams) -> float:
    """D_n(t; a) = R_n(a) - R_n(t, a)."""
    return root_homogeneous(a, p) - root(t, a, p)


def noise_sensitivity(t: float, a: float, p: QabParams) -> float:
    """B_n(t, a) = (n-1) * h(a) / h(R_n(t, a))."""
    n = p.domain.n
    R = root(t, a, p)
    return (n - 1) * profile(a, p) / profile(R, p)


def selectivity(t: float, a: float, p: QabParams) -> float:
    """Gamma_n(t, a) = h(t) / ((n-1) * h(a))  — exact analytic formula."""
    n   = p.domain.n
    u_t = -np.log(t)
    u_a = -np.log(a)
    ratio = (a / t) ** p.alpha * np.exp(
        p.beta * (A_func(p.u0, p.tau, u_t) - A_func(p.u0, p.tau, u_a))
    )
    return ratio / (n - 1)


# ---------------------------------------------------------------------------
# Analytic exact formulas for selectivity floor and regularizer
# ---------------------------------------------------------------------------

def J_selectivity(p: QabParams) -> float:
    """
    J_Gamma = inf_{Omega^+} Gamma_n(t, a).
    Since Gamma depends on (t/a) and on A(-ln t) - A(-ln a),
    and Gamma is decreasing in a/t (fewer incident-to-background distance),
    the infimum is at t = eta*a, a = a_max  (worst case: t as large, a as large).
    """
    d  = p.domain
    t  = d.eta * d.a_max
    a  = d.a_max
    return selectivity(t, a, p)


def J_regularity(p: QabParams) -> float:
    """
    J_reg = alpha * ln(1/t_min) + beta * A(u0, tau, -ln t_min).
    Exact analytic formula from paper.
    """
    d = p.domain
    return p.alpha * np.log(1.0 / d.t_min) + p.beta * A_func(p.u0, p.tau, -np.log(d.t_min))


# ---------------------------------------------------------------------------
# Sandwich (analytic) bounds  (slide 27-28 of presentation)
# ---------------------------------------------------------------------------

def K_constant(p: QabParams) -> float:
    """K = exp(beta * A_inf)."""
    return np.exp(p.beta * A_inf(p.u0, p.tau))


def _S_alpha(t: float, a: float, alpha: float, n: int) -> float:
    """S_alpha(t, a) = (t^{1-alpha} - 1) + (n-1)*(a^{1-alpha} - 1)."""
    return (t ** (1 - alpha) - 1) + (n - 1) * (a ** (1 - alpha) - 1)


def root_lower_bound(t: float, a: float, p: QabParams) -> float:
    """Lower bound on R_n(t, a) from sandwich estimate."""
    K = K_constant(p)
    S = _S_alpha(t, a, p.alpha, p.domain.n)
    return (1 + K * S) ** (-1.0 / (p.alpha - 1))


def root_upper_bound(t: float, a: float, p: QabParams) -> float:
    """Upper bound on R_n(t, a) from sandwich estimate."""
    K = K_constant(p)
    S = _S_alpha(t, a, p.alpha, p.domain.n)
    return (1 + S / K) ** (-1.0 / (p.alpha - 1))


def detection_lower_bound(t: float, a: float, p: QabParams) -> float:
    """Analytic lower bound on D_n(t; a)."""
    K = K_constant(p)
    alpha = p.alpha
    R_lo = root_lower_bound(t, a, p)
    S_diff = (t ** (1 - alpha) - a ** (1 - alpha))  # positive since t < a
    return S_diff / ((alpha - 1) * K) * R_lo ** alpha


def noise_upper_bound(t: float, a: float, p: QabParams) -> float:
    """Analytic upper bound on B_n(t, a)."""
    K = K_constant(p)
    n = p.domain.n
    R_hi = root_upper_bound(t, a, p)
    return (n - 1) * K * (R_hi / a) ** p.alpha


# ---------------------------------------------------------------------------
# Worst-case functionals over Omega^+  (grid search)
# ---------------------------------------------------------------------------

def compute_functionals(p: QabParams,
                        n_a: int = 30,
                        n_t: int = 30,
                        verbose: bool = False):
    """
    Compute worst-case detection J_D = inf D_n  and  noise J_noise = sup B_n
    over the working domain Omega^+ by grid search.

    Also returns exact analytic J_Gamma, J_reg, and analytic bound statistics.

    Parameters
    ----------
    p      : QabParams
    n_a    : grid points for a
    n_t    : grid points for t  (for each a, t in [t_min, eta*a])
    verbose: print progress

    Returns
    -------
    dict with keys: J_D, J_noise, J_Gamma, J_reg,
                    J_D_lb (analytic lower bound on J_D),
                    J_noise_ub (analytic upper bound on J_noise).
    """
    d = p.domain
    a_grid = np.linspace(d.a_min, d.a_max, n_a)

    min_D   = np.inf
    max_B   = -np.inf
    min_D_lb = np.inf
    max_B_ub = -np.inf

    for a in a_grid:
        t_hi = d.eta * a
        t_lo = d.t_min
        if t_lo >= t_hi:
            continue
        t_grid = np.linspace(t_lo, t_hi, n_t)

        for t in t_grid:
            if t >= a:
                continue

            D = detection(t, a, p)
            B = noise_sensitivity(t, a, p)

            if D < min_D:
                min_D = D
            if B > max_B:
                max_B = B

            # analytic bounds (require alpha > 1)
            if p.alpha > 1:
                D_lb = detection_lower_bound(t, a, p)
                B_ub = noise_upper_bound(t, a, p)
                if D_lb < min_D_lb:
                    min_D_lb = D_lb
                if B_ub > max_B_ub:
                    max_B_ub = B_ub

            if verbose:
                print(f"  a={a:.3f}  t={t:.3f}  D={D:.4f}  B={B:.4f}")

    return {
        "J_D"        : min_D,          # inf detection  (exact numerical)
        "J_noise"    : max_B,          # sup noise      (exact numerical)
        "J_Gamma"    : J_selectivity(p),   # inf selectivity (analytic)
        "J_reg"      : J_regularity(p),    # regularizer     (analytic)
        "J_D_lb"     : min_D_lb,       # analytic lower bound on J_D
        "J_noise_ub" : max_B_ub,       # analytic upper bound on J_noise
    }


# ---------------------------------------------------------------------------
# Convenience: scan over alpha or beta
# ---------------------------------------------------------------------------

def scan_alpha(alphas, beta, domain, **kw):
    """Return a list of result dicts for varying alpha at fixed beta."""
    results = []
    for alpha in alphas:
        p = QabParams(alpha, beta, domain)
        res = compute_functionals(p, **kw)
        res["alpha"] = alpha
        res["beta"]  = beta
        results.append(res)
    return results


def scan_beta(alpha, betas, domain, **kw):
    """Return a list of result dicts for varying beta at fixed alpha."""
    results = []
    for beta in betas:
        p = QabParams(alpha, beta, domain)
        res = compute_functionals(p, **kw)
        res["alpha"] = alpha
        res["beta"]  = beta
        results.append(res)
    return results


# ---------------------------------------------------------------------------
# Demo / smoke test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    domain = DomainParams(t_min=0.05, a_min=0.3, a_max=0.85, eta=0.5, n=5)
    p = QabParams(alpha=1.5, beta=0.5, domain=domain)

    print("=== q_{alpha,beta} numerics ===")
    print(f"alpha={p.alpha}, beta={p.beta}, n={domain.n}")
    print(f"u0={p.u0:.4f}, tau={p.tau:.4f}")
    print(f"K = {K_constant(p):.6f}")
    print()

    # Spot-check a single point
    t0, a0 = 0.1, 0.6
    print(f"Spot check at t={t0}, a={a0}:")
    print(f"  g(t)            = {generator(t0, p):.6f}")
    print(f"  g(a)            = {generator(a0, p):.6f}")
    print(f"  R_n(a) [homog]  = {root_homogeneous(a0, p):.6f}")
    print(f"  R_n(t,a)        = {root(t0, a0, p):.6f}")
    print(f"  D_n(t;a)        = {detection(t0, a0, p):.6f}")
    print(f"  B_n(t,a)        = {noise_sensitivity(t0, a0, p):.6f}")
    print(f"  Gamma_n(t,a)    = {selectivity(t0, a0, p):.6f}")
    print()

    # Analytic bounds check
    print(f"Analytic bounds at t={t0}, a={a0}:")
    print(f"  R lower bound   = {root_lower_bound(t0, a0, p):.6f}")
    print(f"  R_n(t,a) exact  = {root(t0, a0, p):.6f}")
    print(f"  R upper bound   = {root_upper_bound(t0, a0, p):.6f}")
    print(f"  D lower bound   = {detection_lower_bound(t0, a0, p):.6f}")
    print(f"  D exact         = {detection(t0, a0, p):.6f}")
    print(f"  B exact         = {noise_sensitivity(t0, a0, p):.6f}")
    print(f"  B upper bound   = {noise_upper_bound(t0, a0, p):.6f}")
    print()

    # Worst-case functionals over Omega^+
    print("Computing worst-case functionals over Omega+ (grid 20x20)...")
    res = compute_functionals(p, n_a=20, n_t=20)
    print(f"  J_D   (inf detection, exact)        = {res['J_D']:.6f}")
    print(f"  J_D   (analytic lower bound)        = {res['J_D_lb']:.6f}")
    print(f"  J_noise (sup noise, exact)          = {res['J_noise']:.6f}")
    print(f"  J_noise (analytic upper bound)      = {res['J_noise_ub']:.6f}")
    print(f"  J_Gamma (selectivity floor, exact)  = {res['J_Gamma']:.6f}")
    print(f"  J_reg   (regularizer, analytic)     = {res['J_reg']:.6f}")
    print()

    # Quick scan over alpha
    print("Scanning alpha in [1.1, 3.0] at beta=0.5:")
    print(f"{'alpha':>6}  {'J_D':>10}  {'J_noise':>10}  {'J_Gamma':>10}  {'J_reg':>8}")
    for r in scan_alpha(np.linspace(1.1, 3.0, 6), beta=0.5, domain=domain, n_a=15, n_t=15):
        print(f"{r['alpha']:6.2f}  {r['J_D']:10.5f}  {r['J_noise']:10.5f}  "
              f"{r['J_Gamma']:10.5f}  {r['J_reg']:8.4f}")
