"""
qab_plot  —  two-dimensional (alpha, beta) heat-maps of all q_{alpha,beta} functionals.

Usage
-----
    python qab_plot.py                     # saves qab_functionals.png
    python qab_plot.py output.png          # saves to custom path

    from qab_plot import qab_plot
    fig, data = qab_plot(n_alpha=30, n_beta=30)
    fig.savefig("out.png", dpi=150)

The six functionals plotted
---------------------------
  J_D        : inf_{Omega+} D_n(t;a)           — worst-case detection (numerical)
  J_noise    : sup_{Omega+} B_n(t,a)           — worst-case noise sensitivity (numerical)
  J_Gamma    : inf_{Omega+} Gamma_n(t,a)       — selectivity floor (analytic)
  J_reg      : regularizer                      — analytic
  J_D_lb     : analytic lower bound on J_D      (only finite when alpha > 1)
  J_noise_ub : analytic upper bound on J_noise  (only finite when alpha > 1)
"""

import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Import qab-numerics (hyphenated filename requires importlib)
# ---------------------------------------------------------------------------
_mod_path = Path(__file__).parent / "qab-numerics.py"
_spec = importlib.util.spec_from_file_location("qab_numerics", _mod_path)
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)

DomainParams        = _mod.DomainParams
QabParams           = _mod.QabParams
compute_functionals = _mod.compute_functionals

# ---------------------------------------------------------------------------
# Functional registry
# Each entry: (dict_key, latex_label, long_title)
# ---------------------------------------------------------------------------
_FUNCTIONALS = [
    ("J_D",        r"$J_D$",
                   "inf Detection"),
    ("J_noise",    r"$J_\mathrm{noise}$",
                   "sup Noise sensitivity"),
    ("J_Gamma",    r"$J_\Gamma$",
                   "inf Selectivity"),
    ("J_reg",      r"$J_\mathrm{reg}$",
                   "Regularizer"),
    ("J_D_lb",     r"$J_D^{\,\mathrm{lb}}$",
                   "Detection lower bound"),
    ("J_noise_ub", r"$J_\mathrm{noise}^{\,\mathrm{ub}}$",
                   "Noise upper bound"),
]


# ---------------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------------

def qab_plot(
    alpha_range=(1.1, 3.0),
    beta_range=(-1.0, 2.0),
    domain=None,
    n_alpha=20,
    n_beta=20,
    n_a=15,
    n_t=15,
    cmap="plasma",
    n_contour_lines=8,
    figsize=None,
):
    """
    Compute and display 2-D heat-maps of all q_{alpha,beta} functionals over
    the (alpha, beta) parameter space.

    Parameters
    ----------
    alpha_range     : (lo, hi) — alpha axis range; must be > 1 for bound functionals
    beta_range      : (lo, hi) — beta axis range
    domain          : DomainParams or None
                      default: t_min=0.05, a_min=0.3, a_max=0.85, eta=0.5, n=5
    n_alpha, n_beta : int — parameter grid resolution
    n_a, n_t        : int — inner (a, t) grid density inside compute_functionals
    cmap            : matplotlib colormap name
    n_contour_lines : number of white iso-lines overlaid on each panel
    figsize         : (width, height) in inches, or None for auto

    Returns
    -------
    fig  : matplotlib.figure.Figure
    data : dict[str, np.ndarray]  — shape (n_alpha, n_beta), keyed by functional name
           Row index → alpha,  column index → beta
    """
    if domain is None:
        domain = DomainParams(t_min=0.05, a_min=0.3, a_max=0.85, eta=0.5, n=5)

    alphas = np.linspace(*alpha_range, n_alpha)
    betas  = np.linspace(*beta_range,  n_beta)

    keys = [k for k, *_ in _FUNCTIONALS]
    data = {k: np.full((n_alpha, n_beta), np.nan) for k in keys}

    total = n_alpha * n_beta
    done  = 0
    step  = max(1, total // 20)
    print(f"Scanning {n_alpha}×{n_beta} (alpha×beta) grid — {total} points ...", flush=True)

    for i, alpha in enumerate(alphas):
        for j, beta in enumerate(betas):
            p = QabParams(alpha, beta, domain)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                try:
                    res = compute_functionals(p, n_a=n_a, n_t=n_t)
                    for k in keys:
                        v = res.get(k, np.nan)
                        data[k][i, j] = v if np.isfinite(v) else np.nan
                except Exception:
                    pass  # leave as nan on numerical failure
            done += 1
            if done % step == 0:
                pct = 100 * done // total
                print(f"  {done}/{total}  ({pct}%)", flush=True)

    print("Scan done.  Rendering figure ...", flush=True)

    # ---- layout -----------------------------------------------------------
    ncols = 3
    nrows = (len(_FUNCTIONALS) + ncols - 1) // ncols
    if figsize is None:
        figsize = (ncols * 4.8, nrows * 4.0)

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)

    # meshgrid: BA and AA both have shape (n_alpha, n_beta)
    BA, AA = np.meshgrid(betas, alphas)

    for idx, (key, label, title) in enumerate(_FUNCTIONALS):
        row, col = divmod(idx, ncols)
        ax = axes[row][col]
        Z  = data[key]

        finite_mask = np.isfinite(Z)
        if finite_mask.any():
            vmin = float(np.nanmin(Z))
            vmax = float(np.nanmax(Z))
            if vmin == vmax:          # constant plane: widen for colorbar legibility
                vmin -= 0.5 * abs(vmin) + 1e-9
                vmax += 0.5 * abs(vmax) + 1e-9

            pc = ax.pcolormesh(BA, AA, Z,
                               cmap=cmap, vmin=vmin, vmax=vmax,
                               shading="auto")

            if n_contour_lines > 0:
                levels = np.linspace(vmin, vmax, n_contour_lines + 2)[1:-1]
                try:
                    cs = ax.contour(BA, AA, Z,
                                    levels=levels,
                                    colors="white", linewidths=0.5, alpha=0.55)
                    ax.clabel(cs, fmt="%.2g", fontsize=6, inline=True)
                except Exception:
                    pass

            cb = fig.colorbar(pc, ax=ax, fraction=0.046, pad=0.04)
            cb.ax.tick_params(labelsize=7)
        else:
            ax.set_facecolor("#222222")
            ax.text(0.5, 0.5, "N/A\n(all NaN)", ha="center", va="center",
                    transform=ax.transAxes, color="white", fontsize=10)

        ax.set_xlabel(r"$\beta$", fontsize=11)
        ax.set_ylabel(r"$\alpha$", fontsize=11)
        ax.set_title(f"{label}  —  {title}", fontsize=10)
        ax.tick_params(labelsize=8)

    # hide unused panels
    for idx in range(len(_FUNCTIONALS), nrows * ncols):
        row, col = divmod(idx, ncols)
        axes[row][col].set_visible(False)

    domain_str = (
        rf"$n={domain.n}$, "
        rf"$t_{{\min}}={domain.t_min}$, "
        rf"$a\in[{domain.a_min},\,{domain.a_max}]$, "
        rf"$\eta={domain.eta}$"
    )
    fig.suptitle(
        rf"$q_{{\alpha,\beta}}$ functionals over $(\alpha,\beta)$  —  {domain_str}",
        fontsize=12,
        y=1.01,
    )
    fig.tight_layout()
    return fig, data


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "qab_functionals.png"
    fig, _ = qab_plot(n_alpha=20, n_beta=20, n_a=12, n_t=12)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved → {out}")
    plt.show()
