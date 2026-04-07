#!/usr/bin/env python3
"""Interactive t-norm explorer: choose a family, tune parameters, see g, h=-g', T(x,y)."""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons
import warnings
warnings.filterwarnings("ignore")

x = np.linspace(1e-6, 1, 400)
X, Y = np.meshgrid(np.linspace(1e-6, 1, 80), np.linspace(1e-6, 1, 80))

# --- generators and t-norms ---

FAMILIES = {
    "Łukasiewicz": {
        "g":    lambda x, p: 1 - x,
        "T":    lambda x, y, p: np.maximum(0, x + y - 1),
        "h":    lambda x, p: np.ones_like(x),
        "param": None,
    },
    "Product": {
        "g":    lambda x, p: -np.log(x),
        "T":    lambda x, y, p: x * y,
        "h":    lambda x, p: 1 / x,
        "param": None,
    },
    "Schweizer–Sklar": {
        "g":    lambda x, p: np.where(p > 0, 1 - x**p, x**p - 1),
        "T":    lambda x, y, p: np.maximum(0, (x**p + y**p - 1))**(1/p) if p != 0 else x*y,
        "h":    lambda x, p: np.abs(p) * x**(p - 1),
        "param": ("λ", -3.0, 3.0, 0.5, 0.5),
    },
    "Yager": {
        "g":    lambda x, p: (1 - x)**p,
        "T":    lambda x, y, p: np.maximum(0, 1 - ((1-x)**p + (1-y)**p)**(1/p)),
        "h":    lambda x, p: p * (1 - x)**(p - 1),
        "param": ("λ≥1", 1.0, 5.0, 1.0, 2.0),
    },
    "Aczel–Alsina": {
        "g":    lambda x, p: (-np.log(np.maximum(x, 1e-300)))**p,
        "T":    lambda x, y, p: np.exp(-((-np.log(np.maximum(x,1e-300)))**p + (-np.log(np.maximum(y,1e-300)))**p)**(1/p)),
        "h":    lambda x, p: p * (-np.log(np.maximum(x, 1e-300)))**(p-1) / np.maximum(x, 1e-300),
        "param": ("λ≥1", 1.0, 5.0, 1.0, 2.0),
    },
    "Dombi": {
        "g":    lambda x, p: ((1 - x) / np.maximum(x, 1e-300))**p,
        "T":    lambda x, y, p: 1 / (1 + (((1-x)/np.maximum(x,1e-300))**p + ((1-y)/np.maximum(y,1e-300))**p)**(1/p)),
        "h":    lambda x, p: p * (1-x)**(p-1) / np.maximum(x, 1e-300)**(p+1),
        "param": ("λ>0", 0.1, 5.0, 0.1, 1.5),
    },
    "Frank": {
        "g":    lambda x, p: np.log((p - 1) / np.maximum(p**x - 1, 1e-300)),
        "T":    lambda x, y, p: np.log(1 + (p**x - 1)*(p**y - 1) / (p - 1)) / np.log(p),
        "h":    lambda x, p: p**x * np.log(p) / np.maximum(p**x - 1, 1e-300),
        "param": ("s>1", 1.01, 20.0, 0.1, 2.0),
    },
    "Hamacher": {
        "g":    lambda x, p: np.log((p + (1 - p) * x) / np.maximum(x, 1e-300)),
        "T":    lambda x, y, p: x*y / (p + (1-p)*(x + y - x*y)),
        "h":    lambda x, p: 1 / (np.maximum(x, 1e-300) * (p + (1-p)*np.maximum(x, 1e-300))),
        "param": ("γ>0", 0.01, 5.0, 0.1, 1.0),
    },
    "Sugeno–Weber": {
        "g":    lambda x, p: np.log((1 + p) / np.maximum(1 + p*x, 1e-300)),
        "T":    lambda x, y, p: np.maximum(0, (x + y - 1 + p*x*y) / (1 + p)),
        "h":    lambda x, p: p / np.maximum(1 + p*x, 1e-300),
        "param": ("λ>−1", -0.99, 10.0, 0.1, 1.0),
    },
}

FAMILY_NAMES = list(FAMILIES.keys())

# --- figure layout ---
fig = plt.figure(figsize=(14, 9))
fig.patch.set_facecolor("#1a1a2e")

ax_radio = fig.add_axes([0.01, 0.3, 0.13, 0.6], facecolor="#16213e")
ax_slider = fig.add_axes([0.01, 0.18, 0.13, 0.06], facecolor="#16213e")

ax_g   = fig.add_axes([0.18, 0.55, 0.25, 0.38])
ax_h   = fig.add_axes([0.50, 0.55, 0.25, 0.38])
ax_T   = fig.add_axes([0.18, 0.08, 0.25, 0.38])
ax_sel = fig.add_axes([0.50, 0.08, 0.25, 0.38])
ax_info = fig.add_axes([0.78, 0.08, 0.20, 0.85], facecolor="#16213e")
ax_info.axis("off")

for ax in [ax_g, ax_h, ax_T, ax_sel]:
    ax.set_facecolor("#0f3460")
    for sp in ax.spines.values():
        sp.set_color("#e94560")

radio = RadioButtons(ax_radio, FAMILY_NAMES, active=1,
                     label_props={"color": ["#e0e0e0"]*len(FAMILY_NAMES), "fontsize": [8]*len(FAMILY_NAMES)})
try:
    for circle in radio.circles:
        circle.set_facecolor("#e94560")
except AttributeError:
    pass

slider_ax_obj = Slider(ax_slider, "param", 0.0, 1.0, valinit=0.5,
                       color="#e94560", track_color="#16213e")
slider_ax_obj.label.set_color("white")
slider_ax_obj.valtext.set_color("white")

state = {"family": "Product", "param": 1.0}

INFO = {
    "Łukasiewicz": "g(x) = 1−x\nNilpotent\ng(0)=1 finite\n\nT(x,y)=max(0,x+y−1)\nh(x)=1 (constant)\n\nProfile flat → equal\nsensitivity everywhere",
    "Product":     "g(x) = −ln x\nStrict\ng(0⁺)=+∞\n\nT(x,y)=xy\nh(x)=1/x\n\nDecreasing profile\n→ weak inputs amplified",
    "Schweizer–Sklar": "g(x)=1−xᵅ (λ>0)\ng(x)=xᵅ−1 (λ<0)\n\nλ=1 → Łukasiewicz\nλ=0 → Product\nλ=−1 → Hamacher γ=0\nλ→−∞ → Min\n\nUnifies nilpotent & strict",
    "Yager":       "g(x)=(1−x)ᵅ\nNilpotent\n\nλ=1 → Łukasiewicz\nλ→∞ → Min\n\nh decreasing iff λ≥1\n→ robustness condition",
    "Aczel–Alsina":"g(x)=(−ln x)ᵅ\nStrict\n\nλ=1 → Product\nλ→∞ → Min\n\nCorresponds to\nGumbel–Hougaard copula",
    "Dombi":       "g(x)=((1−x)/x)ᵅ\nStrict\n\nλ=1 → Hamacher γ=0\nλ→∞ → Min\nλ→0⁺ → Drastic product\n\nh decr. iff λ≥1",
    "Frank":       "g(x)=ln((s−1)/(sˣ−1))\nStrict (s>1)\n\nOnly family where\nT+S=x+y (Frank 1979)\n\ns→0 → Łukasiewicz\ns=1 → Product\ns→∞ → Min",
    "Hamacher":    "g(x)=ln((γ+(1−γ)x)/x)\nStrict (γ>0)\n\nγ=1 → Product\nγ=2 → Einstein product\nγ→0 → xy/(x+y−xy)\n\n≡ AMH copula (θ=1−γ)",
    "Sugeno–Weber":"g(x)=ln((1+λ)/(1+λx))\nNilpotent\n\nλ=0 → Łukasiewicz\nλ→∞ → Product\n\nh(x)=λ/(1+λx)\ndecreasing ✓",
}

def compute_and_plot(family_name, param_val):
    fam = FAMILIES[family_name]
    p = param_val

    # generator
    try:
        gv = fam["g"](x, p)
        gv = np.clip(gv, -1e4, 1e4)
    except Exception:
        gv = np.zeros_like(x)

    # profile h
    try:
        hv = fam["h"](x[1:-1], p)
        hv = np.clip(hv, 0, 50)
    except Exception:
        hv = np.zeros(len(x)-2)

    # T surface
    try:
        Tv = fam["T"](X, Y, p)
        Tv = np.clip(Tv, 0, 1)
    except Exception:
        Tv = np.zeros_like(X)

    # selectivity: gamma(t,a) = h(t)/h(a) for t=0.2, varying a
    a_vals = np.linspace(0.21, 0.99, 200)
    t_fixed = 0.2
    try:
        ht = float(fam["h"](np.array([t_fixed]), p)[0])
        ha = fam["h"](a_vals, p)
        gamma = np.where(ha > 1e-10, ht / ha, np.nan)
        gamma = np.clip(gamma, 0, 20)
    except Exception:
        gamma = np.full(200, np.nan)

    for ax in [ax_g, ax_h, ax_T, ax_sel]:
        ax.cla()
        ax.set_facecolor("#0f3460")
        for sp in ax.spines.values():
            sp.set_color("#e94560")
        ax.tick_params(colors="white", labelsize=7)

    # plot g
    ax_g.plot(x, gv, color="#e94560", lw=2)
    ax_g.set_title(f"Generator g(x)", color="white", fontsize=9)
    ax_g.set_xlabel("x", color="white", fontsize=8)
    ax_g.axhline(0, color="gray", lw=0.5)
    ax_g.grid(True, alpha=0.2)

    # plot h
    ax_h.plot(x[1:-1], hv, color="#f5a623", lw=2)
    ax_h.set_title("Profile h(x) = −g′(x)", color="white", fontsize=9)
    ax_h.set_xlabel("x", color="white", fontsize=8)
    ax_h.grid(True, alpha=0.2)

    # plot T heatmap
    im = ax_T.contourf(X, Y, Tv, levels=20, cmap="plasma")
    ax_T.set_title("T(x, y) heatmap", color="white", fontsize=9)
    ax_T.set_xlabel("x", color="white", fontsize=8)
    ax_T.set_ylabel("y", color="white", fontsize=8)

    # plot selectivity
    ax_sel.plot(a_vals, gamma, color="#50fa7b", lw=2)
    ax_sel.axhline(1, color="gray", lw=0.8, ls="--")
    ax_sel.set_title(f"Selectivity Γ = h(t)/h(a), t={t_fixed}", color="white", fontsize=9)
    ax_sel.set_xlabel("a (background)", color="white", fontsize=8)
    ax_sel.set_ylabel("Γ", color="white", fontsize=8)
    ax_sel.grid(True, alpha=0.2)

    # info panel
    ax_info.cla()
    ax_info.axis("off")
    pstr = f"  param = {param_val:.3f}" if fam["param"] else "  (no param)"
    text = f"{family_name}\n{pstr}\n\n{INFO.get(family_name,'')}"
    ax_info.text(0.05, 0.95, text, transform=ax_info.transAxes,
                 color="white", fontsize=8, va="top", family="monospace",
                 wrap=True)

    fig.suptitle(f"T-norm Explorer  —  {family_name}", color="white", fontsize=12, y=0.99)
    fig.canvas.draw_idle()

def update_family(label):
    state["family"] = label
    fam = FAMILIES[label]
    if fam["param"]:
        pname, vmin, vmax, vstep, vinit = fam["param"]
        slider_ax_obj.label.set_text(pname)
        slider_ax_obj.valmin = vmin
        slider_ax_obj.valmax = vmax
        slider_ax_obj.val = vinit
        slider_ax_obj.ax.set_xlim(vmin, vmax)
        slider_ax_obj.set_val(vinit)
        state["param"] = vinit
        slider_ax_obj.ax.set_visible(True)
    else:
        slider_ax_obj.ax.set_visible(False)
        state["param"] = 1.0
    compute_and_plot(label, state["param"])

def update_slider(val):
    state["param"] = slider_ax_obj.val
    compute_and_plot(state["family"], state["param"])

radio.on_clicked(update_family)
slider_ax_obj.on_changed(update_slider)

update_family("Product")

plt.show()
