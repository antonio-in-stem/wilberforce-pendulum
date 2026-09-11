"""Shared matplotlib style and colour palette for the Wilberforce figures.

A single import point so every figure in the paper shares one visual identity:
Computer-Modern-like serif type, a restrained palette, clean spines, and a
consistent grid. Import `apply_style()` once, then use the `C` colour dict.
"""
import matplotlib as mpl
import matplotlib.pyplot as plt

# --- palette ------------------------------------------------------------------
C = {
    "long":     "#2b6cb0",   # longitudinal (z)      — blue
    "long_lt":  "#9ec5e8",   # longitudinal, light
    "ang":      "#7c3aed",   # angular (theta)       — violet
    "ang_lt":   "#c4a7f0",   # angular, light
    "coup":     "#d97706",   # coupling              — amber
    "total":    "#0f9d78",   # total energy          — teal
    "accent":   "#c026d3",   # phase / highlight      — magenta
    "num":      "#1a1a2e",   # numeric markers        — near-black
    "ref":      "#e11d48",   # reference / rose
    "grid":     "#c9ccd6",
    "ink":      "#1a1a2e",
}


def apply_style():
    mpl.rcParams.update({
        "figure.dpi": 130,
        "savefig.dpi": 220,
        "savefig.bbox": "tight",
        "font.family": "serif",
        "font.serif": ["DejaVu Serif", "CMU Serif", "Times New Roman"],
        "mathtext.fontset": "cm",
        "font.size": 12,
        "axes.titlesize": 13,
        "axes.labelsize": 12.5,
        "axes.edgecolor": "#3a3a48",
        "axes.linewidth": 1.0,
        "axes.grid": True,
        "axes.axisbelow": True,
        "grid.color": C["grid"],
        "grid.linewidth": 0.6,
        "grid.alpha": 0.7,
        "legend.frameon": True,
        "legend.framealpha": 0.92,
        "legend.edgecolor": "#d0d0d8",
        "legend.fontsize": 10.5,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "lines.antialiased": True,
    })


def finish(ax):
    """Trim top/right spines for a cleaner look."""
    if hasattr(ax, "spines"):
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
