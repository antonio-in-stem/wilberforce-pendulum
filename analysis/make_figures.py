#!/usr/bin/env python
"""
make_figures.py — recompute figures without overwriting historical paper assets.

Figures are written to ../generated/figures/ as both PDF (vector, for LaTeX) and PNG
(for quick inspection). Computation uses the validated double-precision reference
in wilberforce.py; where the compiled Fortran binary is available it is run so
that the validation figure overlays the *actual* Fortran output.

    python analysis/make_figures.py            # all figures
    python analysis/make_figures.py fig_beat   # one figure by name
"""
from __future__ import annotations

import os
import sys
import subprocess
import tempfile
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch, Arc
from matplotlib.collections import LineCollection

sys.path.insert(0, os.path.dirname(__file__))
import wilberforce as w
from figstyle import apply_style, finish, C

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.normpath(os.path.join(HERE, "..", "generated", "figures"))
os.makedirs(FIGDIR, exist_ok=True)
apply_style()
_TEMP = tempfile.TemporaryDirectory(prefix="wilberforce-figures-")


# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------
def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGDIR, f"{name}.{ext}"))
    plt.close(fig)
    print(f"  wrote {name}.pdf / .png")


def run_preset(pid, tf=None, h=None):
    """Integrate a preset with the reference RK4; return (t, U, params, meta)."""
    p, d = w.from_preset(pid)
    tf = d["tf"] if tf is None else tf
    h = d["h"] if h is None else h
    u0 = [d["z0"], 0.0, d["theta0"], 0.0]
    n = int(round(tf / h))
    t, U = w.integrate_rk4(p, u0, 0.0, h, n)
    return t, U, p, d


def _run_fortran(pid, outdir):
    """Run the compiled Fortran binary for a preset into <outdir>; return the
    data directory path, or None if the binary is unavailable."""
    for exe in (("build/wilberforce.exe",) if os.name == "nt" else ("build/wilberforce",)):
        exe_path = os.path.join(HERE, "..", exe)
        if os.path.exists(exe_path):
            od = os.path.join(_TEMP.name, str(pid))
            try:
                subprocess.run([os.path.abspath(exe_path), str(pid), "--outdir", od],
                               cwd=os.path.join(HERE, ".."), capture_output=True, timeout=60, check=True)
                return os.path.join(od, "data")
            except (OSError, subprocess.SubprocessError) as exc:
                raise RuntimeError("Available Fortran solver failed; no silent fallback") from exc
    return None


def try_fortran(pid, outdir):
    """Return (t, z_num, th_num) from the Fortran core, or None."""
    d = _run_fortran(pid, outdir)
    if d is None:
        return None
    try:
        dz = np.loadtxt(os.path.join(d, "posicion_longitudinal.txt"))
        dth = np.loadtxt(os.path.join(d, "posicion_angular.txt"))
        return dz[:, 0], dz[:, 1], dth[:, 1]
    except (OSError, ValueError) as exc:
        raise RuntimeError("Invalid Fortran output; no silent fallback") from exc


def fortran_energy(pid, outdir="_figtmp"):
    """Return the Fortran core's energy time series as columns
    (t, E_long, E_ang, E_coup, E_total), or None if the binary is unavailable.
    The original stream starts at h. Prepend the actual initial energy so an
    E(t)-E(0) plot does not silently use E(h) as its reference. These are current
    results, not replacements for the historical paper's conservation numbers."""
    d = _run_fortran(pid, outdir)
    if d is None:
        print("  (Fortran binary not found -> energy figures use the Python "
              "reference; build with `make` for the quad-precision numbers)")
        return None
    try:
        data = np.loadtxt(os.path.join(d, "energias.txt"))
        p, meta = w.from_preset(pid)
        e0 = w.energy(np.array([[meta["z0"], 0.0, meta["theta0"], 0.0]]), p)
        first = [0.0, e0["long"][0], e0["ang"][0], e0["coupling"][0], e0["total"][0]]
        return np.vstack((first, data))
    except (OSError, ValueError) as exc:
        raise RuntimeError("Invalid Fortran output; no silent fallback") from exc


def _period_from_crossings(t, x):
    """Fundamental period from interpolated upward zero-crossings of x(t)."""
    x = np.asarray(x, float) - np.mean(x)
    cr = []
    for i in range(len(x) - 1):
        if x[i] <= 0.0 < x[i + 1]:
            cr.append(t[i] - x[i] * (t[i + 1] - t[i]) / (x[i + 1] - x[i]))
    cr = np.asarray(cr)
    return float(np.mean(np.diff(cr))) if len(cr) >= 2 else float("nan")


# ----------------------------------------------------------------------------
# Fig 1 — schematic of the Wilberforce pendulum
# ----------------------------------------------------------------------------
def fig_schematic():
    fig, ax = plt.subplots(figsize=(5.4, 6.6))
    ax.set_xlim(-3.2, 3.6); ax.set_ylim(-4.6, 4.2); ax.set_aspect("equal")
    ax.axis("off")

    # ceiling
    ax.add_patch(Rectangle((-2.2, 3.6), 4.4, 0.35, facecolor="#4a4954", edgecolor="none"))
    for xx in np.linspace(-2.1, 2.0, 14):
        ax.plot([xx, xx - 0.28], [3.6, 3.95], color="#4a4954", lw=1.2)

    # helical spring (side view: sinusoidal coil)
    n_coils = 9
    tt = np.linspace(0, n_coils * 2 * np.pi, 700)
    yspring = np.linspace(3.6, 0.6, tt.size)
    xspring = 0.55 * np.sin(tt)
    ax.plot(xspring, yspring, color=C["long"], lw=2.2, solid_capstyle="round")
    ax.plot([0, 0], [0.6, 0.2], color=C["long"], lw=2.2)

    # cylindrical bob
    bob = Rectangle((-1.3, -1.8), 2.6, 2.0, facecolor="#b0adc7", edgecolor="#4a4954", lw=1.6)
    ax.add_patch(bob)
    ax.add_patch(plt.matplotlib.patches.Ellipse((0, 0.2), 2.6, 0.5,
                 facecolor="#c9c6dd", edgecolor="#4a4954", lw=1.6))
    ax.add_patch(plt.matplotlib.patches.Ellipse((0, -1.8), 2.6, 0.5,
                 facecolor="#9a97b4", edgecolor="#4a4954", lw=1.6))

    # rod with two slider masses (m_d)
    ax.plot([-2.9, 2.9], [-0.8, -0.8], color="#4a4954", lw=2.6)
    for sx in (-2.35, 2.35):
        ax.add_patch(Rectangle((sx - 0.28, -1.12), 0.56, 0.64,
                     facecolor="#7c3aed", edgecolor="#3a2560", lw=1.4))
    ax.text(2.62, -1.55, r"$m_d$", color="#3a2560", fontsize=12)

    # z coordinate arrow
    ax.add_patch(FancyArrowPatch((2.9, 0.2), (2.9, -1.8), arrowstyle="<->",
                 mutation_scale=15, color=C["long"], lw=1.8))
    ax.text(3.05, -0.8, r"$z$", color=C["long"], fontsize=15, va="center")

    # theta (torsion) curved arrow
    ax.add_patch(Arc((0, -2.35), 2.4, 0.9, angle=0, theta1=200, theta2=340,
                 color=C["ang"], lw=2.0))
    ax.add_patch(FancyArrowPatch((1.05, -2.55), (1.18, -2.28), arrowstyle="-|>",
                 mutation_scale=15, color=C["ang"], lw=2.0))
    ax.text(0, -3.15, r"$\theta$", color=C["ang"], fontsize=15, ha="center")

    # axis line
    ax.plot([0, 0], [-2.0, -3.6], color="#8a8a96", lw=1.0, ls=(0, (4, 3)))
    ax.text(-2.95, 3.15, r"$k,\ \delta$", fontsize=13, color="#3a3a48")
    ax.text(-1.9, -0.8, r"$m,\ I$", fontsize=13, color="#2a2a38")
    ax.set_title("Wilberforce pendulum", pad=6)
    save(fig, "fig_schematic")


# ----------------------------------------------------------------------------
# Fig 2 — the beat: coupled longitudinal & torsional motion (preset 14)
# ----------------------------------------------------------------------------
def fig_beat():
    t, U, p, d = run_preset(14)
    z, th = U[:, 0], U[:, 2]
    za, va, tha, oma = w.analytic(p, t)

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.4, 5.4), sharex=True)
    a1.plot(t, za, color=C["long"], lw=1.4, label="analytic")
    a1.plot(t[::12], z[::12], ls="none", marker="o", ms=2.6, color=C["num"],
            alpha=0.7, label="RK4")
    a1.set_ylabel(r"$z(t)\ \mathrm{[m]}$")
    a1.legend(loc="upper right", ncol=2)
    a1.set_title("Beat: energy shuttles between the longitudinal and torsional modes")

    a2.plot(t, tha, color=C["ang"], lw=1.4, label="analytic")
    a2.plot(t[::12], th[::12], ls="none", marker="o", ms=2.6, color=C["num"],
            alpha=0.7, label="RK4")
    a2.set_ylabel(r"$\theta(t)\ \mathrm{[rad]}$")
    a2.set_xlabel(r"$t\ \mathrm{[s]}$")
    a2.legend(loc="upper right", ncol=2)

    for a in (a1, a2):
        finish(a)
        a.axvline(p.beat_period / 2, color=C["accent"], lw=1.0, ls="--", alpha=0.7)
    a1.text(p.beat_period / 2 + 0.4, a1.get_ylim()[1] * 0.82,
            r"$T_{\rm beat}/2$", color=C["accent"], fontsize=10)
    save(fig, "fig_beat")


# ----------------------------------------------------------------------------
# Fig 3 — energy exchange and conservation (preset 14)
# ----------------------------------------------------------------------------
def fig_energy_exchange():
    fe = fortran_energy(14)
    if fe is not None:
        t, Elong, Eang, Ecoup, Etot = fe[:, 0], fe[:, 1], fe[:, 2], fe[:, 3], fe[:, 4]
    else:
        t, U, p, d = run_preset(14)
        E = w.energy(U, p)
        Elong, Eang, Ecoup, Etot = E["long"], E["ang"], E["coupling"], E["total"]
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    ax.plot(t, Elong, color=C["long"], lw=1.7, label="longitudinal  " + r"$E_z$")
    ax.plot(t, Eang, color=C["ang"], lw=1.7, label="torsional  " + r"$E_\theta$")
    ax.plot(t, Ecoup, color=C["coup"], lw=1.3, label="coupling  " + r"$E_c$")
    ax.plot(t, Etot, color=C["total"], lw=2.2, label="total  " + r"$E$")
    ax.set_xlabel(r"$t\ \mathrm{[s]}$")
    ax.set_ylabel(r"energy $\mathrm{[J]}$")
    ax.set_title("Energy exchange between modes; total energy conserved")
    ax.legend(loc="center right", ncol=1)
    finish(ax)
    save(fig, "fig_energy_exchange")


# ----------------------------------------------------------------------------
# Fig 4 — phase portrait (theta vs z)
# ----------------------------------------------------------------------------
def fig_phase():
    t, U, p, d = run_preset(14)
    z, th = U[:, 0], U[:, 2]
    pts = np.array([th, z]).T.reshape(-1, 1, 2)
    segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
    lc = LineCollection(segs, cmap="viridis", array=t, lw=1.2)
    fig, ax = plt.subplots(figsize=(6.0, 5.2))
    ax.add_collection(lc)
    ax.set_xlim(th.min() * 1.05, th.max() * 1.05)
    ax.set_ylim(z.min() * 1.1, z.max() * 1.1)
    ax.set_xlabel(r"$\theta\ \mathrm{[rad]}$")
    ax.set_ylabel(r"$z\ \mathrm{[m]}$")
    ax.set_title("Configuration-space trajectory")
    cb = fig.colorbar(lc, ax=ax, pad=0.02)
    cb.set_label(r"$t\ \mathrm{[s]}$")
    finish(ax)
    save(fig, "fig_phase")


# ----------------------------------------------------------------------------
# Fig 5 — normal-mode spectrum (FFT reveals the two frequencies)
# ----------------------------------------------------------------------------
def fig_spectrum():
    p, d = w.from_preset(14)
    h, T = 0.01, 400.0
    u0 = [d["z0"], 0.0, d["theta0"], 0.0]
    n = int(T / h)
    t, U = w.integrate_rk4(p, u0, 0.0, h, n)
    th = U[:, 2] - U[:, 2].mean()
    win = np.hanning(th.size)
    F = np.abs(np.fft.rfft(th * win))
    f = np.fft.rfftfreq(th.size, h) * 2 * np.pi  # angular frequency
    F /= F.max()
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.plot(f, F, color=C["ang"], lw=1.6)
    ax.set_xlim(1.8, 2.9)
    for wv, name, col in ((p.w1, r"$\omega_1$", C["ref"]),
                          (p.w2, r"$\omega_2$", C["ref"]),
                          (p.w0, r"$\omega_0$", "#8a8a96")):
        ax.axvline(wv, color=col, lw=1.1, ls="--", alpha=0.8)
        ax.text(wv, 1.02, name, color=col, ha="center", fontsize=11)
    ax.set_xlabel(r"angular frequency $\omega\ \mathrm{[rad/s]}$")
    ax.set_ylabel("normalised spectral power")
    ax.set_title("Two normal modes split symmetrically about " + r"$\omega_0$")
    finish(ax)
    save(fig, "fig_spectrum")


# ----------------------------------------------------------------------------
# Fig 6 — validation: RK4 & DOP853 residuals against the analytic solution
# ----------------------------------------------------------------------------
def fig_validation():
    p, d = w.from_preset(14)
    u0 = [d["z0"], 0.0, d["theta0"], 0.0]
    t, U = w.integrate_rk4(p, u0, 0.0, d["h"], int(d["tf"] / d["h"]))
    za, va, tha, oma = w.analytic(p, t)
    Uref = w.integrate_ref(p, u0, t, t0=0.0)
    err_ref = np.abs(Uref[:, 2] - tha)

    fig, ax = plt.subplots(figsize=(7.4, 4.3))
    # "our solver": prefer the compiled Fortran core; fall back to Python RK4.
    ft = try_fortran(14, "_figtmp")
    if ft is not None:
        tf_, zf_, thf_ = ft
        _, _, thaf, _ = w.analytic(p, tf_)
        ax.semilogy(tf_, np.maximum(np.abs(thf_ - thaf), 1e-16), color=C["long"],
                    lw=1.3, label=r"our RK4 solver ($h{=}0.01$) $-$ analytic")
    else:
        ax.semilogy(t, np.maximum(np.abs(U[:, 2] - tha), 1e-16), color=C["long"],
                    lw=1.3, label=r"our RK4 solver ($h{=}0.01$) $-$ analytic")
    ax.semilogy(t, np.maximum(err_ref, 1e-16), color=C["ref"], lw=1.3,
                label=r"independent DOP853 $-$ analytic")
    ax.set_xlabel(r"$t\ \mathrm{[s]}$")
    ax.set_ylabel(r"$|\theta_{\rm num}-\theta_{\rm ana}|\ \mathrm{[rad]}$")
    ax.set_title("Cross-validation against the closed-form solution")
    ax.legend(loc="lower right")
    ax.set_ylim(1e-13, 1e-4)
    ax.grid(True, which="both", alpha=0.4)
    finish(ax)
    save(fig, "fig_validation")


# ----------------------------------------------------------------------------
# Fig 7 — RK4 fourth-order convergence
# ----------------------------------------------------------------------------
def fig_convergence():
    p, d = w.from_preset(14)
    u0 = [d["z0"], 0.0, d["theta0"], 0.0]
    T = 20.0
    hs = np.array([0.16, 0.08, 0.04, 0.02, 0.01, 0.005, 0.0025])
    errs = []
    for h in hs:
        n = int(round(T / h))
        t, U = w.integrate_rk4(p, u0, 0.0, h, n)
        _, _, tha, _ = w.analytic(p, t)
        errs.append(np.max(np.abs(U[:, 2] - tha)))
    errs = np.array(errs)
    # reference slope-4 line
    c = errs[2] / hs[2] ** 4
    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    ax.loglog(hs, errs, "o-", color=C["long"], lw=1.6, ms=6, label="measured max error")
    ax.loglog(hs, c * hs ** 4, "--", color=C["ref"], lw=1.4, label=r"$\propto h^4$ (theory)")
    # empirical order
    order = np.polyfit(np.log(hs[:-1]), np.log(errs[:-1]), 1)[0]
    ax.set_xlabel(r"step size $h\ \mathrm{[s]}$")
    ax.set_ylabel(r"global error in $\theta$  $\mathrm{[rad]}$")
    ax.set_title(f"Fourth-order convergence of RK4 (fitted order $= {order:.2f}$)")
    ax.legend(loc="upper left")
    ax.grid(True, which="both", alpha=0.4)
    finish(ax)
    save(fig, "fig_convergence")
    return order


# ----------------------------------------------------------------------------
# Fig 8 — coupling-constant sweep: beat period vs epsilon
# ----------------------------------------------------------------------------
def fig_coupling_sweep():
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.6, 4.5))

    # (a) longitudinal energy fraction over time for a few epsilon (same
    #     delta, m, I as presets 6-10 but chosen for clean, well-separated beats)
    base = w.Params(delta=w._s(2.3), eps=0.0, m=w._s(0.1), I=w._s(0.8),
                    z0=w._s(0.08), theta0=0.0)
    T = 130.0
    for eps, col in ((0.20, C["ref"]), (0.10, C["coup"]), (0.05, C["ang"])):
        pe = w.Params(delta=base.delta, eps=eps, m=base.m, I=base.I,
                      z0=base.z0, theta0=base.theta0)
        u0 = [pe.z0, 0.0, pe.theta0, 0.0]
        tt, U = w.integrate_rk4(pe, u0, 0.0, 0.01, int(T / 0.01))
        E = w.energy(U, pe)
        a1.plot(tt, E["long"] / E["total"], color=col, lw=1.6,
                label=fr"$\varepsilon={eps:.2g}$  ($T_{{\rm beat}}={pe.beat_period:.0f}$ s)")
    a1.set_xlabel(r"$t\ \mathrm{[s]}$")
    a1.set_ylabel(r"$E_z / E_{\rm total}$")
    a1.set_title("Slower beats as coupling weakens")
    a1.legend(loc="upper right")
    finish(a1)

    # (b) beat period vs epsilon: exact vs small-eps approximation
    eps = np.linspace(0.01, 1.8, 200)
    p0, _ = w.from_preset(6)  # delta=2.3, m=0.1, I=0.8
    w0 = p0.w0
    sqrtmI = np.sqrt(p0.m * p0.I)
    w1 = np.sqrt(w0**2 + eps / (2 * sqrtmI))
    w2 = np.sqrt(np.clip(w0**2 - eps / (2 * sqrtmI), 1e-9, None))
    Tb_exact = 2 * np.pi / (w1 - w2)
    Tb_approx = 4 * np.pi * w0 * sqrtmI / eps
    a2.plot(eps, Tb_exact, color=C["long"], lw=1.8, label=r"exact $2\pi/(\omega_1-\omega_2)$")
    a2.plot(eps, Tb_approx, "--", color=C["ref"], lw=1.5, label=r"$4\pi\omega_0\sqrt{mI}/\varepsilon$")
    # markers at the preset values
    for pid in (6, 7, 8, 9, 10):
        pp, _ = w.from_preset(pid)
        a2.plot(pp.eps, pp.beat_period, "o", color=C["num"], ms=6, zorder=5)
    a2.set_yscale("log")
    a2.set_xlabel(r"coupling constant $\varepsilon$")
    a2.set_ylabel(r"beat period $T_{\rm beat}\ \mathrm{[s]}$")
    a2.set_title("Beat period diverges as " + r"$\varepsilon\to0$")
    a2.legend(loc="upper right")
    finish(a2)
    save(fig, "fig_coupling_sweep")


# ----------------------------------------------------------------------------
# Fig 9 — measuring epsilon via the periodicity condition
# ----------------------------------------------------------------------------
def fig_coupling_measurement():
    t11, U11, p11, _ = run_preset(11, tf=100.0)   # theta0 = 0  -> non-periodic
    t12, U12, p12, _ = run_preset(12, tf=100.0)   # theta0 = z0 sqrt(m/I) -> periodic

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.6, 4.4))
    a1.plot(t11, U11[:, 0], color=C["ang"], lw=1.3, alpha=0.9,
            label=r"$\theta_0=0$ (two modes)")
    a1.plot(t12, U12[:, 0], color=C["long"], lw=1.5,
            label=r"$\theta_0=z_0\sqrt{m/I}$ (one mode)")
    a1.set_xlim(0, 40)
    a1.set_xlabel(r"$t\ \mathrm{[s]}$"); a1.set_ylabel(r"$z(t)\ \mathrm{[m]}$")
    a1.set_title("Periodicity condition kills the second mode")
    a1.legend(loc="upper right")
    finish(a1)

    # spectra — integrate long (450 s) so the two close modes are resolved
    Tspec = 450.0
    for pid, col, lab in ((11, C["ang"], "two modes"), (12, C["long"], "one mode")):
        ts, Us, ps, _ = run_preset(pid, tf=Tspec)
        s = Us[:, 0] - Us[:, 0].mean()
        win = np.hanning(s.size)
        F = np.abs(np.fft.rfft(s * win)); F /= F.max()
        f = np.fft.rfftfreq(s.size, 0.01) * 2 * np.pi
        a2.plot(f, F, color=col, lw=1.5, label=lab)
    a2.axvline(p11.w1, color=C["ref"], ls="--", lw=1.0)
    a2.axvline(p11.w2, color=C["ref"], ls="--", lw=1.0)
    a2.text(p11.w1, 1.02, r"$\omega_1$", color=C["ref"], ha="center", fontsize=10)
    a2.text(p11.w2, 1.02, r"$\omega_2$", color=C["ref"], ha="center", fontsize=10)
    a2.set_xlim(1.55, 2.05)
    a2.set_xlabel(r"$\omega\ \mathrm{[rad/s]}$"); a2.set_ylabel("spectral power")
    a2.set_title("Spectrum collapses to a single line")
    a2.legend(loc="upper right")
    finish(a2)

    # recovered epsilon: measure the period FROM the single-mode trajectory
    # (not from the exact frequency) so this is a genuine end-to-end measurement.
    T = _period_from_crossings(t12, U12[:, 0])
    eps_rec = 2 * np.sqrt(p12.m * p12.I) * ((2 * np.pi / T) ** 2 - p12.delta / p12.I)
    a1.text(0.03, 0.05, fr"$T={T:.3f}\,$s:$\ \varepsilon_{{\rm true}}={p12.eps:.4f},\ "
            fr"\varepsilon_{{\rm meas}}={eps_rec:.4f}$",
            transform=a1.transAxes, fontsize=10.0,
            bbox=dict(boxstyle="round", fc="white", ec="#d0d0d8"))
    save(fig, "fig_coupling_measurement")


# ----------------------------------------------------------------------------
# Fig 10 — 3D trajectory of the bob rim
# ----------------------------------------------------------------------------
def fig_trajectory3d():
    from mpl_toolkits.mplot3d import Axes3D  # noqa
    t, U, p, d = run_preset(14)
    x, y, z = w.bob_path_3d(U, radius=1.0, z_scale=6.0)
    fig = plt.figure(figsize=(6.6, 5.6))
    ax = fig.add_subplot(111, projection="3d")
    pts = np.array([x, y, z]).T.reshape(-1, 1, 3)
    segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
    from mpl_toolkits.mplot3d.art3d import Line3DCollection
    lc = Line3DCollection(segs, cmap="plasma", array=t, lw=1.1)
    ax.add_collection3d(lc)
    ax.set_xlim(-1.2, 1.2); ax.set_ylim(-1.2, 1.2); ax.set_zlim(z.min(), z.max())
    ax.set_xlabel(r"$x=r\cos\theta$"); ax.set_ylabel(r"$y=r\sin\theta$")
    ax.set_zlabel(r"$z\ \mathrm{(scaled)}$")
    ax.set_title("Motion of a point on the rim of the bob")
    ax.view_init(elev=22, azim=-58)
    cb = fig.colorbar(lc, ax=ax, pad=0.1, shrink=0.6)
    cb.set_label(r"$t\ \mathrm{[s]}$")
    save(fig, "fig_trajectory3d")


# ----------------------------------------------------------------------------
# Fig 11 — energy drift (conservation quality of RK4)
# ----------------------------------------------------------------------------
def fig_energy_drift():
    fe = fortran_energy(14)
    if fe is not None:
        t, Etot = fe[:, 0], fe[:, 4]
    else:
        t, U, p, d = run_preset(14)
        Etot = w.energy(U, p)["total"]
    drift = np.abs(Etot - Etot[0])
    rel = drift.max() / abs(Etot[0])
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    ax.semilogy(t[1:], drift[1:], color=C["total"], lw=1.5)
    ax.set_xlabel(r"$t\ \mathrm{[s]}$")
    ax.set_ylabel(r"$|E(t)-E(0)|\ \mathrm{[J]}$")
    ax.set_ylim(1e-13, 3e-9)
    ax.set_title(fr"Energy drift bounded at $\sim{drift.max():.1e}$ J "
                 fr"(${rel:.1e}$ relative)")
    ax.grid(True, which="both", alpha=0.4)
    finish(ax)
    save(fig, "fig_energy_drift")


ALL = {
    "fig_schematic": fig_schematic,
    "fig_beat": fig_beat,
    "fig_energy_exchange": fig_energy_exchange,
    "fig_phase": fig_phase,
    "fig_spectrum": fig_spectrum,
    "fig_validation": fig_validation,
    "fig_convergence": fig_convergence,
    "fig_coupling_sweep": fig_coupling_sweep,
    "fig_coupling_measurement": fig_coupling_measurement,
    "fig_trajectory3d": fig_trajectory3d,
    "fig_energy_drift": fig_energy_drift,
}


def main():
    targets = sys.argv[1:] or list(ALL)
    print(f"Generating {len(targets)} figure(s) -> {FIGDIR}")
    for name in targets:
        if name not in ALL:
            raise SystemExit(f"Unknown figure: {name}")
        ALL[name]()
    _TEMP.cleanup()


if __name__ == "__main__":
    main()
