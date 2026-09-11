#!/usr/bin/env python
"""
render_animation.py — head-less 3D renderer for the Wilberforce pendulum.

Produces reproducible animation assets from the *same* physics as the Fortran
core and the interactive web simulator:

    python simulator3d/render_animation.py --preset 14 --out anim.gif
    python simulator3d/render_animation.py --preset 1  --filmstrip frames.pdf

  --preset N     built-in configuration 1..14 (default 14)
  --out FILE     write an animated GIF (or .mp4 if ffmpeg is available)
  --filmstrip F  write a static multi-panel PDF/PNG (good for the paper)
  --frames K     number of frames (default 160)
  --fps F        frames per second for the GIF/MP4 (default 25)

The bob's vertical position follows z(t); the cylinder twists by theta(t); the
helical spring stretches between the ceiling and the bob.
"""
from __future__ import annotations
import os, sys, argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "analysis"))
import wilberforce as w

# --- visual constants ---------------------------------------------------------
CEIL = 3.0
BOB_H = 1.05          # cylinder height
R_BOB = 0.85          # cylinder radius
R_SPR = 0.5           # spring radius
COILS = 9
VZ = 1.25             # visual bob travel amplitude


def draw_frame(ax, z_disp, theta, zamp):
    ax.clear()
    ax.set_axis_off()
    ax.set_xlim(-1.6, 1.6); ax.set_ylim(-1.6, 1.6); ax.set_zlim(-2.6, CEIL + 0.4)
    ax.set_box_aspect((1, 1, 1.9))

    drop = np.clip(z_disp / zamp, -1, 1) * VZ
    bob_c = -0.2 - drop
    bob_top = bob_c + BOB_H / 2
    bob_bot = bob_c - BOB_H / 2

    # ceiling plate
    gx, gy = np.meshgrid([-1.3, 1.3], [-1.3, 1.3])
    ax.plot_surface(gx, gy, np.full_like(gx, CEIL), color="#4a4a5c", alpha=0.6,
                    shade=False, zorder=0)

    # helical spring
    u = np.linspace(0, 1, 500)
    ang = u * COILS * 2 * np.pi + 0.4 * theta
    zc = CEIL - 0.1 - u * ((CEIL - 0.1) - bob_top)
    ax.plot(R_SPR * np.cos(ang), R_SPR * np.sin(ang), zc, color="#3f86e0", lw=2.2)

    # cylinder bob
    phi = np.linspace(0, 2 * np.pi, 48) + theta
    Zc = np.array([[bob_bot], [bob_top]])
    Xc = (R_BOB * np.cos(phi))[None, :] * np.ones((2, 1))
    Yc = (R_BOB * np.sin(phi))[None, :] * np.ones((2, 1))
    Zg = np.repeat(Zc, phi.size, axis=1)
    ax.plot_surface(Xc, Yc, Zg, color="#b9bcd8", alpha=0.97, shade=True, zorder=5)
    # caps
    ax.plot(R_BOB * np.cos(phi), R_BOB * np.sin(phi), np.full_like(phi, bob_top),
            color="#e6e8f6", lw=1.5)
    ax.plot(R_BOB * np.cos(phi), R_BOB * np.sin(phi), np.full_like(phi, bob_bot),
            color="#8f8caf", lw=1.2)
    # reference stripe (shows the twist)
    ax.plot([R_BOB * np.cos(theta)] * 2, [R_BOB * np.sin(theta)] * 2,
            [bob_bot, bob_top], color="#ff5c8a", lw=3)
    # rod + sliders
    rod = 1.55
    ax.plot([rod * np.cos(theta), -rod * np.cos(theta)],
            [rod * np.sin(theta), -rod * np.sin(theta)], [bob_c, bob_c],
            color="#59617e", lw=2.5)
    for s in (1, -1):
        ax.scatter([s * rod * np.cos(theta)], [s * rod * np.sin(theta)], [bob_c],
                   s=90, color="#b57cff", edgecolor="#2a1a4a", depthshade=False)
    ax.view_init(elev=16, azim=-60)


def simulate(pid, frames):
    p, d = w.from_preset(pid)
    u0 = [d["z0"], 0.0, d["theta0"], 0.0]
    tf = min(d["tf"], 2.2 * p.beat_period if np.isfinite(p.beat_period) else d["tf"])
    n = 4000
    t, U = w.integrate_rk4(p, u0, 0.0, tf / n, n)
    idx = np.linspace(0, n - 1, frames).astype(int)
    zamp = max(np.max(np.abs(U[:, 0])), 1e-3)
    return t[idx], U[idx], p, zamp


def render_gif(pid, out, frames, fps):
    from matplotlib.animation import FuncAnimation
    t, U, p, zamp = simulate(pid, frames)
    fig = plt.figure(figsize=(4.2, 6.0)); fig.patch.set_facecolor("#0b0e17")
    ax = fig.add_subplot(111, projection="3d"); ax.set_facecolor("#0b0e17")

    def upd(i):
        draw_frame(ax, U[i, 0], U[i, 2], zamp)
        ax.set_title(f"t = {t[i]:5.1f} s", color="#c9d2e8", fontsize=11)

    anim = FuncAnimation(fig, upd, frames=frames, interval=1000 / fps)
    if out.lower().endswith(".mp4"):
        try:
            anim.save(out, writer="ffmpeg", fps=fps, dpi=110,
                      savefig_kwargs={"facecolor": "#0b0e17"})
        except Exception as e:
            out = out[:-4] + ".gif"
            print(f"  ffmpeg unavailable ({e}); writing {out}")
            anim.save(out, writer="pillow", fps=fps)
    else:
        anim.save(out, writer="pillow", fps=fps)
    plt.close(fig)
    print(f"  wrote {out}  ({frames} frames @ {fps} fps)")


def render_filmstrip(pid, out, npanels=6):
    p, d = w.from_preset(pid)
    u0 = [d["z0"], 0.0, d["theta0"], 0.0]
    tf = min(d["tf"], 1.05 * p.beat_period if np.isfinite(p.beat_period) else d["tf"])
    n = 4000
    t, U = w.integrate_rk4(p, u0, 0.0, tf / n, n)
    zamp = max(np.max(np.abs(U[:, 0])), 1e-3)
    idx = np.linspace(0, n - 1, npanels).astype(int)
    fig = plt.figure(figsize=(2.15 * npanels, 3.6)); fig.patch.set_facecolor("white")
    for k, i in enumerate(idx):
        ax = fig.add_subplot(1, npanels, k + 1, projection="3d")
        ax.set_facecolor("white")
        draw_frame(ax, U[i, 0], U[i, 2], zamp)
        ax.set_title(f"t = {t[i]:.1f} s", fontsize=10)
    fig.suptitle(f"Wilberforce pendulum over one beat  (preset {pid})", y=0.98, fontsize=12)
    fig.tight_layout()
    for ext in ({os.path.splitext(out)[1].lstrip(".") or "pdf"} | {"png"}):
        f = os.path.splitext(out)[0] + "." + ext
        fig.savefig(f, dpi=200, facecolor="white", bbox_inches="tight")
        print(f"  wrote {f}")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--preset", type=int, default=14)
    ap.add_argument("--out", default=None)
    ap.add_argument("--filmstrip", default=None)
    ap.add_argument("--frames", type=int, default=160)
    ap.add_argument("--fps", type=int, default=25)
    a = ap.parse_args()
    if not a.out and not a.filmstrip:
        a.filmstrip = os.path.join(os.path.dirname(__file__), "..", "generated",
                                   "figures", "fig_pendulum_frames.pdf")
        os.makedirs(os.path.dirname(a.filmstrip), exist_ok=True)
    if a.out:
        render_gif(a.preset, a.out, a.frames, a.fps)
    if a.filmstrip:
        render_filmstrip(a.preset, a.filmstrip)


if __name__ == "__main__":
    main()
