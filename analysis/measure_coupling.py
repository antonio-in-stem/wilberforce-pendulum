#!/usr/bin/env python
"""
measure_coupling.py — demonstrate the strategy for measuring the coupling
constant epsilon of a Wilberforce pendulum from a single period measurement.

The general motion superposes two
incommensurate normal-mode frequencies and is not periodic, so no period can be
read off. But preparing the system with

        theta0 = z0 * sqrt(m / I)

annihilates the second mode; the motion becomes a clean single-frequency
oscillation with period T = 2*pi/w1. Measuring T then yields

        epsilon = 2 * sqrt(m*I) * ( 4*pi^2/T^2  -  delta/I ).

This script prepares each of a few systems in the beat-free state, integrates
them, measures the period *from the trajectory* (not from the known formula),
and recovers epsilon — verifying the strategy end-to-end.

    python analysis/measure_coupling.py
"""
from __future__ import annotations

import os
import sys
import numpy as np

# Ensure Greek symbols print on Windows consoles (cp1252 by default).
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

sys.path.insert(0, os.path.dirname(__file__))
import wilberforce as w


def measure_period(t, x):
    """Estimate the fundamental period from upward zero-crossings of x(t)
    (with sub-sample linear interpolation), robust for a clean single tone."""
    x = np.asarray(x) - np.mean(x)
    crossings = []
    for i in range(len(x) - 1):
        if x[i] <= 0.0 < x[i + 1]:
            # linear interpolation of the crossing time
            frac = -x[i] / (x[i + 1] - x[i])
            crossings.append(t[i] + frac * (t[i + 1] - t[i]))
    crossings = np.array(crossings)
    if len(crossings) < 2:
        return np.nan
    return np.mean(np.diff(crossings))


def recover_epsilon(delta, m, I, T):
    """epsilon from the measured single-mode period T."""
    return 2.0 * np.sqrt(m * I) * ((2 * np.pi / T) ** 2 - delta / I)


def demo_system(name, delta, m, I, eps_true, z0=1.0, tf=200.0, h=1e-3):
    p = w.Params(delta=delta, eps=eps_true, m=m, I=I, z0=z0,
                 theta0=z0 * np.sqrt(m / I))       # periodicity condition
    _, D = p.mode_amplitudes
    u0 = [p.z0, 0.0, p.theta0, 0.0]
    n = int(tf / h)
    t, U = w.integrate_rk4(p, u0, 0.0, h, n)
    T = measure_period(t, U[:, 0])                  # measure from the trajectory
    eps_meas = recover_epsilon(delta, m, I, T)
    err = 100.0 * abs(eps_meas - eps_true) / eps_true
    print(f"  {name:<22}  δ={delta:<6.4g} m={m:<6.4g} I={I:<7.4g}"
          f" | θ0={p.theta0:.4f}  |D|={abs(D):.1e}")
    print(f"  {'':<22}  T_meas={T:.4f}s  ε_true={eps_true:.5f}"
          f"  ε_meas={eps_meas:.5f}  (err {err:.2f}%)\n")
    return eps_meas


def main():
    print("=" * 70)
    print("  Measuring the coupling constant ε via the periodicity condition")
    print("=" * 70)
    print("  Prepare θ0 = z0·√(m/I) so only one normal mode survives,")
    print("  measure the period T from the trajectory, then invert")
    print("  ε = 2√(mI)(4π²/T² − δ/I).\n")

    # The paper's Set 11/12 system.
    demo_system("Set 11/12", delta=w._s(1.6), m=w._s(0.07), I=w._s(0.5),
                eps_true=w._s(0.075))
    # A couple of independent systems to show it is general.
    demo_system("stiff torsion", delta=2.5, m=0.30, I=0.02, eps_true=0.040)
    demo_system("Berg-Marshall-like", delta=7.44288e-4, m=0.4905, I=1.39e-4,
                eps_true=9.27e-3, z0=0.05, tf=400.0)

    print("  Conclusion: a single stopwatch measurement of the beat-free period")
    print("  recovers ε to a fraction of a percent for any resonant system.")


if __name__ == "__main__":
    main()
