"""
wilberforce.py — an independent, double-precision reference implementation of
the Wilberforce-pendulum model used to cross-check the Fortran solver and to
generate publication figures.

The physics is summarized in ../README.md and ../paper/en.pdf. Symbols:

    z, theta      longitudinal and angular displacements
    k, delta, eps longitudinal / torsional / coupling constants
    m, I          mass and moment of inertia

Equations of motion (resonance k = m*delta/I enforced):

    z''    = -(k/m) z - (eps/2m) theta
    theta''= -(delta/I) theta - (eps/2I) z

This module provides:
    * PRESETS            the 14 built-in configurations (matching the Fortran)
    * Params             a small parameter container with derived frequencies
    * analytic(...)      the closed-form normal-mode solution
    * integrate_rk4(...) a fixed-step RK4 integrator (same scheme as Fortran)
    * integrate_ref(...) an independent adaptive solver (SciPy DOP853)
    * energy(...)        the energy budget
"""
from __future__ import annotations

import numpy as np
from dataclasses import dataclass


def _s(x: float) -> float:
    """Round a literal through IEEE single precision, reproducing the original
    Fortran quirk of assigning bare real literals (e.g. `2.2`) to double
    variables. Keeps this reference numerically aligned with the Fortran core."""
    return float(np.float32(x))


# --- the fourteen presets, mirroring src/wilberforce_presets.f90 --------------
# Physical constants pass through _s() exactly as the Fortran single->double
# promotion; the step h uses a genuine double literal.
_PI_S = _s(4.0 * np.arctan(1.0))  # single-precision pi, as 4*atan(1.0) in Fortran

PRESETS = {
    1:  dict(delta=_s(2.2),  eps=_s(0.05), m=_s(0.5),  I=_s(0.15), z0=_s(0.1),  theta0=_PI_S,        h=0.01,  tf=100.0, label="Set 1"),
    2:  dict(delta=_s(7.0),  eps=_s(1.2),  m=_s(10.0), I=_s(0.9),  z0=_s(0.01), theta0=_s(0.5*_PI_S),h=0.01,  tf=100.0, label="Set 2"),
    3:  dict(delta=_s(5.4),  eps=_s(0.3),  m=_s(0.4),  I=_s(0.2),  z0=_s(0.15), theta0=_s(0.75*_PI_S),h=0.01, tf=10.0,  label="Set 3"),
    4:  dict(delta=_s(0.4),  eps=_s(0.03), m=_s(0.1),  I=_s(0.02), z0=_s(0.1),  theta0=_s(0.25*_PI_S),h=0.01, tf=100.0, label="Set 4"),
    5:  dict(delta=_s(0.9),  eps=_s(0.6),  m=_s(1.0),  I=_s(0.1),  z0=_s(0.2),  theta0=0.0,          h=0.01,  tf=100.0, label="Set 5"),
    6:  dict(delta=_s(2.3),  eps=_s(1.6),  m=_s(0.1),  I=_s(0.8),  z0=_s(0.08), theta0=0.0,          h=0.01,  tf=100.0, label="eps=1.6 (strong)"),
    7:  dict(delta=_s(2.3),  eps=_s(0.8),  m=_s(0.1),  I=_s(0.8),  z0=_s(0.08), theta0=0.0,          h=0.01,  tf=100.0, label="eps=0.8"),
    8:  dict(delta=_s(2.3),  eps=_s(0.2),  m=_s(0.1),  I=_s(0.8),  z0=_s(0.08), theta0=0.0,          h=0.01,  tf=100.0, label="eps=0.2"),
    9:  dict(delta=_s(2.3),  eps=_s(0.05), m=_s(0.1),  I=_s(0.8),  z0=_s(0.08), theta0=0.0,          h=0.01,  tf=100.0, label="eps=0.05"),
    10: dict(delta=_s(2.3),  eps=_s(0.01), m=_s(0.1),  I=_s(0.8),  z0=_s(0.08), theta0=0.0,          h=0.01,  tf=100.0, label="eps=0.01 (weak)"),
    11: dict(delta=_s(1.6),  eps=_s(0.075),m=_s(0.07), I=_s(0.5),  z0=_s(1.0),  theta0=0.0,          h=0.01,  tf=100.0, label="non-periodic"),
    12: dict(delta=_s(1.6),  eps=_s(0.075),m=_s(0.07), I=_s(0.5),  z0=_s(1.0),  theta0=_s(0.374165738),h=0.01,tf=100.0, label="periodicity condition"),
    13: dict(delta=_s(1.6),  eps=_s(0.075),m=_s(0.07), I=_s(0.5),  z0=_s(1.0),  theta0=_s(0.374165738),h=0.001,tf=5.0, label="fine step"),
    14: dict(delta=_s(7.44288e-4), eps=_s(9.27e-3), m=_s(0.4905), I=_s(1.39e-4), z0=0.0, theta0=_s(2*_PI_S), h=0.01, tf=40.0, label="Berg-Marshall / Hill"),
}


@dataclass
class Params:
    delta: float
    eps: float
    m: float
    I: float
    z0: float = 0.0
    theta0: float = 0.0

    def validate(self):
        values = [self.delta, self.eps, self.m, self.I, self.z0, self.theta0]
        if not np.all(np.isfinite(values)) or min(self.delta, self.m, self.I) <= 0:
            raise ValueError("Parameters must be finite; delta, m and I must be positive")
        critical = 2 * self.delta * np.sqrt(self.m / self.I)
        if abs(self.eps) >= critical:
            raise ValueError("The oscillatory model requires eps^2 < 4*k*delta")

    @property
    def k(self) -> float:
        """Longitudinal constant fixed by the resonance condition k = m*delta/I."""
        return self.m * self.delta / self.I

    @property
    def w0(self) -> float:
        """Natural (uncoupled) angular frequency, sqrt(delta/I)."""
        return np.sqrt(self.delta / self.I)

    @property
    def w1(self) -> float:
        """Upper normal-mode frequency."""
        return np.sqrt(self.w0**2 + abs(self.eps) / (2.0 * np.sqrt(self.m * self.I)))

    @property
    def w2(self) -> float:
        """Lower normal-mode frequency."""
        return np.sqrt(self.w0**2 - abs(self.eps) / (2.0 * np.sqrt(self.m * self.I)))

    @property
    def beat_period(self) -> float:
        """Beat period T_beat = 2*pi/(w1 - w2)."""
        return np.inf if self.eps == 0 else 2.0 * np.pi / (self.w1 - self.w2)

    @property
    def mode_amplitudes(self):
        """(B, D) amplitudes of the two normal modes for release from rest."""
        self.validate()
        scale = (-1.0 if self.eps < 0 else 1.0) * np.sqrt(self.I / self.m)
        return ((self.theta0 + self.z0 / scale) / 2,
                (self.theta0 - self.z0 / scale) / 2)


def from_preset(pid: int) -> tuple[Params, dict]:
    d = PRESETS[pid]
    p = Params(delta=d["delta"], eps=d["eps"], m=d["m"], I=d["I"],
               z0=d["z0"], theta0=d["theta0"])
    return p, d


# --- closed-form solution -----------------------------------------------------
def analytic(p: Params, t: np.ndarray):
    """Return z, v, theta, omega from the normal-mode solution at times t."""
    p.validate()
    t = np.asarray(t, dtype=float)
    if not np.all(np.isfinite(t)):
        raise ValueError("Evaluation times must be finite")
    w1, w2 = p.w1, p.w2
    B, D = p.mode_amplitudes
    scale = (-1.0 if p.eps < 0 else 1.0) * np.sqrt(p.I / p.m)
    q1, q2 = B * np.cos(w1*t), D * np.cos(w2*t)
    v1, v2 = -w1*B*np.sin(w1*t), -w2*D*np.sin(w2*t)
    return scale*(q1-q2), scale*(v1-v2), q1+q2, v1+v2


# --- equations of motion & integrators ---------------------------------------
def rhs(u, p: Params):
    """du/dt for the state u = [z, v, theta, omega]:
        dz/dt = v,  dv/dt = z'',  dtheta/dt = omega,  domega/dt = theta''."""
    z, v, th, om = u
    az = -(p.k / p.m) * z - (p.eps / (2 * p.m)) * th
    ath = -(p.delta / p.I) * th - (p.eps / (2 * p.I)) * z
    return np.array([v, az, om, ath])


def integrate_rk4(p: Params, u0, t0: float, h: float, n: int):
    """Fixed-step classical RK4, the same scheme as the Fortran core.
    Returns (t, U) with U of shape (n+1, 4)."""
    p.validate()
    if not isinstance(n, (int, np.integer)) or n < 0:
        raise ValueError("n must be a nonnegative integer")
    if not np.isfinite(t0) or not np.isfinite(h) or h <= 0:
        raise ValueError("t0 and h must be finite; h must be positive")
    if h * p.w1 > np.sqrt(8):
        raise ValueError("RK4 timestep outside oscillatory stability region")
    u0 = np.asarray(u0, dtype=float)
    if u0.shape != (4,) or not np.all(np.isfinite(u0)):
        raise ValueError("u0 must have four finite components")
    t = t0 + h * np.arange(n + 1)
    U = np.empty((n + 1, 4))
    U[0] = u0
    u = np.array(u0, dtype=float)
    for i in range(n):
        k1 = rhs(u, p)
        k2 = rhs(u + 0.5 * h * k1, p)
        k3 = rhs(u + 0.5 * h * k2, p)
        k4 = rhs(u + h * k3, p)
        u = u + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        if not np.all(np.isfinite(u)):
            raise FloatingPointError("Non-finite integrated state")
        U[i + 1] = u
    return t, U


def integrate_ref(p: Params, u0, t_eval, t0: float = 0.0):
    """Independent high-accuracy solver (SciPy DOP853) for cross-validation.
    The initial condition u0 is applied at t0; t_eval must lie within
    [t0, t_eval[-1]]."""
    p.validate()
    t_eval = np.asarray(t_eval, dtype=float)
    u0 = np.asarray(u0, dtype=float)
    if (t_eval.ndim != 1 or not len(t_eval) or not np.isfinite(t0)
            or not np.all(np.isfinite(t_eval)) or t_eval[0] < t0
            or np.any(np.diff(t_eval) <= 0)):
        raise ValueError("t_eval must be finite, nonempty and strictly increasing from t0")
    if u0.shape != (4,) or not np.all(np.isfinite(u0)):
        raise ValueError("u0 must have four finite components")
    if t_eval[-1] == t0:
        return u0.reshape(1, 4).copy()
    from scipy.integrate import solve_ivp
    sol = solve_ivp(lambda t, u: rhs(u, p), (t0, t_eval[-1]), u0,
                    t_eval=t_eval, method="DOP853", rtol=1e-12, atol=1e-14)
    if not sol.success or sol.y.shape != (4, len(t_eval)):
        raise RuntimeError(f"DOP853 failed: {sol.message}")
    return sol.y.T


# --- energy -------------------------------------------------------------------
def energy(U, p: Params):
    """Energy contributions for a trajectory U of shape (N, 4).
    Returns a dict of arrays: kin_long, pot_long, kin_ang, pot_ang, coupling,
    long (=kin+pot longitudinal), ang, total."""
    z, v, th, om = U[:, 0], U[:, 1], U[:, 2], U[:, 3]
    kin_long = 0.5 * p.m * v**2
    pot_long = 0.5 * p.k * z**2
    kin_ang = 0.5 * p.I * om**2
    pot_ang = 0.5 * p.delta * th**2
    coupling = 0.5 * p.eps * z * th
    return dict(
        kin_long=kin_long, pot_long=pot_long, kin_ang=kin_ang, pot_ang=pot_ang,
        coupling=coupling, long=kin_long + pot_long, ang=kin_ang + pot_ang,
        total=kin_long + pot_long + kin_ang + pot_ang + coupling,
    )


def bob_path_3d(U, radius=1.0, z_scale=1.0):
    """Map (z, theta) to a 3D point on the rim of the bob:
    x = r cos(theta), y = r sin(theta), z = z. Useful for 3D trajectory plots."""
    z, th = U[:, 0], U[:, 2]
    x = radius * np.cos(th)
    y = radius * np.sin(th)
    return x, y, z_scale * z


if __name__ == "__main__":
    # Smoke test: integrate preset 14 and report energy conservation.
    p, d = from_preset(14)
    u0 = [d["z0"], 0.0, d["theta0"], 0.0]
    n = int((d["tf"] - 0.0) / d["h"])
    t, U = integrate_rk4(p, u0, 0.0, d["h"], n)
    E = energy(U, p)
    drift = (E["total"].max() - E["total"].min()) / abs(E["total"][0])
    print(f"preset 14: w0={p.w0:.5f} w1={p.w1:.5f} w2={p.w2:.5f} "
          f"T_beat={p.beat_period:.3f}s  energy drift={drift:.2e}")
