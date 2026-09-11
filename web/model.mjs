/** Maintained browser adapter of the original resonant Wilberforce model.
 * State: [z (m), v (m/s), theta (rad), omega (rad/s)]. SI parameters.
 * Fortran is the reference; test/fixtures/fortran.json establishes parity.
 * No dependency, DOM, rendering, or wall-clock state belongs in this module.
 */
export function validate(p) {
  if (![p.m, p.I, p.delta, p.eps, p.z0 ?? 0, p.theta0 ?? 0].every(Number.isFinite)
      || Math.min(p.m, p.I, p.delta) <= 0) {
    throw new RangeError('Finite parameters and positive m, I, delta are required.');
  }
  const k = p.k ?? p.m * p.delta / p.I;
  const omega0Squared = p.delta / p.I;
  const splitting = Math.abs(p.eps) / (2 * Math.sqrt(p.m * p.I));
  if (![k, omega0Squared, splitting].every(Number.isFinite) || k <= 0
      || omega0Squared <= splitting
      || Math.abs(k / p.m - omega0Squared) > 64 * Number.EPSILON * omega0Squared) {
    throw new RangeError('The browser model requires stable resonance: eps² < 4 k delta and k/m = delta/I.');
  }
  return { k, omega0Squared, splitting };
}

function finiteState(u) {
  if (!Array.isArray(u) || u.length !== 4 || !u.every(Number.isFinite)) {
    throw new RangeError('State must contain four finite numbers.');
  }
}

export function normalModes(p) {
  const { omega0Squared, splitting } = validate(p);
  return { plus: Math.sqrt(omega0Squared + splitting),
           minus: Math.sqrt(omega0Squared - splitting), natural: Math.sqrt(omega0Squared) };
}

/** Elapsed time from the supplied initial state, not an absolute clock time. */
export function analyticState(p, elapsed, u0 = [p.z0 ?? 0, 0, p.theta0 ?? 0, 0]) {
  if (!Number.isFinite(elapsed)) throw new RangeError('Time must be finite.');
  finiteState(u0);
  const { plus, minus } = normalModes(p);
  // A signed mass scaling avoids division by eps and cancellation of close frequencies.
  const scale = (p.eps < 0 ? -1 : 1) * Math.sqrt(p.I / p.m);
  const B = (u0[2] + u0[0] / scale) / 2, D = (u0[2] - u0[0] / scale) / 2;
  const Bv = (u0[3] + u0[1] / scale) / 2, Dv = (u0[3] - u0[1] / scale) / 2;
  const q1 = B * Math.cos(plus * elapsed) + Bv * Math.sin(plus * elapsed) / plus;
  const q2 = D * Math.cos(minus * elapsed) + Dv * Math.sin(minus * elapsed) / minus;
  const v1 = -plus * B * Math.sin(plus * elapsed) + Bv * Math.cos(plus * elapsed);
  const v2 = -minus * D * Math.sin(minus * elapsed) + Dv * Math.cos(minus * elapsed);
  const result = [scale * (q1 - q2), scale * (v1 - v2), q1 + q2, v1 + v2];
  finiteState(result);
  return result;
}

export function derivative(p, u) {
  const k = p.k ?? p.m * p.delta / p.I;
  return [u[1], -k / p.m * u[0] - p.eps / (2 * p.m) * u[2],
          u[3], -p.delta / p.I * u[2] - p.eps / (2 * p.I) * u[0]];
}

export function rk4Step(p, u, h) {
  finiteState(u);
  const { plus } = normalModes(p);
  if (!Number.isFinite(h) || h <= 0 || h * plus > Math.sqrt(8)) {
    throw new RangeError('Positive finite h must satisfy h*omega_max <= sqrt(8).');
  }
  const add = (a, b, factor) => a.map((x, i) => x + factor * b[i]);
  const a = derivative(p, u), b = derivative(p, add(u, a, h / 2));
  const c = derivative(p, add(u, b, h / 2)), d = derivative(p, add(u, c, h));
  const result = u.map((x, i) => x + h / 6 * (a[i] + 2 * b[i] + 2 * c[i] + d[i]));
  finiteState(result);
  return result;
}

export function energyOf(p, u) {
  const { k } = validate(p); finiteState(u);
  const verticalEnergy = p.m * u[1] ** 2 / 2 + k * u[0] ** 2 / 2;
  const twistingEnergy = p.I * u[3] ** 2 / 2 + p.delta * u[2] ** 2 / 2;
  const couplingEnergy = p.eps * u[0] * u[2] / 2;
  const totalEnergy = verticalEnergy + twistingEnergy + couplingEnergy;
  if (!Number.isFinite(totalEnergy)) throw new RangeError('Energy overflow.');
  return { verticalEnergy, twistingEnergy, couplingEnergy, totalEnergy };
}

/** Fixed-step accumulator: zero speed never advances; sub-step time is retained.
 * Work per frame is bounded. A capped backlog is intentionally discarded to avoid
 * an unresponsive browser, rather than claiming real-time operation under overload.
 */
export function advanceClock(remainder, elapsed, speed, h, maxSteps = 400) {
  if (![remainder, elapsed, speed, h].every(Number.isFinite)
      || Math.min(remainder, elapsed, speed) < 0 || h <= 0
      || !Number.isInteger(maxSteps) || maxSteps <= 0) throw new RangeError('Invalid clock arguments.');
  if (speed === 0) return { steps: 0, remainder };
  const available = remainder + elapsed * speed;
  const wanted = Math.floor((available + h * 1e-10) / h);
  const steps = Math.min(maxSteps, wanted);
  return { steps, remainder: wanted > maxSteps ? 0 : Math.max(0, available - steps * h) };
}
