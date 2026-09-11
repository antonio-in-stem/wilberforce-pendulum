# Wilberforce pendulum

This repository contains a modular computational study of the Wilberforce pendulum: a helical spring couples the vertical translation and axial rotation of a suspended mass. The project solves the resonant two-coordinate model numerically, compares it with a closed-form normal-mode solution, studies energy exchange and coupling strength, and includes a standalone browser simulator.

The accompanying English paper is available as [`paper/en.pdf`](paper/en.pdf); the Spanish edition is [`paper/es.pdf`](paper/es.pdf).

## Model

The coordinates are the vertical displacement `z` and the angular displacement `theta`. The equations are

```text
m z''     + k z     + (epsilon/2) theta = 0
I theta'' + delta theta + (epsilon/2) z = 0
```

The supplied configurations enforce resonance, `k/m = delta/I`. At resonance, the two normal-mode frequencies split around the uncoupled frequency. Their superposition produces the characteristic beating and transfers energy between longitudinal and torsional motion.

The complete conserved energy includes the coupling term:

```text
E = 1/2 m z'^2 + 1/2 I theta'^2
  + 1/2 k z^2 + 1/2 delta theta^2 + 1/2 epsilon z theta
```

The model is linear, conservative and resonant. It does not include damping, external driving, nonlinear spring behavior, transverse motion or experimental calibration.

## Repository layout

```text
src/          Fortran numerical core: model, RK4, normal modes, energy and I/O
app/          Fortran command-line driver
analysis/     Independent Python reference, diagnostics and figure generators
web/          Dependency-free JavaScript model and the 14 exported configurations
test/         Fortran physics/regression tests and JavaScript model tests
simulator3d/  Standalone Canvas/WebGL simulator and animation renderer
paper/        English and Spanish papers
```

The Fortran package is organized around separate modules for precision and types, equations of motion, the analytical solution, integration, energy diagnostics, presets and output. The state is advanced with classical fourth-order Runge-Kutta; the Fortran implementation uses quadruple-precision state arithmetic.

## Build and run

With GNU Fortran and GNU Make:

```bash
make
./build/wilberforce 14
make test
```

The command-line driver accepts one of the 14 configurations and supports `--help`, `--outdir DIR` and `--plot`. Plotting requires an installed Gnuplot executable.

The project also provides a CMake build:

```bash
cmake -S . -B build-cmake -DCMAKE_Fortran_COMPILER=gfortran
cmake --build build-cmake
ctest --test-dir build-cmake --output-on-failure
```

Fortran output is written to `generated/`, which is ignored by Git. Build directories and compiler artifacts are also ignored.

## Tests and analysis

The Fortran test suite checks resonance, energy conservation, agreement with the analytical solution, the single-mode condition, zero and unstable coupling, initial velocities and fourth-order convergence.

Run the browser-model tests with Node.js:

```bash
node --test test/browser-model.test.mjs
```

The Python analysis layer uses NumPy, SciPy, Matplotlib and ImageIO:

```bash
python -m pip install -r analysis/requirements.txt
python analysis/measure_coupling.py
python analysis/make_figures.py
```

The English paper reports fourth-order convergence, agreement with the analytical solution at approximately `10^-6` relative error for the reference step, and relative total-energy drift of approximately `5.5 x 10^-8` over 4000 steps. It also demonstrates how preparing the initial amplitudes as `theta0 = z0 sqrt(m/I)` suppresses one normal mode and allows the coupling constant `epsilon` to be recovered from a measured period.

## Standalone simulator

The `simulator3d/` directory contains a self-contained browser simulator. Its source templates, shared model, presets, application code and pinned Three.js r128 files are included in the repository. `vendor/LICENSE-three.txt` preserves the Three.js MIT notice.

Generate the two standalone HTML files with:

```bash
python simulator3d/build.py
python -m http.server 8770 --directory simulator3d
```

Then open `http://localhost:8770/index_lite.html` for the Canvas version or `http://localhost:8770/index.html` for the WebGL version. The generated HTML files are build outputs and are intentionally excluded from Git; rebuild them when needed.

The separate renderer can produce a GIF or a static filmstrip using the same model:

```bash
python simulator3d/render_animation.py
python simulator3d/render_animation.py --preset 14 --out generated/animation.gif
```

## Scope and attribution

This project is a numerical and educational model of the Wilberforce pendulum. Its results are simulations, not a substitute for laboratory calibration or a claim about nonlinear behavior in a physical apparatus.

The papers credit José Antonio Martínez Torres, Aarón Marcelo Sandoval Vázquez and Diego Jorge Gómez of the Facultad de Ciencias Físico Matemáticas, Universidad Autónoma de Nuevo León. See [`LICENSE`](LICENSE) and [`CITATION.cff`](CITATION.cff) for the project license and citation metadata.
