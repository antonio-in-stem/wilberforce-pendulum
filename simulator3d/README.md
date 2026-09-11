## Running

```bash
python simulator3d/build.py
python -m http.server 8770 --directory simulator3d
```

Open `http://localhost:8770/index_lite.html` for Canvas 2D, or `http://localhost:8770/index.html` for Three.js/WebGL. The HTML files are generated from the templates, the supplied code, and the shared `web/` model; the build does not require network access. The generated HTML files are not distributed as if they were separate source files.

| Source | Purpose |
|---|---|
| `app.js` | Three.js scene, interaction, and adapter for the shared model |
| `index.template.html` | Original structure and styles for the WebGL interface |
| `index_lite.template.html` | Source for the Canvas alternative, formerly `index_lite.html` |
| `build.py` | Inlines the model/presets and vendors; generates both HTML files |
| `vendor/` | Supplied Three.js r128 and OrbitControls, with no silent upgrade |
| `render_animation.py` | Separate Matplotlib renderer using the Python reference implementation |

The following issues were corrected: loss of the small inertia value in preset 14, the missing preset 13, initial state/energy mixing when parameters change, and the minimum per-frame advance that made the simulation progress even at zero speed. An accumulator maintains the fixed step without rounding it independently on every frame. Under load, the amount of work is capped and accumulated lag is discarded; real-time performance is not promised under all loads. The pure-mode button does not clip the initial condition to fit a slider artificially.

The Canvas version was run in Chromium for all 14 presets, parameter changes, and controls. The WebGL version was inspected, built, and syntax-checked; **its graphical execution was not validated**, because the audit browser could not create a WebGL context. Its `preview.png` and `wilberforce_demo.gif` assets are historical, not evidence from this execution.

## Python rendering

```bash
python simulator3d/render_animation.py
python simulator3d/render_animation.py --preset 14 --out generated/animation.gif
```

The first command generates the frame strip in `generated/figures/`. It does not overwrite the historical figure. Matplotlib does not use the Three.js renderer: it preserves the same physics, not the same pixels.

The simulator is an educational companion to the numerical model described in paper/en.pdf.
