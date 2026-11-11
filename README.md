# PCB Fan Placement: where should the fan go?

**FM116 Advanced Calculus · Plaksha University · Instructor: Dr. Suresh Kumar**

A printed circuit board has three heat-generating chips and one movable, localised fan. The steady temperature field is the solution of a heat equation. Calculus is then used to find the best fan position, to explain *why* it is best, and to check the whole simulation against exact results.

## The model in one line

Heat spreads sideways through the board and leaks into the air from both faces. In steady state the two balance at every point:

```
k·t·∇²T  −  2·h(x,y)·(T − Tₐ)  +  q(x,y)  =  0          insulated edges: ∇T·n = 0

h(x,y) = h_nat + Δh · exp( −[(x − x_f)² + (y − y_f)²] / 2σ² )
```

The fan is a Gaussian patch of extra convection centred at (x_f, y_f). Moving the fan changes only the coefficient *h*, so the question is how the solution depends on (x_f, y_f).

## Headline results

| | |
|---|---|
| Hottest point, no fan | **69.0 °C** (at the CPU) |
| Hottest point, fan at its best position, about (37.7, 50.3) mm | **50.4 °C** |
| Observed order of accuracy of the solver, against an exact solution | **2.00** |
| Relative error of the global energy balance (heat in = heat out) | **≤ 10⁻¹³** |
| Agreement of heat flux through loops with the Divergence theorem | **within 0.1 %** |

**Main finding.** The best position is *not* on the hottest chip. T_max is the maximum over several peaks, so at the optimum the CPU and regulator peaks are equal and their gradients point in opposite directions. That is a Lagrange (KKT) condition, and it is verified numerically: the weighted gradient sum is zero to 10⁻³ K/mm.

## What is on the site

| Page | Contents |
|---|---|
| **Overview** (`index.html`) | Headline results, the equation, the main finding |
| **Live lab** (`lab.html`) | The solver running in your browser: drag the fan, see the field, the gradient of T_max with respect to fan position and the heat-flow arrows update. A live optimiser draws its path to the best position. Scan every fan position and download the CSV. Switch the fan to a Taylor-polynomial model and watch the cooling drop. |
| **Theory** (`theory.html`) | Derivation of the model, the five-point stencil and its Taylor-series error, conservation and the Divergence theorem, the adjoint gradient, KKT optimality, constrained problems, and a table mapping every FM116 course outcome to where it is used |
| **Analysis** (`analysis.html`) | Validation, critical points (second-derivative test), the optimum, the rail and airflow-budget Lagrange problems, flux, circulation, Green and Divergence checks, curvature of a heat-flow line, and a live Taylor-series demo |
| **Code & data** (`downloads.html`) | The Python scripts, the 2091-position scan as CSV, and how to reproduce everything |

## The calculus in the project

- **Partial derivatives, gradients, directional derivatives**: the five-point stencil, the exact gradient of the peak temperature with respect to the fan position (one adjoint solve per chip), steepest-descent optimisation, heat flux −kt∇T
- **Extrema**: three maxima and two saddle points of the temperature field, classified with D = T_xx·T_yy − T_xy²
- **Lagrange multipliers**: a fan constrained to a rail, a fixed airflow budget 2πσ²Δh = C (a polar double integral), and the KKT condition for tied peaks; multipliers interpreted as shadow prices
- **Multiple integrals**: energy balance over the board, the fan's total effect in polar coordinates, areas of hot regions
- **Vector fields and theorems**: the flux is a gradient field (conservative); outward flux through loops against the Divergence theorem; isotherm area by Green's theorem
- **Vector-valued functions**: arc length and curvature of heat-flow lines, with torsion zero for planar curves
- **Taylor series**: the O(ε²) stencil error, and Taylor polynomials of the fan profile with error bounds
- **Differential equations**: a one-dimensional strip reduces to an ODE with an exact solution used to test the solver

## Repository layout

```
index.html  lab.html  theory.html  analysis.html  downloads.html   the pages
assets/
  style.css            shared styles (light and dark themes)
  lab.js               in-browser finite-volume solver, optimiser, scan
  taylor_demo.js       live Taylor-series plot
  figures/             figures produced by the Python scripts
  data/                fan_scan.csv, summary.csv
  code/                the Python scripts behind every number
build.py               regenerates the HTML pages
.nojekyll              tells GitHub Pages to serve the files as they are
```

## Run it locally

The site is static, so any file server works:

```bash
python -m http.server 8000
# open http://localhost:8000
```

The pages load KaTeX (for equations) and Google Fonts from a CDN, so an internet connection is needed for them to look right.

## Reproduce the numbers

Python 3.9+ with `numpy`, `scipy`, `matplotlib` and `sympy`:

```bash
pip install numpy scipy matplotlib sympy
cd assets/code
python simulation.py --test     # built-in self-test
python simulation.py 2          # scan every fan position, 2 mm steps (about 11 s)
python Exact_checks.py          # stencil, 1D fin and 2D convergence studies
python Gradient.py              # adjoint gradients, KKT tie condition, descent
python Lagrange.py              # rail and airflow-budget problems
python VectorCalculus.py        # flux, circulation, Divergence and Green checks
python Taylor.py                # Taylor analysis of the fan profile
```

Each script writes its figure and CSV files into a `results/` folder next to it.

## How the solver works

The board is divided into 2 mm square cells (50 × 40 = 2000 unknowns). A heat balance on each cell gives one linear equation, so the problem becomes **A·θ = b** with θ = T − Tₐ and **A = A_c + 2ε²·diag(h)**. The conduction matrix A_c is symmetric and does not depend on the fan, so only the diagonal changes when the fan moves; that is why a scan of 2091 fan positions takes seconds. The system is solved with a preconditioned conjugate gradient method in the browser and with a sparse direct solver in Python. Both use the same discretisation and agree with each other.

## Edit the site

Text and layout are generated by `build.py`. Set `AUTHORS` and `REPO_URL` at the top of that file, edit the page text, and run:

```bash
python build.py
```

If you change a model parameter in `simulation.py`, change it in `assets/lab.js` too, and rerun the scripts to refresh `assets/figures/` and `assets/data/`.

## Publish on GitHub Pages

This folder is the website root. Upload its contents directly to the top level of the GitHub repository. The repository root should contain `index.html`, `lab.html`, `theory.html`, `analysis.html`, `downloads.html`, `assets/`, and `.nojekyll`—do not upload the parent `Math_project` folder.

For a project site:

1. Create a GitHub repository named `PCB-Fan-Placement` under `RaghavBagga` (or use your preferred repository name).
2. Upload the contents of this folder to the repository root and commit them to the `main` branch.
3. Open **Settings → Pages**. Under **Build and deployment**, choose **Deploy from a branch**, select `main` and `/(root)`, then save.
4. GitHub will publish the site at `https://RaghavBagga.github.io/PCB-Fan-Placement/` (replace the final path with your repository name if it differs). The first deployment can take a few minutes.

For an account homepage, name the repository exactly `RaghavBagga.github.io`; the site URL will be `https://RaghavBagga.github.io/`.

The pages use relative paths, so they work both at the account root and under a project-site path. KaTeX and Google Fonts load from CDNs, so equations and custom fonts need an internet connection; the rest of the site and the live solver are served from this repository.

## Assumptions and limitations

- The chip powers and positions, the effective conductivity (15 W/m·K), the convection coefficients and the Gaussian fan profile are **assumed for illustration**. Nothing has been compared with measurements, so the checks show that the equations are solved correctly, not that the model describes a real board.
- The model is steady state and two-dimensional, with one effective conductivity instead of a layered board, insulated edges, constant natural convection and no radiation.
- A real fan produces a directional flow pattern; the Gaussian ignores this.
- At the 2 mm grid the peak temperature carries a discretisation error of about 0.2 °C, which does not change which fan position is best.
- The airflow-budget optimisation improves on the default fan by only 0.14 °C; the robust finding is the shape of the trade-off, not the exact optimum.

## Reference

Heat-transfer framework adapted from Zhang, Y. & Chen, L. (2024), *A modeling method for thermal steady-state simulation of the four-layer printed circuit board*, PLOS ONE 19(9), e0310237.

Course textbook: Thomas' Calculus, 15th edition.
