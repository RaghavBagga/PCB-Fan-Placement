#!/usr/bin/env python3
"""Generate the static site (index, lab, theory, analysis, downloads). Run:  python build.py"""
from pathlib import Path
ROOT = Path(__file__).parent

# ---- edit these ----------------------------------------------------------------
COURSE = "FM116 Advanced Calculus · Plaksha University"
INSTRUCTOR = "Dr. Suresh Kumar"
AUTHORS = "Raghav Bagga (model, equation, simulation and analysis)"
REPO_URL = "https://github.com/RaghavBagga/PCB-Fan-Placement"  # adds a link on the downloads page
# ---------------------------------------------------------------------------------

PAGES = [("index.html", "Overview"), ("lab.html", "Live lab"), ("theory.html", "Theory"),
         ("analysis.html", "Analysis"), ("downloads.html", "Code & data")]

HEAD = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1,viewport-fit=cover">
<title>{title}</title>
<meta name="description" content="{desc}">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Bricolage+Grotesque:opsz,wght@12..96,600;12..96,700&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap">
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.16.9/katex.min.css">
<link rel="stylesheet" href="assets/style.css">
</head><body>
<a class="skip-link" href="#main">Skip to content</a>
<nav class="top" aria-label="Main navigation"><div class="in">
<a class="brand" href="index.html">PCB Fan Placement</a>
{links}
<button id="themeBtn" aria-label="Toggle colour theme">Theme</button>
</div></nav>
<main class="wrap" id="main">
"""

FOOT = """
</main>
<footer>
{course}. Instructor: {instr}. {authors}.<br>
An educational steady-state approximation, not a validated CFD or product-design tool. All parameters and the Gaussian fan profile are assumptions; nothing has been compared with measurements.
Heat-transfer framework adapted from Zhang, Y. &amp; Chen, L. (2024), <i>A modeling method for thermal steady-state simulation of the four-layer printed circuit board</i>, PLOS ONE 19(9), e0310237.
</footer>
<script src="https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.16.9/katex.min.js"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.16.9/contrib/auto-render.min.js"></script>
<script>
(function(){
  var root=document.documentElement;
  try{var t=localStorage.getItem("theme");if(t)root.setAttribute("data-theme",t);}catch(e){}
  document.getElementById("themeBtn").addEventListener("click",function(){
    var dark=root.getAttribute("data-theme")==="dark"||(!root.getAttribute("data-theme")&&window.matchMedia&&matchMedia("(prefers-color-scheme: dark)").matches);
    var n=dark?"light":"dark";root.setAttribute("data-theme",n);try{localStorage.setItem("theme",n);}catch(e){}
    document.dispatchEvent(new Event("themechange"));
  });
  if(window.renderMathInElement)renderMathInElement(document.body,{delimiters:[{left:"$$",right:"$$",display:true},{left:"$",right:"$",display:false}],throwOnError:false});
})();
</script>
{scripts}
</body></html>
"""


def expand(body):
    import re
    body = re.sub(r"@@FIG:([^|]+)\|([^@]+)@@", lambda m: fig(m.group(1), m.group(2)), body)
    for k, v in EXTRA.items():
        body = body.replace("@@" + k + "@@", v)
    return body


EXTRA = {}


def page(fname, title, desc, body, scripts=""):
    body = expand(body)
    cur = ' aria-current="page"'
    links = "\n".join(
        f'<a class="l" href="{f}"{cur if f == fname else ""}>{n}</a>' for f, n in PAGES)
    def sub(t, **kw):
        for k, v in kw.items():
            t = t.replace("{" + k + "}", v)
        return t
    html = sub(HEAD, title=title, desc=desc, links=links) + body + sub(
        FOOT, course=COURSE, instr=INSTRUCTOR, authors=AUTHORS, scripts=scripts)
    (ROOT / fname).write_text(html, encoding="utf-8")
    print("wrote", fname, len(html) // 1024, "kB")


def fig(src, cap):
    return f'<figure><img src="assets/figures/{src}" alt="{cap}" loading="lazy"><figcaption>{cap}</figcaption></figure>'


# ====================================================================== index
index = """
<header class="hero">
<div class="eyebrow">@@COURSE@@</div>
<h1>Where should the fan go?</h1>
<p>A printed circuit board with three heat-generating chips and one movable, localised fan. The temperature field is the solution of a steady heat equation. Calculus then finds the best fan position, explains <em>why</em> it is best, and checks the whole simulation against exact results.</p>
<div class="btns"><a class="btn" style="background:var(--accent);color:var(--accent-ink);border-color:var(--accent);font-weight:600" href="lab.html">Open the live lab</a><a class="btn" href="analysis.html">Read the analysis</a></div>
</header>

<section class="tiles" aria-label="Headline results">
<div class="tile"><span class="n">69.0<small>°C</small></span><span class="t">Hottest point with no fan (at the CPU)</span></div>
<div class="tile"><span class="n">50.4<small>°C</small></span><span class="t">Hottest point with the fan at its best position, about (37.7, 50.3) mm</span></div>
<div class="tile"><span class="n">2.00</span><span class="t">Measured order of accuracy of the solver, against an exact solution</span></div>
<div class="tile"><span class="n">10<sup style="font-size:.5em">−13</sup></span><span class="t">Relative error of the global energy balance (heat in = heat out)</span></div>
</section>

<section class="grid2">
<div class="card prose">
<h2>The model in one line</h2>
<p>Heat spreads sideways through the board and leaks into the air from both faces. In steady state the two balance at every point:</p>
<div class="eq">$$kt\\,\\nabla^2 T\\;-\\;2h(x,y)\\,(T-T_a)\\;+\\;q(x,y)=0,\\qquad h=h_{\\rm nat}+\\Delta h\\,e^{-r^2/2\\sigma^2}$$</div>
<p>The fan is a Gaussian patch of extra convection centred at $(x_f,y_f)$. Moving it changes only the coefficient $h$, so the question is how the solution $T$ depends on $(x_f,y_f)$.</p>
</div>
<div class="card prose">
<h2>The main finding</h2>
<p>The best position is <strong>not</strong> on the hottest chip. $T_{\\max}$ is the maximum over several peaks, so at the optimum the CPU and regulator peaks are equal and their gradients point in opposite directions. That is a Lagrange (KKT) condition, and the analysis page verifies it numerically: the weighted gradient sum is zero to $10^{-3}$ K/mm.</p>
<p>Along the way the project uses every strand of the course: gradients and directional derivatives, second-derivative tests, Lagrange multipliers, polar double integrals, flux, circulation, the Divergence and Green theorems, curvature, Taylor series and an ODE with an exact solution.</p>
</div>
</section>

<section>
<h2>Explore</h2>
<div class="grid3">
<a class="navcard" href="lab.html"><b>Live lab</b><span>Drag the fan, watch the field and the gradient of T<sub>max</sub> update. Run an optimiser, scan all positions, switch to a Taylor-polynomial fan, show heat-flow arrows.</span></a>
<a class="navcard" href="theory.html"><b>Theory</b><span>The equation, the discretisation, the adjoint gradient and the KKT condition, with a map from each course outcome to where it is used.</span></a>
<a class="navcard" href="analysis.html"><b>Analysis</b><span>Validation, critical points, the optimum, two Lagrange problems, the Divergence theorem checks, and a live Taylor-series demo.</span></a>
<a class="navcard" href="downloads.html"><b>Code &amp; data</b><span>The Python scripts that produced every figure and number, the 2091-position scan as CSV, and how to reproduce it all.</span></a>
</div>
</section>

<div class="callout honest prose"><strong>What this is and is not.</strong> Everything on this site is computed from the stated model. The chip powers, conductivity, convection coefficients and the Gaussian fan shape are assumptions chosen for illustration, and nothing was validated against a physical board. The checks show that the equations are <em>solved correctly</em>, not that the model describes a real product.</div>
"""
page("index.html", "PCB Fan Placement", "Calculus-based thermal model of a PCB with a movable fan: live simulation, optimiser, theory and full analysis.", index.replace("@@COURSE@@", COURSE))

# ====================================================================== lab
lab = """
<header class="hero">
<div class="eyebrow">Live simulation · runs in your browser</div>
<h1>PCB Fan Placement Lab</h1>
<p>Drag the fan across the board. The map is the numerical solution of the 2D heat equation (2 mm cells, 2000 unknowns), recomputed on every move with the same method as the Python code.</p>
</header>

<section class="bench">
<div class="card">
<div class="board-head"><h2>Board temperature</h2><span class="hint">Click or drag to move the fan. Arrow keys also work.</span></div>
<div class="board" id="board">
<canvas id="heat"></canvas><div id="chips"></div>
<div class="fan" id="fan" tabindex="0" role="slider" aria-label="Fan position, use arrow keys to move"><div class="core"><svg width="20" height="20" viewBox="-10 -10 20 20" aria-hidden="true"><g fill="#6fe3ff"><ellipse cx="0" cy="-5" rx="2.6" ry="5"/><ellipse cx="0" cy="5" rx="2.6" ry="5"/><ellipse cx="-5" cy="0" rx="5" ry="2.6"/><ellipse cx="5" cy="0" rx="5" ry="2.6"/></g><circle r="2" fill="#06303a"/></svg></div></div>
</div>
<div class="legend"><span id="lo">25</span><div class="bar" id="bar"></div><span id="hi">70 °C</span></div>
<div class="legend"><div class="key">
<span><i></i>Contour at temperature limit</span><span><b></b>Hottest point</span><span><i class="dash"></i>Fan radius (σ)</span><span><i class="arr"></i>Steepest descent of T<sub>max</sub> / heat flow</span></div></div>
</div>

<div class="side">
<div class="card">
<h2>Readout</h2>
<div class="readout">
<div class="stat"><span class="l">Max temperature</span><span class="v" id="mTmax">–</span><span class="d" id="dTmax"></span></div>
<div class="stat"><span class="l">Area above limit</span><span class="v" id="mArea">–</span><span class="d" id="dArea"></span></div>
<div class="stat"><span class="l">Hottest point</span><span class="v" id="mHot" style="font-size:21px">–</span><span class="d" id="mHotWho"></span></div>
<div class="stat"><span class="l">Mean board temp</span><span class="v" id="mMean" style="font-size:21px">–</span><span class="d" id="dMean"></span></div>
</div>
<div class="balance" id="balance"></div>
</div>

<div class="card ctl">
<h2>Model controls</h2>
<label class="toggle"><input type="checkbox" id="fanOn" checked> Fan running</label>
<label class="toggle"><input type="checkbox" id="fluxOn"> Show heat-flow arrows (F = −kt∇T)</label>
<div class="row"><label for="model">Fan profile</label><span></span>
<select id="model" style="grid-column:1/-1"><option value="gauss">Gaussian (exact)</option><option value="1">Taylor, order 2 in r (paraboloid)</option><option value="2">Taylor, order 4 in r</option><option value="3">Taylor, order 6 in r</option><option value="4">Taylor, order 8 in r</option></select></div>
<div class="row"><label for="dh">Fan strength, Δh</label><output id="dhO"></output><input type="range" id="dh" min="0" max="150" step="5" value="60"></div>
<div class="row"><label for="sig">Fan radius, σ</label><output id="sigO"></output><input type="range" id="sig" min="4" max="30" step="1" value="12"></div>
<div class="row"><label for="lim">Temperature limit</label><output id="limO"></output><input type="range" id="lim" min="35" max="85" step="1" value="60"></div>
<div class="sub">Component power</div>
<div class="row"><label for="p0">CPU</label><output id="p0O"></output><input type="range" id="p0" min="0" max="4" step="0.1" value="2"></div>
<div class="row"><label for="p1">Regulator</label><output id="p1O"></output><input type="range" id="p1" min="0" max="4" step="0.1" value="1"></div>
<div class="row"><label for="p2">Driver</label><output id="p2O"></output><input type="range" id="p2" min="0" max="4" step="0.1" value="0.8"></div>
<div class="btns"><button id="reset">Reset all</button></div>
</div>
</div>
</section>

<section class="card">
<h2>Find the best position</h2>
<p class="hint">The optimiser is a compass search on T<sub>max</sub>(x<sub>f</sub>, y<sub>f</sub>): it tries steps in eight directions, moves when T<sub>max</sub> falls and halves the step otherwise. It converges to a point where two chip peaks are equal, which is exactly the tie condition of the theory page. The dashed cyan path on the board is its trail.</p>
<div class="btns"><button class="primary" id="optBtn">Optimise fan position</button></div>
<div class="scanread" id="optOut" style="margin-top:10px"></div>
</section>

<section class="card">
<h2>Scan every fan position</h2>
<p class="hint" style="margin:0 0 14px">Moves the fan across the whole board in 4 mm steps, solves the equation at each stop, and maps the maximum temperature. Click the map to place the fan there. (The Python scan in the downloads uses 2 mm steps.)</p>
<div class="scan">
<div><div class="scanbox" id="scanbox"><canvas id="scanCv"></canvas><div class="empty" id="scanEmpty">Run a scan to see how the maximum temperature depends on where the fan sits.</div></div>
<div class="scanread" id="scanRead"></div><div class="scanread" id="scanLegend"></div></div>
<div class="best">
<div class="btns"><button class="primary" id="scanBtn">Scan whole board</button><button id="goBest" disabled>Move fan to best</button><button id="csvBtn" disabled>Download CSV</button></div>
<progress id="prog" value="0" max="1" hidden></progress>
<div class="objective" role="radiogroup" aria-label="Objective">
<label><input type="radio" name="obj" value="tmax" checked> Lowest maximum temperature</label>
<label><input type="radio" name="obj" value="area"> Smallest area above the limit</label></div>
<div id="bestOut" class="hint">No scan yet.</div>
</div></div>
</section>

<section class="card prose">
<h2>Things to try</h2>
<ul>
<li><strong>Park the fan on the CPU.</strong> The CPU gets cooler but the regulator becomes the hottest part. The yellow gradient arrow points toward the regulator side, showing why the best spot is a compromise.</li>
<li><strong>Press Optimise.</strong> It stops near (37.7, 50.3) mm for the default fan, where the readout's hottest point jumps between the CPU and the regulator as you nudge the fan.</li>
<li><strong>Switch the fan profile to the order-2 Taylor model.</strong> The paraboloid is the local Taylor model of the fan (each polynomial is used only where it is within 2% of the Gaussian, and set to zero beyond). It recovers well under half of the cooling because the real fan reaches far beyond that region. Higher orders recover more. Expect about 61 °C instead of 50 °C at the optimum.</li>
<li><strong>Show heat-flow arrows.</strong> They are the flux field, which is a gradient field, so it has no circulation. Arrows run from each chip to the fan and toward the edges.</li>
<li><strong>Raise the limit slider to 69 °C</strong> and the no-fan hot region disappears. Lower it to 45 °C and watch the area metric.</li>
</ul>
</section>
"""
page("lab.html", "PCB Fan Placement Lab", "Interactive 2D heat-equation solver with a movable fan, live optimiser, board scan and gradient display.", lab,
     '<script src="assets/lab.js"></script>')

# ====================================================================== theory
rows = [
    ("Surfaces and coordinates", "The temperature field $T(x,y)$ is a surface over the board; the fan profile is written in polar coordinates about its centre.", "Lab, Theory"),
    ("Partial derivatives", "The five-point stencil is a combination of second partial differences; $\\partial T/\\partial x_f$ is the sensitivity of the field to the fan.", "Theory §2, §4"),
    ("Gradient, directional derivative", "Gradient of $T_{\\max}$ with respect to fan position (exact, by an adjoint solve); descent direction; heat flux $-kt\\nabla T$.", "Lab, Analysis §3"),
    ("Extrema, second-derivative test", "Maxima at chip hot spots, saddles in between, classified with $D=T_{xx}T_{yy}-T_{xy}^2$.", "Analysis §2"),
    ("Lagrange multipliers", "Fan on a rail; fixed airflow budget $2\\pi\\sigma^2\\Delta h=C$; the KKT tie condition; multipliers as shadow prices.", "Analysis §4"),
    ("Multiple integrals", "Energy balance over the board; $\\iint \\Delta h\\,e^{-r^2/2\\sigma^2}dA=2\\pi\\sigma^2\\Delta h$ in polar coordinates; areas of hot regions.", "Theory §3, Analysis §4"),
    ("Vector-valued functions", "Heat-flow lines $\\mathbf r(s)$: arc length, curvature $\\kappa$ from the Hessian, torsion $=0$ for a planar curve.", "Analysis §5"),
    ("Vector fields, conservative fields", "$\\mathbf F=-kt\\nabla T$ is a gradient field, so it is conservative with zero curl.", "Analysis §5"),
    ("Line integrals, flux, circulation", "Outward flux through loops around the chips; zero circulation of $\\mathbf F$.", "Analysis §5"),
    ("Green, Stokes and Divergence theorems", "$\\oint\\mathbf F\\cdot\\mathbf n\\,ds=\\iint(q-2h\\theta)\\,dA$ checked on six loops; area $=\\tfrac12\\oint(x\\,dy-y\\,dx)$ of an isotherm.", "Analysis §5"),
    ("Taylor series", "Stencil error $\\frac{\\varepsilon^2}{12}(u_{xxxx}+u_{yyyy})$; Taylor polynomials of the fan profile and their error bounds.", "Theory §2, Analysis §6"),
    ("Introductory ODEs", "The 1D strip reduces to $\\theta''-\\theta/\\ell^2=-(q'/kt)\\delta(x-x_0)$ with an exact solution used to test the solver.", "Theory §1, Analysis §1"),
]
trows = "\n".join(f"<tr><td>{a}</td><td>{b}</td><td class='m'>{c}</td></tr>" for a, b, c in rows)

EXTRA["TROWS"] = trows
theory = """
<header class="hero">
<div class="eyebrow">Theory</div>
<h1>From physics to calculus</h1>
<p>How the heat equation arises, how it becomes a linear system, and how derivatives of the solution with respect to the fan position are obtained exactly.</p>
</header>

<section class="card prose">
<h2>1. The heat equation and its length scale</h2>
<p>For a thin board of thickness $t$ and effective in-plane conductivity $k$, the net heat flowing into a patch by conduction, $kt\\nabla^2T$, balances loss through both faces, $2h(T-T_a)$, plus the generation $q$:</p>
$$kt\\,\\nabla^2T-2h(x,y)(T-T_a)+q(x,y)=0,\\qquad \\nabla T\\cdot\\mathbf n=0\\ \\text{on the edges.}$$
<p>Away from sources and with constant $h$ this becomes $\\nabla^2\\theta=\\theta/\\ell^2$ for $\\theta=T-T_a$, with</p>
$$\\ell=\\sqrt{\\frac{kt}{2h}}=\\sqrt{\\frac{0.024}{20}}\\ \\text{m}\\approx 34.6\\ \\text{mm}.$$
<p>Heat from a source decays like $e^{-r/\\ell}$, so a chip warms material within roughly $\\ell$ of it. In one dimension the equation is the ODE $\\theta''-\\theta/\\ell^2=-(q'/kt)\\,\\delta(x-x_0)$, whose exact solution (with insulated ends) is used as a test of the solver.</p>
</section>

<section class="card prose">
<h2>2. Discretisation and the Taylor-series error</h2>
<p>On a square grid of spacing $\\varepsilon$, adding the Taylor expansions of $u(x\\pm\\varepsilon,y)$ and $u(x,y\\pm\\varepsilon)$ gives</p>
$$\\frac{u_E+u_W+u_N+u_S-4u}{\\varepsilon^2}=\\nabla^2u+\\frac{\\varepsilon^2}{12}\\bigl(u_{xxxx}+u_{yyyy}\\bigr)+O(\\varepsilon^4),$$
<p>so the Laplacian is approximated with error $O(\\varepsilon^2)$. Writing a heat balance for each square cell gives one equation per cell,</p>
$$kt\\sum_{\\rm nbrs}(\\theta_{\\rm nb}-\\theta_c)-2h_c\\varepsilon^2\\theta_c+q_c\\varepsilon^2=0,$$
<p>that is the linear system $A\\boldsymbol\\theta=\\mathbf b$ with $A=A_c+2\\varepsilon^2\\,\\mathrm{diag}(h)$. The conduction matrix $A_c$ is symmetric and independent of the fan, so $A$ is symmetric positive definite and <em>only its diagonal changes when the fan moves</em>. That is why the scan of 2091 positions takes about ten seconds in Python, and why the lab responds instantly. The solver is a preconditioned conjugate gradient method.</p>
</section>

<section class="card prose">
<h2>3. Conservation and the Divergence theorem</h2>
<p>Integrating the equation over the board, the conduction term becomes a boundary integral by the Divergence theorem, and vanishes because the edges are insulated:</p>
$$\\sum_kP_k=\\iint_\\Omega 2h(T-T_a)\\,dA.$$
<p>Heat generated equals heat removed, for every fan position. The local form is $\\operatorname{div}\\mathbf F=q-2h\\theta$ with flux $\\mathbf F=-kt\\nabla T$. The readout in the lab prints both sides of the balance on every update.</p>
<p>The fan's total added convection is a polar integral:</p>
$$\\iint\\Delta h\\,e^{-r^2/2\\sigma^2}dA=\\Delta h\\int_0^{2\\pi}\\!\\!\\int_0^\\infty e^{-r^2/2\\sigma^2}r\\,dr\\,d\\phi=2\\pi\\sigma^2\\Delta h.$$
<p>This becomes the natural "airflow budget" constraint in the Lagrange problem.</p>
</section>

<section class="card prose">
<h2>4. The gradient with respect to the fan position, by an adjoint solve</h2>
<p>Only $h$ depends on the fan position. Differentiating $A\\boldsymbol\\theta=\\mathbf b$ with respect to $x_f$:</p>
$$A\\,\\frac{\\partial\\boldsymbol\\theta}{\\partial x_f}=-\\frac{\\partial A}{\\partial x_f}\\boldsymbol\\theta,\\qquad \\frac{\\partial h}{\\partial x_f}=\\Delta h\\,e^{-r^2/2\\sigma^2}\\frac{x-x_f}{\\sigma^2}.$$
<p>For the temperature in one cell $i$, define $\\boldsymbol\\lambda=A^{-1}\\mathbf e_i$ (one extra solve, using symmetry of $A$). Then</p>
$$\\frac{\\partial T_i}{\\partial x_f}=-\\sum_m\\lambda_m\\,2\\varepsilon^2\\,\\frac{\\partial h_m}{\\partial x_f}\\,\\theta_m,$$
<p>and the same $\\boldsymbol\\lambda$ gives every other parameter derivative. One solve per chip peak yields the whole gradient, instead of two solves per parameter by finite differences. The result agrees with central differences to $1.6\\times10^{-7}$ (relative, worst of 30 partial derivatives).</p>
</section>

<section class="card prose">
<h2>5. Optimality when peaks tie: KKT</h2>
<p>$T_{\\max}=\\max_kf_k$ is smooth except where peaks tie. At a tie of two peaks, the directional derivative of $T_{\\max}$ along $\\mathbf d$ is $\\max(\\nabla f_1\\!\\cdot\\mathbf d,\\nabla f_2\\!\\cdot\\mathbf d)$, so the point is a local minimum exactly when some weight $w\\in[0,1]$ satisfies</p>
$$w\\,\\nabla f_1+(1-w)\\,\\nabla f_2=\\mathbf 0.$$
<p>In the plane this means the two gradients are opposite. This is the Lagrange condition for "minimise $s$ subject to $f_1\\le s,\\ f_2\\le s$", with $w$ as multiplier. Gradient descent uses the minimum-norm vector in the convex hull of the tied gradients.</p>
</section>

<section class="card prose">
<h2>6. Constrained placement problems</h2>
<p><strong>Rail.</strong> If the fan must lie on $y_f=c$, the condition $\\nabla T_{\\max}=\\lambda\\nabla g$ gives $\\partial T_{\\max}/\\partial x_f=0$ along the rail and $\\lambda=\\partial T_{\\max}/\\partial y_f$. The multiplier is the <em>shadow price</em> of the rail: $dT^*_{\\max}/dc=\\lambda$.</p>
<p><strong>Airflow budget.</strong> Hold $g=2\\pi\\sigma^2\\Delta h-C=0$. With $F=\\sum_kw_kf_k$ the weighted peak, the Lagrange condition reduces to an elasticity condition,</p>
$$\\sigma\\,\\frac{\\partial F}{\\partial\\sigma}=2\\Delta h\\,\\frac{\\partial F}{\\partial\\Delta h},$$
<p>meaning the next unit of budget is equally effective whether spent widening or strengthening the fan.</p>
</section>

<section class="card">
<h2>Map from the course outcomes to this project</h2>
<div class="table"><table>
<tr><th>Course outcome</th><th>Where it appears</th><th>Section</th></tr>
@@TROWS@@
</table></div>
</section>
"""
page("theory.html", "Theory · PCB Fan Placement", "Derivation of the PCB heat model, its discretisation, the adjoint gradient and KKT optimality, with a course-outcome map.", theory)

# ====================================================================== analysis
analysis = """
<header class="hero">
<div class="eyebrow">Analysis</div>
<h1>What the calculus shows</h1>
<p>Every number below was produced by the Python scripts in the downloads, and is checked against an exact or independent result wherever one exists.</p>
</header>

<section class="card prose">
<h2>1. Validation of the solver</h2>
<div class="table"><table>
<tr><th>Check</th><th>Result</th></tr>
<tr><td>Five-point stencil on a smooth function, $h=0.2\\dots0.025$</td><td>error ratios 3.97, 3.99, 4.00 (theory: 4)</td></tr>
<tr><td>1D fin with exact cosh solution, $\\varepsilon=4,2,1,0.5$ mm</td><td>relative errors $1.9\\times10^{-3}\\to3.0\\times10^{-5}$, observed order 2.00</td></tr>
<tr><td>2D peak temperature, $\\varepsilon=4\\dots0.25$ mm</td><td>observed order 1.87, 2.03, 2.00; Richardson limit 50.568 °C</td></tr>
<tr><td>Uniform source and uniform $h$ (exact $T=T_a+q/2h$)</td><td>maximum error $2.4\\times10^{-12}$</td></tr>
<tr><td>Energy balance, three fan positions</td><td>relative error $\\le10^{-13}$</td></tr>
<tr><td>Adjoint gradient against central differences</td><td>worst relative error $1.6\\times10^{-7}$</td></tr>
</table></div>
@@FIG:exact_checks.png|Convergence of the stencil, the 1D fin and the 2D peak temperature as the grid is refined, with slope-2 reference lines.@@
<p class="hint">At the lab's 2 mm grid there are only four to six cells across a chip, so the peak temperature carries a discretisation error of order 0.2 °C. That does not change which fan position is best.</p>
</section>

<section class="card prose">
<h2>2. Temperature field, fan scan and critical points</h2>
<p>Without a fan the maximum is 69.0 °C at the CPU. Scanning 2091 fan positions (2 mm steps):</p>
<div class="table"><table>
<tr><th>Case</th><th>Fan (x, y) mm</th><th>T<sub>max</sub> (°C)</th><th>Area &gt; 45 °C (mm²)</th></tr>
<tr><td>No fan</td><td class="num">–</td><td class="num">69.04</td><td class="num">6060</td></tr>
<tr><td>Best for T<sub>max</sub></td><td class="num">(38, 50)</td><td class="num">50.50</td><td class="num">352</td></tr>
<tr><td>Best for area &gt; 45 °C</td><td class="num">(38, 48)</td><td class="num">50.65</td><td class="num">348</td></tr>
<tr><td>Worst</td><td class="num">(100, 0)</td><td class="num">68.21</td><td class="num">4808</td></tr>
</table></div>
@@FIG:fan_scan.png|Maximum board temperature as a function of fan position. The minimum sits between the CPU and the regulator, not on the CPU.@@
<p>With the fan at (38, 50) mm the field has three maxima and two saddles (second-derivative test, $D=T_{xx}T_{yy}-T_{xy}^2$):</p>
<div class="table"><table>
<tr><th>Location (mm)</th><th>T (°C)</th><th>T<sub>xx</sub></th><th>D</th><th>Type</th></tr>
<tr><td class="num">(29.5, 50.0)</td><td class="num">50.5</td><td class="num">−0.229</td><td class="num">+0.0521</td><td>maximum (CPU)</td></tr>
<tr><td class="num">(70.2, 54.9)</td><td class="num">50.4</td><td class="num">−0.315</td><td class="num">+0.0956</td><td>maximum (regulator)</td></tr>
<tr><td class="num">(65.3, 20.0)</td><td class="num">48.0</td><td class="num">−0.188</td><td class="num">+0.0673</td><td>maximum (driver)</td></tr>
<tr><td class="num">(69.9, 36.1)</td><td class="num">41.1</td><td class="num">−0.025</td><td class="num">−0.0011</td><td>saddle (regulator–driver)</td></tr>
<tr><td class="num">(48.5, 50.1)</td><td class="num">37.7</td><td class="num">+0.061</td><td class="num">−0.0005</td><td>saddle (CPU–regulator)</td></tr>
</table></div>
</section>

<section class="card prose">
<h2>3. The optimum, by gradient descent and the KKT condition</h2>
<p>Gradient descent from seven starting points converges to the same place, about $(37.7,\\,50.25)$ mm, with $T_{\\max}\\approx50.40$ °C, matching the scan. There the CPU and regulator peaks are 50.414 and 50.385 °C (the driver is at 47.60 °C). The gradients of the two tied peaks are $(+0.3077,+0.0252)$ and $(-0.2175,-0.0191)$ K/mm, an angle of 179.7° apart, and</p>
$$0.414\\,\\nabla f_{\\rm CPU}+0.586\\,\\nabla f_{\\rm reg}=(+0.00006,\\,-0.00076)\\ \\text{K/mm},$$
<p>which is zero to the accuracy of the tie tolerance. Along a line through the optimum the CPU peak alone is lowest (about 49.3 °C) with the fan almost on the CPU, but then the regulator is the hottest part: the optimum gives up about 1.1 °C on the CPU to help the regulator.</p>
@@FIG:gradient_field.png|Gradients of the individual chip peaks with respect to fan position, and descent paths converging to the tie point.@@
</section>

<section class="card prose">
<h2>4. Constrained problems</h2>
<p><strong>Rail $y=15$ mm.</strong> The best position on the rail is $x=34.09$ mm with $T_{\\max}=62.88$ °C, 12.5 °C worse than the free optimum. The multiplier is $\\lambda=-0.2814$ K/mm. Re-optimising with the rail at 14 mm and 16 mm gives $dT^*_{\\max}/dc=-0.2816$ K/mm, in agreement to 0.1%.</p>
<p><strong>Fixed airflow budget</strong> $C=2\\pi\\sigma^2\\Delta h=0.0543$ W/K (the default fan), with the position re-optimised for each radius:</p>
<div class="table"><table>
<tr><th>σ (mm)</th><th>Δh (W/m²K)</th><th>best fan (x, y) mm</th><th>T<sub>max</sub> (°C)</th></tr>
<tr><td class="num">4</td><td class="num">540</td><td class="num">(39.6, 51.1)</td><td class="num">51.71</td></tr>
<tr><td class="num">6</td><td class="num">240</td><td class="num">(40.6, 51.0)</td><td class="num">50.98</td></tr>
<tr><td class="num">8</td><td class="num">135</td><td class="num">(40.7, 50.8)</td><td class="num">50.48</td></tr>
<tr><td class="num">10</td><td class="num">86.4</td><td class="num">(39.7, 50.6)</td><td class="num">50.27</td></tr>
<tr><td class="num">12 (default)</td><td class="num">60</td><td class="num">(37.7, 50.2)</td><td class="num">50.40</td></tr>
<tr><td class="num">16</td><td class="num">33.8</td><td class="num">(31.9, 50.3)</td><td class="num">52.09</td></tr>
<tr><td class="num">20</td><td class="num">21.6</td><td class="num">(32.7, 50.4)</td><td class="num">54.30</td></tr>
<tr><td class="num">30</td><td class="num">9.6</td><td class="num">(34.3, 49.4)</td><td class="num">58.40</td></tr>
<tr><td class="num">40</td><td class="num">5.4</td><td class="num">(36.1, 49.0)</td><td class="num">61.20</td></tr>
</table></div>
<p>The joint optimum is the fan at (39.45, 50.53) mm with $\\sigma=10.3$ mm and $\\Delta h=81.5$ W/(m²K), where $T_{\\max}=50.264$ °C. The multiplier from the $\\sigma$-component is $-106.89$ K per (W/K) and from the $\\Delta h$-component $-107.00$; the elasticity condition reads $-11.605$ K against $-11.618$ K. The radius optimisation improves on the default fan by only <strong>0.14 °C</strong>, so the robust finding is the shape of the trade-off: wide, weak fans lose steeply and very narrow jets lose gently. The shadow price from re-optimising with the budget changed by ±2% is $-102.7$ K per (W/K), within 4% of $\\lambda$, limited by the optimiser's tie tolerance.</p>
@@FIG:lagrange.png|Left: best T_max along the rail and its multiplier. Right: best T_max against fan radius at fixed airflow budget.@@
</section>

<section class="card prose">
<h2>5. Heat flux, Green's and the Divergence theorem</h2>
<p>The flux $\\mathbf F=-kt\\nabla T$ is a gradient field, hence conservative (zero curl). Its divergence is $q-2h\\theta$, so for any loop</p>
$$\\oint_C\\mathbf F\\cdot\\mathbf n\\,ds=\\iint_D(q-2h\\theta)\\,dA.$$
<p>The left side uses only $\\nabla T$ on the loop; the right side uses only chip power and cooling inside. They are computed independently:</p>
<div class="table"><table>
<tr><th>Loop</th><th>Line integral (W)</th><th>Area integral (W)</th></tr>
<tr><td>Circle R = 15 mm around the CPU</td><td class="num">+0.8615</td><td class="num">+0.8618</td></tr>
<tr><td>Circle R = 12 mm around the regulator</td><td class="num">+0.7862</td><td class="num">+0.7863</td></tr>
<tr><td>Circle R = 10 mm around the driver</td><td class="num">+0.6776</td><td class="num">+0.6777</td></tr>
<tr><td>Circle R = 6 mm in the cooled gap (no chip)</td><td class="num">−0.0571</td><td class="num">−0.0569</td></tr>
<tr><td>Circle R = 6 mm at the saddle point</td><td class="num">−0.1424</td><td class="num">−0.1427</td></tr>
<tr><td>Rectangle 15 ≤ x ≤ 85, 38 ≤ y ≤ 68 mm</td><td class="num">+1.0196</td><td class="num">+1.0204</td></tr>
</table></div>
<p>All six agree to within 0.1%. Negative values mean the fan removes more heat inside the loop than the loop generates, so heat flows in. For the 60 °C isotherm without a fan, the arc length is 54.18 mm and Green's theorem gives an enclosed area of 232.8 mm² against 235.0 mm² from counting cells (1%); the isoperimetric ratio $P^2/4\\pi A=1.004$ says the hot region is almost circular. A heat-flow line from (65, 56) mm is 39.6 mm long with peak curvature $\\kappa=0.227$ mm⁻¹; the Hessian formula $\\kappa=|(I-\\tau\\tau^T)H\\tau|/|\\nabla T|$ and the turning rate of the computed path agree to 0.5%. The line is planar, so its torsion is zero.</p>
<div class="callout honest"><strong>Honest note.</strong> The circulation of $\\mathbf F$ being about $10^{-13}$ W is not an independent test: a discrete gradient field has zero circulation by construction. The flux (Divergence theorem) checks are the independent ones.</div>
@@FIG:vector_calculus.png|Loops, the 60 °C isotherm, and a heat-flow line used for the flux, area, arc-length and curvature checks.@@
</section>

<section class="card prose">
<h2>6. Taylor series of the fan profile (live)</h2>
<p>Write the fan as $e^{-u}$, $u=r^2/2\\sigma^2$. The degree-$m$ Taylor polynomial in $u$ is $P_m(u)=\\sum_{k=0}^m(-u)^k/k!$ with error at most $u^{m+1}/(m+1)!$ (alternating series). Since $u\\propto r^2$, degree $m$ in $u$ is order $2m$ in $r$. The lowest model, $m=1$ (order 2 in $r$), is the paraboloid $h\\approx h_{\\rm nat}+\\Delta h(1-r^2/2\\sigma^2)$. It vanishes at $r=\\sqrt2\\,\\sigma$, and its polar integral over that support is $\\pi\\sigma^2\\Delta h$, <em>exactly half</em> of the Gaussian's $2\\pi\\sigma^2\\Delta h$. Move the sliders:</p>
<div class="grid2" style="align-items:center">
<div><canvas class="plot" id="tplot" aria-label="Plot of exp(-u) and its Taylor polynomial against r over sigma"></canvas>
<div class="hint">Black: exact fan profile. Orange: Taylor polynomial. Shaded: where they agree to within 0.02.</div></div>
<div class="ctl" style="display:flex;flex-direction:column;gap:12px">
<div class="row"><label for="tm">Degree m in u (order 2m in r)</label><output id="tmO"></output><input type="range" id="tm" min="0" max="12" step="1" value="1"></div>
<div class="row"><label for="tr">Probe radius r/σ</label><output id="trO"></output><input type="range" id="tr" min="0" max="4" step="0.05" value="1"></div>
<div class="scanread" id="tread" style="font-size:13px"></div>
</div></div>
<p>Using each polynomial only where it is within 2% and zero beyond (high-degree polynomials grow without bound far away), the full board solution gives:</p>
<div class="table"><table>
<tr><th>Fan model</th><th>Valid out to</th><th>Fan effect captured</th><th>T<sub>max</sub> (°C)</th></tr>
<tr><td>Order 2 in r (paraboloid)</td><td class="num">0.64 σ</td><td class="num">18.6%</td><td class="num">61.20</td></tr>
<tr><td>Order 4 in r</td><td class="num">1.01 σ</td><td class="num">40.5%</td><td class="num">55.76</td></tr>
<tr><td>Order 6 in r</td><td class="num">1.32 σ</td><td class="num">57.6%</td><td class="num">53.20</td></tr>
<tr><td>Order 8 in r</td><td class="num">1.57 σ</td><td class="num">71.4%</td><td class="num">52.21</td></tr>
<tr><td>Exact Gaussian</td><td class="num">all r</td><td class="num">100%</td><td class="num">50.40</td></tr>
<tr><td>No fan</td><td class="num">–</td><td class="num">0%</td><td class="num">69.04</td></tr>
</table></div>
<p>The series needs 5 terms for accuracy $10^{-3}$ at $r=\\sigma$, 10 at $2\\sigma$, 17 at $3\\sigma$ and 26 at $4\\sigma$. A local polynomial model cannot represent a fan whose influence extends over several radii; the lab's "Fan profile" menu lets you see this directly. Used only where accurate, the paraboloid recovers 42% of the exact fan's cooling (the 18.6% column is the share of the fan's total convection it carries). On its full support $r\\le\\sqrt2\\sigma$ it carries exactly half.</p>
@@FIG:taylor.png|Taylor polynomials against the exact fan profile, their error bounds, and the resulting board temperatures.@@
</section>

<section class="card prose">
<h2>Assumptions and limits</h2>
<ul>
<li>Chip powers and positions, conductivity, convection coefficients and the Gaussian fan profile are assumed for illustration. A real fan produces a directional flow pattern.</li>
<li>The model is steady state and two-dimensional, with one effective conductivity instead of a layered stack-up, insulated edges, constant $h_{\\rm nat}$ and no radiation. Zhang and Chen include layered conduction, vias, Joule heating and temperature-dependent boundaries; none of that is modelled here.</li>
<li>No experimental validation was performed. The checks show the equations are solved correctly, not that the model matches a physical board.</li>
<li>$T_{\\max}$ is not smooth where peaks tie; the optimiser handles that with the convex-hull treatment of the theory page.</li>
</ul>
</section>
"""
page("analysis.html", "Analysis · PCB Fan Placement", "Validation, critical points, optimum and KKT condition, Lagrange problems, flux theorems and Taylor series for the PCB fan model.", analysis,
     '<script src="assets/taylor_demo.js"></script>')

# ====================================================================== downloads
repo = f'<a class="btn" href="{REPO_URL}">View repository on GitHub</a>' if REPO_URL else ""
files = [
    ("simulation.py", "Core solver: finite-volume build, cached conduction matrix, fan scan, metrics, self-test. Run <code>python simulation.py</code>."),
    ("Exact_checks.py", "Stencil, 1D fin, 2D convergence study with Richardson extrapolation, critical points. Needs <code>sympy</code>."),
    ("Gradient.py", "Adjoint gradients of each chip's peak, KKT tie condition, minimum-norm descent."),
    ("Lagrange.py", "Rail constraint and fixed airflow budget, with multipliers checked by re-optimisation."),
    ("VectorCalculus.py", "Flux, circulation, Divergence and Green checks, isotherm arc length, flow-line curvature."),
    ("Taylor.py", "Taylor polynomials of the fan profile and board solves with the truncated fan."),
]
frows = "\n".join(f'<tr><td class="m"><a href="assets/code/{n}" download>{n}</a></td><td>{d}</td></tr>' for n, d in files)
EXTRA["REPO"] = repo
EXTRA["FROWS"] = frows
downloads = """
<header class="hero">
<div class="eyebrow">Reproducibility</div>
<h1>Code and data</h1>
<p>Everything on this site comes from these files. Put them in one folder and run them; each script writes its figure into a <code>results/</code> folder next to it.</p>
<div class="btns">@@REPO@@</div>
</header>

<section class="card">
<h2>Python scripts</h2>
<div class="table"><table><tr><th>File</th><th>What it does</th></tr>@@FROWS@@</table></div>
<p class="hint" style="margin-top:12px">Requirements: Python 3.9+, <code>numpy</code>, <code>scipy</code>, <code>matplotlib</code>, and <code>sympy</code> for the exact checks. Install with <code>pip install numpy scipy matplotlib sympy</code>.</p>
</section>

<section class="card">
<h2>Data</h2>
<div class="table"><table>
<tr><td class="m"><a href="assets/data/fan_scan.csv" download>fan_scan.csv</a></td><td>T<sub>max</sub>, hottest point and hot area for all 2091 fan positions (2 mm steps, default fan).</td></tr>
<tr><td class="m"><a href="assets/data/summary.csv" download>summary.csv</a></td><td>No fan, best and worst positions.</td></tr>
</table></div>
</section>

<section class="card prose">
<h2>Reproduce the main results</h2>
<pre class="eq mono">python simulation.py --test     # built-in self-test
python simulation.py 2          # scan, 2 mm steps (~11 s)
python Exact_checks.py          # validation figure
python Gradient.py              # adjoint gradients and optimum
python Lagrange.py              # constrained problems
python VectorCalculus.py        # flux and theorem checks
python Taylor.py                # Taylor analysis</pre>
<p class="hint">The browser lab uses the same discretisation (2 mm cells, 5-point stencil, preconditioned conjugate gradients) as <code>simulation.py</code>, so the two agree; the lab's scan uses 4 mm steps for speed.</p>
</section>
"""
page("downloads.html", "Code & data · PCB Fan Placement", "Python scripts and CSV data behind the PCB fan placement analysis.", downloads)
