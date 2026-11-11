"""
Exact checks for the PCB thermal model.
  1. Taylor-series derivation of the 5-point stencil and its O(dx^2) error (symbolic + numeric)
  2. 1D fin ODE: exact solution vs numerics, convergence order
  3. 2D grid convergence with Richardson extrapolation
  4. Critical points of T(x,y): maxima, saddles, classified by the second-derivative test
"""
import numpy as np, sympy as sp_, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline
from scipy.optimize import fsolve
from scipy.ndimage import minimum_filter
from simulation import *   # core model lives in simulation.py

# ------------------------------------------------------------------ 1. Taylor
print("=== 1. Taylor series -> five-point stencil ===")
x, y, h = sp_.symbols("x y h")
u = sp_.Function("u")
def taylor(dx_, dy_, order=6):
    s = 0
    for n in range(order + 1):
        for m in range(n + 1):
            s += sp_.diff(u(x, y), x, n - m, y, m) * (dx_**(n - m)) * (dy_**m) / (sp_.factorial(n - m) * sp_.factorial(m)) if (n - m or m) else u(x, y)
    return s
# generic Taylor expansion of u(x+dx,y+dy) built term by term, then stencil
stencil = (taylor(h, 0) + taylor(-h, 0) + taylor(0, h) + taylor(0, -h) - 4*u(x, y)) / h**2
stencil = sp_.expand(stencil)
print("  (u_E+u_W+u_N+u_S-4u)/h^2 =", sp_.simplify(stencil.subs(h, 0)), "+ h^2 * [", sp_.simplify(stencil.coeff(h, 2)), "] + O(h^4)")
lap = sp_.diff(u(x, y), x, 2) + sp_.diff(u(x, y), y, 2)
print("  leading term equals Laplacian:", sp_.simplify(stencil.subs(h, 0) - lap) == 0)
print("  -> truncation error = (h^2/12)(u_xxxx + u_yyyy)")

f = lambda X, Y: np.sin(2*X) * np.cos(3*Y)             # smooth test function, exact Laplacian = -13 f
x0, y0 = 0.7, 0.4; exact = -13 * f(x0, y0)
errs = []
for hh in [0.2, 0.1, 0.05, 0.025]:
    st = (f(x0+hh, y0) + f(x0-hh, y0) + f(x0, y0+hh) + f(x0, y0-hh) - 4*f(x0, y0)) / hh**2
    errs.append(abs(st - exact))
print("  stencil error for h = 0.2,0.1,0.05,0.025:", ["%.2e" % e for e in errs])
print("  error ratios (expect ~4):", ["%.2f" % (errs[i]/errs[i+1]) for i in range(3)])

# ------------------------------------------------------------------ 2. 1D ODE
print("\n=== 2. 1D fin ODE:  theta'' - theta/L^2 = -(q'/kt) delta(x-x0),  insulated ends ===")
kt = 15 * 1.6e-3; hh_ = 10.0; L = np.sqrt(kt / (2*hh_)); qline = 1.0; Lx = 0.40
def exact1d(xx, x0, w):
    xl, xg = np.minimum(xx, x0), np.maximum(xx, x0)
    G = (qline * L / kt) * np.cosh(xl / L) * np.cosh((Lx - xg) / L) / np.sinh(Lx / L)
    return G * np.sinh(w / (2*L)) / (w / (2*L))             # finite strip of width w, valid outside the strip
res1 = []
for dx in [4e-3, 2e-3, 1e-3, 0.5e-3]:
    p = Params(Lx=Lx, Ly=4*dx, sources=[], dx=dx)
    nx, ny, X, Y = grid(p)
    ic = int(round(0.2 / dx - 0.5)); x0c = (ic + 0.5) * dx            # source column = one cell
    q = np.zeros_like(X); q[ic, :] = qline / dx
    T = solve(p, q=q, h=np.full_like(X, hh_))
    pts = [int(round((x0c + d) / dx - 0.5)) for d in (-0.06, -0.02, 0.02, 0.06)]
    err = max(abs(T[i, 0] - p.Ta - exact1d(X[i, 0], x0c, dx)) / exact1d(X[i, 0], x0c, dx) for i in pts)
    res1.append((dx, err))
    print(f"  dx = {dx*1e3:4.2f} mm   max relative error at 20 and 60 mm from source = {err:.2e}")
o1 = [np.log2(res1[i][1] / res1[i+1][1]) for i in range(3)]
print("  observed order of accuracy:", ["%.2f" % o for o in o1], "(theory: 2)")

# ------------------------------------------------------------------ 3. 2D convergence
print("\n=== 3. 2D grid convergence (fan at 38,50; chips aligned to a 4 mm grid) ===")
chips = [Source("CPU", .030, .050, .012, .012, 2.0), Source("Reg", .068, .056, .008, .008, 1.0), Source("Drv", .064, .020, .008, .008, 0.8)]
dxs = [4e-3, 2e-3, 1e-3, 0.5e-3, 0.25e-3]
def peak_subcell(T):
    """Max of T located between cell centres: vertex of the parabola through the peak cell and its neighbours."""
    i, j = np.unravel_index(T.argmax(), T.shape); T0 = T[i, j]; add = 0.0
    for a, b in [(T[i-1, j], T[i+1, j]), (T[i, j-1], T[i, j+1])]:
        add += -(b - a)**2 / (8 * (a - 2*T0 + b))
    return T0 + add
def bilinear(T, p, x, y):
    fx, fy = x/p.dx - 0.5, y/p.dx - 0.5; i, j = int(np.floor(fx)), int(np.floor(fy)); u, v = fx - i, fy - j
    return (1-u)*(1-v)*T[i, j] + u*(1-v)*T[i+1, j] + (1-u)*v*T[i, j+1] + u*v*T[i+1, j+1]
rows = []
for dx in dxs:
    p = Params(dx=dx, sources=chips); T = solve(p, (0.038, 0.050))
    rows.append((dx, peak_subcell(T), bilinear(T, p, .030, .050), T.mean() + 0.0))
    print(f"  dx = {dx*1e3:5.2f} mm   raw max = {T.max():.5f}   sub-cell peak = {rows[-1][1]:.5f}   T(30,50) = {rows[-1][2]:.5f}   mean = {rows[-1][3]:.5f} C")
print("  (the raw max over cell centres wobbles by ~|T''|dx^2/8 because the true peak lies between centres)")
plt_dx = plt_err = None
for name, col in [("peak T (sub-cell)", 1), ("T at (30,50)", 2), ("mean T", 3)]:
    v = [r[col] for r in rows]
    d = [v[i] - v[i+1] for i in range(len(v) - 1)]
    order = [np.log2(abs(d[i] / d[i+1])) for i in range(len(d) - 1)]
    extrap = v[-1] - d[-1] / (2**order[-1] - 1)
    print(f"  {name:18s} observed order {['%.2f' % o for o in order]}   Richardson limit = {extrap:.5f}")
    if col == 1: plt_dx = np.array(dxs[:-1]) * 1e3; plt_err = np.abs(np.array(v[:-1]) - extrap)

fig, ax = plt.subplots(1, 2, figsize=(10, 4))
ax[0].loglog([r[0]*1e3 for r in res1], [r[1] for r in res1], "o-", label="1D fin vs exact ODE solution")
ax[0].loglog(plt_dx, plt_err, "s-", label="2D board: |peak T - Richardson limit|")
ref = np.array([0.25, 4.0]); ax[0].loglog(ref, 3e-3 * (ref/4.0)**2, "k--", label="slope 2 (theory)")
ax[0].set_xlabel("grid spacing dx [mm]"); ax[0].set_ylabel("error"); ax[0].legend(fontsize=8); ax[0].set_title("Convergence: error falls as dx²")

# ------------------------------------------------------------------ 4. critical points
print("\n=== 4. Critical points of T(x,y), fan at (38,50) (default chip layout) ===")
p = Params(dx=1e-3); nx, ny, X, Y = grid(p)
T = solve(p, (0.038, 0.050))
xs, ys = X[:, 0] * 1e3, Y[0, :] * 1e3
S = RectBivariateSpline(xs, ys, T, kx=3, ky=3)
grad = lambda v: [S(v[0], v[1], dx=1, grid=False), S(v[0], v[1], dy=1, grid=False)]
gm = np.hypot(np.gradient(T, xs, axis=0), np.gradient(T, ys, axis=1))
from scipy.ndimage import maximum_filter
seeds = np.argwhere(((T == maximum_filter(T, size=15)) | (T == minimum_filter(T, size=15)) |
                     ((gm == minimum_filter(gm, size=15)) & (gm < 0.3))))
crit = []
for i, j in seeds:
    if i < 3 or j < 3 or i > nx-4 or j > ny-4: continue
    sol, info, ok, _ = fsolve(grad, [xs[i], ys[j]], full_output=True)
    if ok != 1 or not (2 < sol[0] < 98 and 2 < sol[1] < 78): continue
    if any(np.hypot(sol[0]-c[0], sol[1]-c[1]) < 1.5 for c in crit): continue
    fxx = S(*sol, dx=2, grid=False); fyy = S(*sol, dy=2, grid=False); fxy = S(*sol, dx=1, dy=1, grid=False)
    D = fxx*fyy - fxy**2
    kind = "local max" if D > 0 and fxx < 0 else "local min" if D > 0 else "saddle"
    crit.append((sol[0], sol[1], float(S(*sol, grid=False)), fxx, D, kind))
for c in sorted(crit, key=lambda c: -c[2]):
    print(f"  ({c[0]:5.1f}, {c[1]:5.1f}) mm  T = {c[2]:5.1f} C   Txx = {c[3]:+.3f}  D = Txx*Tyy - Txy^2 = {c[4]:+.4f}   -> {c[5]}")

ext = [0, 100, 0, 80]
im = ax[1].imshow(T.T, origin="lower", extent=ext, cmap="inferno")
ax[1].contour(xs, ys, T.T, levels=np.arange(30, 70, 5), colors="w", linewidths=.5, alpha=.5)
mk = {"local max": ("c^", "max"), "saddle": ("gs", "saddle"), "local min": ("bv", "min")}
seen = set()
for c in crit:
    s, lab = mk[c[5]]; ax[1].plot(c[0], c[1], s, ms=9, label=lab if lab not in seen else None); seen.add(lab)
ax[1].legend(fontsize=8); ax[1].set_title("Critical points of T (second-derivative test)")
ax[1].set_xlabel("x [mm]"); ax[1].set_ylabel("y [mm]"); plt.colorbar(im, ax=ax[1], label="T [°C]")
plt.tight_layout(); RESULTS.mkdir(exist_ok=True); plt.savefig(RESULTS / "exact_checks.png", dpi=120)
