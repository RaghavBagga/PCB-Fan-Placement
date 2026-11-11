"""
VectorCalculus.py  -  heat-flux field, Divergence/Green theorems, arc length, curvature.

Flux field:      F = -k*t*grad(T)   [W/m]           (heat flows from hot to cold)
Governing law:   div F = q - 2*h*theta,  theta = T - Ta     (this is the PDE:  k t Lap T = 2h theta - q)
Divergence thm:  closed curve C bounding D:   oint_C F.n ds  =  double-integral_D (q - 2 h theta) dA
                 LEFT uses only the gradient of T on the curve; RIGHT uses only the sources and cooling inside.
Conservative:    F = -grad(k t T) is a gradient field, so curl F = 0 and the circulation around any loop is 0.
"""
import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from simulation import *

p = Params(dx=1e-3)                                  # 1 mm grid: chip edges fall exactly on cell edges
KT = p.k_eff * p.t
FAN = (0.0377, 0.05025)                              # optimum found in Gradient.py
nx, ny, X, Y = grid(p)
xs, ys = X[:, 0] * 1e3, Y[0, :] * 1e3                # mm
T = solve(p, FAN)
S = RectBivariateSpline(xs, ys, T, kx=3, ky=3)       # smooth interpolant of the temperature field (mm coordinates)

def gradT(x, y):                                     # K/m, x and y in mm
    return 1e3 * S(x, y, dx=1, grid=False), 1e3 * S(x, y, dy=1, grid=False)
def Fvec(x, y):                                      # W/m
    gx, gy = gradT(x, y); return -KT * gx, -KT * gy
def h_at(x, y):                                      # W/(m^2 K); x, y in mm
    r2 = ((x - FAN[0]*1e3)**2 + (y - FAN[1]*1e3)**2) * 1e-6
    return p.h_nat + p.dh * np.exp(-r2 / (2 * p.sigma**2))
def q_at(x, y):                                      # W/m^2, exact chip rectangles
    out = np.zeros_like(x, dtype=float)
    for s in p.sources:
        m = (np.abs(x - s.x*1e3) <= s.w*500) & (np.abs(y - s.y*1e3) <= s.d*500)
        out += m * s.P / (s.w * s.d)
    return out

def flux_through_circle(cx, cy, R, N=2880):          # all in mm; returns W
    th = (np.arange(N) + .5) * 2 * np.pi / N
    x, y = cx + R*np.cos(th), cy + R*np.sin(th)
    Fx, Fy = Fvec(x, y)
    return np.sum(Fx*np.cos(th) + Fy*np.sin(th)) * (R*1e-3) * 2*np.pi/N
def circulation_on_circle(cx, cy, R, N=2880):
    th = (np.arange(N) + .5) * 2 * np.pi / N
    x, y = cx + R*np.cos(th), cy + R*np.sin(th)
    Fx, Fy = Fvec(x, y)
    return np.sum(-Fx*np.sin(th) + Fy*np.cos(th)) * (R*1e-3) * 2*np.pi/N
def area_integral(cx, cy, R, Nr=700, Nt=1440):       # integral over the disc of (q - 2 h theta) dA   [W], polar midpoint rule
    r = (np.arange(Nr) + .5) * R / Nr; th = (np.arange(Nt) + .5) * 2*np.pi / Nt
    rr, tt = np.meshgrid(r, th, indexing="ij")
    x, y = cx + rr*np.cos(tt), cy + rr*np.sin(tt)
    f = q_at(x, y) - 2 * h_at(x, y) * (S(x, y, grid=False) - p.Ta)
    return np.sum(f * rr) * (R/Nr) * (2*np.pi/Nt) * 1e-6

def flux_through_rect(x0, x1, y0, y1, N=3000):       # outward flux through a rectangle, mm; returns W
    xm = x0 + (np.arange(N) + .5) * (x1 - x0) / N; ym = y0 + (np.arange(N) + .5) * (y1 - y0) / N
    dxs_, dys_ = (x1 - x0) / N * 1e-3, (y1 - y0) / N * 1e-3
    fR = Fvec(np.full(N, x1), ym)[0].sum() * dys_ - Fvec(np.full(N, x0), ym)[0].sum() * dys_
    fT = Fvec(xm, np.full(N, y1))[1].sum() * dxs_ - Fvec(xm, np.full(N, y0))[1].sum() * dxs_
    return fR + fT
def circulation_rect(x0, x1, y0, y1, N=3000):        # counter-clockwise circulation of F around a rectangle [W]
    xm = x0 + (np.arange(N) + .5) * (x1 - x0) / N; ym = y0 + (np.arange(N) + .5) * (y1 - y0) / N
    dxs_, dys_ = (x1 - x0) / N * 1e-3, (y1 - y0) / N * 1e-3
    return (Fvec(xm, np.full(N, y0))[0].sum() * dxs_ + Fvec(np.full(N, x1), ym)[1].sum() * dys_
            - Fvec(xm, np.full(N, y1))[0].sum() * dxs_ - Fvec(np.full(N, x0), ym)[1].sum() * dys_)
def area_rect(x0, x1, y0, y1, Nx=1400, Ny=600):
    xm = x0 + (np.arange(Nx) + .5) * (x1 - x0) / Nx; ym = y0 + (np.arange(Ny) + .5) * (y1 - y0) / Ny
    xx, yy = np.meshgrid(xm, ym, indexing="ij")
    f = q_at(xx, yy) - 2 * h_at(xx, yy) * (S(xx, yy, grid=False) - p.Ta)
    return f.sum() * ((x1 - x0) / Nx) * ((y1 - y0) / Ny) * 1e-6

# --------------------------------------------------------------- 1. Divergence theorem
print("=== 1. Divergence theorem on closed loops:  oint F.n ds  vs  double integral of (q - 2 h theta) dA ===")
print("  loop                                    flux through curve [W]   area integral [W]   difference")
loops = [("circle R=15 around CPU (30,50)", 30, 50, 15), ("circle R=12 around regulator (70,55)", 70, 55, 12),
         ("circle R=10 around driver (65,20)", 65, 20, 10), ("circle R=6 in cooled gap (50,30), no chip", 50, 30, 6), ("circle R=6 at the saddle (48.5,50.1)", 48.5, 50.1, 6)]
rows1 = []
for name, cx, cy, R in loops:
    a, b = flux_through_circle(cx, cy, R), area_integral(cx, cy, R); rows1.append((a, b))
    print(f"  {name:42s} {a:+14.6f}         {b:+14.6f}      {abs(a-b):.1e}")
RECT = (15, 85, 38, 68)
a, b = flux_through_rect(*RECT), area_rect(*RECT)
print(f"  {'rectangle x 15-85, y 38-68 (CPU+regulator)':42s} {a:+14.6f}         {b:+14.6f}      {abs(a-b):.1e}")
print("  (positive = net heat leaves through the curve; negative = heat flows in, the fan removes it)")

# --------------------------------------------------------------- 2. Circulation (Green / Stokes in the plane)
print("\n=== 2. F is conservative: circulation  oint F.dr  around closed loops (analytically 0, since F = -grad(k t T)) ===")
for name, cx, cy, R in loops[:3]:
    print(f"  {name:42s} circulation = {circulation_on_circle(cx, cy, R):+.2e} W   (typical |F|R = {np.hypot(*Fvec(cx+R, cy))*R*1e-3:.2e} W)")
print(f"  {'rectangle x 15-85, y 38-68':42s} circulation = {circulation_rect(*RECT):+.2e} W")

# --------------------------------------------------------------- 3. Isotherm: arc length and Green's area
print("\n=== 3. The T = %.0f C isotherm without a fan: arc length (line integral) and Green's-theorem area ===" % p.T_limit)
T0 = solve(p)
fig0, ax0 = plt.subplots(); cs = ax0.contour(xs, ys, T0.T, [p.T_limit]); plt.close(fig0)
seg = max(cs.allsegs[0], key=len)                    # the closed loop around the CPU (mm)
x, y = seg[:, 0], seg[:, 1]
perim = np.sum(np.hypot(np.diff(x), np.diff(y)))
A_green = 0.5 * abs(np.sum(x[:-1]*y[1:] - x[1:]*y[:-1]))        # (1/2) oint (x dy - y dx)
A_cells = float((T0 > p.T_limit).sum() * p.dx**2 * 1e6)
print(f"  arc length (perimeter)         = {perim:.2f} mm")
print(f"  area by Green's theorem        = {A_green:.1f} mm^2   (1/2 * contour integral of x dy - y dx)")
print(f"  area by counting grid cells    = {A_cells:.1f} mm^2   (double integral of the indicator)")
print(f"  isoperimetric ratio P^2/(4 pi A) = {perim**2/(4*np.pi*A_green):.3f}   (>= 1; equals 1 for a perfect circle)")
mid = np.column_stack([(x[:-1]+x[1:])/2, (y[:-1]+y[1:])/2]); tx, ty = np.diff(x), np.diff(y); tn = np.hypot(tx, ty)
S0 = RectBivariateSpline(xs, ys, T0, kx=3, ky=3)
gx0, gy0 = S0(mid[:, 0], mid[:, 1], dx=1, grid=False), S0(mid[:, 0], mid[:, 1], dy=1, grid=False)
cosang = np.abs((tx*gx0 + ty*gy0) / (tn * np.hypot(gx0, gy0)))
print(f"  heat flux is perpendicular to the isotherm: |cos(angle between F and tangent)| median {np.median(cosang):.3f}, max {cosang.max():.3f}  (max = {np.degrees(np.arcsin(cosang.max())):.1f} deg off perpendicular)")

# --------------------------------------------------------------- 4. A heat-flow line as a space curve: curvature
print("\n=== 4. A heat-flow line r(s) (unit speed, dr/ds = F/|F|): curvature from the Hessian vs from the path ===")
def rhs(s, u):
    fx, fy = Fvec(u[0], u[1]); n = np.hypot(fx, fy); return [fx/n, fy/n]
start = [65.0, 56.0]
ev = lambda s, u: np.hypot(u[0] - FAN[0]*1e3, u[1] - FAN[1]*1e3) - 3.0              # stop within 3 mm of the fan centre
ev.terminal = True
edge = lambda s, u: min(u[0] - 1, 99 - u[0], u[1] - 1, 79 - u[1])            # stop before leaving the board
edge.terminal = True
sol = solve_ivp(rhs, [0, 80], start, max_step=0.1, events=[ev, edge], rtol=1e-9, atol=1e-9, dense_output=False)
px, py, ps = sol.y[0], sol.y[1], sol.t
Hxx = S(px, py, dx=2, grid=False); Hyy = S(px, py, dy=2, grid=False); Hxy = S(px, py, dx=1, dy=1, grid=False)
gx, gy = S(px, py, dx=1, grid=False), S(px, py, dy=1, grid=False)
gn = np.hypot(gx, gy); tau = np.column_stack([-gx, -gy]) / gn[:, None]              # unit tangent (direction of F)
Ht = np.column_stack([Hxx*tau[:, 0] + Hxy*tau[:, 1], Hxy*tau[:, 0] + Hyy*tau[:, 1]])
perp = Ht - np.sum(Ht*tau, axis=1)[:, None]*tau
kappa_H = np.hypot(perp[:, 0], perp[:, 1]) / gn                                        # 1/mm
ang = np.unwrap(np.arctan2(tau[:, 1], tau[:, 0])); kappa_P = np.abs(np.gradient(ang, ps))   # from the path itself
sel = slice(5, -5)
print(f"  path length = {ps[-1]:.1f} mm, starts ({start[0]:.0f}, {start[1]:.1f}) mm, ends ({px[-1]:.1f}, {py[-1]:.1f}) mm")
print(f"  curvature from Hessian: max {kappa_H.max():.4f} 1/mm,   from the path: max {kappa_P[sel].max():.4f} 1/mm")
print(f"  max |difference| = {np.abs(kappa_H - kappa_P)[sel].max():.1e} 1/mm   (torsion = 0: the curve lies in a plane)")

# --------------------------------------------------------------- figure
fig, ax = plt.subplots(1, 3, figsize=(18, 5.2))
im = ax[0].imshow(T.T, origin="lower", extent=[0, 100, 0, 80], cmap="inferno")
gx_, gy_ = np.meshgrid(np.linspace(1, 99, 120), np.linspace(1, 79, 96), indexing="ij")
Fx_, Fy_ = Fvec(gx_, gy_); mag = np.hypot(Fx_, Fy_)
ax[0].streamplot(gx_[:, 0], gy_[0, :], Fx_.T, Fy_.T, color="w", density=1.6, linewidth=0.5 + 1.2*mag.T/mag.max(), arrowsize=0.8)
ax[0].add_patch(plt.Rectangle((RECT[0], RECT[2]), RECT[1]-RECT[0], RECT[3]-RECT[2], fill=False, ec="#6fe3ff", lw=1.5, ls="--"))
for name, cx, cy, R in loops[:3]: ax[0].add_patch(plt.Circle((cx, cy), R, fill=False, ec="#6fe3ff", lw=1.5, ls="--"))
ax[0].plot(px, py, color="#6fe3ff", lw=2.2); ax[0].plot(FAN[0]*1e3, FAN[1]*1e3, "c*", ms=14, mec="k")
ax[0].set_xlim(0, 100); ax[0].set_ylim(0, 80); ax[0].set_title("Heat-flux field F = −k·t·∇T (streamlines), fan at optimum")
ax[0].set_xlabel("x [mm]"); ax[0].set_ylabel("y [mm]"); plt.colorbar(im, ax=ax[0], label="T [°C]")

Rs = np.arange(1, 34.1, 1.0)
fl = [flux_through_circle(30, 50, R) for R in Rs]; ar = [area_integral(30, 50, R) for R in Rs]
ax[1].plot(Rs, ar, "-", lw=5, color="#9ecae1", label="area integral of (q − 2hθ)")
ax[1].plot(Rs, fl, "k.", ms=6, label="flux line integral ∮F·n ds")
ax[1].axhline(2.0, color="grey", lw=.7, ls=":"); ax[1].text(1.5, 2.04, "CPU power 2.0 W", fontsize=8, color="grey")
ax[1].axvline(6*np.sqrt(2), color="grey", lw=.7, ls=":"); ax[1].set_xlabel("radius R of circle centred on the CPU [mm]"); ax[1].set_ylabel("net heat leaving the circle [W]")
ax[1].set_title("Divergence theorem: line integral = area integral"); ax[1].legend(fontsize=8)

ax[2].plot(ps, kappa_H, "-", lw=4, color="#9ecae1", label="from Hessian (analytic formula)")
ax[2].plot(ps[sel], kappa_P[sel], "k.", ms=3, label="from the path itself (dθ/ds)")
ax[2].set_xlabel("arc length s along the heat-flow line [mm]"); ax[2].set_ylabel("curvature κ [1/mm]"); ax[2].set_title("Curvature of a heat-flow line"); ax[2].legend(fontsize=8)
plt.tight_layout(); RESULTS.mkdir(exist_ok=True); plt.savefig(RESULTS / "vector_calculus.png", dpi=115)
print("\nSaved", RESULTS / "vector_calculus.png")
