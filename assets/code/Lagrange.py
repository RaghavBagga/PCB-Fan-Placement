"""
Lagrange.py  -  constrained optimisation of the fan, using the exact gradients (adjoint solve) from Gradient.py.

Problem A (rail):    fan mounted on a rail y = c.     minimise Tmax(xf, yf)  subject to  g = yf - c = 0.
Problem B (budget):  fixed total fan effect.          minimise Tmax(xf, yf, sigma, dh)
                     subject to  g = 2*pi*sigma^2*dh - C = 0,   C = integral of dh*exp(-r^2/2sigma^2) dA (polar coordinates).
Tmax = max over the chip peaks f_k, so the multiplier condition is  sum_k w_k grad f_k = lambda * grad g,
with weights w_k >= 0, sum w_k = 1, nonzero only for the peaks that tie (KKT / Clarke).
"""
import itertools
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import minimize_scalar
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from simulation import *
from simulation import _conduction

p = Params()
nx, ny, X, Y = grid(p)
q = source_field(p, X, Y)
A0, DPOS, DEG = _conduction(nx, ny, p.k_eff * p.t)
CELLS = [np.flatnonzero(((np.abs(X - s.x) <= s.w / 2) & (np.abs(Y - s.y) <= s.d / 2)).ravel()) for s in p.sources]
NAMES = [s.name for s in p.sources]
C_BUDGET = 2 * np.pi * p.sigma**2 * p.dh            # W/K : the default fan's total effect


def peaks(xf, yf, sig, dh, want_grad=True):
    """Chip peaks [C] and, if wanted, gradients [K per mm / K per (W/m2K)] w.r.t. (xf, yf, sigma, dh)."""
    r2 = (X - xf)**2 + (Y - yf)**2
    e = np.exp(-r2 / (2 * sig**2)); g = dh * e
    loss = 2 * (p.h_nat + g).ravel() * p.dx**2
    data = A0.data.copy(); data[DPOS] = DEG + loss
    lu = spla.splu(sp.csc_matrix((data, A0.indices, A0.indptr), shape=A0.shape))
    th = lu.solve(q.ravel() * p.dx**2)
    vals = np.array([th[c].max() + p.Ta for c in CELLS])
    if not want_grad:
        return vals
    dH = np.stack([g * (X - xf) / sig**2, g * (Y - yf) / sig**2, g * r2 / sig**3, e]).reshape(4, -1)
    scale = np.array([1e-3, 1e-3, 1e-3, 1.0])            # per metre -> per mm for the three lengths
    G = []
    for c in CELLS:
        i = c[np.argmax(th[c])]; u = np.zeros(nx * ny); u[i] = 1.0
        w = lu.solve(u) * 2 * p.dx**2 * th
        G.append(-(dH @ w) * scale)
    return vals, np.array(G)


def min_norm_hull(G):
    """Minimum-norm point in the convex hull of the rows of G; returns (vector, weights over rows)."""
    best, bw = None, None
    for r in range(1, len(G) + 1):
        for S in itertools.combinations(range(len(G)), r):
            Gs = G[list(S)]; M = np.zeros((r + 1, r + 1)); M[:r, :r] = Gs @ Gs.T; M[:r, r] = 1; M[r, :r] = 1
            rhs = np.zeros(r + 1); rhs[r] = 1
            lam = np.linalg.lstsq(M, rhs, rcond=None)[0][:r]
            if np.all(lam >= -1e-12):
                d = lam @ Gs
                if best is None or np.linalg.norm(d) < np.linalg.norm(best) - 1e-15:
                    best = d; bw = np.zeros(len(G)); bw[list(S)] = lam
    return best, bw


# =========================================================== Problem A: fan on a rail y = c
print("=== A. Fan mounted on a rail y = c  (constraint g = yf - c = 0, grad g = (0, 1)) ===")
def rail_opt(c_mm):
    f = lambda x: peaks(x * 1e-3, c_mm * 1e-3, p.sigma, p.dh, False).max()
    xs = np.arange(0, 100.01, 1.0); fx = [f(x) for x in xs]; k = int(np.argmin(fx))
    r = minimize_scalar(f, bounds=(max(0, xs[k] - 2), min(100, xs[k] + 2)), method="bounded", options=dict(xatol=1e-4))
    return r.x, r.fun
c = 15.0
xr, Tr = rail_opt(c)
v, G = peaks(xr * 1e-3, c * 1e-3, p.sigma, p.dh)
act = np.flatnonzero(v >= v.max() - 0.02)
_, w = min_norm_hull(G[act][:, :1])                              # weights that cancel the x-components
w_full = np.zeros(3); w_full[act] = w
lam_A = float(w_full @ G[:, 1])
print(f"  rail y = {c:.0f} mm: best x = {xr:.2f} mm, Tmax = {Tr:.3f} C   (active peaks: {[NAMES[i] for i in act]})")
print(f"  weights w = {np.round(w_full, 3)}   sum w*dTx = {w_full @ G[:, 0]:+.5f} K/mm (should be ~0)")
print(f"  multiplier lambda = sum w*dTy = {lam_A:+.4f} K per mm of rail height")
_, Tm = rail_opt(c - 1.0); _, Tp = rail_opt(c + 1.0)
print(f"  check (envelope theorem): d(Tmax*)/dc by re-optimising at c = 14 and 16 mm: {(Tp - Tm) / 2:+.4f} K/mm")

# =========================================================== Problem B: airflow budget
print("\n=== B. Fixed total fan effect: 2*pi*sigma^2*dh = C = %.4f W/K (the default fan) ===" % C_BUDGET)
def reduced(z, C):
    """Values and reduced gradients in z = (xf, yf, sigma) [mm], with dh = C/(2 pi sigma^2) substituted."""
    sig = z[2] * 1e-3; dh = C / (2 * np.pi * sig**2)
    v, G = peaks(z[0] * 1e-3, z[1] * 1e-3, sig, dh)
    Gr = np.column_stack([G[:, 0], G[:, 1], G[:, 2] + G[:, 3] * (-2 * dh / z[2])])    # chain rule: ddh/dsigma = -2 dh/sigma
    return v, Gr, G, dh
def Tmax_z(z, C):
    sig = z[2] * 1e-3
    return peaks(z[0] * 1e-3, z[1] * 1e-3, sig, C / (2 * np.pi * sig**2), False).max()
LO = np.array([0, 0, 4.0]); HI = np.array([100, 80, 40.0])
def descend(z0, C, eps=0.03, maxit=300, fix_sigma=None):
    z = np.array(z0, float)
    for it in range(maxit):
        v, Gr, _, _ = reduced(z, C); Tm = v.max()
        act = np.flatnonzero(v >= Tm - eps)
        Ga = Gr[act].copy()
        if fix_sigma is not None: Ga[:, 2] = 0
        d, _ = min_norm_hull(Ga); gn = np.linalg.norm(d)
        if gn < 2e-3: break
        step = -d / gn; s = 4.0
        while s > 1e-6:
            zn = np.clip(z + s * step, LO, HI)
            if fix_sigma is not None: zn[2] = fix_sigma
            if Tmax_z(zn, C) < Tm - 1e-4 * s * gn: break
            s /= 2
        if s <= 1e-6: break
        z = zn
    return z

sig_list = np.array([4, 5, 6, 8, 10, 12, 14, 16, 20, 25, 30, 40], float)
curve = []
for sg in sig_list:
    z = descend([37.7, 50.2, sg], C_BUDGET, fix_sigma=sg)
    curve.append((sg, Tmax_z(z, C_BUDGET), z[0], z[1], C_BUDGET / (2 * np.pi * (sg * 1e-3)**2)))
print("  sigma [mm]  dh [W/m2K]   best fan (x, y) [mm]    Tmax [C]   (position re-optimised for each sigma)")
for sg, T_, x_, y_, d_ in curve:
    print(f"  {sg:8.1f}   {d_:9.1f}   ({x_:5.1f}, {y_:5.1f})        {T_:8.3f}")

z_star = descend([37.7, 50.2, 12.0], C_BUDGET)
Tstar = Tmax_z(z_star, C_BUDGET)
v, Gr, G, dh_star = reduced(z_star, C_BUDGET)
act = np.flatnonzero(v >= v.max() - 0.05)
_, w = min_norm_hull(Gr[act]); wf = np.zeros(3); wf[act] = w
print(f"\n  joint optimum over (xf, yf, sigma): fan ({z_star[0]:.2f}, {z_star[1]:.2f}) mm, sigma = {z_star[2]:.2f} mm, dh = {dh_star:.1f} W/m2K")
print(f"  Tmax = {Tstar:.3f} C   (default fan, re-optimised position only: {Tmax_z(descend([37.7,50.2,12.0], C_BUDGET, fix_sigma=12.0), C_BUDGET):.3f} C)")
print(f"  active peaks: {[NAMES[i] for i in act]} with weights {np.round(wf[act], 3)}")
Fg = wf @ G                                                      # weighted gradient w.r.t. (xf, yf, sigma, dh) in the ORIGINAL variables
sg_m = z_star[2]
# constraint g = 2*pi*sigma^2*dh - C:  dg/dsigma = 4*pi*sigma*dh,  dg/ddh = 2*pi*sigma^2   (sigma in metres here)
sig_m = sg_m * 1e-3
dg_dsig_mm = 4 * np.pi * sig_m * dh_star * 1e-3
dg_ddh = 2 * np.pi * sig_m**2
lam_s = Fg[2] / dg_dsig_mm; lam_d = Fg[3] / dg_ddh
print(f"  Lagrange condition  sum_k w_k grad f_k = lambda * grad g  in the original variables (sigma, dh):")
print(f"    from the sigma component : lambda = {lam_s:+.2f} K per (W/K)")
print(f"    from the dh component    : lambda = {lam_d:+.2f} K per (W/K)    (the two must agree)")
print(f"    position components: sum w*dT/dxf = {Fg[0]:+.5f}, sum w*dT/dyf = {Fg[1]:+.5f} K/mm (must vanish: position is unconstrained)")
print(f"  elasticity form: sigma*dT/dsigma = {sg_m * Fg[2]:+.4f} K,  2*dh*dT/ddh = {2 * dh_star * Fg[3]:+.4f} K   (equal at the optimum)")
Tp_ = Tmax_z(descend(z_star, C_BUDGET * 1.02), C_BUDGET * 1.02); Tm_ = Tmax_z(descend(z_star, C_BUDGET * 0.98), C_BUDGET * 0.98)
print(f"  shadow price check: dTmax*/dC from re-optimising at C*0.98 and C*1.02 = {(Tp_ - Tm_) / (0.04 * C_BUDGET):+.2f} K per (W/K); lambda = {lam_s:+.2f}")

# =========================================================== figure
fig, ax = plt.subplots(1, 3, figsize=(18, 5))
xs_ = np.arange(0, 101, 4.0); ys_ = np.arange(0, 81, 4.0)
TM = np.array([[peaks(a * 1e-3, b * 1e-3, p.sigma, p.dh, False).max() for b in ys_] for a in xs_])
im = ax[0].contourf(xs_, ys_, TM.T, 24, cmap="viridis"); plt.colorbar(im, ax=ax[0], label="Tmax [°C]")
ax[0].contour(xs_, ys_, TM.T, [Tr], colors="w", linewidths=1.0, linestyles="--")
ax[0].axhline(c, color="#ff6b3d", lw=2, label=f"rail y = {c:.0f} mm  (g = 0)")
ax[0].plot(xr, c, "o", color="#ff6b3d", mec="w", ms=10, label=f"constrained optimum ({xr:.1f}, {c:.0f})")
gx, gy = G[:, :2][act].T
ax[0].annotate("", xy=(xr, c + 11 * np.sign(-lam_A) * 1), xytext=(xr, c), arrowprops=dict(arrowstyle="->", color="w", lw=2))
for s in p.sources: ax[0].plot(s.x * 1e3, s.y * 1e3, "w+", ms=10)
ax[0].set_title("A: level curve of Tmax touches the rail (tangency)"); ax[0].set_xlabel("fan x [mm]"); ax[0].set_ylabel("fan y [mm]"); ax[0].legend(fontsize=8, loc="upper right")

cs = np.array(curve)
ax[1].plot(cs[:, 0], cs[:, 1], "o-", label="best Tmax for each σ (position re-optimised)")
ax[1].plot(z_star[2], Tstar, "r*", ms=15, label=f"joint optimum σ = {z_star[2]:.1f} mm")
ax[1].axvline(12, color="grey", lw=.8, ls=":"); ax[1].text(12.3, cs[:, 1].max() * .995, "default σ = 12", fontsize=8, color="grey", va="top")
ax[1].set_xlabel("fan radius σ [mm]   (Δh = C / 2πσ²)"); ax[1].set_ylabel("Tmax [°C]"); ax[1].set_title("B: strong narrow jet vs weak wide jet, same total effect"); ax[1].legend(fontsize=8)

sg_g = np.linspace(5, 30, 26); dh_g = np.linspace(15, 260, 26)
TT = np.array([[peaks(z_star[0] * 1e-3, z_star[1] * 1e-3, s_ * 1e-3, d_, False).max() for d_ in dh_g] for s_ in sg_g])
cf = ax[2].contourf(sg_g, dh_g, TT.T, 22, cmap="magma_r"); plt.colorbar(cf, ax=ax[2], label="Tmax [°C]")
sc = np.linspace(5, 30, 200); ax[2].plot(sc, C_BUDGET / (2 * np.pi * (sc * 1e-3)**2), "w-", lw=2.2, label="constraint 2πσ²Δh = C")
ax[2].plot(z_star[2], dh_star, "c*", ms=15, mec="k", label="constrained optimum")
ax[2].set_ylim(15, 260); ax[2].set_xlabel("σ [mm]"); ax[2].set_ylabel("Δh [W/(m²·K)]"); ax[2].set_title("B: constraint curve vs level sets (fan fixed at optimum)"); ax[2].legend(fontsize=8)
plt.tight_layout(); RESULTS.mkdir(exist_ok=True); plt.savefig(RESULTS / "lagrange.png", dpi=115)
print("\nSaved", RESULTS / "lagrange.png")
