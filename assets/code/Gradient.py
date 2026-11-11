"""
Gradient.py  -  gradient of the peak temperature with respect to the fan position.

Board:  A(h) theta = b,  theta = T - Ta,  A = A_cond + diag(2 h dx^2),  b = q dx^2 (independent of the fan).
Differentiate with respect to the fan coordinate xf (only h depends on it):
        A dtheta/dxf = -(dA/dxf) theta = -diag(2 dx^2 dh/dxf) theta
For the temperature at one cell i, A is symmetric, so with  A lam = e_i  (one extra solve)
        dT_i/dxf = - sum_m lam_m * 2 dx^2 * (dh_m/dxf) * theta_m ,   same for yf.
dh/dxf = dh*exp(-r^2/2s^2) * (x - xf)/s^2 .

Each chip has a peak f_k(xf,yf) (hottest cell under the chip). Tmax = max_k f_k is smooth except where two peaks
tie. At the optimum two peaks tie and their gradients point in opposite directions (KKT / Lagrange condition).
"""
import itertools
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from simulation import *
from simulation import _conduction

p = Params()
nx, ny, X, Y = grid(p)
q = source_field(p, X, Y)
SIG2 = p.sigma**2


def footprint_cells(s):
    m = (np.abs(X - s.x) <= s.w / 2) & (np.abs(Y - s.y) <= s.d / 2)
    return np.flatnonzero(m.ravel())
CELLS = [footprint_cells(s) for s in p.sources]


def peaks_and_grads(fan, want_grad=True):
    """Return Tmax-per-chip values [C] and gradients [K/mm] (array n_chips x 2) at fan position (m)."""
    xf, yf = fan
    g = p.dh * np.exp(-((X - xf)**2 + (Y - yf)**2) / (2 * SIG2))
    h = p.h_nat + g
    A0, dpos, deg = _conduction(nx, ny, p.k_eff * p.t)
    loss = 2 * h.ravel() * p.dx**2
    data = A0.data.copy(); data[dpos] = deg + loss
    A = sp.csc_matrix((data, A0.indices, A0.indptr), shape=A0.shape)
    lu = spla.splu(A)
    theta = lu.solve(q.ravel() * p.dx**2)                       # theta = T - Ta   (b = q dx^2)
    vals, grads = [], []
    dhx = (g * (X - xf) / SIG2).ravel(); dhy = (g * (Y - yf) / SIG2).ravel()
    for cells in CELLS:
        i = cells[np.argmax(theta[cells])]
        vals.append(theta[i] + p.Ta)
        if want_grad:
            e = np.zeros(nx * ny); e[i] = 1.0
            lam = lu.solve(e)
            w = lam * 2 * p.dx**2 * theta
            grads.append([-np.dot(w, dhx) * 1e-3, -np.dot(w, dhy) * 1e-3])     # K/m -> K/mm
    return np.array(vals), np.array(grads)


def Tmax_at(fan):
    return float(solve(p, fan).max())


# ------------------------------------------------------------- 1. check against finite differences
print("=== 1. Analytic gradient vs central finite differences (K per mm) ===")
print("  fan (mm)    chip        analytic d/dx   FD d/dx      analytic d/dy   FD d/dy      rel. error")
worst = 0
for fan_mm in [(38, 50), (30, 45), (60, 30), (72, 58), (20, 20)]:
    fan = np.array(fan_mm) * 1e-3
    v, G = peaks_and_grads(fan)
    for k, s in enumerate(p.sources):
        fd = []
        for ax in (0, 1):
            e = np.zeros(2); e[ax] = 1e-5
            vp, _ = peaks_and_grads(fan + e, False); vm, _ = peaks_and_grads(fan - e, False)
            fd.append((vp[k] - vm[k]) / 2e-5 * 1e-3)
        err = np.hypot(*(G[k] - fd)) / max(np.hypot(*fd), 1e-12); worst = max(worst, err)
        print(f"  {str(fan_mm):10s} {s.name:10s} {G[k,0]:+13.6f} {fd[0]:+12.6f}   {G[k,1]:+13.6f} {fd[1]:+12.6f}   {err:.1e}")
print(f"  worst relative error: {worst:.1e}")


# ------------------------------------------------------------- 2. minimax descent (min-norm element of the active gradients)
def min_norm_hull(G):
    """Point of minimum norm in the convex hull of the rows of G (exact for a few rows: try every support)."""
    best, bl = None, None
    for r in range(1, len(G) + 1):
        for S in itertools.combinations(range(len(G)), r):
            Gs = G[list(S)]; M = np.zeros((r + 1, r + 1)); M[:r, :r] = Gs @ Gs.T; M[:r, r] = 1; M[r, :r] = 1
            rhs = np.zeros(r + 1); rhs[r] = 1
            sol = np.linalg.lstsq(M, rhs, rcond=None)[0]; lam = sol[:r]
            if np.all(lam >= -1e-12):
                d = lam @ Gs
                if best is None or np.linalg.norm(d) < np.linalg.norm(best) - 1e-15:
                    best, bl = d, (S, lam)
    return best, bl


def descend(start_mm, eps=0.03, maxit=400):
    x = np.array(start_mm, float) * 1e-3; path = [x.copy()]
    for it in range(maxit):
        v, G = peaks_and_grads(x); Tm = v.max()
        act = np.flatnonzero(v >= Tm - eps)
        d, info = min_norm_hull(G[act]); gn = np.linalg.norm(d)       # K/mm
        if gn < 2e-3: break
        dirn = -d / gn; s = 4e-3
        while s > 1e-7:
            xn = np.clip(x + s * dirn, [0, 0], [p.Lx, p.Ly])
            if Tmax_at(xn) < Tm - 1e-4 * (s * 1e3) * gn: break
            s /= 2
        if s <= 1e-7: break
        x = xn; path.append(x.copy())
    return np.array(path) * 1e3, it


print("\n=== 2. Minimax gradient descent from different starting points ===")
starts = [(10, 10), (90, 70), (60, 40), (20, 70), (85, 15), (50, 5), (5, 45)]
paths, finals = [], []
for st in starts:
    P, its = descend(st); paths.append(P)
    f = P[-1] * 1e-3; finals.append(f)
    print(f"  start {str(st):9s} -> ({P[-1,0]:5.2f}, {P[-1,1]:5.2f}) mm after {len(P)-1:3d} steps,  Tmax = {Tmax_at(f):.3f} C")

# ------------------------------------------------------------- 3. KKT condition at the optimum
opt = finals[0]
v, G = peaks_and_grads(opt)
order = np.argsort(-v)
print("\n=== 3. Optimality (KKT) condition at the optimum ===")
print(f"  optimum fan position: ({opt[0]*1e3:.2f}, {opt[1]*1e3:.2f}) mm")
for k in order:
    print(f"  peak under {p.sources[k].name:10s}: {v[k]:.3f} C   gradient = ({G[k,0]:+.4f}, {G[k,1]:+.4f}) K/mm   |g| = {np.linalg.norm(G[k]):.4f}")
a, b = order[0], order[1]
cosang = G[a] @ G[b] / (np.linalg.norm(G[a]) * np.linalg.norm(G[b]))
lam = (G[b] @ (G[b] - G[a])) / np.sum((G[a] - G[b])**2)
print(f"  two active peaks: {p.sources[a].name} and {p.sources[b].name}; value difference {abs(v[a]-v[b]):.4f} C")
print(f"  angle between their gradients: {np.degrees(np.arccos(np.clip(cosang, -1, 1))):.1f} degrees (180 = exactly opposite)")
print(f"  multiplier: {lam:.3f}*grad({p.sources[a].name}) + {1-lam:.3f}*grad({p.sources[b].name}) = "
      f"({lam*G[a,0]+(1-lam)*G[b,0]:+.5f}, {lam*G[a,1]+(1-lam)*G[b,1]:+.5f}) K/mm  (should be ~0)")

# ------------------------------------------------------------- 4. figure
xs = np.arange(0, 101, 5) * 1e-3; ys = np.arange(0, 81, 5) * 1e-3
TM = np.zeros((len(xs), len(ys))); U = np.zeros_like(TM); V = np.zeros_like(TM)
for i, xf in enumerate(xs):
    for j, yf in enumerate(ys):
        v_, G_ = peaks_and_grads((xf, yf)); k = np.argmax(v_)
        TM[i, j] = v_[k]; U[i, j], V[i, j] = -G_[k]                # descent direction
fig, ax = plt.subplots(1, 3, figsize=(18, 4.9))
im = ax[0].imshow(TM.T, origin="lower", extent=[-2.5, 102.5, -2.5, 82.5], cmap="viridis", interpolation="bilinear")
mag = np.hypot(U, V); ax[0].quiver(xs * 1e3, ys * 1e3, (U / (mag + 1e-12)).T, (V / (mag + 1e-12)).T, color="w", scale=30, width=.0035)
for P in paths: ax[0].plot(P[:, 0], P[:, 1], "-", color="#ff6b3d", lw=1.4); ax[0].plot(P[0, 0], P[0, 1], "o", color="#ff6b3d", ms=4)
ax[0].plot(opt[0]*1e3, opt[1]*1e3, "r*", ms=15, mec="w")
for s in p.sources: ax[0].plot(s.x*1e3, s.y*1e3, "w+", ms=9)
ax[0].set_title("−∇(Tmax) over fan positions, with descent paths"); ax[0].set_xlabel("fan x [mm]"); ax[0].set_ylabel("fan y [mm]")
plt.colorbar(im, ax=ax[0], label="Tmax [°C]")

xline = np.arange(15, 65.1, 0.5) * 1e-3
F = np.array([peaks_and_grads((x_, opt[1]), False)[0] for x_ in xline])
for k, s in enumerate(p.sources): ax[1].plot(xline * 1e3, F[:, k], label=f"peak under {s.name}")
ax[1].plot(xline * 1e3, F.max(axis=1), "k--", lw=1.2, label="Tmax = max of peaks")
ax[1].axvline(opt[0] * 1e3, color="r", lw=.8); ax[1].set_xlabel(f"fan x [mm]  (y = {opt[1]*1e3:.1f} mm)"); ax[1].set_ylabel("temperature [°C]")
ax[1].set_title("The optimum is where two peaks tie"); ax[1].legend(fontsize=8)

dd = []
for fan_mm in [(a_, b_) for a_ in (10, 30, 50, 70, 90) for b_ in (10, 30, 50, 70)]:
    fan = np.array(fan_mm) * 1e-3; _, G_ = peaks_and_grads(fan); v0, _ = peaks_and_grads(fan, False)
    for k in range(3):
        for ax_i in (0, 1):
            e = np.zeros(2); e[ax_i] = 1e-5
            fd = (peaks_and_grads(fan + e, False)[0][k] - peaks_and_grads(fan - e, False)[0][k]) / 2e-5 * 1e-3
            dd.append((G_[k, ax_i], fd))
dd = np.array(dd); lim = np.abs(dd).max() * 1.05
ax[2].plot([-lim, lim], [-lim, lim], "k-", lw=.8); ax[2].plot(dd[:, 0], dd[:, 1], "o", ms=4)
ax[2].set_xlabel("analytic derivative [K/mm]"); ax[2].set_ylabel("finite-difference derivative [K/mm]"); ax[2].set_title("Gradient check (120 derivatives)")
ax[2].set_aspect("equal")
plt.tight_layout(); RESULTS.mkdir(exist_ok=True); plt.savefig(RESULTS / "gradient_field.png", dpi=115)
print("\nSaved", RESULTS / "gradient_field.png")
