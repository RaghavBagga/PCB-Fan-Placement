"""
simulation.py  -  core model + fan-position scan for the PCB thermal project.

Model (SI units, T in deg C):
    k*t*Laplacian(T) - 2*h(x,y)*(T - Ta) + q(x,y) = 0,   insulated board edges
    h(x,y) = h_nat + dh*exp(-((x-xf)^2 + (y-yf)^2) / (2*sigma^2))

Use as a library (Exact_checks.py does this):   from simulation import *
Run as a script:
    python simulation.py            scan the whole board in 2 mm steps -> results/
    python simulation.py 5          same, 5 mm steps (faster)
    python simulation.py --test     quick self-test of the solver
"""
import csv, sys, time
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"


# ----------------------------------------------------------------- parameters
@dataclass
class Source:
    name: str
    x: float      # centre x [m]
    y: float      # centre y [m]
    w: float      # width  [m]
    d: float      # height [m]
    P: float      # power  [W]


@dataclass
class Params:
    Lx: float = 0.100          # board length [m]
    Ly: float = 0.080          # board width  [m]
    t: float = 1.6e-3          # thickness [m]
    k_eff: float = 15.0        # effective in-plane conductivity [W/(m K)]
    h_nat: float = 10.0        # natural convection [W/(m^2 K)]
    dh: float = 60.0           # peak extra convection under the fan [W/(m^2 K)]
    sigma: float = 0.012       # fan influence radius [m]
    Ta: float = 25.0           # ambient [deg C]
    T_limit: float = 60.0      # overheating threshold [deg C]
    dx: float = 2.0e-3         # grid spacing [m]
    sources: list = field(default_factory=lambda: [
        Source("CPU",       0.030, 0.050, 0.012, 0.012, 2.0),
        Source("Regulator", 0.070, 0.055, 0.008, 0.008, 1.0),
        Source("Driver",    0.065, 0.020, 0.010, 0.006, 0.8),
    ])


# ----------------------------------------------------------------- grid and fields
def grid(p):
    nx, ny = int(round(p.Lx / p.dx)), int(round(p.Ly / p.dx))
    X, Y = np.meshgrid((np.arange(nx) + .5) * p.dx, (np.arange(ny) + .5) * p.dx, indexing="ij")
    return nx, ny, X, Y


def source_field(p, X, Y):
    """Heat input per unit area q(x,y) [W/m^2]; each chip spreads its power P over its footprint."""
    q = np.zeros_like(X)
    for s in p.sources:
        m = (np.abs(X - s.x) <= s.w / 2) & (np.abs(Y - s.y) <= s.d / 2)
        if m.any():
            q[m] += s.P / (m.sum() * p.dx**2)
    return q


def h_field(p, X, Y, fan=None):
    h = np.full_like(X, p.h_nat)
    if fan is not None:
        h = h + p.dh * np.exp(-((X - fan[0])**2 + (Y - fan[1])**2) / (2 * p.sigma**2))
    return h


# ----------------------------------------------------------------- solver
@lru_cache(maxsize=8)
def _conduction(nx, ny, kt):
    """Conduction part of the matrix. It never changes when the fan moves, so build it once per grid."""
    idx = np.arange(nx * ny).reshape(nx, ny)
    a = np.concatenate([idx[:-1, :].ravel(), idx[:, :-1].ravel()])
    b = np.concatenate([idx[1:, :].ravel(), idx[:, 1:].ravel()])
    n = nx * ny
    off = sp.coo_matrix((np.full(2 * a.size, -kt), (np.r_[a, b], np.r_[b, a])), shape=(n, n))
    deg = -np.asarray(off.sum(axis=1)).ravel()                   # kt * (number of neighbours)
    A = (off + sp.diags(deg)).tocsc()
    A.sort_indices()
    cols = np.repeat(np.arange(n), np.diff(A.indptr))
    diag_pos = np.flatnonzero(A.indices == cols)                 # where each diagonal entry lives in A.data
    return A, diag_pos, deg


def solve(p, fan=None, q=None, h=None):
    """Steady-state temperature field T [deg C], shape (nx, ny). Only the diagonal changes between fan positions."""
    nx, ny, X, Y = grid(p)
    q = source_field(p, X, Y) if q is None else q
    h = h_field(p, X, Y, fan) if h is None else h
    A0, diag_pos, deg = _conduction(nx, ny, p.k_eff * p.t)
    loss = 2.0 * h.ravel() * p.dx**2
    data = A0.data.copy()
    data[diag_pos] = deg + loss
    A = sp.csc_matrix((data, A0.indices, A0.indptr), shape=A0.shape)
    rhs = q.ravel() * p.dx**2 + loss * p.Ta
    return spla.spsolve(A, rhs).reshape(nx, ny)


# ----------------------------------------------------------------- measurements
def peak_subcell(T):
    """Peak of T between cell centres: vertex of the parabola through the hottest cell and its neighbours."""
    i, j = np.unravel_index(np.argmax(T), T.shape)
    T0, add = T[i, j], 0.0
    if 0 < i < T.shape[0] - 1 and 0 < j < T.shape[1] - 1:
        for a, b in [(T[i-1, j], T[i+1, j]), (T[i, j-1], T[i, j+1])]:
            den = a - 2 * T0 + b
            if den < 0:
                add += -(b - a)**2 / (8 * den)
    return T0 + add


def metrics(p, T):
    i, j = np.unravel_index(np.argmax(T), T.shape)
    _, _, X, Y = grid(p)
    return dict(Tmax=float(T[i, j]), Tpeak=float(peak_subcell(T)),
                x_hot=float(X[i, j]), y_hot=float(Y[i, j]),
                area_over=float((T > p.T_limit).sum() * p.dx**2),     # [m^2]
                Tmean=float(T.mean()))


def heat_balance(p, T, fan=None):
    """Return (heat in, heat out) in watts. They must agree: this is the divergence theorem on the board."""
    _, _, X, Y = grid(p)
    inn = source_field(p, X, Y).sum() * p.dx**2
    out = (2 * h_field(p, X, Y, fan) * (T - p.Ta)).sum() * p.dx**2
    return inn, out


# ----------------------------------------------------------------- scan
def scan_fan(p, step=2e-3, limits=(60.0, 45.0)):
    """Solve for every fan position on a regular grid. Returns a list of dict rows (millimetres in the output)."""
    xs = np.arange(0, p.Lx + 1e-9, step)
    ys = np.arange(0, p.Ly + 1e-9, step)
    rows = []
    for xf in xs:
        for yf in ys:
            T = solve(p, (xf, yf))
            m = metrics(p, T)
            row = dict(x_fan_mm=xf * 1e3, y_fan_mm=yf * 1e3, Tmax_C=m["Tmax"], Tpeak_C=m["Tpeak"],
                       x_hot_mm=m["x_hot"] * 1e3, y_hot_mm=m["y_hot"] * 1e3, Tmean_C=m["Tmean"])
            for L in limits:
                row[f"area_over_{L:g}C_mm2"] = float((T > L).sum() * p.dx**2 * 1e6)
            rows.append(row)
    return xs, ys, rows


def save_scan(p, xs, ys, rows, base, limits=(60.0, 45.0), outdir=RESULTS):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    outdir.mkdir(exist_ok=True)
    fields = list(rows[0].keys())
    with open(outdir / "fan_scan.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: round(v, 4) for k, v in r.items()})
    a1, a2 = f"area_over_{limits[0]:g}C_mm2", f"area_over_{limits[1]:g}C_mm2"
    best = min(rows, key=lambda r: r["Tmax_C"])
    best_area = min(rows, key=lambda r: (r[a2], r["Tmax_C"]))
    worst = max(rows, key=lambda r: r["Tmax_C"])
    with open(outdir / "summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["case", "x_fan_mm", "y_fan_mm", "Tmax_C", a1, a2])
        w.writerow(["no fan", "", "", round(base["Tmax_C"], 3), round(base[a1], 1), round(base[a2], 1)])
        for name, r in [("best Tmax", best), (f"best area > {limits[1]:g} C", best_area), ("worst Tmax", worst)]:
            w.writerow([name, r["x_fan_mm"], r["y_fan_mm"], round(r["Tmax_C"], 3), round(r[a1], 1), round(r[a2], 1)])
    shape = (len(xs), len(ys))
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for a, key, title in [(ax[0], "Tmax_C", "Max temperature vs fan position [C]"),
                          (ax[1], a2, f"Area above {limits[1]:g} C vs fan position [mm²]")]:
        im = a.imshow(np.array([r[key] for r in rows]).reshape(shape).T, origin="lower",
                      extent=[0, p.Lx * 1e3, 0, p.Ly * 1e3], cmap="viridis")
        a.plot(best["x_fan_mm"], best["y_fan_mm"], "r*", ms=14)
        for s in p.sources:
            a.plot(s.x * 1e3, s.y * 1e3, "w+")
        a.set_title(title); a.set_xlabel("fan x [mm]"); a.set_ylabel("fan y [mm]")
        plt.colorbar(im, ax=a)
    plt.tight_layout(); plt.savefig(outdir / "fan_scan.png", dpi=110); plt.close(fig)
    return best, best_area, worst


# ----------------------------------------------------------------- self-test
def selftest():
    p = Params(sources=[], dx=2e-3)
    _, _, X, Y = grid(p)
    T = solve(p, q=np.full_like(X, 500.0), h=np.full_like(X, 10.0))
    print(f"uniform source test  : max |T - exact| = {np.abs(T - (p.Ta + 25.0)).max():.1e}   (expect ~1e-12)")
    p = Params()
    for fan in (None, (0.03, 0.05), (0.08, 0.02)):
        i, o = heat_balance(p, solve(p, fan), fan)
        print(f"energy balance fan={fan}: in {i:.4f} W, out {o:.4f} W, relative error {abs(o - i) / i:.1e}")


def main():
    if "--test" in sys.argv:
        selftest(); return
    step = (float(sys.argv[1]) if len(sys.argv) > 1 else 2.0) * 1e-3
    p, limits = Params(), (60.0, 45.0)
    T0 = solve(p); m0 = metrics(p, T0)
    base = dict(Tmax_C=m0["Tmax"], **{f"area_over_{L:g}C_mm2": float((T0 > L).sum() * p.dx**2 * 1e6) for L in limits})
    print(f"No fan: Tmax = {m0['Tmax']:.2f} C at ({m0['x_hot']*1e3:.0f}, {m0['y_hot']*1e3:.0f}) mm")
    t0 = time.time()
    xs, ys, rows = scan_fan(p, step, limits)
    print(f"Scanned {len(rows)} fan positions in {time.time() - t0:.1f} s")
    best, best_area, worst = save_scan(p, xs, ys, rows, base, limits)
    a2 = f"area_over_{limits[1]:g}C_mm2"
    for name, r in [("best Tmax", best), ("best area", best_area), ("worst", worst)]:
        print(f"{name:10s} fan = ({r['x_fan_mm']:.0f}, {r['y_fan_mm']:.0f}) mm   Tmax = {r['Tmax_C']:.2f} C   area > {limits[1]:g} C = {r[a2]:.0f} mm²")
    print("Saved to", RESULTS)


if __name__ == "__main__":
    main()
