"""
Taylor.py  -  Taylor series of the Gaussian fan profile.

h(r) = h_nat + dh * G(r),   G = exp(-u),  u = r^2 / (2 sigma^2).
Taylor polynomial in u of degree m (order 2m in r):  P_m(u) = sum_{k=0}^{m} (-u)^k / k!.
Alternating-series bound:  |exp(-u) - P_m(u)| <= u^(m+1) / (m+1)!.
Questions:  (1) how far from the fan centre is the polynomial a good model?
            (2) how many terms does the series need at distance r?
            (3) does a truncated Taylor model of the fan change the board temperature?
"""
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from simulation import *

p = Params()
FAN = (0.0377, 0.05025)
nx, ny, X, Y = grid(p)

def P(u, m):
    return sum((-u)**k / math.factorial(k) for k in range(m + 1))

orders = [(1, "order 2 (paraboloid)"), (2, "order 4"), (3, "order 6"), (4, "order 8")]
s = np.linspace(0, 4, 801)                           # r / sigma
U = s**2 / 2

# ---------------------------------------------------------------- 1. radius of validity
print("=== 1. Where is the Taylor polynomial of the fan profile accurate? (error as a fraction of the peak dh) ===")
print("  polynomial        r/sigma where |error| first exceeds 1%    ... 5%     zero crossing of P")
for m, name in orders:
    err = np.abs(P(U, m) - np.exp(-U))
    r1 = s[np.argmax(err > 0.01)]; r5 = s[np.argmax(err > 0.05)]
    zc = s[np.argmax(P(U, m) < 0)] if (P(U, m) < 0).any() else float("nan")
    print(f"  {name:20s}  {r1:8.2f}                              {r5:6.2f}      {zc:6.2f}")
print("  (r/sigma = 1 is one fan radius; the board-scale fan has sigma = 12 mm)")

# ---------------------------------------------------------------- 2. terms needed
print("\n=== 2. Number of terms of the series needed at distance r (error < 1e-3 of the peak) ===")
need = []
for ss in s[1:]:
    u = ss**2 / 2; m = 0
    while abs(P(u, m) - math.exp(-u)) >= 1e-3: m += 1
    need.append(m + 1)
need = np.array(need)
for rs in (0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4):
    print(f"  r = {rs:3.1f} sigma ({rs*12:5.1f} mm):  {need[np.argmin(np.abs(s[1:] - rs))]:2d} terms")
print("  The series converges for every r (exp is entire), but the number of terms grows like u = r^2/(2 sigma^2).")

# ---------------------------------------------------------------- 3. effect of the fan: integral and temperature
print("\n=== 3. How much of the fan's total effect does the truncated polynomial capture? ===")
print("  exact:      double integral of dh*exp(-r^2/2s^2) dA over the plane = 2*pi*s^2*dh")
rr = np.linspace(0, 40, 400001)
exact_I = 1.0                                         # in units of 2 pi s^2 dh
print("  paraboloid on its own support (r <= sqrt(2) s):  integral = pi s^2 dh  ->  exactly 50% of the fan effect")
Ic = {}
for m, name in orders:
    u = rr**2 / 2
    f = np.maximum(P(u, m), 0.0)
    keep = np.abs(P(u, m) - np.exp(-u)) <= 0.02       # use the polynomial only where it is within 2% of dh
    cut = rr[np.argmin(keep)]
    f = np.where(rr <= cut, f, 0.0)
    I = np.trapezoid(f * rr, rr) / 1.0                # integral of G * r dr  (sigma = 1) ; exact = 1
    Ic[m] = (cut, I)
    print(f"  {name:20s} used for r <= {cut:4.2f} s : captures {100*I:5.1f}% of the exact fan effect")

print("\n=== 4. Board temperature with the Taylor-truncated fan (polynomial used only for r <= its 2% radius, zero beyond) ===")
def solve_with(m):
    r2 = (X - FAN[0])**2 + (Y - FAN[1])**2
    u = r2 / (2 * p.sigma**2)
    if m is None:
        g = np.exp(-u)
    else:
        cut = Ic[m][0]
        g = np.where(r2 <= (cut * p.sigma)**2, np.maximum(P(u, m), 0.0), 0.0)
    h = p.h_nat + p.dh * g
    return solve(p, h=h)
Tex = solve_with(None); mex = metrics(p, Tex)
print(f"  exact Gaussian: Tmax = {mex['Tmax']:.3f} C, mean = {mex['Tmean']:.3f} C")
res = []
for m, name in orders:
    T = solve_with(m); mm = metrics(p, T)
    res.append((name, mm["Tmax"], mm["Tmean"]))
    print(f"  {name:20s} Tmax = {mm['Tmax']:.3f} C ({mm['Tmax']-mex['Tmax']:+.3f}),  mean = {mm['Tmean']:.3f} C ({mm['Tmean']-mex['Tmean']:+.3f})")
Tno = metrics(p, solve(p))
print(f"  no fan:         Tmax = {Tno['Tmax']:.3f} C   -> the paraboloid alone recovers {(Tno['Tmax']-res[0][1])/(Tno['Tmax']-mex['Tmax'])*100:.0f}% of the cooling")

# ---------------------------------------------------------------- figure
fig, ax = plt.subplots(1, 3, figsize=(18, 5))
ax[0].plot(s, np.exp(-U), "k", lw=3, label="exact  exp(−r²/2σ²)")
for (m, name), c in zip(orders, ["#4c78a8", "#f58518", "#54a24b", "#b279a2"]):
    yv = P(U, m); yv = np.where((yv > -1.3) & (yv < 2.3), yv, np.nan); ax[0].plot(s, yv, color=c, lw=1.6, label=name)
ax[0].axhline(0, color="grey", lw=.6); ax[0].set_ylim(-1.3, 2.3); ax[0].set_xlim(0, 4)
ax[0].set_xlabel("distance from fan centre r / σ"); ax[0].set_ylabel("fan profile  (h − h_nat) / Δh"); ax[0].set_title("Taylor polynomials of the Gaussian fan profile"); ax[0].legend(fontsize=8)

for (m, name), c in zip(orders, ["#4c78a8", "#f58518", "#54a24b", "#b279a2"]):
    err = np.abs(P(U, m) - np.exp(-U)); ax[1].semilogy(s, np.maximum(err, 1e-12), color=c, lw=1.6, label=name)
    ax[1].semilogy(s, U**(m + 1) / math.factorial(m + 1), color=c, lw=1, ls=":")
ax[1].axhline(0.01, color="grey", lw=.7, ls="--"); ax[1].text(0.05, 0.012, "1% of peak", fontsize=8, color="grey")
ax[1].set_ylim(1e-8, 1e3); ax[1].set_xlim(0, 4); ax[1].set_xlabel("r / σ"); ax[1].set_ylabel("absolute error (fraction of Δh)")
ax[1].set_title("Error grows with distance; dotted = bound uᵐ⁺¹/(m+1)!"); ax[1].legend(fontsize=8)

ax[2].step(s[1:], need, where="post", color="k", lw=1.8)
ax[2].set_xlabel("r / σ"); ax[2].set_ylabel("terms needed for error < 10⁻³"); ax[2].set_title("The series always converges, but needs more terms far away")
ax[2].set_xlim(0, 4)
plt.tight_layout(); RESULTS.mkdir(exist_ok=True); plt.savefig(RESULTS / "taylor.png", dpi=115)
print("\nSaved", RESULTS / "taylor.png")
