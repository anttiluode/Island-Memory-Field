#!/usr/bin/env python3
"""
rh_blind_zone.py
================
Reproduces every claim in blind_zone_of_pointwise_convexity.md:

  1. The xi-only identity  G(t) = Re(xi'/xi)'(1/2+it) = Hadamard zero sum,
     verified against direct high-precision computation of xi.
  2. The detection radius of the pointwise convexity criterion:
     d* ~ (one local mean zero gap) ~ 2*pi/log T.
  3. The blind-zone counterexample at T = 1e9 (off-line zero at beta = 0.92,
     allowed by the classical zero-free region, invisible to convexity at
     every t) and the demonstration that the Lagarias global criterion
     Re xi'/xi > 0 on sigma > 1/2 catches it.
  4. The figure.

pip install numpy mpmath matplotlib
Antti Luode (PerceptionLab, Helsinki) + Claude (Anthropic), June 2026.
Do not hype. Do not lie. Just show.
"""
import numpy as np
import mpmath as mp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

mp.mp.dps = 25

# ----------------------------------------------------------------------
# 1. Identity check: G from xi directly vs G from the zero multiset
# ----------------------------------------------------------------------

def xi(s):
    s = mp.mpc(s)
    return 0.5*s*(s-1)*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)

def G_direct(t):
    s0 = mp.mpc(0.5, t)
    x0 = xi(s0); x1 = mp.diff(xi, s0); x2 = mp.diff(xi, s0, 2)
    return float(mp.re((x2*x0 - x1*x1)/(x0*x0)))

print("="*72)
print("1. IDENTITY: G(t) is a function of the zero multiset alone")
print("="*72)
print("fetching first 200 zeros (takes ~a minute) ...")
gam = np.array([float(mp.zetazero(n).imag) for n in range(1, 201)])

def density_tail(t, X):
    f = lambda x: (mp.log(x/(2*mp.pi))/(2*mp.pi)) * (1/(x-t)**2 + 1/(x+t)**2)
    return float(mp.quad(f, [X, mp.inf]))

def G_zero_sum(t):
    X = float(gam[-1]) + 0.5*float(gam[-1]-gam[-2])
    return float(np.sum(1.0/(t-gam)**2 + 1.0/(t+gam)**2)) + density_tail(t, X)

for t in [5.0, 10.0, 14.1347, 17.578, 25.0, 33.0, 49.7]:
    gd, gz = G_direct(t), G_zero_sum(t)
    print(f"  t={t:8.4f}   direct={gd:+.6e}   zero-sum={gz:+.6e}   reldiff={abs(gd-gz)/abs(gd):.2e}")

t = 21.5; eps = mp.mpf('1e-6')
F = lambda sig: abs(xi(mp.mpc(sig, t)))**2
d2F = float((F(mp.mpf('0.5')+eps) - 2*F(mp.mpf('0.5')) + F(mp.mpf('0.5')-eps))/eps**2)
rhs = 2*float(abs(xi(mp.mpc(0.5,t)))**2)*G_direct(t)
print(f"  lemma F''(1/2)=2|xi|^2 G at t=21.5: reldiff = {abs(d2F-rhs)/abs(rhs):.2e}")

# ----------------------------------------------------------------------
# 2. Detection radius
# ----------------------------------------------------------------------

print()
print("="*72)
print("2. DETECTION RADIUS of pointwise convexity")
print("="*72)

print("\n  A. Low height (REAL zeros), quadruple at mid-gap gamma0 = 17.578:")
t0 = 0.5*(14.134725 + 21.022040)
S = G_zero_sum(t0)          # all 200 zeros are on the line -> this IS S_on
print(f"     S_on = {S:.4f},  d* = sqrt(2/S_on) = {np.sqrt(2/S):.2f}  (> strip half-width 0.5)")
for d in [0.499, 0.25, 0.05]:
    Q = -2/d**2 + 2*((2*t0)**2-d**2)/(((2*t0)**2+d**2)**2)
    print(f"     d = {d:5.3f}: G(gamma0) = {S+Q:+9.3f}  -> {'DETECTED' if S+Q<0 else 'invisible'}")
print("     => at low heights every off-line zero is detected (why April looked clean)")

print("\n  B. High height (synthetic on-line zeros at correct local density):")
def S_on_synth(delta, conserve, K=200000):
    k = np.arange(1 if conserve else 0, K) + 0.5
    return np.sum(2.0/((k*delta)**2)) + 2.0/(delta*(K*delta))

R_ZF = 5.5587
table = []
print(f"     {'T':>8} {'delta':>7} {'d* (cons)':>10} {'zero-free d_max':>16}   blind zone")
for T in [1e3, 1e6, 1e7, 1e9, 1e12, 1e15, 1e18]:
    delta = 2*np.pi/np.log(T/(2*np.pi))
    dstar = min(np.sqrt(2/S_on_synth(delta, True)), 0.5)
    dmax  = 0.5 - 1/(R_ZF*np.log(T))
    table.append((T, delta, dstar, dmax))
    bz = f"d in ({dstar:.3f}, {dmax:.3f})" if dmax > dstar else "none"
    print(f"     {T:8.0e} {delta:7.4f} {dstar:10.4f} {dmax:16.4f}   {bz}")
print(f"     theory: d* = 1.034*delta (count conserved), 0.450*delta (not)")

# ----------------------------------------------------------------------
# 3. The counterexample at T = 1e9 and the Lagarias repair
# ----------------------------------------------------------------------

print()
print("="*72)
print("3. COUNTEREXAMPLE (T=1e9, beta=0.92) and the Lagarias island")
print("="*72)

T = 1e9
delta = 2*np.pi/np.log(T/(2*np.pi))
K = 200000
off = (np.arange(1, K) + 0.5)*delta          # conserved count

def G_synth(t_rel, d):
    v = np.concatenate([t_rel - off, t_rel + off])
    s = np.sum(1.0/v**2) + 2.0/(delta*(K*delta))
    return s + 2*(t_rel**2 - d**2)/((t_rel**2 + d**2)**2)

ts = np.linspace(-1.5*delta, 1.5*delta, 4001)[1:-1]   # avoid exact zero hits
G_blind   = np.array([G_synth(t, 0.42) for t in ts])
G_control = np.array([G_synth(t, 0.25) for t in ts])
print(f"  d = 0.42 (beta = 0.92, allowed by zero-free region): min G = {G_blind.min():+.3f}")
print(f"     -> convexity holds at EVERY t. RH false, criterion silent.")
print(f"  d = 0.25 control:                                    min G = {G_control.min():+.3f}  (detected)")

def Phi_synth(sigma, d, t_rel=0.0):
    u = sigma - 0.5
    s = np.sum(u/(u**2 + (t_rel-off)**2)) + np.sum(u/(u**2 + (t_rel+off)**2))
    s += (sigma-(0.5+d))/((sigma-(0.5+d))**2 + t_rel**2)
    s += (sigma-(0.5-d))/((sigma-(0.5-d))**2 + t_rel**2)
    return s

sigs = np.linspace(0.505, 0.918, 400)
phis = np.array([Phi_synth(s, 0.42) for s in sigs])
print(f"  Lagarias Phi(sigma, gamma0): Phi(0.55) = {Phi_synth(0.55,0.42):+.3f} (boundary blind), "
      f"min over sigma = {phis.min():+.1f} (global criterion catches it)")

# ----------------------------------------------------------------------
# 4. Figure
# ----------------------------------------------------------------------

fig, ax = plt.subplots(1, 3, figsize=(17, 4.6))
fig.suptitle('The blind zone of pointwise log-convexity — stress test of the RH criterion (ii)\n'
             'Antti Luode (PerceptionLab) + Claude (Anthropic), June 2026', fontsize=11)

Ts  = [r[0] for r in table]; dc = [r[2] for r in table]; dz = [r[3] for r in table]
ax[0].plot(Ts, dc, 'o-', color='#1a1a2e', lw=2, label='convexity detection radius d*')
ax[0].plot(Ts, dz, 's--', color='#c0392b', lw=1.5, label='classical zero-free bound on d')
ax[0].axhline(0.5, color='gray', lw=1, ls=':', label='strip half-width')
ax[0].fill_between(Ts, dc, np.minimum(dz, 0.5), where=np.array(dz) > np.array(dc),
                   color='#f39c12', alpha=0.35, label='BLIND ZONE')
ax[0].set_xscale('log'); ax[0].set_xlabel('height T'); ax[0].set_ylabel('d = |β − 1/2|')
ax[0].set_title('Convexity sees d ≲ δ(T) ~ 2π/log T;\nzero-free regions see d ≈ 1/2. Middle: nobody.')
ax[0].legend(fontsize=7, loc='center left'); ax[0].set_ylim(0, 0.55); ax[0].grid(alpha=0.3)

ax[1].plot(ts, G_blind, color='#1a1a2e', lw=1.8, label='d = 0.42 (blind: G > 0 everywhere)')
ax[1].plot(ts, G_control, color='#c0392b', lw=1.8, label='d = 0.25 (detected: G < 0)')
ax[1].axhline(0, color='gray', lw=1); ax[1].set_ylim(-20, 60)
ax[1].set_xlabel('t − γ₀'); ax[1].set_ylabel('G(t)')
ax[1].set_title('T = 10⁹: convexity holds at every t\ndespite a zero at β = 0.92')
ax[1].legend(fontsize=7); ax[1].grid(alpha=0.3)

ax[2].plot(sigs, phis, color='#1a1a2e', lw=1.8)
ax[2].axhline(0, color='gray', lw=1)
ax[2].axvline(0.92, color='#c0392b', ls='--', lw=1, label='off-line zero β = 0.92')
ax[2].set_ylim(-60, 10); ax[2].set_xlabel('σ'); ax[2].set_ylabel("Φ(σ, γ₀) = Re ξ'/ξ")
ax[2].set_title('Same configuration, Lagarias criterion:\nΦ > 0 near the line, negative island deeper in')
ax[2].legend(fontsize=7); ax[2].grid(alpha=0.3)

plt.tight_layout(rect=[0, 0, 1, 0.93])
plt.savefig('blind_zone.png', dpi=140)
print("\n[saved] blind_zone.png")
