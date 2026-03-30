#!/usr/bin/env python3
"""
GPU MG Poisson performance model as a function of CPU-GPU bandwidth
and OMP threads per GPU.

Measured data: Cosmo1024 test, 3 coarse steps from restart.
  n4 r=8 H100 NVL: CPU MG AMR = 1093.6 s, GPU MG AMR = 665.1 s

Two-stage model:
  (a,b) Bandwidth model at fixed r=8:
     T_mg(B) = T_cpu_fixed + T_gpu_kernel + V_pcie / B
  (c)  r-dependent model calibrated from measurements:
     T_gpu(r) = a * T_cpu(r) + b(card)
     T_cpu(r) = T_s * (1/r + alpha)   [empirical OMP scaling]
     S(r)     = T_cpu(r) / T_gpu(r)

  Card-dependent b decomposition:
     b(card) = T_kernel(card) + V_pcie / B_pcie(card)
     T_kernel ∝ 1/N_SM  (latency-bound AMR stencil, SM-count scaling)
  Calibrated from H100 & A40 data at n4 (r=4, r=8).
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ================================================================
# BANDWIDTH MODEL (panels a, b) — fixed r=8, n4
# ================================================================
r_bench = 8
T_cpu_mg = 1093.6   # CPU MG AMR (s) at r=8, n4

# Original bandwidth-model parameters (for panels a, b only)
f_cpu_fixed = 0.35
gpu_kernel_speedup = 4.0
T_cpu_fixed  = f_cpu_fixed * T_cpu_mg          # 382.8 s
T_smooth_cpu = (1 - f_cpu_fixed) * T_cpu_mg    # 710.8 s
T_gpu_kernel = T_smooth_cpu / gpu_kernel_speedup  # 177.7 s
B_h100_measured = 22.0

V_pcie_bw = (643.4 - T_cpu_fixed - T_gpu_kernel) * B_h100_measured  # 1824 GB

def T_mg_gpu_B(B):
    return T_cpu_fixed + T_gpu_kernel + V_pcie_bw / B

max_speedup_bw = T_cpu_mg / (T_cpu_fixed + T_gpu_kernel)

# ================================================================
# IMPROVED r-DEPENDENT MODEL (panel c) — calibrated from H100 & A40
# ================================================================
# T_gpu(r) = a * T_cpu(r) + b
# Measured n4 data:
#   H100: (r=4, Tc=1925.75, Tg=963.68), (r=8, Tc=1093.61, Tg=665.08)
#   A40:  (r=4, Tc=1924.14, Tg=1048.84), (r=8, Tc=1094.72, Tg=748.06)

# CPU-fixed fraction (slope of T_gpu vs T_cpu):
a_fixed = (963.683 - 665.075) / (1925.750 - 1093.609)  # 0.3588

# OMP scaling: T_cpu(r) = T_s * (1/r + alpha)
_alpha = (1925.750/1093.609 * 0.125 - 0.25) / (1.0 - 1925.750/1093.609)
_T_serial = 1093.609 / (0.125 + _alpha)

# b values from direct fit:
b_h100_fit = 665.075 - a_fixed * 1093.609   # 272.6
b_a40_fit  = 748.058 - a_fixed * 1094.720   # 355.2

# Decompose b = Tk + V4/B using SM scaling (Tk ∝ 1/N_SM)
# and cross-card constraint:
#   Tk_h + V4/22 = 272.6
#   1.571*Tk_h + V4/B_a = 355.2    (SM ratio 132/84 = 1.571)
#   dV/22 = 226.8                  (from n2-n4 H100 data)
# Solution: Tk_h=176.0, V4=2125, B_a40=27.0
N_SM_h100 = 132
Tk_h100 = 176.0
V4_pcie = 2125.0
B_a40_eff = 27.0  # fitted from A40 data

def T_cpu_r(r):
    return _T_serial * (1.0/r + _alpha)

def b_for_card(n_sm, B_pcie):
    Tk = Tk_h100 * (N_SM_h100 / n_sm)
    return Tk + V4_pcie / B_pcie

def S_r_model(r, b_val):
    tc = T_cpu_r(r)
    return tc / (a_fixed * tc + b_val)

# GPU card specs: (name, B_pcie, N_SM, color, ls, lw)
gpu_cards = [
    ('V100',      12,   80, 'grey',      '--', 1.2),
    ('A40',       B_a40_eff,   84, '#8c564b',   '-',  1.5),
    ('A100',      25,  108, '#2ca02c',   '--', 1.2),
    ('H100 NVL',  22,  132, '#d62728',   '-',  2.0),
    ('GH200',    450,  132, '#9467bd',   '--', 1.2),
]

# Verify
print(f"r-model: a={a_fixed:.4f}, T_s={_T_serial:.0f}, alpha={_alpha:.4f}")
print(f"Tk(H100)={Tk_h100:.1f}s, V4={V4_pcie:.0f}GB, B(A40)={B_a40_eff:.1f}GB/s")
print(f"Check b(H100)={b_for_card(132,22):.1f} (fit: {b_h100_fit:.1f})")
print(f"Check b(A40) ={b_for_card(84,B_a40_eff):.1f} (fit: {b_a40_fit:.1f})")

# ================================================================
# GPU/interconnect specs for panels (a), (b)
# ================================================================
gpus_ab = [
    ('PCIe Gen3',           12,  'grey',    's', 9,  False,
     ( 0, 14, 'center'),  ( 0,-16, 'left')),
    ('H100 NVL\n(measured)',22,  '#d62728', '*', 14, True,
     ( 0,-16, 'center'),  ( 0,-44, 'left')),
    ('PCIe Gen4\n(A40, A100)',    25,  '#2ca02c', 'D', 10, False,
     ( 8, 12, 'left'),    ( 0,-16, 'left')),
    ('PCIe Gen5\n(theoretical)',  50,  '#1f77b4', 'D', 10, False,
     ( 0,-16, 'center'),  ( 0,-16, 'left')),
    ('GH200\n(NVLink-C2C)', 450, '#9467bd', 'p', 12, False,
     (-8, 12, 'right'),   ( 0, 12, 'right')),
]

# ================================================================
# FIGURE
# ================================================================
B_range = np.logspace(np.log10(5), np.log10(600), 500)
T_model_B = T_mg_gpu_B(B_range)
speedup_B = T_cpu_mg / T_model_B

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 5.5))

# --- (a) MG AMR wall time vs bandwidth ---
ax1.plot(B_range, T_model_B, 'k-', lw=2, label='GPU MG model')
ax1.axhline(T_cpu_mg, color='steelblue', ls='--', lw=1.5,
            label=f'CPU-only ({T_cpu_mg:.0f} s)')
ax1.axhline(T_cpu_fixed + T_gpu_kernel, color='gray', ls=':', lw=1,
            label=f'$B\\to\\infty$ limit ({T_cpu_fixed + T_gpu_kernel:.0f} s)')

for name, bw, color, marker, ms, measured, off_a, off_b in gpus_ab:
    t = T_mg_gpu_B(bw)
    ax1.plot(bw, t, marker, color=color, ms=ms, zorder=5,
             markeredgecolor='k', markeredgewidth=0.5)
    label = name.replace('\n', ' ')
    va = 'top' if off_a[1] < 0 else 'bottom'
    ax1.annotate(label, (bw, t), textcoords='offset points',
                 xytext=(off_a[0], off_a[1]), fontsize=9.5, color=color,
                 fontweight='bold', ha=off_a[2], va=va)

ax1.set_xlabel('Effective CPU$-$GPU Bandwidth (GB/s)', fontsize=13)
ax1.set_ylabel('MG AMR Wall Time (s)', fontsize=13)
ax1.set_xscale('log')
ax1.set_xlim(5, 600)
ax1.set_ylim(0, 1300)
ax1.tick_params(labelsize=11)
ax1.text(0.03, 0.97, '(a) Multigrid Poisson Solver',
         transform=ax1.transAxes, fontsize=13, fontweight='bold',
         va='top', ha='left')
ax1.legend(fontsize=9.5, loc='lower center',
           bbox_to_anchor=(0.55, 0.12))
ax1.grid(True, alpha=0.3)
ax1.fill_between([5, 20], 0, 1300, alpha=0.04, color='red')
ax1.fill_between([100, 600], 0, 1300, alpha=0.04, color='blue')
ax1.text(8, 150, 'transfer-\ndominated', fontsize=8.5, color='red', alpha=0.6)
ax1.text(200, 150, 'compute-\nlimited', fontsize=8.5, color='blue', alpha=0.6)

# --- (b) Speedup vs bandwidth ---
ax2.plot(B_range, speedup_B, 'k-', lw=2, label='MG AMR speedup')
ax2.axhline(1.0, color='steelblue', ls='--', lw=1, alpha=0.5, label='Break-even')
ax2.axhline(max_speedup_bw, color='gray', ls=':', lw=1,
            label=f'Asymptotic ({max_speedup_bw:.2f}$\\times$)')

for name, bw, color, marker, ms, measured, off_a, off_b in gpus_ab:
    s = T_cpu_mg / T_mg_gpu_B(bw)
    ax2.plot(bw, s, marker, color=color, ms=ms, zorder=5,
             markeredgecolor='k', markeredgewidth=0.5)
    label = name.replace('\n', ' ')
    if 'H100' in name:
        ax2.annotate(f'{label}\n({s:.2f}$\\times$)', (bw, s),
                     textcoords='data', xytext=(bw, 1.53),
                     fontsize=9.5, color=color, fontweight='bold',
                     ha='left', va='top')
    else:
        va = 'top' if off_b[1] < 0 else 'bottom'
        ax2.annotate(f'{label}\n({s:.2f}$\\times$)', (bw, s),
                     textcoords='offset points', xytext=(off_b[0], off_b[1]),
                     fontsize=9.5, color=color, fontweight='bold',
                     ha=off_b[2], va=va)

ax2.set_xlabel('Effective CPU$-$GPU Bandwidth (GB/s)', fontsize=13)
ax2.set_ylabel('Speedup vs CPU-only', fontsize=13)
ax2.set_xscale('log')
ax2.set_xlim(5, 600)
ax2.set_ylim(0.8, 2.2)
ax2.tick_params(labelsize=11)
ax2.text(0.03, 0.97, '(b) GPU Speedup vs Bandwidth',
         transform=ax2.transAxes, fontsize=13, fontweight='bold',
         va='top', ha='left')
ax2.legend(fontsize=9.5, loc='lower right')
ax2.grid(True, alpha=0.3)

# --- (c) Speedup vs r (improved model) ---
# GH200: 72 Grace ARM cores / 1 GPU = 72 → extend x-axis slightly beyond
r_gh200 = 72   # natural operating point for GH200
r_range = np.linspace(0.5, 80, 500)

# Overlap model: T = max(a*T_cpu, Tk, V/B)  (multi-stream pipeline)
# Since Tk > V/B for all cards, simplifies to max(a*T_cpu, Tk)
def S_r_overlap(r, n_sm):
    tc = T_cpu_r(r)
    Tk = Tk_h100 * (N_SM_h100 / n_sm)
    return tc / max(a_fixed * tc, Tk)

S_ceiling = 1.0 / a_fixed   # theoretical ceiling = 2.79x

print(f"\n{'Card':<12s} {'B':>6s} {'SM':>4s} {'Tk':>7s} {'b':>7s} "
      f"{'S(8)':>6s} {'S_ov(8)':>8s} {'S_ov(72)':>9s}")
for gname, gbw, gsm, gcol, gls, glw in gpu_cards:
    bc = b_for_card(gsm, gbw)
    Tk = Tk_h100 * N_SM_h100 / gsm
    # Serial model (solid)
    sr = np.array([S_r_model(r, bc) for r in r_range])
    ax3.plot(r_range, sr, color=gcol, ls='-', lw=glw,
             label=f'{gname} ({gsm} SM)')
    # Overlap model (dashed, same color)
    sr_ov = np.array([S_r_overlap(r, gsm) for r in r_range])
    ax3.plot(r_range, sr_ov, color=gcol, ls='--', lw=max(glw-0.3, 0.8), alpha=0.7)
    print(f"{gname:<12s} {gbw:>6.1f} {gsm:>4d} {Tk:>7.1f} {bc:>7.1f} "
          f"{S_r_model(8,bc):>6.3f} {S_r_overlap(8,gsm):>8.3f} {S_r_overlap(72,gsm):>9.3f}")

# Break-even and ceiling
ax3.axhline(1.0, color='steelblue', ls='-', lw=0.8, alpha=0.4)
ax3.axhline(S_ceiling, color='k', ls=':', lw=1.0, alpha=0.4)
ax3.text(1.5, S_ceiling + 0.04, f'$1/a = {S_ceiling:.2f}\\times$',
         fontsize=9, color='k', alpha=0.6, va='bottom')

# Legend entries for serial vs overlap
ax3.plot([], [], 'k-', lw=1.5, label='serial (single stream)')
ax3.plot([], [], 'k--', lw=1.0, alpha=0.7, label='overlap (multi-stream)')

# Mark GH200 natural operating point (r=72)
ax3.axvline(r_gh200, color='#9467bd', ls=':', lw=1.0, alpha=0.4)
S_gh200_serial = S_r_model(r_gh200, b_for_card(132, 450))
S_gh200_overlap = S_r_overlap(r_gh200, 132)
ax3.plot(r_gh200, S_gh200_overlap, 'p', color='#9467bd', ms=11, zorder=7,
         markeredgecolor='k', markeredgewidth=0.8)
ax3.annotate(f'GH200 $r={r_gh200}$\n'
             f'serial {S_gh200_serial:.2f}$\\times$\n'
             f'overlap {S_gh200_overlap:.2f}$\\times$',
             (r_gh200, S_gh200_overlap), textcoords='offset points',
             xytext=(-12, 8), fontsize=8.5, color='#9467bd',
             fontweight='bold', ha='right', va='bottom')

# Mark benchmark (r=8) on H100 curve
S_bench = S_r_model(r_bench, b_h100_fit)
ax3.plot(r_bench, S_bench, 'h', color='#d62728', ms=11, zorder=5,
         markeredgecolor='k', markeredgewidth=0.8)

# === Measured data points (n4 only) ===
measured_n4 = [
    ('H100', 4, 1925.750/963.683, '#d62728', 'h'),
    ('H100', 8, 1093.609/665.075, '#d62728', 'h'),
    ('A40',  4, 1924.141/1048.841, '#8c564b', 's'),
    ('A40',  8, 1094.720/748.058,  '#8c564b', 's'),
]

for mname, mr, ms_val, mcol, mmk in measured_n4:
    ax3.plot(mr, ms_val, mmk, color=mcol, ms=10, zorder=6,
             markeredgecolor='k', markeredgewidth=1.0)

ax3.plot([], [], 'h', color='#d62728', ms=9, markeredgecolor='k',
         markeredgewidth=1.0, label='H100 measured')
ax3.plot([], [], 's', color='#8c564b', ms=8, markeredgecolor='k',
         markeredgewidth=1.0, label='A40 measured')

ax3.fill_between(r_range, 0, 1.0, alpha=0.03, color='steelblue')
ax3.text(65, 0.55, 'CPU faster', fontsize=10, color='steelblue',
         fontstyle='italic', ha='center', va='top')

ax3.set_xlabel('OMP Threads per GPU ($r$)', fontsize=13)
ax3.set_ylabel('Speedup vs CPU-only', fontsize=13)
ax3.set_xlim(0, 80)
ax3.set_ylim(0.0, 3.1)
ax3.tick_params(labelsize=11)
ax3.text(0.03, 0.97, '(c) Speedup vs $r$',
         transform=ax3.transAxes, fontsize=13, fontweight='bold',
         va='top', ha='left')
ax3.legend(fontsize=6.5, loc='upper right', ncol=1)
ax3.grid(True, alpha=0.3)

fig.tight_layout()
fig.savefig('/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/fig_gpu_mg_model.pdf',
            bbox_inches='tight', dpi=150)
fig.savefig('/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/fig_gpu_mg_model.png',
            bbox_inches='tight', dpi=150)
print("\nSaved fig_gpu_mg_model.pdf/png")

# === Validation ===
print(f"\n{'Card':<6s} {'r':>3s} {'S_meas':>8s} {'S_model':>8s} {'err':>6s}")
for mname, mr, ms_val, _, _ in measured_n4:
    bv = b_h100_fit if mname == 'H100' else b_a40_fit
    sm = S_r_model(mr, bv)
    err = (sm - ms_val) / ms_val * 100
    print(f"{mname:<6s} {mr:>3d} {ms_val:>8.3f} {sm:>8.3f} {err:>+5.1f}%")
