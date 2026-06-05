"""Extend the alpha-stability trend to alpha=0.85 and 0.90.

Apples-to-apples with the only CONVERGED trend point (alpha=0.8, u0=4, w=4, 32^3:
residual 0.13, Floquet rho=3.59). Same TREND_OPTS, same size, only alpha changes.
Decisive question: does the Floquet radius cross 1 (linearly stable) as alpha rises
from 0.8 toward 0.9? Uses the analytic band_top override (no dense eigensystem).

Run detached:  PYTHONPATH=<repo> python3 run_alpha_hi.py
"""
import csv
import json

from run_breather_sweep import run_one, TREND_OPTS, TREND_CSV_PATH, OUT_DIR
from branesim.solver.breather import omega_band_top_analytic

BAND_TOP = omega_band_top_analytic(1.0, 1.0, 3)  # k_s=1, mass=1, dim=3 -> sqrt(12)

CONFIGS = [
    ("trend_a0p85_u4_w4_32", 0.85, 4.0, 4.0, 32),
    ("trend_a0p9_u4_w4_32",  0.90, 4.0, 4.0, 32),
]

rows = []
for label, alpha, u0, w, grid_n in CONFIGS:
    row = run_one(label, alpha, u0, w, grid_n, TREND_OPTS, band_top_override=BAND_TOP)
    rows.append(row)
    (OUT_DIR / f"{label}_result.json").write_text(json.dumps(row, indent=2, default=str))
    print(f"[RESULT] {label}: converged={row['converged']} "
          f"resid_final={row['residual_final']:.3g} "
          f"rho={row['floquet_spectral_radius']:.4g} stable={row['floquet_stable']} "
          f"R_h={row['R_h']:.3g} walltime={row['walltime_s']:.0f}s")

# Append to the existing trend CSV (header already present).
with open(TREND_CSV_PATH, "a", newline="") as f:
    wr = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    for row in rows:
        wr.writerow(row)

print("\n[ALPHA-HI DONE]")
for row in rows:
    print(f"  alpha={row['alpha']}: rho={row['floquet_spectral_radius']:.4g} "
          f"(resid {row['residual_final']:.3g}, stable={row['floquet_stable']})")
