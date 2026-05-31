"""Standalone local validation for the 3D topological-breather baryon candidate.

This is a VALIDATION DRIVER, not a pytest gate.  It may not converge.
An honest partial result (including "phases unlock") is a high-value finding.

Config: grid 11^3 (or 13^3 if time allows), OPEN boundary, m_ambient=4,
alpha=0.5, u0=0.2, w=1.5*a, P=16, profile power2, k_s=m=a=1.

GO/NO-GO metrics (checked in order, stops at first hard failure):
  #0  convergence:   residual_norm < 1e-8
  #1  phase-lock:    lateral carrier phase within ~0.1 rad of X4 carrier
  #2  winding:       B(l) ≈ 1 per slice (lattice degree)
  #3  confinement:   spread_ratio << 1
  #4  stability:     Floquet spectral radius <= 1.05 (expensive, only if #0-#3 pass)

Usage:
    python run_topological_breather_validation.py

Reports to stdout.  No files written.
"""

from __future__ import annotations

import math
import time
import sys

import numpy as np

from branesim.core.conventions import ActionParams, LatticeParams
from branesim.core.lattice import SpacelikeLattice
from branesim.diagnostics.confinement import confinement_summary
from branesim.solver.breather import (
    BreatherOpts,
    breather_seed_skyrmion,
    floquet_multipliers,
    solve_breather,
)

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

GRID_N = 11          # 11^3 = 1331 nodes
M_AMBIENT = 4
ALPHA = 0.5
U0 = 0.2
A = 1.0              # lattice spacing
W = 1.5 * A
P = 16               # temporal slices per period
PROFILE = "power2"
K_S = 1.0
RHO = 1.0
MASS = RHO * A ** 3  # = 1.0 for a=1


# ---------------------------------------------------------------------------
# GO/NO-GO thresholds
# ---------------------------------------------------------------------------

CONV_TOL = 1e-8
PHASE_LOCK_TOL = 0.1    # radians
WINDING_TOL = 0.10      # |B - 1| < 0.10
SPREAD_RATIO_THRESHOLD = 0.8  # spread_ratio < this = "confined"
FLOQUET_PASS_THRESHOLD = 1.05
X4_WEIGHT_PASS_THRESHOLD = 0.30  # dominant eigenvector >= 30% X4

VERDICT_PASS = "PASS"
VERDICT_FAIL = "FAIL"
VERDICT_PARTIAL = "PARTIAL"


def sep(title: str = "") -> None:
    width = 72
    if title:
        pad = max(0, width - len(title) - 4)
        print(f"\n--- {title} {'-' * pad}")
    else:
        print("-" * width)


def header() -> None:
    print("=" * 72)
    print("  BraneSim  |  Topological Breather Validation  |  3D Skyrmion")
    print(f"  Grid: {GRID_N}^3 ({GRID_N**3} nodes)  m={M_AMBIENT}  alpha={ALPHA}")
    print(f"  u0={U0}  w={W}  P={P}  profile={PROFILE}  k_s=m=a={A}")
    print("=" * 72)


# ---------------------------------------------------------------------------
# Helper: phase-lock measurement
# ---------------------------------------------------------------------------


def measure_phase_lock(slices: np.ndarray, peak_node: int, dim: int) -> dict:
    """Measure per-component carrier phase relative to X4 at the peak node.

    For each component c (lateral 0..dim-1 and X4 = dim), extract the
    time series at the peak node and fit the carrier phase as:
        phi_c = argmax_l A_c * cos(2π l/P + phi_c)
    via the discrete Fourier transform at frequency 1 (one cycle over P).

    Returns phase_x4, phases_lateral, and max|phase_lateral - phase_x4|.
    """
    P, n_nodes, m_ambient = slices.shape
    x4_comp = dim

    # DFT at frequency 1 (one full cycle over P slices)
    l_arr = np.arange(P)
    freqs = np.exp(-2j * np.pi * l_arr / P)

    phases = {}
    amplitudes = {}
    for c in range(m_ambient):
        signal = slices[:, peak_node, c]
        coeff = np.sum(signal * freqs)
        phases[c] = float(np.angle(coeff))
        amplitudes[c] = float(np.abs(coeff))

    phase_x4 = phases[x4_comp]
    lateral_phases = [phases[c] for c in range(dim)]
    lateral_amps = [amplitudes[c] for c in range(dim)]

    # Phase differences (wrapped to [-pi, pi])
    phase_diffs = []
    for c in range(dim):
        diff = phases[c] - phase_x4
        # wrap to [-pi, pi]
        diff = (diff + math.pi) % (2 * math.pi) - math.pi
        phase_diffs.append(diff)

    max_diff = max(abs(d) for d in phase_diffs) if phase_diffs else float("nan")

    return {
        "phase_x4": phase_x4,
        "lateral_phases": lateral_phases,
        "lateral_amps": lateral_amps,
        "amp_x4": amplitudes[x4_comp],
        "phase_diffs_rad": phase_diffs,
        "max_phase_diff_rad": max_diff,
        "locked": max_diff < PHASE_LOCK_TOL,
    }


# ---------------------------------------------------------------------------
# Helper: approximate topological winding number per slice
# ---------------------------------------------------------------------------


def approximate_winding_per_slice(slices: np.ndarray, ref: np.ndarray, lattice: SpacelikeLattice, dim: int) -> np.ndarray:
    """Approximate topological degree B per temporal slice.

    Uses the south-to-north pole sweep of the X4 component: the degree is
    estimated by checking the fraction of nodes mapping near each pole versus
    a uniform baseline.  This is a fast APPROXIMATE check sufficient for a
    go/no-go decision; it is not the rigorous 24pi^2 lattice degree integral.

    For the Skyrme map:
      - South pole (r≈0): X4 disp ≈ -u0   →  B=1: present
      - North pole (r>>w): X4 disp ≈ +u0  →  B=1: present

    The signed "polar asymmetry" is a proxy: if it tracks |x4_min/u0| * |x4_max/u0| ≈ 1.
    A random (B=0) map would give near-zero min or a non-pole max.

    Returns an array of shape (P,) with the approximate degree at each slice.
    """
    P, n_nodes, m_ambient = slices.shape
    x4_comp = dim

    degrees = np.zeros(P)
    for l_idx in range(P):
        disp_x4 = slices[l_idx, :, x4_comp] - ref[:, x4_comp]
        x4_min = float(np.min(disp_x4))
        x4_max = float(np.max(disp_x4))
        # Carrier at this slice (should = cos(2pi*l/P) * u0 * cos(pi) = -u0*carrier)
        # The degree depends only on the direction map, not the S^3 radius.
        # As long as min<0 and max>0, the map sweeps the sphere.
        if abs(x4_max) > 1e-6 and abs(x4_min) > 1e-6 and x4_min < 0 < x4_max:
            # Estimate relative coverage: both poles present => degree ≈ 1
            # More precisely: |x4_min/u0_effective| * sign is the cos(F(0)) proxy
            carrier = math.cos(2.0 * math.pi * l_idx / P)
            if abs(carrier) > 0.1:
                effective_u0 = abs(carrier) * U0
                south_ok = abs(x4_min / effective_u0) > 0.3
                north_ok = abs(x4_max / effective_u0) > 0.3
                degrees[l_idx] = 1.0 if (south_ok and north_ok) else 0.5
            else:
                # carrier ~ 0 near l=P/4 and l=3P/4: spatial structure collapses
                degrees[l_idx] = float("nan")
        else:
            degrees[l_idx] = 0.0

    return degrees


# ---------------------------------------------------------------------------
# Floquet: eigenvector X4 weight
# ---------------------------------------------------------------------------


def floquet_x4_weight(fl: dict, n_nodes: int, m_ambient: int) -> float:
    """Fraction of energy in the X4 channel for the dominant eigenvector.

    The monodromy eigenvectors are 2N-dim state vectors [dR^0; dR^{-1}],
    each dR has shape (n_nodes, m_ambient).  We report the fraction of the
    L2 norm in the X4 component (index m_ambient-1) for the largest-|rho|
    eigenvector.

    This is only available from the dense monodromy path.  If the Arnoldi
    path was used (fl["dense"]=False), return nan.
    """
    if not fl.get("dense", False):
        return float("nan")

    # fl["multipliers"] are eigenvalues but we need eigenvectors.
    # The monodromy was assembled in floquet_multipliers as M[:, j] = M_apply(e_j).
    # We don't have eigenvectors directly — they weren't returned.
    # Return nan here; the caller can report "eigenvector X4 weight: not available"
    return float("nan")


# ---------------------------------------------------------------------------
# Main validation
# ---------------------------------------------------------------------------


def run_validation() -> str:
    """Run the validation and return overall verdict string."""
    header()

    # Build lattice
    lp = LatticeParams(
        grid_shape=(GRID_N, GRID_N, GRID_N),
        spacing=A,
        periodic_axes=(False, False, False),  # OPEN boundary
    )
    lattice = SpacelikeLattice(lp)
    params = ActionParams(
        k_s=K_S, alpha=ALPHA, rho=RHO, dt=0.1, n_slices=1, m_ambient=M_AMBIENT,
    )
    n_nodes = lattice.n_nodes
    dim = lattice.dim
    ref = lattice.reference_positions(M_AMBIENT)

    print(f"\nLattice: {GRID_N}^3 = {n_nodes} nodes,  dim={dim},  m_ambient={M_AMBIENT}")
    print(f"DOF per slice: {n_nodes * M_AMBIENT},  total unknowns: {P * n_nodes * M_AMBIENT + 1}")

    # Build seed
    sep("Seed construction")
    t_seed_start = time.perf_counter()
    slices_seed, T_seed = breather_seed_skyrmion(
        lattice, m=M_AMBIENT, u0=U0, w=W, P=P, profile=PROFILE
    )
    t_seed = time.perf_counter() - t_seed_start
    print(f"  Seed built in {t_seed:.2f}s")
    print(f"  T_seed = {T_seed:.4f}  (omega_seed = {2*math.pi/T_seed:.4f})")
    print(f"  X4 at peak node l=0: {slices_seed[0, lattice.n_nodes//2, 3]:.4f}  "
          f"(expected ~ {-U0:.4f})")

    # ---------------------------------------------------------------------------
    # METRIC #0: Convergence
    # ---------------------------------------------------------------------------
    sep("METRIC #0: Newton-Krylov convergence")

    opts = BreatherOpts(
        tol=1e-8,
        max_iter=200,
        inner_maxiter=3000,
        method="lgmres",
        verbose=True,
    )

    print(f"\n  Starting Newton-Krylov (max_iter={opts.max_iter}, "
          f"inner_maxiter={opts.inner_maxiter}, tol={opts.tol:.0e})...")
    print(f"  Grid: {GRID_N}^3,  P={P},  mode=topological")
    t_solve_start = time.perf_counter()

    result = solve_breather(
        lattice, params, MASS,
        P=P,
        amplitude=U0,
        seed=slices_seed,
        T_seed=T_seed,
        mode="topological",
        skyrme_profile=PROFILE,
        skyrme_w=W,
        opts=opts,
    )

    walltime = result["walltime_s"]
    res_norm = result["residual_norm"]
    converged = result["converged"]
    slices_sol = result["slices"]
    T_sol = result["T"]
    omega_sol = result["omega"]
    peak_node = result["peak_node"]

    print(f"\n  Residual norm:  {res_norm:.3e}")
    print(f"  Converged:      {converged}")
    print(f"  Walltime:       {walltime:.1f}s")
    print(f"  T_converged:    {T_sol:.4f}  (omega = {omega_sol:.4f})")

    if not converged:
        print(f"\n  METRIC #0: FAIL — residual {res_norm:.3e} > tol {CONV_TOL:.0e}")
        print(f"  (Honest partial result.  Reporting best-iterate metrics below.)")
        verdict_0 = VERDICT_FAIL
    else:
        print(f"\n  METRIC #0: PASS — residual {res_norm:.3e} < {CONV_TOL:.0e}")
        verdict_0 = VERDICT_PASS

    # ---------------------------------------------------------------------------
    # METRIC #1: Phase-lock (cheap, decisive)
    # ---------------------------------------------------------------------------
    sep("METRIC #1: Phase-lock (carrier phase coherence)")

    pl = measure_phase_lock(slices_sol, peak_node, dim)
    print(f"  Phase X4:         {pl['phase_x4']:.4f} rad  (amp={pl['amp_x4']:.4f})")
    for i, (ph, amp) in enumerate(zip(pl["lateral_phases"], pl["lateral_amps"])):
        diff = pl["phase_diffs_rad"][i]
        print(f"  Phase lateral[{i}]: {ph:.4f} rad  (amp={amp:.4f})  |diff|={abs(diff):.4f} rad")
    print(f"  Max phase diff:   {pl['max_phase_diff_rad']:.4f} rad  (threshold = {PHASE_LOCK_TOL:.2f} rad)")

    if pl["locked"]:
        print(f"  METRIC #1: PASS — phases locked within {PHASE_LOCK_TOL} rad")
        verdict_1 = VERDICT_PASS
    else:
        print(f"  METRIC #1: FAIL — phases NOT locked: max diff = {pl['max_phase_diff_rad']:.4f} rad")
        print(f"  This REFUTES the topological-stabilization hypothesis before Floquet.")
        verdict_1 = VERDICT_FAIL

    if verdict_0 == VERDICT_FAIL and not converged:
        # Still check phase lock on the best iterate (informative even if not converged)
        pass

    # ---------------------------------------------------------------------------
    # METRIC #2: Winding number
    # ---------------------------------------------------------------------------
    sep("METRIC #2: Winding number B(l) per slice")

    degrees = approximate_winding_per_slice(slices_sol, ref, lattice, dim)
    # Only check non-NaN (carrier amplitude > 0.1*u0) slices
    valid_mask = ~np.isnan(degrees)
    if valid_mask.sum() > 0:
        B_mean = float(np.nanmean(degrees[valid_mask]))
        B_min = float(np.nanmin(degrees[valid_mask]))
        B_max = float(np.nanmax(degrees[valid_mask]))
        print(f"  B per slice (valid slices={valid_mask.sum()}/{P}):")
        for l_idx in range(P):
            if valid_mask[l_idx]:
                print(f"    l={l_idx:2d}: B ~ {degrees[l_idx]:.2f}")
            else:
                print(f"    l={l_idx:2d}: carrier ~ 0  (skip)")
        print(f"  B mean={B_mean:.3f}  min={B_min:.3f}  max={B_max:.3f}")
        if abs(B_mean - 1.0) < WINDING_TOL:
            print(f"  METRIC #2: PASS — B ≈ {B_mean:.3f} (within {WINDING_TOL} of 1.0)")
            verdict_2 = VERDICT_PASS
        else:
            print(f"  METRIC #2: FAIL — B ≈ {B_mean:.3f} (expected 1.0 ± {WINDING_TOL})")
            verdict_2 = VERDICT_FAIL
    else:
        print(f"  No valid slices found.  METRIC #2: PARTIAL")
        verdict_2 = VERDICT_PARTIAL

    # ---------------------------------------------------------------------------
    # METRIC #3: Confinement
    # ---------------------------------------------------------------------------
    sep("METRIC #3: Confinement (spread_ratio)")

    conf = confinement_summary(slices_sol, ref, dim=dim)
    spread_mean = float(np.mean(conf["spread_ratio"]))
    spread_final = float(conf["final"]["spread_ratio"])
    w_over_a = W / A
    print(f"  box_fill_radius:  {conf['box_fill_radius']:.3f}")
    print(f"  spread_ratio mean: {spread_mean:.4f}  final: {spread_final:.4f}")
    print(f"  w/a = {w_over_a:.2f}  (object scale relative to lattice spacing)")
    if spread_final < SPREAD_RATIO_THRESHOLD:
        print(f"  METRIC #3: PASS — spread_ratio {spread_final:.4f} << 1 (confined)")
        verdict_3 = VERDICT_PASS
    else:
        print(f"  METRIC #3: FAIL — spread_ratio {spread_final:.4f} >= {SPREAD_RATIO_THRESHOLD} "
              f"(not well-confined)")
        verdict_3 = VERDICT_FAIL

    # ---------------------------------------------------------------------------
    # METRIC #4: Floquet stability (only if #0-#3 all PASS)
    # ---------------------------------------------------------------------------
    verdict_4 = "SKIPPED"
    if verdict_0 == VERDICT_PASS and verdict_1 == VERDICT_PASS and verdict_3 == VERDICT_PASS:
        sep("METRIC #4: Floquet spectral radius (stability)")
        n_state = 2 * n_nodes * M_AMBIENT
        print(f"  2N state dimension: {n_state}  (dense threshold = 400)")

        if n_state > 400:
            print(f"  Using Arnoldi (matrix-free) path, k=6 multipliers...")

        fl_opts = {
            "n_multipliers": 6,
            "dense_threshold": 400,  # will use Arnoldi for 11^3 (n_state >> 400)
            "stability_tol": 0.05,
        }

        t_fl_start = time.perf_counter()
        try:
            fl = floquet_multipliers(
                slices_sol, T_sol, lattice, params, MASS, **fl_opts
            )
            t_fl = time.perf_counter() - t_fl_start

            rho_max = fl["spectral_radius"]
            n_unstable = fl["n_unstable"]
            stable = fl["stable"]
            mags = np.abs(fl["multipliers"])

            print(f"  Spectral radius rho_max: {rho_max:.4f}")
            print(f"  n_unstable:              {n_unstable}")
            print(f"  stable (|rho|<=1.05):    {stable}")
            print(f"  |rho| values: {np.round(mags, 4)}")
            print(f"  Floquet walltime: {t_fl:.1f}s")

            if stable and rho_max <= FLOQUET_PASS_THRESHOLD:
                print(f"  METRIC #4: PASS — rho_max = {rho_max:.4f} <= {FLOQUET_PASS_THRESHOLD}")
                verdict_4 = VERDICT_PASS
            else:
                print(f"  METRIC #4: FAIL — rho_max = {rho_max:.4f} > {FLOQUET_PASS_THRESHOLD}")
                print(f"  (bare transverse breather has rho_max ~ 2.0 in 1D)")
                verdict_4 = VERDICT_FAIL

        except Exception as exc:
            print(f"  Floquet failed with exception: {exc}")
            verdict_4 = f"ERROR: {exc}"
    else:
        sep("METRIC #4: Floquet (SKIPPED — prerequisites not met)")
        failed = [f"#{i}" for i, v in enumerate([verdict_0, verdict_1, verdict_3])
                  if v != VERDICT_PASS]
        print(f"  Skipping Floquet because metrics {failed} did not pass.")

    # ---------------------------------------------------------------------------
    # Summary
    # ---------------------------------------------------------------------------
    sep("SUMMARY")
    print(f"\n  #0  Convergence:    {verdict_0}  (residual = {res_norm:.3e}, wall = {walltime:.1f}s)")
    print(f"  #1  Phase-lock:     {verdict_1}  (max diff = {pl['max_phase_diff_rad']:.4f} rad)")
    print(f"  #2  Winding:        {verdict_2}")
    print(f"  #3  Confinement:    {verdict_3}  (spread_ratio = {spread_final:.4f})")
    print(f"  #4  Floquet:        {verdict_4}")

    # Overall verdict
    all_verdicts = [verdict_0, verdict_1, verdict_2, verdict_3]
    if all(v == VERDICT_PASS for v in all_verdicts) and verdict_4 == VERDICT_PASS:
        overall = "TOPOLOGICAL BREATHER: FULL PASS"
    elif verdict_0 == VERDICT_FAIL:
        overall = "NO_CONVERGENCE — solve did not reach tolerance"
    elif verdict_1 == VERDICT_FAIL:
        overall = "PHASES UNLOCKED — topological-stabilization hypothesis REFUTED before Floquet"
    elif all(v == VERDICT_PASS for v in all_verdicts) and verdict_4 == VERDICT_FAIL:
        overall = "CONVERGED + LOCKED + CONFINED, but FLOQUET UNSTABLE"
    elif all(v == VERDICT_PASS for v in all_verdicts[:3]) and verdict_4 == "SKIPPED":
        overall = "CONVERGED + LOCKED + WOUND + CONFINED (Floquet not reached)"
    else:
        overall = "PARTIAL — mixed metrics (see above)"

    print(f"\n  OVERALL VERDICT: {overall}")
    sep()
    return overall


if __name__ == "__main__":
    verdict = run_validation()
    # Exit 0 if any PASS was reached; 1 only for total failure
    if "REFUTED" in verdict or "NO_CONVERGENCE" in verdict:
        sys.exit(1)
    sys.exit(0)
