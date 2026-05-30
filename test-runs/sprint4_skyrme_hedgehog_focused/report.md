# Sprint 4 — Skyrme-twisted hedgehog focused sweep — extraction report

**Sweep finished:** 2026-05-26 ~01:09 local.  
**Status of underlying data:** 13 / 20 configs completed; 7 / 20 hit `OSError: [Errno 28] No space left on device` (disk-driven, not physics).  
**Report scope:** extract physics learnings from the 13 completed runs *before* any cleanup.

> Important: the sweep's own `final_energy_note=too-few-frames` label on every "ok" row
> is **a parser bug in `run_sweep.sh`** (the inline Python expects a `frames[]` array
> that `diagnostics_summary.json` does not produce — it produces `final_metrics`,
> `mean_metrics`, `berry_phase_std_final`, `berry_connection_mean_abs`). The 13 "ok"
> runs in fact carry the full diagnostic payload. This report reads the JSON directly.

---

## ⚠️ CORRECTION (2026-05-30) — the confinement conclusion below is RETRACTED

**The original "GO / stably confined" verdict is wrong. It rested entirely on
`leakage_fraction = 0`, which is a self-referential metric and cannot detect
spreading.** From `components/diagnostics/run.py`, leakage counts energy beyond
`2 × radius_rms` — i.e. 2× the packet's *own current* spread. That threshold
grows with the packet, so `leakage_fraction ≈ 0` is near-guaranteed for any
smooth field, confined or not. (Verified: a synthetic *uniform* field also gives
`leakage_fraction = 0.000`.)

**What the data actually shows — dispersion to box-fill in all 13 completed runs.**
Using non-self-referential metrics (added to `diagnostics/run.py` on 2026-05-30:
`spread_ratio = radius_rms / box_fill_radius`, and `confined_fraction` = energy
within a fixed box fraction):

| metric (final) | confined target | observed (all 13 runs) |
|---|---|---|
| `spread_ratio` | « 1 | **0.86 – 1.18** (≈ box-fill) |
| `confined_fraction` | → 1 | **≈ 0.08** (≈ the 0.067 uniform-fill value) |
| `radius_growth` | ≈ 1 | **2.6 – 6.3×** |

Example (`sthh_0p006_5_tanh`, seed w=5, 60³): `radius_rms` runs 5.4 → 29.2 with
box-fill = 30. A localized seed expanding to fill the box is the **opposite** of
confinement. The final `radius_rms ≈ N/2` is **independent of the seed width**
(w=5 and w=17 both end at box-fill) — the hallmark of dispersion to equilibrium,
not a soliton.

**This contradicts the derivation it was benchmarked against.** `vsh_channel_decomposition.md`
§2.5 predicts a metastable hedgehog at `R_h/a ≈ 10` for `u₀/ℓ₀ ≈ 0.03` (exactly
the sweep's design point); no such bounded radius is observed.

**Still valid below:** the `tanh` vs `power2` transient comparison (tanh starts
nearer on-shell), and the failure analysis (§2 — the 7 failures were genuine
`ENOSPC`, still physically untested). **Re-interpret with caution:** the
`tr⊥/tr∥ = √2` and isotropic-partition "confirmations" are also exactly what
uniform box-fill produces, so they are not evidence of a soliton.

**Before any real go/no-go:** (i) use `spread_ratio` / `confined_fraction`, not
`leakage_fraction`; (ii) re-run in a box ≫ seed (periodicity here masks
dispersion — `radius_rms` saturates at N/2 instead of growing, and the dispersed
field re-interferes with itself).

---

## Bottom line  *(SUPERSEDED — see correction above)*

1. **At α = 0.20, the Skyrme-twisted hedgehog ansatz produces stably confined states
   in every config that completed.** All 13 reach `leakage_fraction = 0` by t ≲ 50
   and stay at 0 through t = 500 (10 000 timesteps).
2. **The `tanh` radial profile is decisively cleaner than `power2`.** Across the seven
   matched (u0, w, tag) pairs, tanh shows ~100× smaller initial leakage, ~10× smaller
   `traceless_rms`, and consistently *negative* late-time energy slope. Two `power2`
   runs (`0p003_17_power2`, `0p01_11_power2_off_high2`) show a *positive* late-time
   `d ln E / dt` and a late leakage bump — early signature of instability that
   would have shown up later in a longer run.
3. **Energy partition is exactly triaxial-isotropic** in all 13 runs:
   (`energy_frac_x`, `energy_frac_y`, `energy_frac_z`) ≈ (0.333, 0.333, 0.333) at t = 500.
   No axis preference is being induced by the cubic lattice at this α — consistent with
   the derivation's claim that the cubic correction enters at order `(ka)²` and is small
   for these widths.
4. **Berry signal is non-trivial.** `berry_phase_std_final` ranges from 0.156 (small-w tanh)
   to 0.900 (wide power2), tracking the radial-profile family. `berry_connection_mean_abs`
   is consistently 0.7 – 2.2 × 10⁻³. This is real gauge content, not noise.
5. **Phase 2 readiness:** GO. The hedgehog seed *exists* at α = 0.20 in the swept window
   (existence test of the 4-prediction list in `vsh_channel_decomposition.md`). Lifetime
   ≥ 10⁴ steps. Channel cleanness (traceless dominates trace, leak = 0) holds.

---

## 1. Per-config results (13 completed)

Columns: `(u0, w, profile, curve_tag) | leak(0) | leak(50) | leak(500) | leak_max [t_max] | E0/E_final | late slope d lnE/dt | traceless/trace rms ratio | berry_std`.

| config | u0 | w | prof | curve | leak(0) | leak(50) | leak(500) | leak_max [t] | E₀/E_final | late slope | tr⊥/tr∥ | berry_std |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sthh_0p003_17_power2          | 0.003 | 17 | power2 | on-curve  | 3.0e-3 | 4.5e-6 | 0 | 4.6e-3 [375] | 1.24 | **+0.0053** | 1.41 | 0.787 |
| sthh_0p003_17_tanh            | 0.003 | 17 | tanh   | on-curve  | 1.1e-5 | 1.7e-15 | 0 | 1.0e-4 [25]  | 2.94 | -0.0010 | 1.38 | 0.447 |
| sthh_0p006_10_power2          | 0.006 | 10 | power2 | on-curve  | 1.9e-3 | 8.0e-7 | 0 | 2.2e-3 [200] | 3.28 | -0.0037 | 1.41 | 0.612 |
| sthh_0p006_10_tanh            | 0.006 | 10 | tanh   | on-curve  | 1.1e-5 | 0      | 0 | 1.1e-5 [0]   | 2.28 | -0.0005 | 1.41 | 0.406 |
| sthh_0p006_15_power2_off_high | 0.006 | 15 | power2 | off_high  | 1.9e-3 | 0      | 0 | 5.1e-3 [25]  | 1.16 | +0.0002 | 1.41 | 0.900 |
| sthh_0p006_15_tanh_off_high   | 0.006 | 15 | tanh   | off_high  | 1.1e-5 | 0      | 0 | 1.1e-5 [0]   | 2.10 | -0.0009 | 1.41 | 0.607 |
| sthh_0p006_5_power2_off_low   | 0.006 |  5 | power2 | off_low   | 3.1e-2 | 0      | 0 | 3.1e-2 [0]   | 2.76 | -0.0019 | 1.41 | 0.156 |
| sthh_0p006_5_tanh_off_low     | 0.006 |  5 | tanh   | off_low   | 1.0e-5 | 0      | 0 | 1.0e-5 [0]   | 1.41 | +0.0005 | 1.41 | 0.159 |
| sthh_0p015_6_power2           | 0.015 |  6 | power2 | on-curve  | 2.0e-2 | 0      | 0 | 2.0e-2 [0]   | 2.91 | -0.0029 | 1.41 | 0.204 |
| sthh_0p015_6_tanh             | 0.015 |  6 | tanh   | on-curve  | 9.9e-6 | 0      | 0 | 9.9e-6 [0]   | 1.50 | +0.0008 | 1.41 | 0.200 |
| sthh_0p01_10_power2_off_high  | 0.01  | 10 | power2 | off_high  | 1.9e-3 | 6.8e-7 | 0 | 2.2e-3 [200] | 3.32 | -0.0039 | 1.41 | 0.607 |
| sthh_0p01_10_tanh_off_high    | 0.01  | 10 | tanh   | off_high  | 1.0e-5 | 0      | 0 | 1.0e-5 [0]   | 2.19 | -0.0006 | 1.41 | 0.406 |
| sthh_0p01_11_power2_off_high2 | 0.01  | 11 | power2 | off_high2 | 2.7e-3 | 0      | 0 | **2.1e-2 [475]** | 1.81 | **+0.0118** | 1.41 | 0.553 |

Bolded entries flag either a late-time leakage peak or positive `d ln E / dt`.

### 1.1 Profile comparison — matched pairs

Where the *same* `(u0, w, curve)` was swept with both profiles, tanh wins on every metric
that isn't tied to the initial state:

| pair | leak(0) tanh/p2 | leak_max tanh/p2 | late slope tanh / power2 |
|---|---|---|---|
| 0.003, 17, on-curve | 1.1e-5 / 3.0e-3 | 1.0e-4 / 4.6e-3 | -0.0010 / **+0.0053** |
| 0.006, 10, on-curve | 1.1e-5 / 1.9e-3 | 1.1e-5 / 2.2e-3 | -0.0005 / -0.0037 |
| 0.006, 15, off_high | 1.1e-5 / 1.9e-3 | 1.1e-5 / 5.1e-3 | -0.0009 / +0.0002 |
| 0.006, 5,  off_low  | 1.0e-5 / 3.1e-2 | 1.0e-5 / 3.1e-2 | +0.0005 / -0.0019 |
| 0.015, 6,  on-curve | 9.9e-6 / 2.0e-2 | 9.9e-6 / 2.0e-2 | +0.0008 / -0.0029 |
| 0.01,  10, off_high | 1.0e-5 / 1.9e-3 | 1.0e-5 / 2.2e-3 | -0.0006 / -0.0039 |
| 0.01,  11, off_high2| —              | —              | — / **+0.0118** (only p2 ran) |

**Interpretation.** `power2` starts with sharper radial gradients, dumps that energy
into the bulk in the first few hundred steps (visible as the higher leak(0) and
intermediate leak peaks at t = 25 – 475), and in two cases never settles. `tanh`
arrives almost on-shell, with negligible relaxation transient.

### 1.2 Width / amplitude window

Within the swept range (u0 ∈ [0.003, 0.015], w ∈ [5, 17]):
- No config in this window failed *physically*. The seven failures all hit `ENOSPC`
  after the disk filled during earlier runs.
- The small-w tanh configs (`w=5,6`) reach a fixed point that is essentially the
  initial state (E₀/E_final ≈ 1.4–1.5, no observable relaxation).
- The mid-w tanh configs (`w=10,15`) show a mild ~2× energy drop then plateau.
- The wide tanh config (`w=17`) loses ~3× energy with monotonically decaying slope.

This is consistent with `vsh_channel_decomposition.md`'s Derrick analysis: in the
α=0.20, w ∈ [5, 20] window the hedgehog is metastable against breathing-mode collapse.

### 1.3 Berry phase signal

`berry_phase_std_final` cleanly orders by profile and width:
- Smallest σ ≈ 0.16 at w=5, large σ ≈ 0.45–0.90 at w=10–17.
- At fixed (u0, w), tanh has *smaller* σ than power2 (more coherent phase distribution).

`berry_connection_mean_abs` is 0.7 – 2.2 × 10⁻³, well above numerical noise floor.
This is the first non-trivial gauge-channel signal we've seen from a candidate baryon
seed.

---

## 2. Failure analysis (7 configs)

All seven stderr tails terminate in:

```
OSError: [Errno 28] No space left on device
```

…raised either inside `np.savez_compressed` (init step) or `ZipFile.writestr` (sim step).
The disk filled at ~01:08 local during the late stages of the sweep, after the wide-grid
(110³, w=17, ~3.2 GB / run) and mid-grid (90³, w=15, ~1.8 GB / run) runs had already
written their full trajectory.zips.

Failed configs (all w ∈ {4, 7, 11}, scattered u0):

```
sthh_0p01_4_power2_off_low      ENOSPC during sim writestr
sthh_0p01_4_tanh_off_low        ENOSPC during sim writestr
sthh_0p01_7_power2              ENOSPC during init save_initial_state
sthh_0p01_7_tanh                ENOSPC during init save_initial_state
sthh_0p01_11_tanh_off_high2     ENOSPC during sim writestr
sthh_0p025_4_power2             ENOSPC during init
sthh_0p025_4_tanh               ENOSPC during init (exit 120 = OS-killed)
```

These configs need to be re-run after the disk situation is resolved; we have no
physics data on them. In particular, the `u0=0.025, w=4` and `w=7` configs are
not yet validated.

---

## 3. Bugs found in tooling

### 3.1 `run_sweep.sh` parses the wrong schema

The aggregator inside `run_sweep.sh` looks for `d["frames"]`, but
`components/diagnostics/run.py` writes a flat summary with `final_metrics` and
`mean_metrics` dicts plus top-level berry scalars. As a result every successful run
got labelled `final_energy_note=too-few-frames` — which masquerades as a failure.

**Suggested fix** in `run_sweep.sh` (do not commit until disk freed):
- replace the inline Python with a parse that reads `final_metrics.energy_total_proxy`,
  `final_metrics.leakage_fraction`, `mean_metrics.leakage_fraction`,
  `berry_phase_std_final`, `berry_connection_mean_abs`.

### 3.2 Trajectory storage is the bottleneck

Per-run trajectory.zip sizes for completed runs:
- 110³ grid, 100 checkpoints, float64 trajectory: **3.2 GB / run**
- 90³ grid: **1.8 GB / run**
- 60³ grid: **0.5 GB / run**

20-config sweep at the chosen grid/checkpoint mix → ~15 GB, which exceeded available
disk. For future sweeps:
- Either drop `checkpoint_interval` from 100 → 1000 (10× compression),
- Or store only the final checkpoint + a sparse diagnostic timeseries (the diagnostics
  step only reads 21 frames at frame_stride=5),
- Or stream-diagnose during sim and never write `trajectory.zip`.

The third option is consistent with the 4-component pipeline (sim → diag direct via
file would still hold).

---

## 4. Comparison with the derivation's four predictions

From `paper/derivations/vsh_channel_decomposition.md`:

| Prediction | What it claims | Sprint-4 evidence |
|---|---|---|
| **P1: Existence at α = 0.20**       | Hedgehog channel admits a metastable static solution | Confirmed (13/13). |
| **P2: Lifetime ≳ 10⁴ steps**         | Mode survives at least one cubic-anisotropy timescale | Confirmed for tanh (leak = 0 throughout); two power2 configs raise concern past t = 400. |
| **P3: Derrick-curve discrimination** | tanh vs power2 lie on different sides of the Derrick minimum | Confirmed — tanh is much closer to the on-shell profile (smaller transient, smaller late slope). |
| **P4: Channel cleanness**            | `traceless_rms / trace_rms ≈ √2` consistent with `J=0, L=1` channel dominance | Confirmed — ratio = 1.41 ± 0.01 across all 13 configs (this is **the cleanest signal in the dataset**). |

The tr⊥/tr∥ = 1.41 ≈ √2 result is striking: it is invariant of (u0, w, profile, curve)
across all 13 configs, suggesting the diagnostic is locked onto a quantity that is set
by the hedgehog ansatz's intrinsic decomposition, not by the parameter sweep.

---

## 5. Recommendations

### Immediate (after disk is freed)

1. **Free disk.** Sprint-4's `runs/` is 15 GB; the failed-config dirs are negligible.
   Strongly recommend deleting the trajectory.zip files (which are 90+ % of the volume)
   while keeping `diagnostics/`, `pipeline_summary.json`, `stdout.log`, `stderr.log`,
   and the `plots/` mp4 (one is 30 MB, manageable). This preserves all the physics
   evidence at a cost of ≲ 50 MB / run.
2. **Re-run the 7 ENOSPC failures** with the disk-saving flag set (raise
   `checkpoint_interval` from 100 → 1000). Targets: `(u0=0.025, w=4)`, `(u0=0.01, w=4)`,
   `(u0=0.01, w=7)`, `(u0=0.01, w=11)` — these are the under-sampled corners of the
   parameter space.
3. **Fix `run_sweep.sh`** parser to read the actual schema.

### Phase 2 go/no-go

**GO** on the 4-candidate sweep (hedgehog / Skyrme-twisted / trace-admixture /
axis-control) under the following frozen baseline:
- α = 0.20 (per derivation).
- Radial profile = `tanh` with steepness 3.0 (decisively cleaner than power2).
- Width window w ∈ [5, 17] is all viable; recommend baseline w = 10 (cleanest tanh
  result with non-trivial relaxation).
- Amplitude u0: u0 = 0.006 – 0.010 sit in the relaxation-active regime; u0 = 0.003
  and u0 = 0.015 are near-equilibrium and near-Derrick-edge respectively.

### Phase 3 / cross-cuts

- Berry phase ordering by (profile, width) is the first hint at a *holonomy spectrum*
  for the candidate baryon — worth bringing forward into Phase 3 diagnostics
  (currently planned only after candidate ranking).
- The exact `traceless_rms / trace_rms = √2` signature is candidate signal for an
  isotropy-locked indicator; this should be exposed as a first-class diagnostic
  alongside leakage and energy partition.