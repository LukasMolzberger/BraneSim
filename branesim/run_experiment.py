"""branesim experiment entrypoint (CLI).

Reads a JSON config, builds the spacelike lattice + action parameters + a seed,
runs a block-solve (or the IVP special case), and writes a world-volume zip plus
a summary JSON to an output directory. This is the unit the AWS launcher runs
remotely (see orchestration/aws/ and DEPLOYMENT.md).

Status: runs the VALIDATED paths only —
  - mode="ivp"           : forward Verlet march (rest start), validated.
  - mode="bvp_dirichlet" : JFNK block-solve with Dirichlet two-time BCs derived
                           from an IVP march; validated to recover the march.
  - mode="bvp_chiral"    : Chiral Cauchy BC march from two past slices (R0, R1).
                           Well-posed for all N; κ bounded and N-independent.
                           Verdict (a) implemented 2026-05-30.

Config schema (JSON)::

    {
      "lattice": {"grid_shape": [16,16,16], "spacing": 1.0,
                  "periodic_axes": [true,true,true], "axial_weight": 1.0},
      "action":  {"k_s": 1.0, "alpha": 0.2, "rho": 1.0, "dt": 0.1,
                  "n_slices": 64, "m_ambient": 4, "r_t": 0.0},
      "seed":    {"kind": "plane_wave", "amplitude": 1e-3,
                  "k_index": [1,0,0], "polarization": [0,1,0,0], "rng_seed": 0},
      "solver":  {"mode": "ivp", "tol": 1e-9, "max_iter": 100, "warm_start": true}
    }

Usage::

    branesim-run --config config.json --output-dir ./run-out
    python -m branesim.run_experiment --config config.json --output-dir ./run-out
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from branesim.core.lattice import SpacelikeLattice
from branesim.core.conventions import LatticeParams, ActionParams
from branesim.core.residual import residual_norm
from branesim.solver.ivp import IVPProblem, march
from branesim.solver.bvp import BoundaryProblem, SolveOpts, solve_block
from branesim.solver.boundary import DirichletBC, ChiralBC
from branesim.io.contracts import WorldVolumeWriter


def _build_seed(seed_cfg: dict, lattice: SpacelikeLattice, m: int) -> np.ndarray:
    """Return the past-slice configuration R0 (n_nodes, m) from the seed config."""
    ref = lattice.reference_positions(m)
    kind = seed_cfg.get("kind", "flat")
    amp = float(seed_cfg.get("amplitude", 1e-3))
    dim = lattice.params.dim

    if kind == "flat":
        return ref.copy()

    if kind == "random":
        rng = np.random.default_rng(int(seed_cfg.get("rng_seed", 0)))
        return ref + amp * rng.standard_normal((lattice.n_nodes, m))

    if kind == "plane_wave":
        grid = np.asarray(lattice.params.grid_shape, dtype=np.float64)
        a = lattice.params.spacing
        n_idx = np.asarray(seed_cfg.get("k_index", [1] + [0] * (dim - 1)), dtype=np.float64)
        if n_idx.shape[0] != dim:
            raise ValueError(f"k_index must have length dim={dim}, got {n_idx.shape[0]}")
        kvec = 2.0 * np.pi * n_idx / (grid * a)        # (dim,)
        phase = ref[:, :dim] @ kvec                    # (n_nodes,)
        pol = np.asarray(seed_cfg.get("polarization", [0.0, 1.0] + [0.0] * (m - 2)), dtype=np.float64)
        if pol.shape[0] != m:
            raise ValueError(f"polarization must have length m_ambient={m}, got {pol.shape[0]}")
        return ref + amp * np.cos(phase)[:, None] * pol[None, :]

    if kind in ("hedgehog", "skyrme_twisted", "axis_triplet"):
        from branesim.initialization.seeds import (
            hedgehog as _hedgehog,
            skyrme_twisted_hedgehog as _skyrme,
            axis_triplet as _axis_triplet,
        )
        u0 = float(seed_cfg.get("u0", amp))
        w = float(seed_cfg.get("w", 5.0))
        profile_shape = str(seed_cfg.get("profile_shape", "gaussian"))
        if kind == "hedgehog":
            R0, _ = _hedgehog(lattice, m, u0=u0, w=w, profile_shape=profile_shape)
        elif kind == "skyrme_twisted":
            tanh_steepness = float(seed_cfg.get("tanh_steepness", 3.0))
            R0, _ = _skyrme(
                lattice, m, u0=u0, w=w,
                profile_shape=profile_shape,
                tanh_steepness=tanh_steepness,
            )
        else:  # axis_triplet
            weights = seed_cfg.get("weights", None)
            R0, _ = _axis_triplet(
                lattice, m, u0=u0, w=w,
                weights=weights,
                profile_shape=profile_shape,
            )
        return R0

    raise ValueError(f"Unknown seed kind: {kind!r}")


def _memory_estimate_gb(n_slices: int, n_nodes: int, m: int, krylov_vectors: int = 30) -> float:
    """Rough JFNK working-set estimate in GB (ARCHITECTURE.md / DEPLOYMENT.md).

    One world-volume vector is (n_slices+1)*n_nodes*m*8 bytes; JFNK/GMRES holds
    ~`krylov_vectors` of them. This is the figure that drives instance sizing.
    """
    one_vec = (n_slices + 1) * n_nodes * m * 8.0
    return krylov_vectors * one_vec / 1e9


def run(config: dict, output_dir: Path) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)

    lcfg = config["lattice"]
    acfg = config["action"]
    scfg = config.get("seed", {"kind": "flat"})
    solver_cfg = config.get("solver", {"mode": "ivp"})

    lp = LatticeParams(
        grid_shape=tuple(int(v) for v in lcfg["grid_shape"]),
        spacing=float(lcfg.get("spacing", 1.0)),
        periodic_axes=tuple(bool(v) for v in lcfg.get("periodic_axes", [True] * len(lcfg["grid_shape"]))),
        axial_weight=float(lcfg.get("axial_weight", 1.0)),
    )
    lattice = SpacelikeLattice(lp)
    m = int(acfg.get("m_ambient", lp.dim + 1))
    r_t_cfg = acfg.get("r_t", None)
    ap = ActionParams(
        k_s=float(acfg.get("k_s", 1.0)),
        alpha=float(acfg.get("alpha", 0.2)),
        rho=float(acfg.get("rho", 1.0)),
        dt=float(acfg.get("dt", 0.1)),
        n_slices=int(acfg["n_slices"]),
        m_ambient=m,
        r_t=float(r_t_cfg) if r_t_cfg is not None else None,
    )
    mass = ap.rho * lp.spacing ** lp.dim
    N = ap.n_slices

    mem_gb = _memory_estimate_gb(N, lattice.n_nodes, m)
    print(f"[branesim] lattice dim={lp.dim} grid={lp.grid_shape} n_nodes={lattice.n_nodes} "
          f"m_ambient={m} n_slices={N}")
    print(f"[branesim] JFNK working-set estimate ~{mem_gb:.2f} GB "
          f"(={(N+1)*lattice.n_nodes*m*8/1e9:.3f} GB/vector x ~30)")

    R0 = _build_seed(scfg, lattice, m)
    mode = solver_cfg.get("mode", "ivp")
    t0 = time.perf_counter()

    if mode == "ivp":
        prob = IVPProblem(lattice=lattice, params=ap, mass=mass, R0=R0, R1=R0.copy())
        wv = march(prob)
        report = {"mode": "ivp"}

    elif mode == "bvp_dirichlet":
        # Ground-truth march supplies the two-time Dirichlet data (l=0, l=N).
        gt = march(IVPProblem(lattice=lattice, params=ap, mass=mass, R0=R0, R1=R0.copy()))
        bc = DirichletBC(R0=gt.slices[0].copy(), RN=gt.slices[N].copy())
        opts = SolveOpts(
            tol=float(solver_cfg.get("tol", 1e-9)),
            max_iter=int(solver_cfg.get("max_iter", 100)),
            warm_start=bool(solver_cfg.get("warm_start", True)),
        )
        wv = solve_block(BoundaryProblem(lattice, ap, mass, bc), opts)
        report = dict(wv.solver_report)
        report["recovery_max_abs_vs_march"] = float(np.abs(wv.slices[1:N] - gt.slices[1:N]).max())

    elif mode == "bvp_chiral":
        # Chiral Cauchy BC: build R1 from a one-step IVP march from R0.
        # For a plane-wave seed this gives the exact second eigenmode slice;
        # for other seeds it gives the stationary-start first step (zero vel).
        R1 = R0.copy()
        if scfg.get("kind") == "plane_wave":
            # One Verlet step from R0 with zero velocity to get R1.
            import math
            dim = lp.dim
            n_idx = np.asarray(scfg.get("k_index", [1] + [0] * (dim - 1)), dtype=np.float64)
            kvec = 2.0 * np.pi * n_idx / (np.array(lp.grid_shape, dtype=np.float64) * lp.spacing)
            pol_vec = np.asarray(scfg.get("polarization", [0.0, 1.0] + [0.0] * (m - 2)), dtype=np.float64)
            pol_axis = int(np.argmax(np.abs(pol_vec)))
            from branesim.core.conventions import d_of_k_eigenvalues
            eig = d_of_k_eigenvalues(kvec, ap.alpha, ap.k_s, ap.rho, lp.spacing)
            # arccos-domain guard only (math, not a physics clamp): a sub-CFL
            # propagating mode has 1 - 0.5 dt^2 omega^2 in [-1,1]; a super-CFL
            # (evanescent) mode would hit the clamp -> revisit dt if that happens.
            theta_k = math.acos(max(-1.0, min(1.0, 1.0 - 0.5 * ap.dt ** 2 * eig[pol_axis])))
            ref = lattice.reference_positions(m)
            phase = ref[:, :dim] @ kvec
            amp = float(scfg.get("amplitude", 1e-3))
            R1 = ref + amp * np.cos(phase - theta_k)[:, None] * pol_vec[None, :]
        bc = ChiralBC(R0=R0, R1=R1, chirality=solver_cfg.get("chirality", "forward"))
        wv = solve_block(BoundaryProblem(lattice, ap, mass, bc))
        report = dict(wv.solver_report)
    else:
        raise ValueError(f"Unknown solver mode: {mode!r}")

    walltime = time.perf_counter() - t0

    # Diagnostics summary (lightweight; full diagnostics port is future work).
    disp = wv.slices - lattice.reference_positions(m)[None, :, :]
    res_norm = float(residual_norm(wv.slices, lattice, ap, mass))
    summary = {
        "config": config,
        "mode": mode,
        "n_slices": N,
        "n_nodes": lattice.n_nodes,
        "m_ambient": m,
        "memory_estimate_gb": mem_gb,
        "walltime_s": walltime,
        "max_abs_displacement": float(np.abs(disp).max()),
        "interior_residual_norm": res_norm,
        "solver_report": report,
    }

    wv_path = output_dir / "worldvolume.zip"
    # Not using the context manager: we need to pass solver_report to close().
    # Include lattice metadata so that confinement_from_worldvolume can read dim
    # without needing the original config (contracts.py / confinement.py API).
    manifest_lattice = {
        "grid_shape": list(lp.grid_shape),
        "spacing": lp.spacing,
        "periodic_axes": list(lp.periodic_axes),
        "axial_weight": lp.axial_weight,
        "dim": lp.dim,
    }
    w = WorldVolumeWriter(wv_path, manifest_extra={"mode": mode, "lattice": manifest_lattice})
    for l in range(wv.slices.shape[0]):
        w.write_slice(l, l * ap.dt, wv.slices[l])
    w.write_npy("aux/ref_positions.npy", lattice.reference_positions(m))
    w.close(solver_report=report)

    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"[branesim] wrote {wv_path} and summary.json  (walltime {walltime:.1f}s)")
    return summary


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run a branesim block-solve experiment from a JSON config.")
    p.add_argument("--config", required=True, help="Path to JSON experiment config.")
    p.add_argument("--output-dir", required=True, help="Directory for worldvolume.zip + summary.json.")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    config = json.loads(Path(args.config).read_text(encoding="utf-8"))
    run(config, Path(args.output_dir))


if __name__ == "__main__":
    main()
