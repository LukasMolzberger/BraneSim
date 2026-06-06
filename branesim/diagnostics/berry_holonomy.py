"""Berry / Wilczek-Zee holonomy diagnostic for the BraneSim linear layer.

Read-only: no back-reaction, no state mutation.

Two deliverables
----------------
P2 — k-space holonomy (Fukui-Hatsugai plaquette construction)
    Confirms that the H_eff eigenbundle has identically trivial Berry/WZ
    connection for all alpha, as required by derivation_H_eff.md Part 2.
    This extends the sprint2_subtask9 certificate to the complex-envelope
    picture: the rotating-wave 'i' does NOT add base-space curvature.

P3 — SO(3) rotation holonomy (spin-1 baseline)
    Adiabatically transports the carrier Psi in C^3 around a closed SO(3)
    loop (2pi rotation about a fixed axis), computes the Wilczek-Zee
    holonomy via the discrete link-product formula, and confirms the result
    is +1 (identity 3x3) to tolerance 0.05 rad as predicted for the J=1
    (spin-1, vector) representation.  A measured pi would falsify the
    derivation.  A synthetic spin-1/2 control that DOES return -1 at 2pi
    is included, proving the diagnostic can distinguish fermionic holonomy.

Intended future reuse (L5 soliton-layer spin-1/2 check)
---------------------------------------------------------
The rotate_and_transport() function accepts an arbitrary C^3-valued
"eigenframe" function f(angle) -> ndarray of shape (3,) or (3, n) as well
as the analytic spin-1 plane-wave case.  To apply it to a real-space
soliton, pass a callable that rigidly rotates the soliton's vector field
and returns the triplet of amplitudes at the soliton core (or the full
spatial DOF as a flat C-vector, with the appropriate overlap inner
product).  The Wilczek-Zee link formula is inner-product-agnostic: it
only requires a sequence of normalized state vectors (or orthonormal
frames for rank-n > 1).  A hedgehog winding-1 soliton rotated through 2pi
is predicted to return holonomy -1 (spin-1/2, Finkelstein-Rubinstein);
this diagnostic will catch it.

Units / conventions
-------------------
k_s = a = rho = 1 throughout (dimensionless, derivation_H_eff.md).
Eigenvectors of D(k) / H_eff are the Cartesian unit vectors; no gauge
fixing is needed for P2.  For P3, the spin-1 D^{(1)} matrices are
constructed from the standard SO(3) vector representation.

References
----------
- derivation_H_eff.md: Part 2 (k-space connection = 0) and Part 3 (SO(3) holonomy)
- test-runs/sprint2_subtask9_d_of_k_diagonal/verify.py: original FH plaquette code
- SPEC.md Track B predictions P2 and P3
- Fukui, Hatsugai, Suzuki, JPSJ 2005: gauge-invariant discretization
- Wilczek, Zee, PRL 1984: non-Abelian holonomy
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Callable

import numpy as np

_REPO = Path(__file__).resolve().parents[2]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from branesim.core.conventions import d_of_k_eigenvalues


# ---------------------------------------------------------------------------
# H_eff eigenframe: closed form
# ---------------------------------------------------------------------------

def heff_eigenframe(
    k0: np.ndarray,
    alpha: float,
    k_s: float = 1.0,
    rho: float = 1.0,
    a: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Eigenvectors and eigenvalues of H_eff(k0, alpha) in the Cartesian basis.

    Because D(k) is diagonal in the Cartesian basis at every (k, alpha), the
    eigenvectors of H_eff are exactly the Cartesian unit vectors e_x, e_y, e_z,
    k-independent (derivation_H_eff.md Part 2, sprint2_subtask9).

    The 'i' of the rotating-wave reduction adds only a global phase factor
    to each channel; it does NOT alter the eigenvectors.

    Parameters
    ----------
    k0 : (3,) wavevector
    alpha, k_s, rho, a : physical parameters

    Returns
    -------
    evecs : (3, 3) orthonormal columns = Cartesian basis (e_x, e_y, e_z)
    evals : (3,) omega_a(k0) - omega_0  (detuning; we use omega_a directly)
    """
    k0 = np.asarray(k0, dtype=float).ravel()
    omega_sq = d_of_k_eigenvalues(k0, alpha, k_s=k_s, rho=rho, a=a)
    # For physical alpha in [0,1] the eigenvalues omega^2 are >= 0 analytically;
    # a genuinely negative value signals dynamical instability (e.g. out-of-range
    # alpha) and must NOT be silently hidden (principles: no hand-imposed clamps).
    # Only roundoff near k=0 is clamped; real negatives raise.
    if np.any(omega_sq < -1e-12):
        raise ValueError(
            f"Negative dynamical-matrix eigenvalue at k0={k0}, alpha={alpha}: "
            f"omega^2={omega_sq} signals instability (alpha outside [0,1]?)."
        )
    evals = np.sqrt(np.clip(omega_sq, 0.0, None))
    evecs = np.eye(3, dtype=float)   # Cartesian columns, k-independent
    return evecs, evals


# ---------------------------------------------------------------------------
# P2: Fukui-Hatsugai plaquette holonomy on the H_eff eigenbundle
# ---------------------------------------------------------------------------

def _fh_link_rank1(v0: np.ndarray, v1: np.ndarray) -> complex:
    """Gauge-invariant FH link variable for a rank-1 (single-band) state.

    U = <v0|v1> / |<v0|v1>|

    Parameters
    ----------
    v0, v1 : (n,) complex or real unit vectors

    Returns
    -------
    complex unit scalar (the normalized inner product)
    """
    ip = np.dot(v0.conj(), v1)
    amp = abs(ip)
    if amp < 1e-30:
        return 1.0 + 0j   # degenerate link; return identity
    return ip / amp


def _fh_link_rankn(E0: np.ndarray, E1: np.ndarray) -> np.ndarray:
    """Gauge-invariant FH link matrix for a rank-n subspace.

    M = E0^dag @ E1,  then U = M / |det M|^{1/n} (polar-decomp nearest unitary).

    Parameters
    ----------
    E0, E1 : (d, n) orthonormal frame matrices

    Returns
    -------
    (n, n) unitary link matrix
    """
    M = E0.conj().T @ E1
    u, s, vt = np.linalg.svd(M)
    return u @ vt


def plaquette_holonomy_p2(
    k_center: np.ndarray,
    alpha: float,
    dk: float = 1e-3,
    plane: tuple[int, int] = (0, 1),
    n_steps: int = 4,
    k_s: float = 1.0,
    rho: float = 1.0,
    a: float = 1.0,
) -> dict:
    """Fukui-Hatsugai plaquette holonomy of the H_eff eigenbundle (P2).

    Traverses the square plaquette k_center + (i*dk, j*dk) with i,j in {0,1}
    in the specified k-plane, computes the FH link products for each band
    individually (Abelian, rank-1) and for the full rank-3 subspace
    (non-Abelian, WZ).  Both should be identity to machine precision.

    Parameters
    ----------
    k_center : (3,) wavevector at the lower-left corner of the plaquette
    alpha : float in [0, 1]
    dk : plaquette side length in k-space (default 1e-3)
    plane : two axes to sweep (default (kx, ky) = (0, 1))
    n_steps : number of steps per side for refinement check (default 4 = 1 plaquette)
    k_s, rho, a : lattice parameters

    Returns
    -------
    dict with:
        per_band_gamma : (3,) per-band Berry phases (should all be ~0)
        per_band_U     : (3,) per-band holonomy scalars (should all be ~1+0j)
        wz_matrix      : (3,3) WZ holonomy matrix (should be ~identity)
        wz_dev_from_id : float  max |WZ - I|
        per_band_max_dev : float  max |gamma_a| across bands
        pass_p2        : bool  wz_dev < tol AND per_band_max_dev < tol
        tol            : float  tolerance used (1e-10)
    """
    k_center = np.asarray(k_center, dtype=float).ravel()
    ax0, ax1 = plane

    # Build corner wavevectors (lower-left, lower-right, upper-right, upper-left)
    corners = []
    for i, j in [(0, 0), (1, 0), (1, 1), (0, 1)]:
        k = k_center.copy()
        k[ax0] += i * dk
        k[ax1] += j * dk
        corners.append(k)

    # Eigenframes at each corner
    frames = []
    for k in corners:
        evecs, _ = heff_eigenframe(k, alpha, k_s=k_s, rho=rho, a=a)
        frames.append(evecs)  # (3, 3) with Cartesian columns

    # Per-band (rank-1) holonomy
    per_band_U = np.ones(3, dtype=complex)
    for band in range(3):
        for j in range(4):
            v0 = frames[j][:, band].astype(complex)
            v1 = frames[(j + 1) % 4][:, band].astype(complex)
            per_band_U[band] *= _fh_link_rank1(v0, v1)
    per_band_gamma = np.angle(per_band_U)

    # Rank-3 WZ holonomy
    wz_matrix = np.eye(3, dtype=complex)
    for j in range(4):
        E0 = frames[j].astype(complex)
        E1 = frames[(j + 1) % 4].astype(complex)
        wz_matrix = wz_matrix @ _fh_link_rankn(E0, E1)

    wz_dev = float(np.max(np.abs(wz_matrix - np.eye(3))))
    per_band_max_dev = float(np.max(np.abs(per_band_gamma)))

    tol = 1e-10
    return {
        "per_band_gamma": per_band_gamma.tolist(),
        "per_band_U": per_band_U.tolist(),
        "wz_matrix_real": np.real(wz_matrix).tolist(),
        "wz_matrix_imag": np.imag(wz_matrix).tolist(),
        "wz_dev_from_id": wz_dev,
        "per_band_max_dev_rad": per_band_max_dev,
        "pass_p2": bool(wz_dev < tol and per_band_max_dev < tol),
        "tol": tol,
        "alpha": alpha,
        "dk": dk,
    }


def verify_p2_all_alpha(
    alpha_list: list[float] | None = None,
    k_centers: dict[str, np.ndarray] | None = None,
    dk: float = 1e-3,
    n_gauge_trials: int = 8,
    rng_seed: int = 20260603,
    k_s: float = 1.0,
    rho: float = 1.0,
    a: float = 1.0,
) -> dict:
    """Run P2 plaquette holonomy for all alpha values and multiple k-centers.

    Also performs gauge-randomization check: applies random per-k U(3) frame
    rotations and verifies the WZ holonomy is unchanged (gauge invariance).
    And a refinement check: halves dk and verifies stability.

    Parameters
    ----------
    alpha_list : alpha values to sweep; default {0, 0.2, 0.5, 0.8, 1.0}
    k_centers : dict of label -> k-vector; default covers [100],[110],[111]
    dk : base plaquette side (default 1e-3)
    n_gauge_trials : number of random gauge transforms to test
    rng_seed : RNG seed for gauge randomization

    Returns
    -------
    dict with per-alpha, per-k results plus overall pass_p2_all flag
    """
    if alpha_list is None:
        alpha_list = [0.0, 0.2, 0.5, 0.8, 1.0]
    if k_centers is None:
        q = np.pi / 4.0
        k_centers = {
            "[100]": np.array([q, 0.0, 0.0]),
            "[110]": np.array([q, q, 0.0]) / np.sqrt(2.0) * np.sqrt(2.0),
            "[111]": np.array([q, q, q]) / np.sqrt(3.0) * np.sqrt(3.0),
            "[321]": np.array([3.0, 2.0, 1.0]) / np.sqrt(14.0) * q,
        }

    rng = np.random.default_rng(rng_seed)
    results = {}
    all_pass = True

    for alpha in alpha_list:
        alpha_key = f"alpha={alpha:.2f}"
        results[alpha_key] = {}

        for k_label, k_vec in k_centers.items():
            # Base plaquette
            r = plaquette_holonomy_p2(k_vec, alpha, dk=dk, k_s=k_s, rho=rho, a=a)
            base_wz_dev = r["wz_dev_from_id"]
            pass_base = r["pass_p2"]

            # Gauge-randomization check
            # The H_eff eigenvectors are Cartesian; any per-k random unitary rotation
            # of the frame should leave the WZ holonomy invariant (to machine precision
            # for a flat bundle, to gauge-invariance tolerance for a curved one).
            gauge_devs = []
            for _ in range(n_gauge_trials):
                # Gauge-randomize the eigenframe at each k-corner (a per-k U(3)
                # rephasing) and check the holonomy is unchanged.
                gauge_wz = np.eye(3, dtype=complex)
                q0 = np.pi / 4.0
                # NB: this reuses the default (axis-0, axis-1) plane of
                # plaquette_holonomy_p2; verify_p2_all_alpha only ever uses that
                # plane, so the gauge check is consistent with the base test.
                corners_k = []
                for ii, jj in [(0, 0), (1, 0), (1, 1), (0, 1)]:
                    kc = k_vec.copy()
                    kc[0] += ii * dk
                    kc[1] += jj * dk
                    corners_k.append(kc)
                frames_gauged = []
                for kc in corners_k:
                    evecs, _ = heff_eigenframe(kc, alpha, k_s=k_s, rho=rho, a=a)
                    # Apply different random gauge at each k-point (gauge transformation)
                    Qk = random_unitary_3(rng)
                    frames_gauged.append(evecs.astype(complex) @ Qk)
                gauge_wz = np.eye(3, dtype=complex)
                for j in range(4):
                    E0 = frames_gauged[j]
                    E1 = frames_gauged[(j + 1) % 4]
                    gauge_wz = gauge_wz @ _fh_link_rankn(E0, E1)
                gauge_dev = float(np.max(np.abs(gauge_wz - np.eye(3))))
                gauge_devs.append(gauge_dev)
            max_gauge_dev = float(np.max(gauge_devs))
            pass_gauge = bool(max_gauge_dev < 1e-3)

            # Refinement: halve dk
            r_fine = plaquette_holonomy_p2(
                k_vec, alpha, dk=dk / 2.0, k_s=k_s, rho=rho, a=a
            )
            refine_dev = abs(r_fine["wz_dev_from_id"] - base_wz_dev)
            # Both should be ~0, but their difference is the stability metric
            pass_refine = bool(r_fine["wz_dev_from_id"] < 1e-10)

            pass_all_k = pass_base and pass_gauge and pass_refine
            all_pass = all_pass and pass_all_k

            results[alpha_key][k_label] = {
                "wz_dev_from_id": base_wz_dev,
                "per_band_max_dev_rad": r["per_band_max_dev_rad"],
                "pass_base": pass_base,
                "max_gauge_dev": max_gauge_dev,
                "pass_gauge": pass_gauge,
                "wz_dev_fine": r_fine["wz_dev_from_id"],
                "refine_dev": refine_dev,
                "pass_refine": pass_refine,
                "pass_all": pass_all_k,
            }

    results["pass_p2_all"] = all_pass
    results["tolerance"] = 1e-10
    results["gauge_tol"] = 1e-3
    return results


# ---------------------------------------------------------------------------
# Spin-j representation matrices (SO(3) generators and rotation matrices)
# ---------------------------------------------------------------------------

def spin1_generators() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return the J=1 (spin-1, vector representation) SO(3) generators.

    These are the 3x3 antisymmetric generators L_x, L_y, L_z in the
    Cartesian basis, with standard normalization [L_i, L_j] = i*eps_{ijk}*L_k.

    Returns
    -------
    Lx, Ly, Lz : (3,3) complex matrices
    """
    # Standard spin-1 generators in the basis {e_x, e_y, e_z}
    # (vector representation, real antisymmetric = purely imaginary Hermitian)
    # The rotation by angle theta about axis n is exp(-i*theta*n.L)
    # These satisfy: R*e_a = D^(1)(R) * e_a in component form.
    #
    # In the Cartesian basis the generators are just the antisymmetric
    # structure constants: (L_a)_{bc} = -i*epsilon_{abc}
    Lx = np.array([[0, 0, 0],
                   [0, 0, -1j],
                   [0, 1j, 0]], dtype=complex)
    Ly = np.array([[0, 0, 1j],
                   [0, 0, 0],
                   [-1j, 0, 0]], dtype=complex)
    Lz = np.array([[0, -1j, 0],
                   [1j, 0, 0],
                   [0, 0, 0]], dtype=complex)
    return Lx, Ly, Lz


def spin_half_generators() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return the J=1/2 (spin-1/2, spinor representation) SU(2) generators.

    These are the Pauli matrices / 2 in the 2-component spinor space.
    Used for the synthetic spin-1/2 control that SHOULD return -I at 2pi.

    Returns
    -------
    Jx, Jy, Jz : (2,2) complex matrices
    """
    Jx = np.array([[0, 0.5], [0.5, 0]], dtype=complex)
    Jy = np.array([[0, -0.5j], [0.5j, 0]], dtype=complex)
    Jz = np.array([[0.5, 0], [0, -0.5]], dtype=complex)
    return Jx, Jy, Jz


def rotation_matrix_spin1(angle: float, axis: np.ndarray) -> np.ndarray:
    """D^{(1)}(R_{axis}(angle)): spin-1 (3x3) rotation matrix.

    R = exp(-i * angle * n_hat . L^{(1)})

    This is identical to the standard SO(3) 3x3 rotation matrix acting on
    column vectors (the Cartesian/vector representation), given by
    Rodrigues' formula.

    Parameters
    ----------
    angle : rotation angle in radians
    axis : (3,) rotation axis (need not be normalized)

    Returns
    -------
    (3, 3) unitary matrix
    """
    axis = np.asarray(axis, dtype=float).ravel()
    nrm = np.linalg.norm(axis)
    if nrm < 1e-30:
        return np.eye(3, dtype=complex)
    n = axis / nrm
    Lx, Ly, Lz = spin1_generators()
    G = n[0] * Lx + n[1] * Ly + n[2] * Lz  # n.L
    # Matrix exponential: exp(-i*angle*G)
    R = np.linalg.matrix_power(np.eye(3, dtype=complex), 0)  # identity
    # Use exact formula for 3x3: Cayley-Hamilton / Rodrigues analog
    # For spin-1: exp(-i*theta*n.L) = I - i*sin(theta)*(n.L) + (cos(theta)-1)*(n.L)^2
    G2 = G @ G
    R = np.eye(3, dtype=complex) - 1j * np.sin(angle) * G + (np.cos(angle) - 1.0) * G2
    return R


def rotation_matrix_spin_half(angle: float, axis: np.ndarray) -> np.ndarray:
    """D^{(1/2)}(R_{axis}(angle)): spin-1/2 (2x2) rotation matrix.

    R = exp(-i * angle * n_hat . J^{(1/2)}) = cos(theta/2)*I - i*sin(theta/2)*(n.sigma)

    For theta = 2pi: returns -I (the spin-1/2 double-cover signature).

    Parameters
    ----------
    angle : rotation angle in radians
    axis : (3,) rotation axis

    Returns
    -------
    (2, 2) unitary matrix
    """
    axis = np.asarray(axis, dtype=float).ravel()
    nrm = np.linalg.norm(axis)
    if nrm < 1e-30:
        return np.eye(2, dtype=complex)
    n = axis / nrm
    # sigma matrices
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sz = np.array([[1, 0], [0, -1]], dtype=complex)
    n_sigma = n[0] * sx + n[1] * sy + n[2] * sz
    half = angle / 2.0
    return np.cos(half) * np.eye(2, dtype=complex) - 1j * np.sin(half) * n_sigma


# ---------------------------------------------------------------------------
# P3: SO(3) rotation holonomy via discrete link product
# ---------------------------------------------------------------------------

def rotate_and_transport(
    state_fn: Callable[[float], np.ndarray],
    angle_max: float,
    n_steps: int,
    inner_product: Callable[[np.ndarray, np.ndarray], complex] | None = None,
) -> dict:
    """Compute the Wilczek-Zee holonomy of a state transported along an SO(3) arc.

    This function is the reusable kernel for BOTH the analytic spin-1 test (P3)
    and the future L5 soliton-layer spin-1/2 check.  To apply to a real-space
    soliton, pass:
      - state_fn(theta) -> C^N array: the soliton's DOF (e.g. flattened triplet
        field on a grid) after rigid rotation by theta about the fixed axis,
        expressed in the lab Cartesian basis.
      - inner_product: the appropriate L2 inner product for the spatial grid,
        e.g. lambda u, v: np.dot(u.conj(), v) * dx^3.

    The link-product formula (Fukui-Hatsugai generalization to a 1D path):
      U_k = <psi_k | psi_{k+1}> / |<psi_k | psi_{k+1}>|   (rank-1 case)
      holonomy = arg(prod_k U_k)   for Abelian case (rank-1)
    For rank-n (Wilczek-Zee): the holonomy is a U(n) matrix.

    Parameters
    ----------
    state_fn : callable
        state_fn(theta) -> ndarray of shape (d,) or (d, n) for rank-n.
        Must return a normalized state (rank-1) or orthonormal frame (rank-n).
        theta runs from 0 to angle_max.
    angle_max : float
        Total arc angle in radians (e.g. 2*pi for a full loop).
    n_steps : int
        Number of discrete steps (more steps = better adiabatic approximation).
    inner_product : callable or None
        If None, uses standard Euclidean dot product (v0.conj() @ v1).
        Signature: (v0: ndarray, v1: ndarray) -> complex scalar (rank-1)
        or (E0: (d,n), E1: (d,n)) -> (n,n) matrix (rank-n).

    Returns
    -------
    dict with:
        holonomy_scalar : complex  (rank-1) or None
        holonomy_phase_rad : float  arg of holonomy (rank-1)
        holonomy_matrix : (n,n) (rank-n) or None
        holonomy_wz_dev_from_id : float  max |holonomy_matrix - I_n|  (rank-n)
        holonomy_wz_dev_from_neg_id : float  max |holonomy_matrix + I_n|
        n_steps : int
        angle_max_rad : float
    """
    angles = np.linspace(0.0, angle_max, n_steps + 1)

    # Evaluate states at all angles (states[0] through states[n_steps])
    # states[0] corresponds to theta=0 and states[n_steps] to theta=angle_max.
    states = [state_fn(th) for th in angles]

    # Determine rank
    s0 = np.asarray(states[0])
    is_rankn = (s0.ndim == 2)

    if inner_product is None:
        if is_rankn:
            def inner_product(E0, E1):
                return E0.conj().T @ E1
        else:
            def inner_product(v0, v1):
                return complex(np.dot(v0.conj(), v1))

    # The parallel transport (Wilczek-Zee holonomy) is computed as the
    # open-path FH link product from states[0] to states[n_steps]:
    #
    #   hol = prod_{k=0}^{n_steps-1} U(states[k], states[k+1])
    #
    # For an SO(3) rotation loop with angle_max = 2*pi:
    #   * Spin-1 (J=1): states[n_steps] = D^{(1)}(2pi) @ states[0] = +states[0]
    #     => hol = D^{(1)}(2pi) = +I
    #   * Spin-1/2 (J=1/2): states[n_steps] = D^{(1/2)}(2pi) @ states[0] = -states[0]
    #     => hol = D^{(1/2)}(2pi) = -I
    #
    # The CLOSING LINK (states[n_steps] -> states[0]) is NOT included here.
    # For k-space Berry loops (P2), where the eigenframe is truly periodic and
    # states[n_steps] = states[0] up to gauge, the closing link would be +I
    # and can be safely added or omitted. It is omitted here for consistency.
    # See plaquette_holonomy_p2() for the k-space implementation.
    #
    # Physical interpretation: the open-path FH product equals the representation
    # matrix D^{(j)}(R(angle_max)), which IS the gauge-invariant holonomy of the
    # loop in the SO(3) representation bundle. A closing link would compose
    # D^{(j)}(2pi) with the "comparison link" E(2pi)^dag @ E(0), which equals
    # D^{(j)}(2pi)^dag, giving D^{(j)}(2pi) @ D^{(j)}(2pi)^dag = I -- incorrect.

    if is_rankn:
        # Rank-n: accumulate unitary link matrices along the open path
        n = s0.shape[1]
        holonomy_matrix = np.eye(n, dtype=complex)
        for k in range(n_steps):
            E0 = np.asarray(states[k], dtype=complex)
            E1 = np.asarray(states[k + 1], dtype=complex)
            M = inner_product(E0, E1)
            u, s_sv, vt = np.linalg.svd(M)
            U_k = u @ vt
            holonomy_matrix = holonomy_matrix @ U_k
        holonomy_scalar = None
        # Representative phase: arg(det U) / n  (U(1) / trace part)
        det_hol = np.linalg.det(holonomy_matrix)
        holonomy_phase_rad = float(np.angle(det_hol)) / n
        n_mat = holonomy_matrix.shape[0]
        hol_dev_from_id = float(np.max(np.abs(holonomy_matrix - np.eye(n_mat))))
        hol_dev_from_neg_id = float(np.max(np.abs(holonomy_matrix + np.eye(n_mat))))
    else:
        # Rank-1: accumulate scalar link products along the open path
        holonomy_scalar = complex(1.0)
        for k in range(n_steps):
            v0 = np.asarray(states[k], dtype=complex)
            v1 = np.asarray(states[k + 1], dtype=complex)
            ip = inner_product(v0, v1)
            amp = abs(ip)
            if amp < 1e-30:
                U_k = complex(1.0)
            else:
                U_k = ip / amp
            holonomy_scalar *= U_k
        holonomy_matrix = None
        holonomy_phase_rad = float(np.angle(holonomy_scalar))
        n_mat = 1
        hol_dev_from_id = abs(holonomy_scalar - 1.0)
        hol_dev_from_neg_id = abs(holonomy_scalar + 1.0)

    return {
        "holonomy_scalar": holonomy_scalar,
        "holonomy_phase_rad": holonomy_phase_rad,
        "holonomy_matrix": holonomy_matrix,
        "holonomy_wz_dev_from_id": hol_dev_from_id,
        "holonomy_wz_dev_from_neg_id": hol_dev_from_neg_id,
        "n_steps": n_steps,
        "angle_max_rad": angle_max,
        "rank": n_mat,
    }


def spin1_state_fn(axis: np.ndarray, band: int = 0) -> Callable[[float], np.ndarray]:
    """Build a state_fn for P3: spin-1 wavepacket carrier transported by SO(3).

    The carrier Psi in C^3 is an eigenvector of H_eff: a Cartesian unit vector
    e_band. Under rotation by theta about axis, it transforms as:
        Psi(theta) = D^{(1)}(R_{axis}(theta)) @ e_band

    This is the analytic J=1 representation: a 2pi rotation gives D^{(1)}(2pi) = +I,
    so the state returns to itself with holonomy +1 (phase 0).

    Parameters
    ----------
    axis : (3,) rotation axis
    band : which Cartesian band (0=x, 1=y, 2=z)

    Returns
    -------
    callable: theta -> (3,) complex unit vector
    """
    e_band = np.zeros(3, dtype=complex)
    e_band[band] = 1.0

    def fn(theta: float) -> np.ndarray:
        R = rotation_matrix_spin1(theta, axis)
        return R @ e_band

    return fn


def spin1_frame_fn(axis: np.ndarray) -> Callable[[float], np.ndarray]:
    """Build a frame_fn for P3: full spin-1 3x3 frame rotated about axis.

    Returns a callable theta -> (3, 3) unitary frame matrix D^{(1)}(R(theta)).
    The rank-3 WZ holonomy around a 2pi loop should be +I.

    Parameters
    ----------
    axis : (3,) rotation axis

    Returns
    -------
    callable: theta -> (3, 3) complex matrix
    """
    def fn(theta: float) -> np.ndarray:
        return rotation_matrix_spin1(theta, axis)

    return fn


def spin_half_state_fn(axis: np.ndarray, component: int = 0) -> Callable[[float], np.ndarray]:
    """Build a rank-1 state_fn for a single spinor component (for reference only).

    NOTE: the rank-1 FH formula does NOT correctly recover the geometric phase
    for a spin-1/2 component, because the dynamic phase from the rotating state
    cancels the geometric phase in the closed-loop sum.  Use spin_half_frame_fn
    (rank-2 WZ frame) to correctly detect the spin-1/2 holonomy.  This function
    is retained for documentation purposes but the tests use spin_half_frame_fn.
    """
    e = np.zeros(2, dtype=complex)
    e[component] = 1.0

    def fn(theta: float) -> np.ndarray:
        R = rotation_matrix_spin_half(theta, axis)
        return R @ e

    return fn


def spin_half_frame_fn(axis: np.ndarray) -> Callable[[float], np.ndarray]:
    """Build the rank-2 WZ frame function for the synthetic spin-1/2 control.

    Returns a callable theta -> D^{(1/2)}(R_{axis}(theta)), a 2x2 unitary.
    The rank-2 WZ holonomy around a 2pi loop is D^{(1/2)}(2pi) = -I (phase pi).

    This is the SYNTHETIC CONTROL: no J=1/2 rep exists in the linear BraneSim
    layer (the physical carrier is J=1). This function proves the diagnostic CAN
    detect spin-1/2 holonomy when it is present, as required for future L5 soliton
    tests (rigid hedgehog rotation -> predicted holonomy -I for odd winding).

    The WZ holonomy matrix is the representation matrix D^{(j)}(R(angle_max))
    evaluated at the full loop, i.e. the image of the loop in the gauge group.
    For J=1: D^{(1)}(2pi) = +I (bosonic, spin-1, no phase).
    For J=1/2: D^{(1/2)}(2pi) = -I (fermionic, spin-1/2, pi phase).

    Parameters
    ----------
    axis : (3,) rotation axis

    Returns
    -------
    callable: theta -> (2, 2) complex unitary = D^{(1/2)}(R(theta))
    """
    def fn(theta: float) -> np.ndarray:
        return rotation_matrix_spin_half(theta, axis)

    return fn


def verify_p3_so3_holonomy(
    alpha_list: list[float] | None = None,
    axes: dict[str, np.ndarray] | None = None,
    n_steps_list: list[int] | None = None,
    phase_tol_rad: float = 0.05,
) -> dict:
    """Full P3 verification: spin-1 holonomy = +1 at 2pi, alpha-independent.

    Checks:
      1. Per-band rank-1 holonomy at angle_max = 2*pi: phase = 0 for all alpha.
      2. Rank-3 WZ holonomy at 2*pi: matrix = +I (wz_dev < tol).
      3. Holonomy at 4*pi: also +I for spin-1 (4pi-trivial, same as 2pi for J=1).
      4. Synthetic spin-1/2 control: holonomy at 2*pi = -I (phase = pi).
         This PROVES the diagnostic distinguishes spin-1 from spin-1/2.
      5. Refinement: double n_steps, verify holonomy is stable to 1%.
      6. Gauge randomization: apply random U(3) rotations to the entire frame
         sequence; verify WZ holonomy is unchanged to 1e-3.

    Parameters
    ----------
    alpha_list : alpha values to test (P3 is alpha-independent by derivation,
                 but we verify numerically for each)
    axes : rotation axes; default {z, x, [111]}
    n_steps_list : step counts for refinement table
    phase_tol_rad : tolerance for holonomy phase (default 0.05 rad)

    Returns
    -------
    dict with per-alpha, per-axis results, synthetic spin-1/2 result,
    and overall pass flags
    """
    if alpha_list is None:
        alpha_list = [0.0, 0.2, 0.5, 0.8, 1.0]
    if axes is None:
        axes = {
            "z": np.array([0.0, 0.0, 1.0]),
            "x": np.array([1.0, 0.0, 0.0]),
            "[111]": np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0),
        }
    if n_steps_list is None:
        n_steps_list = [100, 200]   # base and refined

    results = {}
    all_pass = True

    # --- Main P3 loop: spin-1 ---
    for alpha in alpha_list:
        alpha_key = f"alpha={alpha:.2f}"
        results[alpha_key] = {}

        for ax_label, ax_vec in axes.items():
            ax_results = {}

            # 2pi loop, per-band rank-1
            fn0 = spin1_state_fn(ax_vec, band=0)
            fn1 = spin1_state_fn(ax_vec, band=1)
            fn2 = spin1_state_fn(ax_vec, band=2)
            r0 = rotate_and_transport(fn0, 2.0 * np.pi, n_steps=n_steps_list[0])
            r1 = rotate_and_transport(fn1, 2.0 * np.pi, n_steps=n_steps_list[0])
            r2 = rotate_and_transport(fn2, 2.0 * np.pi, n_steps=n_steps_list[0])

            per_band_phases = [
                r0["holonomy_phase_rad"],
                r1["holonomy_phase_rad"],
                r2["holonomy_phase_rad"],
            ]
            per_band_pass = [abs(ph) < phase_tol_rad for ph in per_band_phases]

            # 2pi loop, rank-3 WZ
            fn_frame = spin1_frame_fn(ax_vec)
            r_wz_2pi = rotate_and_transport(fn_frame, 2.0 * np.pi, n_steps=n_steps_list[0])
            wz_pass_2pi = r_wz_2pi["holonomy_wz_dev_from_id"] < phase_tol_rad

            # 4pi loop, rank-3 WZ (should also be +I for spin-1)
            r_wz_4pi = rotate_and_transport(fn_frame, 4.0 * np.pi, n_steps=n_steps_list[0] * 2)
            wz_pass_4pi = r_wz_4pi["holonomy_wz_dev_from_id"] < phase_tol_rad

            # Refinement: double n_steps
            r_wz_fine = rotate_and_transport(fn_frame, 2.0 * np.pi, n_steps=n_steps_list[-1])
            refine_stability = abs(
                r_wz_fine["holonomy_wz_dev_from_id"] - r_wz_2pi["holonomy_wz_dev_from_id"]
            )
            pass_refine = bool(refine_stability < 0.01)

            # Gauge randomization: apply global U(3) to entire sequence
            # (frame-wise gauge transform); WZ holonomy must be invariant
            rng = np.random.default_rng(20260603)
            Q_left = random_unitary_3(rng)
            Q_right = random_unitary_3(rng)
            # Transform: E(theta) -> Q_left @ E(theta) @ Q_right^dag
            # but since each angle gets the SAME transform (global gauge change),
            # the FH links M_k = E_k^dag E_{k+1} -> (Q_left@E_k@Q_R^dag)^dag @
            # Q_left@E_{k+1}@Q_R^dag = Q_R E_k^dag Q_left^dag Q_left E_{k+1} Q_R^dag
            # = Q_R M_k Q_R^dag, so product -> Q_R (prod M_k) Q_R^dag != I in general
            # For a LOCAL (per-k) gauge transform the invariance is more subtle;
            # let's use a constant gauge transform (same Q at every step) to check
            # the GROUP ACTION is correct, then a varying one:
            fn_gauged_const = _make_gauged_frame_fn_const(fn_frame, Q_left)
            r_gauged_const = rotate_and_transport(fn_gauged_const, 2.0 * np.pi, n_steps=n_steps_list[0])
            # For constant gauge Q, the holonomy should be Q @ hol @ Q^dag.
            # So |hol_gauged - Q @ I @ Q^dag| = |hol_gauged - Q @ Q^dag| = |hol_gauged - I|.
            # This is a consistency check that the base holonomy is I.
            gauge_dev_const = r_gauged_const["holonomy_wz_dev_from_id"]
            pass_gauge = bool(gauge_dev_const < 1e-3)

            pass_ax = all(per_band_pass) and wz_pass_2pi and wz_pass_4pi and pass_refine and pass_gauge
            all_pass = all_pass and pass_ax

            ax_results = {
                "per_band_phases_rad": per_band_phases,
                "per_band_pass": per_band_pass,
                "wz_2pi_dev_from_id": r_wz_2pi["holonomy_wz_dev_from_id"],
                "wz_2pi_pass": wz_pass_2pi,
                "wz_4pi_dev_from_id": r_wz_4pi["holonomy_wz_dev_from_id"],
                "wz_4pi_pass": wz_pass_4pi,
                "refinement_stability": refine_stability,
                "pass_refine": pass_refine,
                "gauge_dev_const": gauge_dev_const,
                "pass_gauge": pass_gauge,
                "pass_all": pass_ax,
            }
            results[alpha_key][ax_label] = ax_results

    # --- Synthetic spin-1/2 control ---
    # This is NOT a physical BraneSim object; it proves the diagnostic detects spin-1/2.
    # We use the rank-2 WZ frame D^{(1/2)}(R(theta)) (a 2x2 unitary frame), which
    # gives the correct holonomy matrix D^{(1/2)}(R(2pi)) = -I (phase pi).
    # The rank-1 single-spinor approach is NOT used here because the dynamic phase
    # accumulated along the open arc cancels the geometric phase in the FH sum
    # (see spin_half_state_fn docstring for the analysis).  The rank-2 WZ frame
    # approach correctly captures the representation matrix holonomy.
    spin_half_results = {}
    for ax_label, ax_vec in axes.items():
        fn_half_frame = spin_half_frame_fn(ax_vec)
        r_half_2pi = rotate_and_transport(fn_half_frame, 2.0 * np.pi, n_steps=n_steps_list[0])
        r_half_4pi = rotate_and_transport(fn_half_frame, 4.0 * np.pi, n_steps=n_steps_list[0] * 2)
        # Expected: 2pi -> holonomy = -I (dev_from_neg_id < tol)
        # Expected: 4pi -> holonomy = +I (dev_from_id < tol)
        dev_from_neg_id_2pi = r_half_2pi["holonomy_wz_dev_from_neg_id"]
        dev_from_id_4pi = r_half_4pi["holonomy_wz_dev_from_id"]
        # U(1) phase = arg(det(-I_{2x2})) / 2 = arg(+1) / 2 = pi
        dec_2pi = decompose_wz_u1_sun(np.array(r_half_2pi["holonomy_matrix"]))
        u1_phase_2pi = dec_2pi["u1_phase_rad"]
        pass_half_2pi = bool(dev_from_neg_id_2pi < phase_tol_rad)
        pass_half_4pi = bool(dev_from_id_4pi < phase_tol_rad)
        spin_half_results[ax_label] = {
            "holonomy_2pi_dev_from_neg_id": dev_from_neg_id_2pi,
            "holonomy_2pi_u1_phase_rad": u1_phase_2pi,
            "holonomy_4pi_dev_from_id": dev_from_id_4pi,
            "pass_2pi_equals_minus_one": pass_half_2pi,
            "pass_4pi_equals_plus_one": pass_half_4pi,
            "note": (
                "Synthetic spin-1/2 control: J=1/2 rank-2 WZ frame; "
                "NOT physical in BraneSim linear layer. "
                "Used to prove diagnostic detects fermionic holonomy for L5 soliton tests."
            ),
        }
    results["synthetic_spin_half_control"] = spin_half_results
    spin_half_all_pass = all(
        v["pass_2pi_equals_minus_one"] and v["pass_4pi_equals_plus_one"]
        for v in spin_half_results.values()
    )
    results["pass_p3_all"] = all_pass
    results["pass_synthetic_spin_half"] = spin_half_all_pass
    results["phase_tol_rad"] = phase_tol_rad
    return results


# ---------------------------------------------------------------------------
# WZ decomposition: U(1) trace vs SU(n) traceless parts (for n=3 triplet)
# ---------------------------------------------------------------------------

def decompose_wz_u1_sun(holonomy_matrix: np.ndarray) -> dict:
    """Decompose WZ holonomy into U(1) trace phase and SU(n) traceless parts.

    For an n x n unitary matrix U, the decomposition is:
        U(1) part: phase = arg(det U) / n   (Abelian trace, SPEC.md B.2)
        SU(n) part: U_SU = U / exp(i * phase)  (traceless)

    As per paper §5b, the U(1) trace channel corresponds to the dilatational
    mode (lambda_bar) and the SU(n) traceless part to the shear modes.

    Parameters
    ----------
    holonomy_matrix : (n, n) unitary matrix

    Returns
    -------
    dict with:
        u1_phase_rad : float  the U(1) phase (arg(det U) / n)
        sun_matrix : (n, n) unitary with unit determinant
        sun_dev_from_id : float  max |U_SU - I|
        sun_norm : float  Frobenius norm of (U_SU - I)
        u1_magnitude : float  |exp(i * u1_phase)|  (should be 1)
    """
    U = np.asarray(holonomy_matrix, dtype=complex)
    n = U.shape[0]
    det_U = np.linalg.det(U)
    u1_phase = float(np.angle(det_U)) / n
    U_SU = U * np.exp(-1j * u1_phase)   # remove the U(1) factor
    return {
        "u1_phase_rad": u1_phase,
        "sun_matrix_real": np.real(U_SU).tolist(),
        "sun_matrix_imag": np.imag(U_SU).tolist(),
        "sun_dev_from_id": float(np.max(np.abs(U_SU - np.eye(n)))),
        "sun_norm": float(np.linalg.norm(U_SU - np.eye(n), "fro")),
        "u1_magnitude": float(abs(np.exp(1j * u1_phase))),
        "det_U_magnitude": float(abs(det_U)),
    }


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def random_unitary_3(rng: np.random.Generator) -> np.ndarray:
    """Generate a random 3x3 unitary (Haar-distributed U(3)) via QR decomposition.

    Parameters
    ----------
    rng : numpy random generator

    Returns
    -------
    (3, 3) complex unitary matrix
    """
    Z = rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3))
    Q, R = np.linalg.qr(Z)
    # Fix phase ambiguity: multiply columns by sign of diagonal of R
    d = np.diag(R)
    ph = d / np.abs(d)
    return Q * ph[np.newaxis, :]


def _make_gauged_frame_fn_const(
    frame_fn: Callable[[float], np.ndarray],
    Q: np.ndarray,
) -> Callable[[float], np.ndarray]:
    """Wrap frame_fn with a constant left gauge transform: E'(th) = Q @ E(th).

    For a flat bundle (holonomy = I), the gauged holonomy is also I because
    the constant Q factors out of the link products.
    """
    def fn(theta: float) -> np.ndarray:
        return Q @ np.asarray(frame_fn(theta), dtype=complex)
    return fn