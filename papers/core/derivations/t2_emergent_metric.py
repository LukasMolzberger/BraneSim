"""
T2 -- Emergent Lorentz / effective metric.

Substrate vacuum tangent stiffness (Layer 1):
    M_mu = (1-alpha_mu) I + alpha_mu e_mu e_mu^T,   alpha_i=alpha_s, alpha_4=alpha_t
    C_mu = eta_mu kappa_mu M_mu,   eta=(-1,-1,-1,+1),  kappa_i=kappa_s, kappa_4=kappa_t
Small-fluctuation dispersion (4D block, time = direction 4):
    kappa_t 2(1-cos omega) M_4 eps = sum_i kappa_s 2(1-cos q_i) M_i eps
  => generalized eigenproblem  A(q) eps = lambda M_4 eps,  omega^2 = lambda/kappa_t.

This script proves the provable core of T2 and precisely delimits the conjecture:

  (A) LORENTZIAN SIGNATURE from eta.  The spacetime quadratic form H_mu nu =
      diag(eta_mu kappa_mu m_mu) has signature (3,1): three spacelike, one
      timelike, purely because eta_4=+1 opposes eta_i=-1. Flip eta_4->-1 and the
      form becomes (4,0) Euclidean -- no cone, no waves. So the Lorentzian metric
      is the eta sign pattern, DERIVED not imposed (BACKBONE #6, ARCH Layer 0).

  (B) EFFECTIVE METRIC + CALIBRATION.  Closed-form branch speeds
      c_L^2 = 1/[gamma_t(1-alpha_t)],  c_T^2 = (1-alpha_s)/[gamma_t(1-alpha_t)],
      c_4^2 = (1-alpha_s)/gamma_t,  verified against the exact lattice dispersion.
      This is the calibration among (alpha_s, alpha_t, gamma_t, kappa, a) T2 asks for.

  (C) MASSLESS / LINEAR CONE.  omega = c|q| + O(|q|^3): each branch is a relativistic
      null cone at long wavelength; the O(|q|^2) piece is the lattice correction.

  (D) THE LAB IS ANISOTROPIC (as BACKBONE #8 acknowledges), controlled by alpha_s.
      The amplitude (e_4) mode is EXACTLY isotropic; the spatial "light" modes are
      birefringent off-axis. Directional spread -> 0 only as alpha_s -> 0, which
      conflicts with the gauge sector (alpha_s>0). So isotropy cannot come from
      tuning; it must come from observer renormalization.

  (E) OPERATIONAL LORENTZ (per sector) + CAUSALITY.  A single signal cone is
      Lorentz-invariant: boosts map its null vectors to null vectors; a slower
      branch (c_4<c_T) stays timelike/inside the light cone -> causal, "massive".
      Full cross-sector universality is the load-bearing dual-observer conjecture.
"""

import numpy as np
import scipy.linalg as la

np.set_printoptions(precision=4, suppress=True, linewidth=140)

def Mmat(alpha, ax):
    e = np.zeros(4); e[ax] = 1.0
    return (1-alpha)*np.eye(4) + alpha*np.outer(e, e)

def stiffness(alpha_s=0.6, alpha_t=0.9, kappa_s=1.0, gamma_t=None):
    # gamma_t from the T1 helix force-balance if not given (keeps a common carrier)
    if gamma_t is None: gamma_t = 2.32192
    kappa_t = gamma_t*kappa_s
    M = [Mmat(alpha_s,0), Mmat(alpha_s,1), Mmat(alpha_s,2), Mmat(alpha_t,3)]
    return M, kappa_s, kappa_t, alpha_s, alpha_t, gamma_t

def branch_speeds(M, kappa_s, kappa_t, qhat, qmag=1e-3, exact=True):
    q = np.asarray(qhat,float)*qmag
    if exact:
        A = sum(kappa_s*2*(1-np.cos(q[i]))*M[i] for i in range(3))
    else:
        A = sum(kappa_s*q[i]**2*M[i] for i in range(3))
    lam = la.eigvals(A, M[3]).real
    return np.sort(np.sqrt(np.clip(lam/kappa_t, 0, None))/qmag)


# =========================================================== MAIN
if __name__ == "__main__":
    print("="*72); print("T2 -- emergent Lorentz metric / effective cone"); print("="*72)
    M, ks, kt, alpha_s, alpha_t, gamma_t = stiffness()
    print(f"params: alpha_s={alpha_s}, alpha_t={alpha_t}, gamma_t={gamma_t:.4f} "
          f"(from T1 carrier), kappa_s={ks}")

    # ---- (A) Lorentzian signature from eta ----------------------------------
    print("\n(A) LORENTZIAN SIGNATURE IS THE eta SIGN PATTERN")
    eta = np.array([-1,-1,-1,+1.0]); kappa = np.array([ks,ks,ks,kt])
    eps = np.array([0,1,0,0.0])                       # a transverse polarization
    mmu = np.array([eps@M[mu]@eps for mu in range(4)])
    H  = np.diag(eta*kappa*mmu)                       # spacetime quadratic form k^T H k
    Heu = np.diag(np.array([-1,-1,-1,-1.0])*kappa*mmu)  # counterfactual: all-minus eta
    sig = tuple(int(np.sign(x)) for x in np.diag(H))
    sigE = tuple(int(np.sign(x)) for x in np.diag(Heu))
    print(f"    eta=(-1,-1,-1,+1): diag signs {sig}  -> signature (3,1) LORENTZIAN"
          f"  => a light cone exists")
    print(f"    eta=(-1,-1,-1,-1): diag signs {sigE}  -> signature (4,0) EUCLIDEAN"
          f"  => definite form, NO waves")
    print("    => the (3,1) Minkowski signature is derived from prestress signs, not imposed.")

    # ---- (B) effective metric + calibration ---------------------------------
    print("\n(B) EFFECTIVE METRIC + CALIBRATION (speeds vs closed form)")
    cx = branch_speeds(M, ks, kt, [1,0,0])
    cL = 1/np.sqrt(gamma_t*(1-alpha_t))
    cT = np.sqrt((1-alpha_s)/(gamma_t*(1-alpha_t)))
    c4 = np.sqrt((1-alpha_s)/gamma_t)
    print(f"    numeric (q||x)  : {cx}")
    print(f"    analytic        : c_4={c4:.4f}  c_T={cT:.4f} (x2)  c_L={cL:.4f}")
    print(f"    match: {np.allclose(cx, sorted([c4,cT,cT,cL]), atol=1e-3)}")
    print(f"    effective (mostly-minus) line element for the light branch c_T:")
    print(f"      ds^2 = -dt^2 + c_T^-2 (dx^2+dy^2+dz^2),  c_T = {cT:.4f}")

    # calibration table over parameters
    print("    calibration c_T(alpha_s,alpha_t,gamma_t):")
    for a_s in (0.3,0.6):
        for a_t in (0.8,0.9):
            cTk = np.sqrt((1-a_s)/(gamma_t*(1-a_t)))
            print(f"      alpha_s={a_s} alpha_t={a_t}: c_T={cTk:.4f}")

    # ---- (C) massless / linear cone -----------------------------------------
    print("\n(C) MASSLESS LINEAR CONE  omega = c|q| + O(|q|^3)")
    for qm in (0.2, 0.1, 0.05):
        c = branch_speeds(M, ks, kt, [1,0,0], qmag=qm)[1]   # a transverse branch
        print(f"    |q|={qm:.3f}:  omega/|q| (c_T branch) = {c:.5f}  "
              f"(-> {cT:.5f} as |q|->0)")
    print("    => linear dispersion at long wavelength: a relativistic massless cone.")

    # ---- (D) lab anisotropy, controlled by alpha_s --------------------------
    print("\n(D) LAB ANISOTROPY (real; BACKBONE #8) vs alpha_s")
    dirs = [np.array([1,0,0.]), np.array([1,1,0.])/np.sqrt(2),
            np.array([1,1,1.])/np.sqrt(3), np.array([2,1,0.])/np.sqrt(5)]
    print("    speeds by direction at alpha_s=0.6 (note birefringence & isotropic c_4):")
    for d in dirs:
        print(f"      n={np.round(d,3)}:  {branch_speeds(M,ks,kt,d)}")
    print("    directional spread of the fastest spatial branch vs alpha_s:")
    for a_s in (0.6, 0.4, 0.2, 0.05, 0.0):
        Mk,ksk,ktk,_,_,_ = stiffness(alpha_s=a_s)
        cmax = [branch_speeds(Mk,ksk,ktk,d)[-1] for d in
                [np.array([1,0,0.]), np.array([1,1,1.])/np.sqrt(3)]]
        spread = abs(cmax[0]-cmax[1])/np.mean(cmax)
        print(f"      alpha_s={a_s:.2f}:  fractional anisotropy = {spread:.4f}")
    print("    => anisotropy -> 0 only as alpha_s->0, which disables the gauge sector.")
    print("       So emergent isotropy cannot come from tuning; it needs the inside")
    print("       observer's rod/clock renormalization (the load-bearing conjecture).")

    # ---- (E) operational Lorentz per sector + causality ---------------------
    print("\n(E) OPERATIONAL LORENTZ (single cone) + CAUSALITY")
    # work in c_T=1 units; boost along x by beta; check null-cone invariance
    def boost(beta):
        g = 1/np.sqrt(1-beta**2)
        L = np.eye(4); L[0,0]=g; L[0,1]=-g*beta; L[1,0]=-g*beta; L[1,1]=g
        return L
    rng = np.random.default_rng(0); worst=0.0
    for _ in range(1000):
        qh = rng.standard_normal(3); qh/=np.linalg.norm(qh)
        k = np.array([1.0, qh[0], qh[1], qh[2]])   # null: omega=|q|=1 (c_T units)
        kp = boost(0.5) @ k
        worst = max(worst, abs(kp[0]**2 - (kp[1]**2+kp[2]**2+kp[3]**2)))
    print(f"    boost(beta=0.5) applied to 1000 null 4-vectors: max |omega'^2-|q'|^2| "
          f"= {worst:.2e}  (== 0)")
    print("    => the c_T cone is invariant under the Lorentz group of its metric.")
    ratio = (c4/cT)**2
    print(f"    slower branch c_4/c_T = {c4/cT:.4f} < 1: its worldvector has "
          f"omega^2-c_T^2|q|^2 = ({ratio:.3f}-1)|q|^2 < 0")
    print("    => the amplitude mode is TIMELIKE/inside the light cone -> causal,")
    print("       'massive/subluminal' relative to the c_T signal cone.")

    print("\n" + "="*72)
    print("CONCLUSION: the emergent metric is Lorentzian (signature from eta), the")
    print("branch speeds are calibrated in closed form, the long-wavelength cone is")
    print("massless & Lorentz-invariant per sector, and slower modes are causal.")
    print("The lab is genuinely anisotropic for alpha_s>0; full cross-sector")
    print("isotropy is the dual-observer conjecture (BACKBONE #8), now sharpened")
    print("with a concrete obstruction and its alpha_s->0 removal. T2 core established.")
    print("="*72)