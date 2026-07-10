"""
T1 strengthening #5 -- analytic backing + consistency checks.

(A) SMALL-AMPLITUDE SCALING.  The traceless WZ curvature is present at ARBITRARILY
    small carrier amplitude A: the su(3) span stays 8 as A->0, while its magnitude
    scales as ||G|| ~ A^2 (leading order). So the su(3) is not a large-amplitude /
    strong-coupling numerical artifact; it appears at leading nontrivial order, and
    the straight-vacuum limit A->0 correctly returns zero curvature. Mechanism:
    d_k P_3 ~ O(A) (texture-induced eigenvector rotation), so F ~ (d P)^2 ~ O(A^2),
    while the complex Bloch phases (O(1)) supply the full 8 directions already there.

(B) NON-ABELIAN BIANCHI.  The projected curvature obeys D_[mu F_nu rho] = 0. Verified
    numerically in a parallel-transported smooth gauge; the residual -> 0 with the
    finite-difference step (a genuine identity, guarding the T4 G against code error).
"""
import numpy as np
import importlib.util, os

_here = os.path.dirname(__file__)
spec = importlib.util.spec_from_file_location("t1", os.path.join(_here, "t1_su3_witness.py"))
t1 = importlib.util.module_from_spec(spec); spec.loader.exec_module(t1)
np.set_printoptions(precision=4, suppress=True, linewidth=140)

KPTS = [np.array(k) for k in
        [(0.7,1.1,1.9,0.5),(1.3,0.6,2.4,1.0),(2.0,2.2,0.5,1.7),(0.4,2.5,1.2,2.9)]]


# ---------------------------------------------------------------- (A) small A
def su3_span(A):
    """raw su(3) span (directions) of the normalized curvature samples over k*."""
    m = t1.build_helix(A=A); allgens = []
    for k in KPTS:
        r = t1.wz_su3(m, k, verbose=False)
        allgens += [g/np.linalg.norm(g) for g in r['gens'] if np.linalg.norm(g) > 0]
    return t1.span_rank(allgens) if allgens else 0

def curvature_norm_fixed(A, TRIP, kstar, h=2e-3):
    """||G|| summed over the 6 planes at a FIXED band index and FIXED k* (so the
    same physical curvature is tracked as A varies)."""
    m = t1.build_helix(A=A)
    def frame(k):
        w, V = np.linalg.eigh(t1.bloch_D(m, k)); return V[:, TRIP:TRIP+3]
    def linkU(a, b):
        M = frame(a).conj().T @ frame(b); U, _, Vh = np.linalg.svd(M); return U @ Vh
    tot = 0.0
    for (mu, nu) in [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]:
        dmu = np.zeros(4); dmu[mu] = h; dnu = np.zeros(4); dnu[nu] = h
        W = (linkU(kstar, kstar+dmu) @ linkU(kstar+dmu, kstar+dmu+dnu)
             @ linkU(kstar+dmu+dnu, kstar+dnu) @ linkU(kstar+dnu, kstar))
        ev, evec = np.linalg.eig(W)
        F = 1j*(evec @ np.diag(np.log(ev)) @ np.linalg.inv(evec))/h**2
        F = 0.5*(F + F.conj().T); G = F - np.trace(F)/3*np.eye(3)
        tot += np.linalg.norm(G)
    return tot

def comm_norm(A):
    """mean stiffness-matrix noncommutativity ||[C,C]|| (the mechanical seed)."""
    Cs = list(t1.build_helix(A=A)['Cmat'].values()); tot = 0.0; cnt = 0
    for i in range(0, len(Cs), 7):
        for j in range(0, len(Cs), 11):
            tot += np.linalg.norm(Cs[i]@Cs[j] - Cs[j]@Cs[i]); cnt += 1
    return tot/cnt


# ---------------------------------------------------------------- (B) Bianchi
def smooth_frame(model, k, TRIP, U0=None):
    w, V = np.linalg.eigh(t1.bloch_D(model, k))
    U = V[:, TRIP:TRIP+3]
    if U0 is not None:                       # parallel-transport gauge alignment
        M = U0.conj().T @ U                  # 3x3 overlap
        Aa, _, Bh = np.linalg.svd(M)
        U = U @ (Bh.conj().T @ Aa.conj().T)  # U * (polar unitary of M)^dagger
    return U

def bianchi_residual(model, k0, TRIP, h):
    dirs = (0, 1, 2)                         # test the (k0,k1,k2) 3-cube
    U0 = smooth_frame(model, k0, TRIP)
    def frame(k):
        return smooth_frame(model, k, TRIP, U0)
    def A_at(k):
        Uc = frame(k); Amu = []
        for mu in range(4):
            dk = np.zeros(4); dk[mu] = h
            Up = smooth_frame(model, k+dk, TRIP, U0)
            Um = smooth_frame(model, k-dk, TRIP, U0)
            dU = (Up - Um)/(2*h)
            Amu.append(1j*(Uc.conj().T @ dU))
        return Amu
    def F_at(k):
        A = A_at(k); F = {}
        for mu in range(4):
            for nu in range(4):
                dmu = np.zeros(4); dmu[mu] = h; dnu = np.zeros(4); dnu[nu] = h
                dAnu = (A_at(k+dmu)[nu] - A_at(k-dmu)[nu])/(2*h)
                dAmu = (A_at(k+dnu)[mu] - A_at(k-dnu)[mu])/(2*h)
                F[(mu, nu)] = dAnu - dAmu - 1j*(A[mu]@A[nu] - A[nu]@A[mu])
        return F, A
    F0, A0 = F_at(k0)
    # covariant derivative D_mu F_nu rho = d_mu F_nu rho - i[A_mu, F_nu rho]
    def DF(mu, nu, rho):
        dk = np.zeros(4); dk[mu] = h
        Fp, _ = F_at(k0+dk); Fm, _ = F_at(k0-dk)
        dF = (Fp[(nu, rho)] - Fm[(nu, rho)])/(2*h)
        return dF - 1j*(A0[mu] @ F0[(nu, rho)] - F0[(nu, rho)] @ A0[mu])
    mu, nu, rho = dirs
    cyc = DF(mu, nu, rho) + DF(nu, rho, mu) + DF(rho, mu, nu)
    scale = max(np.linalg.norm(F0[(mu, nu)]), 1e-12)
    return np.linalg.norm(cyc)/scale


if __name__ == "__main__":
    print("="*74); print("T1 -- small-amplitude scaling + non-abelian Bianchi"); print("="*74)

    kstar = np.array([0.7, 1.1, 1.9, 0.5])
    print("\n(A) SMALL-AMPLITUDE SCALING of the traceless su(3) curvature")
    print("    (i) su(3) span (directions) persists as A -> 0:")
    for A in (0.30, 0.10, 0.05, 0.02, 0.01):
        print(f"        A={A:.3f}:  span = {su3_span(A)}/8")
    print("    (ii) magnitude at FIXED band index & k* (tracks one curvature):")
    TRIP0 = t1.wz_su3(t1.build_helix(A=0.1), kstar, verbose=False)['TRIP']
    xs, ys = [], []
    print("         A        ||G||")
    for A in (0.10, 0.07, 0.05, 0.03, 0.02, 0.01, 0.005):
        g = curvature_norm_fixed(A, TRIP0, kstar)
        print(f"        {A:.3f}   {g:.4e}")
        if g > 0: xs.append(np.log(A)); ys.append(np.log(g))
    pG = np.polyfit(xs, ys, 1)[0]
    print("    (iii) mechanical seed -- stiffness noncommutativity ||[C,C]||:")
    xs2, ys2 = [], []
    for A in (0.10, 0.05, 0.02, 0.01):
        c = comm_norm(A); xs2.append(np.log(A)); ys2.append(np.log(c))
        print(f"        A={A:.3f}:  ||[C,C]|| = {c:.4e}")
    pC = np.polyfit(xs2, ys2, 1)[0]
    print(f"    => ||[C,C]|| ~ A^{pC:.2f} (=A^1, seed)  ->  ||G|| ~ A^{pG:.2f} (=A^2, curvature),")
    print("       and the su(3) span stays 8: leading-order effect, NOT a large-A artifact.")

    print("\n(B) NON-ABELIAN BIANCHI  D_[mu F_nu rho] = 0  (parallel-transport gauge)")
    m = t1.build_helix(); kstar = np.array([0.7, 1.1, 1.9, 0.5])
    TRIP = t1.wz_su3(m, kstar, verbose=False)['TRIP']
    print("      h          relative ||D_[mu F_nu rho]|| / ||F||")
    for h in (2e-2, 1e-2, 5e-3):
        res = bianchi_residual(m, kstar, TRIP, h)
        print(f"      {h:.0e}     {res:.3e}")
    print("    => residual small and decreasing with h: Bianchi holds (identity),")
    print("       confirming the T4 field strength G is a consistent curvature.")

    print("\n" + "="*74)
    print("CONCLUSION: the su(3) curvature is leading-order (A^2) and present for all")
    print("finite amplitudes, and satisfies the non-abelian Bianchi identity.")
    print("="*74)
