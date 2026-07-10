"""
T2 test -- is the emergent GAUGE-mode cone more isotropic than the bare phonons?

T2 found the bare acoustic branches are anisotropic (fastest-branch directional spread
~0.25 at alpha_s=0.6). But "light" is the gauge/Berry mode on the nonlinear carrier,
not a bare phonon. The gauge kinetic tensor is the substrate quantum metric
    g_ij(k) = 1/2 Re Tr[ (d_i P_3)(d_j P_3) ]   (spatial i,j),
whose spatial anisotropy sets the emergent photon/gluon cone. We measure its
eigenvalue spread and compare to the bare acoustic anisotropy.

A single helix picks its propagation axis, so its g_ij is trivially axis-anisotropic;
the physical (orientation-averaged) vacuum is modelled by averaging g_ij over the three
propagation-axis choices. We report both the single-orientation and the
orientation-averaged gauge anisotropy.
"""
import numpy as np
import importlib.util, os

_here = os.path.dirname(__file__)
def _load(n):
    s = importlib.util.spec_from_file_location(n, os.path.join(_here, f"{n}.py"))
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
t1 = _load("t1_su3_witness"); gen = _load("t1_genericity_robustness")
np.set_printoptions(precision=4, suppress=True, linewidth=140)

KPTS = [np.array(k) for k in
        [(0.7,1.1,1.9,0.5),(1.3,0.6,2.4,1.0),(2.0,2.2,0.5,1.7),(0.4,2.5,1.2,2.9)]]

def projector(model, k, TRIP):
    w, V = np.linalg.eigh(t1.bloch_D(model, k)); U = V[:, TRIP:TRIP+3]
    return U @ U.conj().T

def local_triplet(model, k):
    w, _ = np.linalg.eigh(t1.bloch_D(model, k)); nb = len(w); best = (None, -1)
    for j in range(1, nb-3):
        trip = w[j:j+3]; gap = min(trip.min()-w[j-1], w[j+3]-trip.max())
        if gap > 0 and gap/max(trip.max()-trip.min(), 1e-9) > best[1]:
            best = (j, gap/max(trip.max()-trip.min(), 1e-9), gap)
    return best

def quantum_metric_spatial(model, kpts=KPTS, h=1e-3, gap_min=0.2):
    """spatial 3x3 quantum-metric tensor g_ij, averaged over isolated-carrier k."""
    G = np.zeros((3, 3)); n = 0
    for k in kpts:
        TRIP, ratio, gap = local_triplet(model, k)
        if gap < gap_min: continue
        dP = []
        for i in range(3):
            dk = np.zeros(4); dk[i] = h
            dP.append((projector(model, k+dk, TRIP) - projector(model, k-dk, TRIP))/(2*h))
        for i in range(3):
            for j in range(3):
                G[i, j] += 0.5*np.real(np.trace(dP[i] @ dP[j]))
        n += 1
    return G/max(n, 1), n

def anisotropy(M):
    ev = np.linalg.eigvalsh(0.5*(M+M.T)); ev = np.abs(ev)
    return (ev.max()-ev.min())/max(ev.mean(), 1e-15), ev


if __name__ == "__main__":
    print("="*74); print("T2 -- gauge-mode cone isotropy vs bare phonons"); print("="*74)

    print("\nbare acoustic anisotropy (T2, alpha_s=0.6): fastest-branch spread ~ 0.254")

    print("\ngauge kinetic tensor g_ij (quantum metric of the carrier):")
    Gsum = np.zeros((3, 3)); tot = 0
    for prop in (0, 1, 2):
        m = gen.build_helix_general(prop_axis=prop, alpha_s=0.6, alpha_t=0.9)
        G, n = quantum_metric_spatial(m)
        an, ev = anisotropy(G)
        Gsum += G; tot += 1
        print(f"    prop_axis={prop}:  eigenvalues {np.round(ev,4)}   anisotropy = {an:.3f}"
              f"   ({n} k-pts)")
    Gavg = Gsum/tot
    an_avg, ev_avg = anisotropy(Gavg)
    print(f"\n  orientation-averaged gauge tensor eigenvalues: {np.round(ev_avg,4)}")
    print(f"  orientation-averaged gauge-cone anisotropy   = {an_avg:.3f}")
    print(f"  bare acoustic anisotropy                     = 0.254")
    if an_avg < 0.254:
        print("  => the gauge-mode cone is MORE isotropic than the bare phonons:")
        print("     the emergent 'light' is closer to Lorentz-invariant than sound.")
    else:
        print("  => the gauge-mode cone is NOT more isotropic than the bare phonons;")
        print("     the T2 anisotropy obstruction persists in the gauge sector.")

    print("\n" + "="*74)
    print("HONEST READ: single-helix carriers pick a propagation axis (large apparent")
    print("anisotropy); the orientation-averaged gauge cone is the fair measure. The")
    print("number above is the result -- it either softens or confirms the T2")
    print("dual-observer obstruction; either way it is now quantified, not assumed.")
    print("="*74)
