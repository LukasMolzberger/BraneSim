"""
T1 -- SU(3) Wilczek-Zee sector: the decisive calculation (OPEN_TASKS §5, core D5).

Substrate (Layer 0):
    S = 1/2 sum_{n,mu} eta_mu kappa_mu (L_{n mu} - r_mu)^2
    L_{n mu} = |R_{n+e_mu} - R_n|            (Pythagorean -- the only nonlinearity)
    eta = (-1,-1,-1,+1),  r_i = alpha_s a, r_4 = alpha_t a,
    kappa_i = kappa_s,    kappa_4 = kappa_t = gamma_t kappa_s.

Chain executed here (no gauge field is ever inserted by hand):
    R_bar  (exact stationary, finite-amplitude, periodic vacuum carrier)
      -> anisotropic stiffness tensor C_{n mu}[R_bar]        (Layer 1)
      -> Bloch fluctuation operator D_Rbar(k)                (Layer 2)
      -> isolated rank-3 carrier projector P_3(k)            (Layer 2/3)
      -> gauge-invariant Wilczek-Zee curvature F_{mu nu}(k)  (Layer 3)
      -> trace/traceless split -> G_{mu nu} in su(3)         (Layer 3/4)
      -> Lie-closure rank test: 8 == genuine su(3)           (settles D5)

Background: a closed-form circularly-polarized helix -- propagation in the (1,4)
plane, polarization in the (2,3) plane -- so every link length is n-independent
and the transverse force balance collapses to ONE scalar equation, solved for
gamma_t. The result is an EXACT critical point of S (verified ||grad S|| ~ 1e-16),
i.e. the universal nonlinear periodic vacuum carrier of Decision M.
"""

import numpy as np

np.set_printoptions(precision=4, suppress=True, linewidth=140)

# fixed Gell-Mann-adapted coordinates of an anti-Hermitian traceless 3x3 (dim 8)
def su3_coords(X):
    return np.array([X[0,1].real, X[0,1].imag, X[0,2].real, X[0,2].imag,
                     X[1,2].real, X[1,2].imag, X[0,0].imag, X[1,1].imag], float)

def span_rank(mats, tol=1e-7):
    if not mats: return 0
    M = np.array([su3_coords(X) for X in mats])
    s = np.linalg.svd(M, compute_uv=False)
    return int((s > tol*s.max()).sum())

def lie_closure(gen, max_iter=12):
    basis = list(gen); rank = span_rank(basis)
    for _ in range(max_iter):
        new = [basis[i]@basis[j]-basis[j]@basis[i]
               for i in range(len(basis)) for j in range(i+1, len(basis))]
        r2 = span_rank(basis+new)
        if r2 == rank: break
        rank, basis = r2, basis+new
    return rank

# so(3) content = span of the 3 antisymmetric directions only (imag off-diagonal)
def so3_coords(X):
    return np.array([X[0,1].imag, X[0,2].imag, X[1,2].imag], float)


def build_helix(a=1.0, alpha_s=0.6, alpha_t=0.9, kappa_s=1.0, A=0.30,
                N1=3, N4=3, w1=1, w4=1,
                pol=(np.array([0,1,0,0.]), np.array([0,0,1,0.]))):
    """construct the exact stationary helical carrier; return everything needed."""
    K1, K4 = 2*np.pi*w1/N1, 2*np.pi*w4/N4
    Kvec = np.array([K1, 0.0, 0.0, K4])
    eta  = np.array([-1.0, -1.0, -1.0, +1.0])
    p, q = pol
    def Rbar(n):
        n = np.asarray(n, float); th = Kvec @ n
        return a*n + A*(np.cos(th)*p + np.sin(th)*q)
    L1 = np.sqrt(a*a + 2*A*A*(1-np.cos(K1)))
    L4 = np.sqrt(a*a + 2*A*A*(1-np.cos(K4)))
    r_s, r_t = alpha_s*a, alpha_t*a
    # transverse force balance -> gamma_t
    lhs = kappa_s*(L1-r_s)/L1*(1-np.cos(K1))
    rhs = (L4-r_t)/L4*(1-np.cos(K4))
    kappa_t = lhs/rhs; gamma_t = kappa_t/kappa_s
    kappa = np.array([kappa_s, kappa_s, kappa_s, kappa_t])
    r     = np.array([r_s, r_s, r_s, r_t])
    N = np.array([N1,1,1,N4])
    sites = [np.array([i,0,0,l]) for i in range(N1) for l in range(N4)]
    site_index = {tuple(t):i for i,t in enumerate(sites)}
    e = np.eye(4, dtype=int)
    def Cmatrix(n, mu):
        Q = Rbar(n+e[mu]) - Rbar(n); L = np.linalg.norm(Q); Qh = Q/L
        return eta[mu]*kappa[mu]*((1-r[mu]/L)*np.eye(4) + (r[mu]/L)*np.outer(Qh,Qh))
    Cmat = {(tuple(t),mu):Cmatrix(t,mu) for t in sites for mu in range(4)}
    # stationarity residual
    F = {tuple(t):np.zeros(4) for t in sites}
    for t in sites:
        for mu in range(4):
            Q = Rbar(t+e[mu])-Rbar(t); L=np.linalg.norm(Q)
            force = eta[mu]*kappa[mu]*(L-r[mu])*Q/L
            F[tuple(t)] += force; F[tuple((t+e[mu])%N)] -= force
    resid = max(np.linalg.norm(v) for v in F.values())
    return dict(N=N, sites=sites, site_index=site_index, e=e, Cmat=Cmat,
                gamma_t=gamma_t, L1=L1, L4=L4, K1=K1, K4=K4, resid=resid, Ns=len(sites))


def bloch_D(model, kvec):
    N, sites, si, e, Cmat, Ns = (model['N'], model['sites'], model['site_index'],
                                 model['e'], model['Cmat'], model['Ns'])
    H = np.zeros((4*Ns, 4*Ns), complex)
    for t in sites:
        i = si[tuple(t)]
        for mu in range(4):
            C = Cmat[(tuple(t), mu)]
            p_int = t+e[mu]; s = p_int % N; delta = p_int - s
            j = si[tuple(s)]; phase = np.exp(1j*(kvec@delta))
            bi, bj = 4*i, 4*j
            H[bi:bi+4, bi:bi+4] += C; H[bj:bj+4, bj:bj+4] += C
            H[bi:bi+4, bj:bj+4] += -C*phase; H[bj:bj+4, bi:bi+4] += -C*np.conj(phase)
    return 0.5*(H+H.conj().T)


def wz_su3(model, kstar, h=2e-3, ball=0.05, nbase=6, seed=0, verbose=True):
    """isolate a rank-3 carrier at kstar, compute WZ curvature, return ranks."""
    def eigs(k):
        w,V = np.linalg.eigh(bloch_D(model,k)); return w,V
    w0,_ = eigs(kstar); nb = len(w0)
    # best locally-isolated interior triplet
    best=(None,-1)
    for j in range(1, nb-3):
        trip=w0[j:j+3]; eps=trip.max()-trip.min()
        gap=min(trip.min()-w0[j-1], w0[j+3]-trip.max())
        if gap>0 and gap/max(eps,1e-9)>best[1]: best=(j, gap/max(eps,1e-9))
    TRIP=best[0]
    def frame(k):
        w,V=eigs(k); return V[:,TRIP:TRIP+3]
    def gapok(k):
        w,_=eigs(k); return w[TRIP]-w[TRIP-1]>1e-3 and w[TRIP+3]-w[TRIP+2]>1e-3
    def linkU(a_,b_):
        M=frame(a_).conj().T@frame(b_); U,_,Vh=np.linalg.svd(M); return U@Vh
    def wz(kb,mu,nu):
        dmu=np.zeros(4);dmu[mu]=h; dnu=np.zeros(4);dnu[nu]=h; k=np.array(kb,float)
        if not all(gapok(c) for c in (k,k+dmu,k+dmu+dnu,k+dnu)): return None
        W=linkU(k,k+dmu)@linkU(k+dmu,k+dmu+dnu)@linkU(k+dmu+dnu,k+dnu)@linkU(k+dnu,k)
        ev,evec=np.linalg.eig(W); logW=evec@np.diag(np.log(ev))@np.linalg.inv(evec)
        F=1j*logW/(h*h); return 0.5*(F+F.conj().T)
    planes=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
    rng=np.random.default_rng(seed)
    bases=[np.array(kstar,float)]+[np.array(kstar,float)+ball*rng.standard_normal(4)
                                   for _ in range(nbase-1)]
    G, trmag, tlmag = [], [], []
    for kb in bases:
        for (mu,nu) in planes:
            F=wz(kb,mu,nu)
            if F is None: continue
            tr=np.trace(F)/3.0; Gm=F-tr*np.eye(3)
            trmag.append(abs(tr)); tlmag.append(np.linalg.norm(Gm))
            if np.linalg.norm(Gm)>1e-6: G.append(1j*Gm)
    # decompose one representative into so(3) / symmetric-traceless(+Cartan)
    so3_dim = span_rank([g for g in G]) if not G else None
    info = dict(TRIP=TRIP, ratio=best[1], nsamp=len(G), gens=G,
                trace=np.mean(trmag) if trmag else 0.0,
                traceless=np.mean(tlmag) if tlmag else 0.0,
                raw_rank=span_rank(G), lie_rank=lie_closure(G),
                so3_only=span_rank([1j*(g/1j) for g in []]))  # placeholder
    # how much of the span is the antisymmetric so(3) subspace vs the full 8
    if G:
        Mso = np.array([so3_coords(g) for g in G])
        info['so3_span'] = int((np.linalg.svd(Mso,compute_uv=False)>1e-7).sum())
        # symmetric-traceless real off-diagonal presence (lambda1,4,6) + Cartan (lambda3,8)
        sym = np.array([[g[0,1].real,g[0,2].real,g[1,2].real,g[0,0].imag,g[1,1].imag]
                        for g in G])
        info['sym_cartan_span'] = int((np.linalg.svd(sym,compute_uv=False)>1e-7).sum())
    return info


def control_real_frame(model, kstar, h=2e-3, ball=0.05, nbase=6, seed=0):
    """FAIR so(3) control: a genuine real transverse frame.
    Use a single site's four (real symmetric, non-commuting) stiffness tensors as
    a 4-band operator M(k)=sum_mu 2(1-cos k_mu) C_mu. Real symmetric at every k =>
    real eigenvectors => WZ connection is so(3)-valued (D5's 'real frame' case)."""
    t0 = model['sites'][len(model['sites'])//2]
    Cs = [model['Cmat'][(tuple(t0),mu)] for mu in range(4)]
    def M(k): return sum(2*(1-np.cos(k[mu]))*Cs[mu] for mu in range(4))
    def frame(k):
        w,V=np.linalg.eigh(M(k)); return V[:,0:3]        # 3 of 4 bands
    def gapok(k):
        w,_=np.linalg.eigh(M(k)); return w[3]-w[2]>1e-4
    def linkU(a_,b_):
        Mo=frame(a_).conj().T@frame(b_); U,_,Vh=np.linalg.svd(Mo); return U@Vh
    def wz(kb,mu,nu):
        dmu=np.zeros(4);dmu[mu]=h; dnu=np.zeros(4);dnu[nu]=h; k=np.array(kb,float)
        if not all(gapok(c) for c in (k,k+dmu,k+dmu+dnu,k+dnu)): return None
        W=linkU(k,k+dmu)@linkU(k+dmu,k+dmu+dnu)@linkU(k+dmu+dnu,k+dnu)@linkU(k+dnu,k)
        ev,evec=np.linalg.eig(W); F=1j*(evec@np.diag(np.log(ev))@np.linalg.inv(evec))/h**2
        return 0.5*(F+F.conj().T)
    planes=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
    rng=np.random.default_rng(seed)
    bases=[np.array(kstar,float)]+[np.array(kstar,float)+ball*rng.standard_normal(4)
                                   for _ in range(nbase-1)]
    G=[]
    for kb in bases:
        for (mu,nu) in planes:
            F=wz(kb,mu,nu)
            if F is None: continue
            Gm=F-np.trace(F)/3*np.eye(3)
            if np.linalg.norm(Gm)>1e-6: G.append(1j*Gm)
    return dict(nsamp=len(G), raw_rank=span_rank(G), lie_rank=lie_closure(G))


# ============================================================ MAIN
if __name__ == "__main__":
    print("="*72)
    print("PRIMARY WITNESS")
    print("="*72)
    m = build_helix()
    print(f"gamma_t (from force balance) = {m['gamma_t']:.5f}   "
          f"L1={m['L1']:.5f} L4={m['L4']:.5f}")
    print(f"exact-stationarity residual  ||grad S|| = {m['resid']:.2e}")
    kstar = np.array([0.7, 1.1, 1.9, 0.5])
    r = wz_su3(m, kstar)
    print(f"\ncarrier = bands [{r['TRIP']},{r['TRIP']+1},{r['TRIP']+2}]  "
          f"local gap/spread = {r['ratio']:.2f}  ({r['nsamp']} curvature samples)")
    print(f"U(1) trace curvature |tr F|      ~ {r['trace']:.3e}")
    print(f"SU(3) traceless curvature ||G||  ~ {r['traceless']:.3e}")
    print(f"so(3) (antisymmetric) span       = {r['so3_span']}  / 3")
    print(f"symmetric+Cartan span            = {r['sym_cartan_span']}  / 5")
    print(f"raw curvature span               = {r['raw_rank']}  / 8")
    print(f"Lie-closure rank                 = {r['lie_rank']}  / 8")
    verdict = "FULL su(3)" if r['lie_rank']==8 else \
              ("only so(3)" if r['lie_rank']==3 else f"rank {r['lie_rank']}")
    print(f"\n   >>> traceless WZ curvature generates: {verdict}")

    print("\n"+"="*72)
    print("FAIR CONTROL -- real transverse frame (D5's so(3) case)")
    print("="*72)
    c = control_real_frame(m, kstar)
    print(f"real-symmetric single-site operator: {c['nsamp']} samples, "
          f"raw span={c['raw_rank']}, Lie rank={c['lie_rank']}  (expected <=3)")

    print("\n"+"="*72)
    print("ROBUSTNESS -- rank across parameters, base points, plaquette size")
    print("="*72)
    cases = [
        dict(alpha_s=0.6, alpha_t=0.9, A=0.30, N1=3, N4=3),
        dict(alpha_s=0.4, alpha_t=0.8, A=0.25, N1=3, N4=3),
        dict(alpha_s=0.7, alpha_t=0.95, A=0.40, N1=4, N4=4),
        dict(alpha_s=0.5, alpha_t=0.85, A=0.20, N1=4, N4=3, w1=1, w4=1),
        dict(alpha_s=0.6, alpha_t=0.9,  A=0.35, N1=5, N4=5, w1=2, w4=1),
    ]
    for cfg in cases:
        mm = build_helix(**cfg)
        ok = True
        for ks, hh in [((0.7,1.1,1.9,0.5),2e-3),((1.3,0.6,2.4,1.0),1e-3),
                       ((2.0,2.2,0.5,1.7),3e-3)]:
            rr = wz_su3(mm, np.array(ks), h=hh, verbose=False)
            ok = ok and (rr['lie_rank']==8)
        print(f"  {cfg}  resid={mm['resid']:.0e}  ->  su(3) rank 8 at all "
              f"(k*,h): {ok}")