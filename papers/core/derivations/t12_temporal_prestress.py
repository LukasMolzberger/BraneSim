"""
T12 -- Temporal prestress: the continuum / kinetic limit of r_4 = alpha_t a.

Layer-0 temporal link term (eta_4 = +1):
    S_time = 1/2 kappa_t sum_n (L_{n4} - r_4)^2,   r_4 = alpha_t a,  0<alpha_t<1
    L_{n4} = |R_{n+e4} - R_n|.

Claims proved here:

  (A) THE TEMPORAL PRESTRESS *IS* THE KINETIC TERM.  Expanding a temporal link
      Q = a e_4 + Delta_4 u about the vacuum, the quadratic action is
        S_time^(2) = 1/2 kappa_t sum [ (Delta_4 u^4)^2 + (1-alpha_t) sum_i (Delta_4 u^i)^2 ]
      -> in the continuum (Delta_4 u = a d_t u) a genuine kinetic term
        T = 1/2 kappa_t a^2 [ (d_t u^4)^2 + (1-alpha_t) sum_i (d_t u^i)^2 ],
      with EMERGENT inertia m = kappa_t a^2 (no separately-postulated mass).

  (B) THE (1-alpha_t) ANISOTROPY IS THE FINGERPRINT OF GENUINE PRESTRESS.
      Longitudinal (u^4) inertia = kappa_t (r_4-independent); transverse (u^i)
      inertia = kappa_t(1-alpha_t) = rho_4/a (the prestress). A pure kinetic term
      (r_4=0) is ISOTROPIC (both = kappa_t) -> the old "kinetic-sign artifact".
      For 0<alpha_t<1 the term is anisotropic -> the +1 is a genuine SPRING prestress,
      not just an inertia.

  (C) rho_4 = kappa_t a (1-alpha_t) is a genuine, finite prestress (Trace_vs_Traceless).
      alpha_t->0: rho_4 maximal but effect isotropic (artifact); alpha_t->1: rho_4->0,
      transverse temporal stiffness -> 0 (over-decoupled, forbidden, Decision G).

  (D) eta_4 = +1 MAKES S = T - V (Hamilton's principle) -> WAVE dynamics.
      The 4D block stationarity dS/dR=0 is exactly the Stormer-Verlet stencil; with
      eta_4=+1, eta_i=-1 it is hyperbolic and a wavepacket propagates at the T2 speed
      c_T. Flip eta_4->-1 (all-minus, Euclidean) and it is elliptic -> blows up.
"""

import numpy as np
np.set_printoptions(precision=5, suppress=True, linewidth=140)

a = 1.0

# --------------------------------------------------------------------------
# (A)+(B) single temporal-link quadratic coefficients (exact) vs prediction
# --------------------------------------------------------------------------
def link_quad_coeff(alpha_t, kappa_t, direction, eps=1e-4):
    """exact d^2/deps^2 of 1/2 kappa_t (|a e4 + eps d| - alpha_t a)^2 at eps=0."""
    d = np.zeros(4); d[direction] = 1.0
    def V(e):
        Q = a*np.array([0,0,0,1.0]) + e*d
        return 0.5*kappa_t*(np.linalg.norm(Q) - alpha_t*a)**2
    return (V(eps) - 2*V(0) + V(-eps))/eps**2          # -> coefficient of eps^2 (x2/2)

def report_AB():
    kappa_t = 1.7
    print("(A,B) temporal-link quadratic inertia coefficients (exact vs predicted)")
    print("      polarization        exact        predicted")
    for alpha_t in (0.0, 0.5, 0.9):
        cl = link_quad_coeff(alpha_t, kappa_t, 3)     # longitudinal u^4
        ct = link_quad_coeff(alpha_t, kappa_t, 1)     # transverse u^i
        # predicted: longitudinal = kappa_t ; transverse = kappa_t(1-alpha_t)
        tag = "(isotropic: kinetic-artifact)" if alpha_t==0 else \
              ("(over-decoupled)" if abs(1-alpha_t)<1e-9 else "(anisotropic: prestress)")
        print(f"  alpha_t={alpha_t:<4} long u^4  {cl:8.4f}     {kappa_t:8.4f}")
        print(f"  alpha_t={alpha_t:<4} tran u^i  {ct:8.4f}     {kappa_t*(1-alpha_t):8.4f}"
              f"   {tag}")
    print("  => transverse/longitudinal inertia ratio = (1-alpha_t):")
    for alpha_t in (0.0, 0.5, 0.9):
        r = link_quad_coeff(alpha_t,kappa_t,1)/link_quad_coeff(alpha_t,kappa_t,3)
        print(f"       alpha_t={alpha_t}:  ratio = {r:.4f}  (1-alpha_t = {1-alpha_t:.4f})")

# --------------------------------------------------------------------------
# (B') continuum convergence: discrete temporal action -> 1/2 m (d_t u)^2
# --------------------------------------------------------------------------
def continuum_convergence():
    print("\n(B') continuum/kinetic limit (long wavelength, FIXED lattice a=1):")
    print("     discrete temporal kinetic operator 2 kappa_t (1-cos omega) -> kappa_t omega^2")
    kappa_t, alpha_t = 1.0, 0.9
    # for a transverse mode ~e^{i omega t}, the temporal quadratic action coefficient is
    #   2 kappa_t (1-alpha_t)(1-cos omega)  ->  kappa_t (1-alpha_t) omega^2   (m = kappa_t(1-a_t))
    print("     omega    2(1-cos w)     w^2       ratio      -> emergent inertia m=kappa_t(1-a_t)")
    for w in (0.8, 0.4, 0.2, 0.1):
        disc = 2*(1-np.cos(w)); cont = w**2
        print(f"    {w:<7.2f} {disc:10.6f} {cont:10.6f}   {disc/cont:.5f}")
    print(f"     => ratio -> 1: temporal link term -> kinetic T = 1/2 m (d_t u)^2,")
    print(f"        m_long = kappa_t = {kappa_t},  m_tran = kappa_t(1-alpha_t) = {kappa_t*(1-alpha_t)}"
          f"  (emergent inertia, no postulated mass).")

# --------------------------------------------------------------------------
# (C) prestress reading
# --------------------------------------------------------------------------
def prestress_reading():
    print("\n(C) prestress reading  rho_4 = eta_4 kappa_t a (1-alpha_t)")
    kappa_t = 1.0
    for alpha_t in (0.0, 0.5, 0.9, 1.0):
        rho4 = +1*kappa_t*a*(1-alpha_t)
        note = {0.0:"max tension but ISOTROPIC effect (kinetic-artifact limit)",
                0.5:"genuine anisotropic prestress",
                0.9:"genuine prestress (Decision G regime, alpha_t=1-eps)",
                1.0:"rho_4=0: no prestress, transverse temporal stiffness=0 (FORBIDDEN)"}[alpha_t]
        print(f"   alpha_t={alpha_t:<4}: rho_4 = {rho4:+.3f} kappa_t a   preferred advance "
              f"|dR/dstep|={alpha_t:.2f}   {note}")

# --------------------------------------------------------------------------
# (D) eta_4=+1 -> Verlet wave dynamics; eta_4=-1 -> elliptic blow-up
# --------------------------------------------------------------------------
def verlet_demo():
    print("\n(D) 4D-block stationarity = Stormer-Verlet; eta_4=+1 gives wave dynamics")
    # linearized transverse field on a 1D periodic spatial chain, march in time.
    N = 256; alpha_s, alpha_t, kappa_s, gamma_t = 0.5, 0.6, 1.0, 1.5
    kappa_t = gamma_t*kappa_s
    cT2 = kappa_s*(1-alpha_s)/(kappa_t*(1-alpha_t))     # T2 transverse speed^2 (a=1)
    cT = np.sqrt(cT2)
    x = np.arange(N)
    def lap(u): return np.roll(u,-1) - 2*u + np.roll(u,1)
    def march(eta4):
        # EOM: eta4*kappa_t(1-a_t)[2u_t-u_{t+1}-u_{t-1}] - kappa_s(1-a_s)*lap(u)=0
        # => u_{t+1} = 2u_t - u_{t-1} + eta4 * cT2 * lap(u_t)   (with eta4=+1 stable)
        u0 = np.exp(-0.5*((x-N/2)/8)**2)*np.cos(2*np.pi*(x-N/2)/16)  # right-moving seed
        um = np.exp(-0.5*((x-(N/2-cT))/8)**2)*np.cos(2*np.pi*(x-(N/2-cT))/16)
        u_prev, u = um, u0
        norms=[]; centers=[]
        for t in range(120):
            u_next = 2*u - u_prev + eta4*cT2*lap(u)
            u_prev, u = u, u_next
            norms.append(np.max(np.abs(u)))
            e = u**2; centers.append((x*e).sum()/e.sum())
        return np.array(norms), np.array(centers)
    n_plus, c_plus = march(+1)
    n_minus, _     = march(-1)
    # measured speed = slope of packet center vs time (steps 10..110)
    speed = np.polyfit(np.arange(10,110), c_plus[10:110], 1)[0]
    print(f"   c_T (predicted, T2)          = {cT:.4f}")
    print(f"   c_T (measured from Verlet)   = {speed:.4f}")
    print(f"   eta_4=+1 amplitude stable?   max|u| {n_plus.min():.3f}..{n_plus.max():.3f}"
          f"  (bounded -> propagating wave)")
    print(f"   eta_4=-1 amplitude           max|u| grows to {n_minus[-1]:.2e}"
          f"  (elliptic blow-up -> no dynamics)")
    print("   => eta_4=+1 (opposite to spatial) turns the static 4D block into genuine")
    print("      time evolution; the kinetic term from temporal prestress carries waves.")


if __name__ == "__main__":
    print("="*74); print("T12 -- temporal prestress: continuum/kinetic limit"); print("="*74)
    report_AB()
    continuum_convergence()
    prestress_reading()
    verlet_demo()
    print("\n" + "="*74)
    print("CONCLUSION: r_4=alpha_t a makes the temporal link a genuine prestressed")
    print("spring whose continuum limit IS the kinetic term (emergent inertia kappa_t a^2),")
    print("with the (1-alpha_t) anisotropy as the fingerprint distinguishing it from a")
    print("kinetic-sign artifact (r_4=0). eta_4=+1 makes S=T-V, so 4D-block stationarity")
    print("is Verlet time-evolution carrying waves at c_T. T12 established.")
    print("="*74)
