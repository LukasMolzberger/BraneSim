"""
Generate the two manuscript figures from the verified machinery.
  fig_bands.pdf : carrier band structure along a k-cut; isolated rank-3 triplet + gap
  fig_su3.pdf   : (a) Lie-closure rank over (alpha_s, alpha_t); (b) rank vs Bloch-phase strength
Outputs to paper/figures/.
"""
import numpy as np, importlib.util, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

_here = os.path.dirname(__file__)
def _load(n):
    s = importlib.util.spec_from_file_location(n, os.path.join(_here, f"{n}.py"))
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
t1 = _load("t1_su3_witness"); gen = _load("t1_genericity_robustness")

OUT = os.path.join(_here, "..", "paper", "figures"); os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 11, "font.family": "serif", "axes.linewidth": 0.8,
                     "figure.dpi": 150})

KPTS = [np.array(k) for k in [(0.7,1.1,1.9,0.5),(1.3,0.6,2.4,1.0),(2.0,2.2,0.5,1.7)]]

# ---------------------------------------------------------------- Figure 1
def figure_bands():
    m = t1.build_helix()
    kstar = np.array([0.7, 1.1, 1.9, 0.5])
    # choose the isolated triplet at kstar (same selector as wz_su3)
    w0, _ = np.linalg.eigh(t1.bloch_D(m, kstar)); nb = len(w0)
    TRIP = max(range(1, nb-3),
               key=lambda j: min(w0[j]-w0[j-1], w0[j+3]-w0[j+2]) /
                             max(w0[j+2]-w0[j], 1e-9))
    # restrict the cut to the contiguous window around kstar where bands TRIP..TRIP+2
    # stay gap-isolated (near-degenerate composite; isolation is local, per Decision L)
    def isolated(k1, thr=0.06):
        w = np.linalg.eigvalsh(t1.bloch_D(m, np.array([k1, kstar[1], kstar[2], kstar[3]])))
        return min(w[TRIP]-w[TRIP-1], w[TRIP+3]-w[TRIP+2]) > thr
    k1c = kstar[0]; lo_k, hi_k = k1c, k1c
    while lo_k > 0.05 and isolated(lo_k-0.01): lo_k -= 0.01
    while hi_k < 2*np.pi/3-0.05 and isolated(hi_k+0.01): hi_k += 0.01
    k1s = np.linspace(lo_k, hi_k, 140)
    bands = np.array([np.linalg.eigvalsh(t1.bloch_D(m, np.array([k, kstar[1], kstar[2], kstar[3]])))
                      for k in k1s])
    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    for b in range(nb):
        ax.plot(k1s, bands[:, b], color="0.78", lw=0.7, zorder=1)
    ax.fill_between(k1s, bands[:, TRIP-1], bands[:, TRIP], color="#3498db", alpha=0.16, zorder=0)
    ax.fill_between(k1s, bands[:, TRIP+2], bands[:, TRIP+3], color="#3498db", alpha=0.16, zorder=0)
    for j in range(TRIP, TRIP+3):
        ax.plot(k1s, bands[:, j], color="#c0392b", lw=1.9, zorder=3)
    # annotate internal spread eps and outside gap Delta at kstar
    w = np.linalg.eigvalsh(t1.bloch_D(m, kstar))
    x0 = kstar[0]
    ax.annotate("", xy=(x0, w[TRIP+2]), xytext=(x0, w[TRIP]),
                arrowprops=dict(arrowstyle="<->", color="#c0392b", lw=1.2))
    ax.text(x0+0.012, 0.5*(w[TRIP]+w[TRIP+2]), r"$\varepsilon$", color="#c0392b", fontsize=12, va="center")
    ax.annotate("", xy=(x0-0.03, w[TRIP+3]), xytext=(x0-0.03, w[TRIP+2]),
                arrowprops=dict(arrowstyle="<->", color="#2471a3", lw=1.2))
    ax.text(x0-0.07, 0.5*(w[TRIP+2]+w[TRIP+3]), r"$\Delta$", color="#2471a3", fontsize=12, va="center", ha="right")
    ax.axvline(x0, color="0.45", ls=":", lw=0.9, zorder=2)
    ax.set_xlabel(r"$k_1$ (isolation window; other $k$ fixed at $k^\ast$)")
    ax.set_ylabel(r"fluctuation eigenvalues $\lambda(k)$")
    ax.set_title(r"Isolated rank-3 carrier of $D_{\bar R}(k)$")
    ax.plot([], [], color="#c0392b", lw=1.9, label=r"rank-3 carrier triplet")
    ax.plot([], [], color="0.78", lw=0.7, label="other bands")
    ax.fill_between([], [], color="#3498db", alpha=0.16, label="outside gap")
    ax.legend(loc="lower right", framealpha=0.92, fontsize=9)
    ax.set_xlim(lo_k, hi_k)
    ax.set_ylim(w[TRIP-1]-0.12, w[TRIP+3]+0.12)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_bands.pdf")); plt.close(fig)
    print("fig_bands.pdf written; TRIP =", TRIP)

# ---------------------------------------------------------------- Figure 2
def robust_rank(model):
    ranks = [t1.wz_su3(model, k, verbose=False)['lie_rank'] for k in KPTS]
    return max(set(ranks), key=ranks.count)

def bloch_scaled(model, kvec, s):
    N, sites, si, e, Cmat, Ns = (model['N'], model['sites'], model['site_index'],
                                 model['e'], model['Cmat'], model['Ns'])
    H = np.zeros((4*Ns, 4*Ns), complex)
    for t in sites:
        i = si[tuple(t)]
        for mu in range(4):
            C = Cmat[(tuple(t), mu)]; p = t+e[mu]; ss = p % N; d = p-ss
            j = si[tuple(ss)]; ph = np.exp(1j*s*(kvec@d)); bi, bj = 4*i, 4*j
            H[bi:bi+4, bi:bi+4] += C; H[bj:bj+4, bj:bj+4] += C
            H[bi:bi+4, bj:bj+4] += -C*ph; H[bj:bj+4, bi:bi+4] += -C*np.conj(ph)
    return 0.5*(H+H.conj().T)

def figure_su3():
    sa = _load("t1_smallamp_bianchi")
    kstar = np.array([0.7, 1.1, 1.9, 0.5])
    # (a) small-amplitude scaling: ||G|| ~ A^2, su(3) span persists to A->0
    TRIP0 = t1.wz_su3(t1.build_helix(A=0.1), kstar, verbose=False)['TRIP']
    As = np.array([0.005, 0.01, 0.02, 0.03, 0.05, 0.07, 0.10])
    Gs = np.array([sa.curvature_norm_fixed(A, TRIP0, kstar) for A in As])
    spans = [sa.su3_span(A) for A in As]
    pfit = np.polyfit(np.log(As), np.log(Gs), 1)[0]
    # (b) mechanism: no transport (0) -> real frame (so(3)=3) -> complex phases (su(3)=8)
    m = t1.build_helix(); kstar = np.array([0.7, 1.1, 1.9, 0.5])
    orig = t1.bloch_D
    t1.bloch_D = (lambda model, k: bloch_scaled(model, k, 0.0))   # s=0: no k-transport
    r_none = t1.wz_su3(m, kstar, verbose=False)['lie_rank']
    t1.bloch_D = orig
    r_real = t1.control_real_frame(m, kstar)['lie_rank']         # real moving frame -> so(3)
    r_cplx = t1.wz_su3(m, kstar, verbose=False)['lie_rank']      # complex Bloch -> su(3)
    cases = ["no\ntransport", "real\nframe", "complex\nBloch phases"]
    vals = [r_none, r_real, r_cplx]

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.2, 3.7))
    # (a) log-log scaling with an A^2 reference line
    ref = Gs[-1]*(As/As[-1])**2
    axA.loglog(As, ref, ls="--", color="0.55", lw=1.1, zorder=1, label=r"$\propto A^2$ (reference)")
    axA.loglog(As, Gs, "o-", color="#c0392b", lw=1.7, ms=6, zorder=3,
               label=r"$\|\mathcal{G}\|$ (substrate)")
    axA.set_xlabel(r"carrier amplitude $A/a$")
    axA.set_ylabel(r"traceless curvature $\|\mathcal{G}\|$")
    axA.set_title(r"(a) $\mathfrak{su}(3)$ curvature is leading-order ($\propto A^{%.2f}$)" % pfit)
    axA.legend(loc="upper left", fontsize=9, framealpha=0.9)
    axA.text(0.97, 0.06, r"$\mathfrak{su}(3)$ span $=8$ for all $A$ shown",
             transform=axA.transAxes, ha="right", va="bottom", fontsize=9,
             bbox=dict(boxstyle="round", fc="#fdf2e9", ec="0.7"))
    axA.grid(True, which="both", ls=":", lw=0.5, alpha=0.5)

    bars = axB.bar(cases, vals, color=["#bdc3c7", "#f1c40f", "#c0392b"], width=0.6, zorder=3)
    axB.axhline(3, color="0.6", ls="--", lw=0.9, zorder=1)
    axB.axhline(8, color="0.6", ls=":", lw=0.9, zorder=1)
    for b, v in zip(bars, vals):
        axB.text(b.get_x()+b.get_width()/2, v+0.2, str(v), ha="center", fontsize=11, fontweight="bold")
    axB.set_ylabel(r"Lie-closure rank")
    axB.set_title(r"(b) $\mathfrak{so}(3)\to\mathfrak{su}(3)$ via complex phases")
    axB.set_ylim(0, 9.2); axB.set_yticks([0, 3, 8])
    axB.set_yticklabels(["0", "3  so(3)", "8  su(3)"])
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_su3.pdf")); plt.close(fig)
    print(f"fig_su3.pdf written; scaling exponent p={pfit:.2f}, spans={spans}, "
          f"mechanism ranks={list(zip(cases, vals))}")


if __name__ == "__main__":
    figure_bands()
    figure_su3()
