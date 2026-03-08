# QCD and Lattice QCD Diagnostics to Refine a Baryon Simulation Roadmap

## Why Lattice QCD is a useful template for your experiment architecture

A key practical reason **lattice QCD** is relevant to your goals is that it already embodies the separation you want: a “clean” simulation core that evolves the system, and an **external** measurement/analysis layer that extracts spectra, matrix elements, and diagnostics from stored outputs. In modern lattice workflows, the *generation of configurations / trajectories* and the *measurement of observables* are deliberately separable so that analysis can be repeated, audited, upgraded, or replaced without changing the simulation history. As a result, methodological reliability is often discussed in terms of (i) what was actually simulated and stored, and (ii) what was later extracted from that stored data. This same separation is emphasized when discussing how energy levels are extracted: the core outputs (field configurations / trajectories) enable correlators, which enable spectral fits. citeturn51view0turn55view0turn55view2

This separation is not only a software preference—it is also tied to **fundamental lattice systematics**. Spectroscopy and scattering require correlator construction and fitting choices (operator bases, generalized eigenvalue problems, fit windows, treatment of finite time extent, etc.), which are *analysis-layer* decisions; the simulation core should not “peek” at whether a fit looks good and react mid-run. Lattice references explicitly warn that choices such as operator sets and Euclidean-time windows can determine whether energy levels are missed or contaminated, and that large Euclidean times reduce excited-state contamination but worsen noise and can invite thermal wrap-around artifacts. citeturn55view2turn53view1turn15view1turn51view1

## Working purely in simulation units and moving dimensional mapping into preparation

Your desire to remove dimensional mapping from experiment scripts is strongly aligned with how lattice QCD is normally conceptualized: the simulation natively produces **dimensionless** quantities (e.g., masses in units of the inverse lattice spacing, or lengths in lattice units), and *only later* do you decide how to map those to physical units by choosing a **scale-setting observable**. citeturn55view0turn55view1turn10view0

A robust pattern (mirroring lattice practice) is:

- **Experiment preparation**: decide a target physical point (or family of points), choose a scale-setting convention, produce a *single* “experiment input package” expressed entirely in simulation units plus provenance metadata.
- **Experiment run**: consume only simulation-unit inputs; produce only simulation-unit outputs.
- **Analysis**: interpret outputs, optionally convert to physical units *after the run* for reporting/comparison—not for in-run logic.

In lattice QCD, several widely used scale-setting choices are based on the **gradient flow**, which provides convenient reference scales often denoted \(t_0\) or \(w_0\). citeturn8view1turn11view0turn55view0  
A crucial lattice lesson is that *different* reasonable scale-setting choices mainly shift **cutoff artifacts** at finite resolution; only in a controlled continuum limit should physical predictions become independent of that choice. Translating that principle to your project: you can allow different preparation-time mappings (for exploratory calibration), but your *claims* should emphasize **dimensionless ratios and invariants** that do not depend on an arbitrary mapping choice. citeturn55view1turn55view0

A direct architectural consequence for your codebase is to keep anything analogous to scale-setting, calibration, or “dimensional mapping” in a **strictly separate preparation namespace**, so that experiment scripts cannot accidentally import it. This mirrors the idea that scale setting is a “meta-step” applied to results, not a dynamical ingredient of the simulation itself. citeturn55view0turn10view0

## QCD-inspired invariants and diagnostics that are especially informative for baryon claims

This section focuses on *what to measure* (offline) and *what to treat as invariants / consistency checks* that can tighten a baryon-oriented experimental setup.

image_group{"layout":"carousel","aspect_ratio":"16:9","query":["lattice QCD effective mass plot baryon correlator","lattice QCD Polyakov loop order parameter deconfinement","center vortex area law diagram lattice gauge theory","lattice QCD baryon spectrum octet decuplet plot"],"num_per_query":1}

### Confinement and screening observables that are hard to fake accidentally

If your theory claims to reproduce baryons as emergent bound states in a strong-interaction-like regime, you want at least one confinement-oriented diagnostic that is **not** “baryon-specific” (so it cannot be tuned only for baryons).

A standard lattice set of confinement diagnostics uses:

- **Area-law falloff of large gauge-invariant loop observables** (a proxy for a nonzero asymptotic string tension in confining regimes).
- A **temporal-loop order parameter** whose expectation value vanishes in the center-symmetric (confining) phase in pure-gauge settings, and becomes nonzero in a deconfined phase.
- Dual loop diagnostics (center-vortex free energy / dual loops) as complementary order parameters. citeturn49view0turn49view1turn49view2

Even if your microscopic ingredients differ from QCD, these suggest a general strategy: define *a* loop-like nonlocal diagnostic tied to your symmetry structure and test for a sharp qualitative change under temperature-like or coupling-like control parameters—*without using baryon measurements as input*.

### Baryon spectroscopy extraction as a “no free parameters” diagnostic layer

Lattice QCD extracts hadron masses and excited levels from the Euclidean-time behavior of correlation functions. The core idea is that correlators admit a spectral decomposition as sums of exponentials, and at sufficiently large Euclidean time the lowest state dominates. citeturn55view0turn51view1turn51view0

Operationally (and in a way that maps cleanly into your “analysis only” requirement):

- Build baryon-like operators (or whatever your theory’s baryon probes are).
- Compute two-point correlators \(C(t)\).
- Extract energies from exponential decay, often using an **effective mass** definition like \(\,m_\text{eff}(t) = -\log(C(t)/C(t-1))\,\) (in lattice-normalized time). citeturn51view1turn55view0

To obtain **multiple levels** (ground + excitations) and reduce operator-basis dependence, lattice spectroscopy frequently uses a correlator matrix and a generalized eigenvalue problem (GEVP). The GEVP framework makes precise statements about how extracted effective energies approach the true energies and how corrections fall exponentially with Euclidean time under suitable conditions. citeturn53view0turn53view1turn55view0

The translation to your project is straightforward: your simulation core should only output sufficient raw information (fields/particle trajectories) to build correlators. Your analysis layer then implements (i) simple effective-mass extraction and (ii) a correlator-matrix / eigenvalue method as a higher-grade diagnostic.

### Dimensionless baryon mass ratios and SU(3)-pattern tests

A particularly persuasive “phenomenological landscape” claim is not just “we got one baryon-like bound state,” but that you recover **structured families** of baryons with **mass relations** and **hierarchies** resembling QCD.

A minimal experimental target set (in nature) includes:

- Nucleon (proton), \(\Delta(1232)\), \(\Lambda\), \(\Sigma\), \(\Xi\), and \(\Omega\) benchmark masses. citeturn56view1turn56view2turn56view3turn56view4turn56view5
- The charged pion mass as the light “Goldstone-like” benchmark scale. citeturn44view0turn44view1

From these, you can define **dimensionless** diagnostics that are mapping-independent once you pick a single scale (or even without picking a scale, if you compare ratios only). Examples using PDG summary values:

- \(m_p/m_{\pi^\pm} \approx 938.27/139.57 \approx 6.72\). citeturn56view0turn44view0  
- \(m_\Omega/m_p \approx 1672.45/938.27 \approx 1.78\). citeturn56view5turn56view0  
- A classic SU(3)-octet mass relation is \(\frac{m_N+m_\Xi}{2} \approx \frac{3m_\Lambda+m_\Sigma}{4}\). Using representative PDG entries (e.g., \(p,\Sigma^+,\Xi^0,\Lambda\)) this holds at the sub‑percent level, making it a **tight pattern** to test in a simulation family scan. citeturn56view2turn56view3turn56view4turn56view0  

Lattice calculations explicitly investigate how small the deviations from this octet mass relation are as a function of quark masses, finding deviations to be small across a range of simulated pion masses—i.e., the relation is not a fragile coincidence of just the physical point, which makes it an excellent robustness diagnostic. citeturn26search10

For your roadmap, the actionable advice is: treat these relations as **diagnostic invariants** that your analysis layer computes automatically from stored spectra, while your simulation layer remains blind to them.

## Edge cases and systematics from lattice QCD that can make baryon evidence look stronger or weaker than it is

This is the area where lattice QCD offers the most “trap avoidance” value: it provides concrete examples where a qualitatively correct theory can look wrong (or vice versa) due to finite-size, contamination, or measurement-layer artifacts.

### Finite-volume sensitivity and the role of the lightest mode

A widely used lattice rule of thumb is that \(m_\pi L\) should be “comfortably large” (often quoted as \(m_\pi L \gtrsim 4\)) to suppress finite-volume distortions. citeturn23view0  
However, a key refinement (motivated by the physical picture that the baryon’s long-range tail lives *outside* the core) is that an even more relevant condition can be \(m_\pi (L-2R) \gg 4\), where \(R\) is an effective confinement/core radius; violating this refined condition can lead to qualitatively misleading behavior in hadronic properties as the pion mass is lowered. citeturn23view0turn23view1

For baryon-relevant observables, this matters because:

- Finite-volume effects can be *especially dramatic* for certain matrix elements; for example, strong finite-volume dependence has been highlighted for axial observables in model studies designed to mimic lattice boundary effects. citeturn23view2turn23view0  
- Finite-volume corrections are routinely expressed and estimated using effective field theory / forward-scattering inputs, and reliable evaluation may require going beyond leading asymptotic terms in regimes where \(m_\pi L\) is not large. citeturn25view1turn25view2turn25view0  

For your experimental design, the translatable rule is: whenever you are using a “lightest excitation” to set the long-distance physics (your analogue of the pion), treat the dimensionless product \((m_\text{light} \times L)\) as a *first-class control parameter* in your initialization scans, and record it in metadata for each run.

### Exponentially worsening baryon signal-to-noise and its consequences

Baryon correlators are notorious in lattice QCD because the signal-to-noise ratio degrades exponentially with Euclidean time. A standard quantified example is that the nucleon signal-to-noise ratio squared behaves like \(\exp\{-(2M_N-3M_\pi)\,t\}\), making long time separations simultaneously desirable (for suppressing excited states) and expensive (for statistics). citeturn15view1

Two practical consequences you can adopt as “edge-case checks”:

- If an extracted baryon mass depends strongly on the chosen analysis time window, you may be seeing exactly the excited-state/noise tradeoff that lattice warns about.
- For matrix elements such as an axial coupling, lattice-oriented discussions note that physically clean source–sink separations can be substantially larger than what is often affordable, which can bias results if not addressed. citeturn15view1turn55view2

### Excited-state contamination, operator incompleteness, and “missing levels”

Lattice spectroscopy texts emphasize that (i) the effective-mass approach assumes the ground state dominates at large Euclidean time, and (ii) practical extractions require good overlap of operators with the state and sufficient statistics. citeturn51view1turn55view0

In practice, to obtain multiple levels and reduce contamination, the GEVP approach is used; it provides explicit formulas for how effective energies approach the true spectrum and highlights that small energy gaps can make corrections large unless time separations are chosen carefully. citeturn53view1turn55view0

A particularly important “gotcha” for resonance/scattering channels is that **operator choice can cause you to miss finite-volume energy levels**, leading to a wrong physical interpretation. Reviews explicitly state that including both single-hadron and multi-hadron operators can be necessary; otherwise energy levels may be absent from the extracted spectrum. citeturn55view2turn55view0

### Finite time extent and thermal wrap-around (“thermal pollution”)

With periodic (or anti-periodic) boundary conditions in the time direction, correlators can receive contributions from states that wrap around the temporal extent, modifying naive single-exponential behavior. This is flagged both in basic effective-mass discussions and in resonance/scattering reviews as a practical source of distortion at large Euclidean times. citeturn51view1turn55view2

For your design, the actionable point is: if your simulation uses periodic time (or long but finite observation windows), treat wrap-around contributions as an analysis-layer systematic and include a diagnostic that checks for it by varying temporal extent.

## A concrete baryon-focused experimental roadmap in simulation units

The roadmap below is designed to (i) stay strictly in simulation units during runs, (ii) avoid any in-run feedback from diagnostics, and (iii) steadily increase the “hard-to-fake” content of your baryon claims.

### Establish a scale-free baseline and long-distance control parameter

Start by ensuring you can define a **lightest excitation scale** (your analogue of the pion mass) and a finite-volume control parameter \(m_\text{light}L\). This is foundational because lattice experience shows that finite-volume systematics are driven by the lightest modes and can affect even heavy states if they couple to the light sector. citeturn23view0turn25view0turn55view1

Your analysis should report, per run, dimensionless “environment” descriptors such as:

- \(m_\text{light}L\) (finite volume),
- \(T_\text{time} \times m_\text{light}\) (finite temporal extent),
- any resolution/cutoff proxy your discretization implies (your analogue of “lattice spacing artifacts”). citeturn55view1turn55view2

### Demonstrate a single-baryon ground state with lattice-style extraction discipline

Next, focus on a single baryon channel and adopt the lattice extraction discipline:

- Define baryon-like operators and compute two-point correlators.
- Use an effective-mass style estimator as a first pass. citeturn51view1turn55view0
- Upgrade to a correlator matrix + GEVP workflow to extract the ground state and at least one excitation in a controlled way, with an explicit study of time-window dependence (excited-state contamination vs noise). citeturn53view0turn53view1turn55view0

**Deliverable-style diagnostic outputs** that strengthen credibility (without feeding back into the run) include: fit stability plots vs fit window, comparisons of effective-mass plateaus vs GEVP energies, and a record of how results change with operator basis size.

### Expand to a baryon family and test SU(3)-pattern invariants

Once a stable baryon extraction exists, expand initialization to support a *family* of baryon-like states (by tuning whatever parameter in your theory plays the role of flavor content or mass splitting). The analysis should then test “pattern invariants” such as the octet mass relation \(\frac{m_N+m_\Xi}{2} \approx \frac{3m_\Lambda+m_\Sigma}{4}\) and quantify deviation as a dimensionless fraction. citeturn26search10turn56view2turn56view3turn56view4turn56view0

For reporting (and for demonstrating contact with real-world phenomenology), it is reasonable to compare a small set of mass ratios to experimental ones quoted by entity["organization","Particle Data Group","particle properties review"] (e.g., \(m_p/m_{\pi^\pm}\), \(m_\Omega/m_p\), etc.), while making clear that *only one* observable was used for any mapping and all other comparisons are predictions. citeturn44view0turn56view5turn56view0

### Add baryon matrix-element probes only after volume systematics are controlled

Matrix elements (axial couplings, form-factor-like observables, etc.) can be much more sensitive to finite volume and excited-state contamination than masses. Model-based studies emphasize large finite-volume dependence for axial-type observables at light masses/long ranges, and lattice discussions stress the tension between needing large separations for clean extraction and the exponential noise problem in baryonic correlators. citeturn23view2turn15view1turn55view2

So, make matrix-element experiments a later stage where you already have:

- at least two volumes for a finite-volume trend check,
- a demonstrated ability to vary Euclidean-time (or time-window) separation and show stability within uncertainties,
- an operator basis that includes multi-hadron-like structures if your channel allows mixing. citeturn55view2turn53view1

### Move to baryon interactions and finite-volume scattering only after the spectrum is robust

A credible demonstration of a “rich landscape” often includes not only bound states but also interactions. Lattice QCD connects finite-volume energy levels to scattering information via finite-volume quantization conditions, while explicitly warning that discretization artifacts and exponential finite-volume effects must be controlled to apply these methods cleanly. citeturn55view1turn55view2

The practical translation: if you plan to infer interaction parameters from energy shifts in a finite box, you must (i) generate a sufficiently dense set of energy levels, (ii) vary volume, and (iii) treat exponential corrections as a first-class systematic.

## Code separation guidance that matches these constraints

A clean implementation pattern that enforces your “no feedback from analysis into the run” rule is to define **one-way data flow**:

- preparation → experiment config (simulation units only),
- experiment run → raw outputs + minimal in-run invariants (energy conservation logs, constraint violations),
- analysis → derived observables, diagnostics, plots.

This parallels the lattice narrative that energy levels come from correlators built from stored trajectories/configurations, and that the quality of extraction depends on analysis choices rather than on altering the run mid-stream. citeturn55view0turn55view2turn51view1

A folder split that enforces this mechanically (imports and CI checks can make it hard to violate) is:

- `prep/`  
  Purpose: scale setting, dimensional mapping, physical-to-simulation conversions, calibration scripts.  
  Output: a sealed “experiment package” (e.g., `experiment_config_sim_units.json` + provenance record).
- `experiments/`  
  Purpose: initialization + simulation evolution only, consuming simulation-unit configs and producing raw outputs.  
  Rule: must not import from `prep/`.
- `analysis/`  
  Purpose: correlators, GEVP, spectrum fits, finite-volume studies, pattern tests (octet relations), plotting.  
  Rule: reads experiment outputs; never writes back into experiment state (except human-authored decisions).

To make the split auditable, store in every run’s metadata:

- the exact simulation-unit parameters used,
- the preparation provenance (hash of prep scripts + mapping constants),
- the finite-volume control parameter \(m_\text{light}L\) (once measured),
- temporal extent and boundary conditions (for thermal pollution checks). citeturn23view0turn51view1turn55view2

Finally, if you want an external “quality stamp” template, the entity["organization","Flavour Lattice Averaging Group","lattice qcd averages"] approach is instructive: it emphasizes judging results by whether major systematics (finite volume, cutoff effects, and unphysical mass extrapolations) are controlled and documented—exactly the kind of structure that turns an impressive plot into a persuasive scientific claim. citeturn55view1turn55view2turn6view1turn6view0