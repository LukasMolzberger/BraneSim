Here’s a “Claude-ready” implementation plan, structured **by file** and with enough detail that it can be applied directly to the LaTeX sources.

---

## 1. `draft-paper.tex` – global structure & appendices

1. **Keep the current top-level section order**, but make sure the `\input` order reflects this:

   ```tex
   \input{abstract}
   \input{introduction}
   \input{conceptual-model}
   \input{reconstructing-physics}
   \input{experimental-setting}
   \input{discussion}
   \input{conclusion}
   \appendix
   \input{appendix-curvature}
   % (new) optional: \input{appendix-electron-tube}
   ```

2. **Do not input the disclaimer any more**:

   * Ensure there is **no** line like `\input{disclaimer}` anywhere in `draft-paper.tex`.
   * The file `disclaimer.tex` can stay in the repo as scratchpad, but is not used.

3. **Add a hook for the new electron appendix (optional but recommended)**:

   * After `\input{appendix-curvature}`, add:

     ```tex
     \input{appendix-electron-tube}
     ```
   * Claude can create `appendix-electron-tube.tex` as described in section 5 below.

---

## 2. `conceptual-model.tex` – reordering, renaming, clarifications

### 2.1 Reorder linear / threshold subsections

Assuming you currently have subsections like:

```tex
\subsection{Continuous brane embedding}
...
\subsection{Elastic Lagrangian and geometric coupling}
...
\subsection{Equations of motion}
...
\subsection{Threshold localization}
...
\subsection{Linear wave regime and isotropy}
...
\subsection{Gravity as an effective metric}
...
```

1. **Rename and move the threshold subsection**:

   * Change the title of `\subsection{Threshold localization}` to:

     ```tex
     \subsection{Geometric nonlinearity and localization threshold}
     \label{subsec:geom-threshold}
     ```
   * Move the *entire* subsection (from this `\subsection` line down to just before the next `\subsection`) so that it appears **after** the “Linear wave regime and isotropy” subsection, i.e. the order becomes:

     ```tex
     \subsection{Continuous brane embedding}
     ...
     \subsection{Elastic Lagrangian and geometric coupling}
     ...
     \subsection{Equations of motion}
     ...
     \subsection{Linear wave regime and isotropy}
     ...
     \subsection{Geometric nonlinearity and localization threshold}
     ...
     \subsection{Gravity as an effective metric}
     ...
     ```

2. In the text where you first talk about “threshold localization”, update wording to “geometric nonlinearity and localization threshold” to match the new title.

### 2.2 Add a short “static gauge / index conventions” paragraph

In the **“Continuous brane embedding”** subsection, after you first introduce the embedding map and the induced metric (g_{ij}), add a short explicit paragraph such as:

```tex
We work throughout in a static gauge. Spatial brane coordinates are denoted by $x^a$ with indices $a,b=1,2,3$, and time is a separate absolute coordinate $t$. In this gauge we set
\begin{equation}
X^a(x,t) = x^a, \qquad X^4(x,t) = \xi(x,t),
\end{equation}
so that all deformations occur in the fourth (lateral) direction. Latin indices $i,j$ are used for spatial indices $1,2,3$, while time derivatives are written explicitly as $\partial_t$.
```

Make sure this paragraph comes **before** the stress–strain and energy density formulas.

### 2.3 Clarify wave-speed identification

In the **“Linear wave regime and isotropy”** subsection:

1. Where you state that small-amplitude waves propagate with speed (c^2 = T/\rho_m), make this explicit. For example, replace the existing sentence around the wave equation with:

   ```tex
   In the small-slope, small-strain limit the transverse displacement $\xi(x,t)$ obeys a linear wave equation of the form
   \begin{equation}
   \partial_t^2 \xi = c^2 \nabla^2 \xi + \dots,
   \label{eq:linear-wave-equation}
   \end{equation}
   where the effective wave speed is
   \begin{equation}
   c^2 = \frac{T}{\rho_m}.
   \end{equation}
   Throughout this work we identify this emergent wave speed with the physical speed of light $c$.
   ```

2. If `eq:linear-wave-equation` already exists, just ensure the wording around it contains the explicit identification and is **not** contradictory.

### 2.4 Link to curvature appendix from threshold subsection

In the **renamed “Geometric nonlinearity and localization threshold”** subsection, add an explicit reference to the Gaussian-bulge curvature appendix:

```tex
As a concrete example of this geometric nonlinearity, Appendix~\ref{app:bulge-curvature} shows that a Gaussian bulge of amplitude $A$ and width $\sigma$ in the fourth direction induces a three-dimensional scalar curvature at its center that scales as
\begin{equation}
R_3(0) \propto \frac{A^2}{\sigma^4}.
\end{equation}
This illustrates how high local curvature (or, equivalently, large amplitude gradients on Compton scales) can act as a geometric threshold for localization.
```

*(Adjust the label `app:bulge-curvature` if needed after step 5.1.)*

---

## 3. `reconstructing-physics.tex` – section names, gravity eq, tubular derivation → appendix

### 3.1 Rename some sections for clarity

Find the section headings:

* `\section{Wave Mechanics and Quantum Behavior}`
* `\section{Charge as Amplitude Deformation}`

and change them to:

```tex
\section{Quantum sector: envelope and Schr\"odinger limit}
\label{sec:quantum-envelope}
```

```tex
\section{Electromagnetic sector: charge from amplitude}
\label{sec:charge-amplitude}
```

If `\label{sec:charge-amplitude}` is already present, keep the same label; just ensure the section title is updated.

### 3.2 Fix the broken reference to `subsec:charge-amplitude`

Search in *all* `.tex` files for `\ref{subsec:charge-amplitude}` and replace each occurrence by:

```tex
\ref{sec:charge-amplitude}
```

Do **not** keep the `subsec:` prefix anywhere.

### 3.3 Re-introduce `eq:PhiG_from_gradX4` in a compact form

To fix the broken references and keep the gravitational story coherent:

1. In `reconstructing-physics.tex`, in the part where you discuss **emergent Newtonian gravity / gravitational potential** (likely in a subsection under an “Emergent gravity / Einstein-like sector” heading), add the following equation and label:

   ```tex
   We define an emergent Newtonian gravitational potential $\Phi_G(x)$ directly from the embedding into the fourth dimension by
   \begin{equation}
   \Phi_G(x) = \zeta_G c^2 \bigl( X^4(x) - X^4_0 \bigr),
   \label{eq:PhiG_from_gradX4}
   \end{equation}
   where $X^4_0$ is a reference level of the brane and $\zeta_G$ is a dimensionless coupling parameter. The corresponding gravitational acceleration seen by internal observers is
   \begin{equation}
   g_i(x) = -\partial_i \Phi_G(x) = -\zeta_G c^2 \,\partial_i X^4(x).
   \end{equation}
   ```

2. Make sure the text around this equation explains briefly that this is a **definition** of the Newtonian potential in terms of the embedding.

3. Do **not** add an extra Poisson equation here if it is already defined elsewhere; if needed, just reference the existing Poisson-like relation from this point.

Now the references to `eq:PhiG_from_gradX4` in other files will become valid.

### 3.4 Strengthen the carrier–envelope / Schrödinger derivation (no re-order, just clarify)

In the **“Quantum sector: envelope and Schrödinger limit”** section:

1. Just before or after the derivation where you write down the carrier–envelope ansatz, add an explicit small-parameter statement. For example:

   ```tex
   We assume a narrow-band excitation around a carrier frequency $\omega_C$ and wavenumber $k_C$, such that the envelope $\Phi(x,t)$ varies slowly in space and time:
   \begin{equation}
   \bigl|\partial_t \Phi\bigr| \ll \omega_C |\Phi|, 
   \qquad
   |\nabla \Phi| \ll k_C |\Phi|.
   \end{equation}
   This hierarchy justifies keeping only first-order derivatives of the envelope in the effective equation of motion.
   ```

2. Where you define or discuss the effective (\hbar), add one clarifying sentence about dimensions. For example:

   ```tex
   Here $\hbar$ is defined by
   \begin{equation}
   \hbar = \mathcal{C}_{\text{brane}} \,\omega_C,
   \end{equation}
   where $\mathcal{C}_{\text{brane}}$ is a brane-dependent constant with dimensions of action times time. This identification is chosen such that the envelope energy per Compton cell matches the standard quantum-mechanical normalization.
   ```

   Adjust the wording and symbol `\mathcal{C}_{\text{brane}}` to match your existing notation.

### 3.5 Move heavy tubular-electron algebra into a new appendix

You said you generally want complex derivations in appendices, so:

1. In `reconstructing-physics.tex`, locate the section(s) and subsection(s) that contain the detailed Frenet–Serret / tubular coordinate derivation around the toroidal electron. Typical headings to move:

   * `\subsection{Tubular coordinates around the toroidal electron}`
   * `\subsubsection{Frenet--Serret frame and tubular coordinates}`
   * `\subsubsection{Metric components and energy density}`
   * Or similar sections where you write down explicit metric components, determinants, curvature terms, and long formulas.

2. **Keep a short summary** in the main text:

   * Replace the bulk of that derivation by ~1–3 paragraphs that:

     * explain the physical picture (closed filament, tubular neighborhood),
     * state the key **scaling results** (e.g. how energy depends on radius, how spin arises),
     * and refer to the new appendix, e.g.:

     ```tex
     A detailed derivation of the tubular metric, energy density, and curvature around such a toroidal filament is given in Appendix~\ref{app:electron-tube}. Here we only summarize the qualitative picture and the key scalings.
     ```

3. **Move the full derivation** (all long equations and accompanying explanations) into a new file `appendix-electron-tube.tex` and wrap it as an appendix section, e.g.:

   ```tex
   \section{Tubular electron soliton: detailed derivation}
   \label{app:electron-tube}

   % (Paste the full Frenet--Serret / tubular coordinates derivation here)
   ```

4. Ensure `draft-paper.tex` includes this new appendix after `appendix-curvature.tex` (see 1.3 above).

---

## 4. `experimental-setting.tex` – nonlinearity cleanup, references, SI mapping

### 4.1 Remove / quarantine explicit “saturation potential” text

Search in `experimental-setting.tex` for any text that:

* talks about a “saturation potential”,
* describes a “hard cutoff” on amplitude or strain as part of the **model**,
* or introduces a separate nonlinear potential term like (V(\epsilon)) that is not purely geometric.

Then:

1. If such a mode is not actually used for reported results, either **delete** those paragraphs entirely or rephrase them as a short historical note. Recommended:

   ```tex
   % Old exploratory mode removed:
   % In early experiments we also implemented an additional material saturation potential...
   ```

   Or, if you want to keep a short mention:

   ```tex
   In some exploratory simulations (not reported here) we experimented with an additional material saturation potential that explicitly clipped large strains. In the present work, however, all simulations and conclusions are based solely on the geometric nonlinearity induced by the brane's embedding; no extra saturation potential is active.
   ```

2. Ensure that nowhere in the main narrative do you suggest that such a saturation potential is part of the **core** model.

### 4.2 Fix references to `eq:PhiG_from_gradX4` and conceptual model label

1. Search for all occurrences of `eq:PhiG_from_gradX4` in `experimental-setting.tex`. These references should now **point to the equation you reintroduced in section 3.3**, so leave the references themselves as they are.

2. Search for `\ref{sec:conceptual-model}` or `\ref{sec:conceptual_model}` etc.; ensure they use the **actual label** defined in `conceptual-model.tex`.

   * If the label in `conceptual-model.tex` is `\label{conceptual-model}`, then in `experimental-setting.tex` all references should be:

     ```tex
     \ref{conceptual-model}
     ```

   * Replace any `\ref{sec:conceptual-model}` or similar variants accordingly.

### 4.3 Add a concise explanation of the SI→dimensionless mapping

In the part of `experimental-setting.tex` where you explain how physical parameters are mapped to dimensionless units (the mapping layer you mentioned), add a short, explicit paragraph like:

```tex
All simulations are carried out in dimensionless units. We choose the Compton wavelength $\lambda_C$ and the corresponding frequency $\omega_C = 2\pi c/\lambda_C$ as basic scales, and measure
\begin{equation}
x' = \frac{x}{\lambda_C}, \qquad t' = \omega_C t,
\end{equation}
so that the dimensionless wave speed is $c' = 1$. The elastic parameters $(T,\rho_m,\kappa)$ are rescaled accordingly. Importantly, this mapping leaves the algebraic structure of the equations of motion unchanged; it serves only to keep numerical values in a tractable range and does not add extra constraints to the model.
```

Adjust variable names / primes to match your existing notation. The key point: explicitly state that the mapping is **transparent**.

### 4.4 Tighten “Interpretation layers” subsection

In the subsection that explains “interpretation layers” (field measurements vs fundamental variables):

1. Compress the text into a bullet list, e.g.:

   ```tex
   \subsection{Interpretation layers (clarity note)}

   For clarity it is useful to distinguish three layers:
   \begin{enumerate}
     \item The fundamental brane dynamics: the embedding $X^A(x,t)$ and its elastic equations of motion.
     \item Derived observables: local energy density, curvature, and other geometric quantities computed directly from $X^A$.
     \item Emergent fields: effective gravitational and electromagnetic fields inferred from suitable averages of the derived observables.
   \end{enumerate}
   Only the first layer is hard-coded into the simulation; the latter two are purely diagnostic and can be modified without changing the underlying dynamics.
   ```

2. Make sure this subsection is referenced once from the Discussion (e.g. “As discussed in Section~\ref{...}, these emergent fields are interpretation-layer objects…”).

---

## 5. `appendix-curvature.tex` – label & minor wording

### 5.1 Ensure the label matches `app:bulge-curvature`

At the top of `appendix-curvature.tex`, where the main appendix section starts, set:

```tex
\section{Curvature induced by a spherical bulge}
\label{app:bulge-curvature}
```

If there is some other label currently (e.g. `\label{app:curvature-bulge}`), change it to `app:bulge-curvature`.

Then search in all `.tex` files for `app:bulge-curvature` and update references if needed (or vice versa).

### 5.2 Emphasize scaling rather than exact coefficient

Where you state the final result for the scalar curvature at the center of the Gaussian bulge, slightly soften the equality, e.g. replace

```tex
R_3(0) = \frac{6 A^2}{\sigma^4}
```

by

```tex
R_3(0) \approx \frac{6 A^2}{\sigma^4},
```

and then add a short sentence like:

```tex
The precise numerical coefficient depends on the details of our approximations, but the scaling $R_3(0) \propto A^2/\sigma^4$ is robust.
```

---

## 6. New `appendix-electron-tube.tex` – full tubular derivation

Create a new file `appendix-electron-tube.tex` with something like:

```tex
\section{Tubular electron soliton: detailed derivation}
\label{app:electron-tube}

% (Paste here the full Frenet--Serret and tubular metric derivation that was
%  removed from the main text in Section~\ref{sec:charge-amplitude} or its
%  surrounding subsections.)
```

* Paste all the detailed Frenet–Serret formulas, metric components, determinants, curvature, and energy integrals that are now removed from `reconstructing-physics.tex`.
* Keep the structure of subsections if needed (`\subsection{Frenet--Serret frame}`, etc.), but this all lives under the appendix section.
* Ensure all equation labels used in the moved text remain unique.

---

## 7. `discussion.tex` and `conclusion.tex` – entanglement & parameter-match clarity

### 7.1 Explicitly mark entanglement as an open problem

In `discussion.tex`, find the paragraph where you mention entanglement / Bell violations / nonlocality.

* Replace any phrasing that sounds like “we might explain entanglement via X” by something firm like:

  ```tex
  We emphasize that the present brane model does not yet provide a detailed account of entanglement, Bell-inequality violations, or the full quantum measurement problem. These remain open questions for future work.
  ```

* Make sure this statement is unambiguous and not hedged.

### 7.2 Clarify robustness of parameter matches

In either `discussion.tex` or `conclusion.tex` (where you summarise the parameter identification), add a short paragraph like:

```tex
Among the various parameter identifications, the mapping between the brane wave speed and the speed of light $c$, and the use of the Compton frequency to define $\hbar$, are relatively robust consequences of the elastic brane dynamics. By contrast, the proposed mappings for the elementary charge $e$ and the gravitational constant $G$ involve more phenomenological assumptions about how internal observers infer effective electromagnetic and gravitational fields. These identifications should therefore be regarded as conjectural and subject to further refinement.
```

This keeps the storyline honest with respect to mathematical status.

---

## 8. Quick implementation checklist for Claude

Here’s a compact checklist Claude can follow:

1. **draft-paper.tex**

   * Ensure correct `\input` order, no `\input{disclaimer}`.
   * Add `\input{appendix-electron-tube}` after `appendix-curvature`.

2. **conceptual-model.tex**

   * Move and rename “Threshold localization” → “Geometric nonlinearity and localization threshold” after “Linear wave regime and isotropy”.
   * Add explicit static gauge / index paragraph.
   * Make explicit (c^2 = T/\rho_m) and identification with physical (c).
   * Reference Appendix `app:bulge-curvature` from the threshold subsection.

3. **reconstructing-physics.tex**

   * Rename sections to “Quantum sector: envelope and Schrödinger limit” and “Electromagnetic sector: charge from amplitude” (keep label `sec:charge-amplitude`).
   * Reintroduce `eq:PhiG_from_gradX4` as definition of (\Phi_G) from (X^4).
   * Add explicit small-parameter conditions for the envelope approximation and clarify the (\hbar) definition.
   * Move the heavy tubular-electron derivation into `appendix-electron-tube.tex` and leave a summarized version with a reference to `\ref{app:electron-tube}`.

4. **experimental-setting.tex**

   * Remove or quarantine any “saturation potential” or “hard cutoff” nonlinearity that is not purely geometric.
   * Ensure references to `eq:PhiG_from_gradX4` and `\ref{conceptual-model}` are correct.
   * Add a short explicit paragraph on the SI→dimensionless mapping (Compton-based scaling).
   * Compress and clarify the “interpretation layers” subsection.

5. **appendix-curvature.tex**

   * Set the main label to `\label{app:bulge-curvature}` and ensure all references match.
   * Soften the central curvature result from `=` to `\approx` and emphasize scaling.

6. **appendix-electron-tube.tex**

   * Create this file and move the full tubular-electron derivation here, with label `\label{app:electron-tube}`.

7. **discussion.tex / conclusion.tex**

   * Make the entanglement / Bell sector clearly marked as an open problem.
   * Add 1–2 sentences explaining which parameter identifications are robust and which are more conjectural.

If you’d like, the next step I can do is draft the actual LaTeX snippets for some of the heavier edits (e.g. the SI→dimensionless mapping paragraph, the entanglement paragraph, etc.) in final form that you can paste directly.
