Here’s a “Claude-code ready” TODO list, grouped by file, with concrete search/replace or edit instructions.

---

## 1. `draft-paper.tex`

### 1.1 Remove the Disclaimer chapter entirely

1. **Remove the inclusion of the disclaimer file**

   Search for:

   ```tex
   \input{disclaimer}
   ```

   and **delete this line**.

2. You can keep `disclaimer.tex` in the repo for your own reference, but it will no longer be part of the compiled paper.

---

## 2. `conceptual-model.tex`

### 2.1 Remove the saturation potential from the Lagrangian

1. **Edit the Lagrangian definition** (around `\label{eq:brane-lagrangian}`):

   Current equation:

   ```tex
   \begin{equation}
   \mathcal{L}
   =\frac{\rho_m}{2}|\partial_t\mathbf{X}|^2
   -\frac{T}{2}\mathrm{tr}(E)
   -\mu|E|^2
   -\frac{\kappa}{2}|b|^2
   -W_\text{sat}(E),
   \label{eq:brane-lagrangian}
   \end{equation}
   ```

   Replace by:

   ```tex
   \begin{equation}
   \mathcal{L}
   =\frac{\rho_m}{2}|\partial_t\mathbf{X}|^2
   -\frac{T}{2}\mathrm{tr}(E)
   -\mu|E|^2
   -\frac{\kappa}{2}|b|^2,
   \label{eq:brane-lagrangian}
   \end{equation}
   ```

2. **Remove the text that says saturation enters via (W_\text{sat})**

   In the paragraph just before the Lagrangian, there is a sentence like:

   > “Nonlinear saturation effects enter via (W_\text{sat}(E)) (defined below).”

   Delete that sentence (or the clause that mentions (W_\text{sat})) so that nonlinearity in this section is attributed **only** to the geometric dependence of (g_{ij}, E_{ij}, b_{ij}) on (X^A).

---

### 2.2 Remove the explicit definition of (W_\text{sat})

1. **Delete the “A convenient convex, smooth saturation is:” block**

   Find the block:

   ```tex
   A convenient convex, smooth saturation is:
   \begin{equation}
   W_\text{sat}(E)=\frac{k\epsilon_\text{cr}^2}{2}
   \ln\Big(1+\frac{\mathrm{tr}(E^2)}{\epsilon_\text{cr}^2}\Big),
   \label{eq:saturation-potential}
   \end{equation}
   ...
   ```

   Delete this entire paragraph and the equation environment, including any text that discusses (k), (\epsilon_\text{cr}) and “stiffness / critical strain scale” **in the context of (W_\text{sat})**.

2. Remove any remaining isolated references to `W_\text{sat}`, `\epsilon_\text{cr}` or “saturation” in this section, unless they are explicitly saying that such terms are **not** used in the present paper.

---

### 2.3 Keep only geometric nonlinearity in “Sources of nonlinearity”

1. In the paragraph:

   ```tex
   \paragraph{Sources of nonlinearity.}
   ```

   you currently have an `enumerate` with two items:

   * Pure geometric nonlinearity
   * Optional material or strong-field nonlinearity (saturation (W_\text{sat}), (\epsilon_\text{cr}), etc.)

2. **Keep the first item and delete the second item.**

   That is, keep:

   ```tex
   \item \emph{Pure geometric nonlinearity.} Even if the elastic law
         is Hookean in the Green strain, the induced metric $g_{ij}$,
         the strain tensor $E_{ij}$ and the curvature tensor $b_{ij}$
         are nonlinear functionals of the embedding $X^A(x^\mu)$.
         Large slopes and large curvatures therefore produce nonlinear
         behavior purely from the geometry of the brane.
         A compact quantitative estimate of when this geometric nonlinearity
         becomes important for photon-like modes is given in
         Sec.~\ref{threshold-localization}, ...
   ```

   and **remove** the item that starts with something like:

   ```tex
   \item \emph{Optional material or strong-field nonlinearity.} On top
         of this geometric effect we may introduce smooth saturation or
         strong-field potentials, schematically written as $W_\text{sat}$,
         ...
   ```

3. At the end of this paragraph, you can optionally add a sentence like:

   ```tex
   In the present paper we restrict attention to case~(1) and do not
   introduce any additional material saturation terms; all nonlinearity
   considered here is purely geometric.
   ```

---

### 2.4 Remove (W_\text{sat}) from the energy density and threshold localization

1. **Edit the energy density equation** (`\label{eq:energy-density}`):

   Current:

   ```tex
   \begin{equation}
   \mathcal{E}=\frac{\rho_m}{2}|\partial_t\mathbf{X}|^2
   +\frac{T}{2}\mathrm{tr}(E)
   +\mu|E|^2
   +\frac{\kappa}{2}|b|^2
   +W_\text{sat}(E).
   \label{eq:energy-density}
   \end{equation}
   ```

   Replace with:

   ```tex
   \begin{equation}
   \mathcal{E}=\frac{\rho_m}{2}|\partial_t\mathbf{X}|^2
   +\frac{T}{2}\mathrm{tr}(E)
   +\mu|E|^2
   +\frac{\kappa}{2}|b|^2.
   \label{eq:energy-density}
   \end{equation}
   ```

2. **Adjust the text right after Eq.~\eqref{eq:threshold-criterion}**

   There is a sentence like:

   > “Above threshold, the saturation (W_\text{sat}) naturally locks in localized high-curvature 4D bulges …”

   Replace this with something that attributes localization purely to geometry, for example:

   ```tex
   Above threshold, the geometric nonlinearity inherent in the full
   strain and curvature tensors can naturally lock in localized
   high-curvature 4D bulges, turning a propagating photon-like mode
   into a self-confined, particle-like soliton.
   ```

   (Adjust wording to your taste, keeping the “purely geometric” emphasis.)

---

## 3. `experimental-setting.tex`

Here we need to remove the discrete saturation potential (\phi(\epsilon)) and keep only linear (Hooke-like) springs.

### 3.1 Replace saturating (\phi(\epsilon)) with a simple quadratic Hooke law

1. In the subsection that describes stretching energy (look for:

   ```tex
   \item \textbf{Stretching energy (saturation built-in):}
   ```

   or similar), you currently have:

   ```tex
   \item \textbf{Stretching energy (saturation built-in):}
     \[
       U_\text{str}=\sum_{(p,q)\in\mathcal E}
       \phi\left(\frac{|\mathbf R_q-\mathbf R_p|-h}{h}\right),
     \]
     with a smooth, convex-near-zero, saturating potential
     $\phi(\epsilon)$ that behaves like $\tfrac12 k\epsilon^2$
     for small strain and saturates for $|\epsilon|\gg\epsilon_\text{cr}$.
     Recommended:
     \[
       \phi(\epsilon)=\tfrac12 k\epsilon_\text{cr}^2\ln\big(1+(\epsilon/\epsilon_\text{cr})^2\big),
     \]
     ...
   ```

2. Replace this entire bullet item by a **purely linear** stretching energy:

   ```tex
   \item \textbf{Stretching energy (Hooke-like):}
     \[
       U_\text{str}=\sum_{(p,q)\in\mathcal E}
       \phi\left(\frac{|\mathbf R_q-\mathbf R_p|-h}{h}\right),
     \]
     with a simple quadratic potential
     \[
       \phi(\epsilon)=\frac{1}{2}k\,\epsilon^2,
     \]
     so that the discrete stretching law reduces locally to a
     Hooke-like response. In the present work we do not include
     any additional material saturation; all nonlinear behavior
     arises from the geometric dependence of distances and curvature
     on the brane embedding.
   ```

3. Delete all references in this file to:

   * `\epsilon_\text{cr}`
   * “saturates for (|\epsilon|\gg\epsilon_\text{cr})”
   * “nonlinearity scale (\epsilon_\text{cr})”
   * “saturation built-in”
   * Any enumeration of “saturation modes” that describes different (\phi) shapes.

   If you want to keep a note about previous experiments using saturation, move that to a short parenthetical remark or footnote and make clear it is **outside the present paper’s focus**.

4. Near the end of the experimental section, you already have sentences like:

   > “No explicit amplitude or strain clamps are applied; any propagation or threshold behavior arises from the elastic dynamics rather than from artificial hard cutoffs.”

   Make sure that, after removing (\phi)-saturation, this statement is still consistent and doesn’t mention saturation as a modeling option.

---

## 4. `reconstructing-physics.tex`

### 4.1 Clean up references to saturation (W_\text{sat})

1. Anywhere in this file where you have phrases of the form:

   * “tension, bending, and saturation”
   * “corrections from (\mu) and (W_\text{sat})”
   * “stretching, bending, and saturation energy functional”

   **remove “and saturation”, “(W_\text{sat})” etc.**

   For example, change:

   ```tex
   ... plus corrections from $\mu$ and $W_\text{sat}$.
   ```

   to:

   ```tex
   ... plus corrections from the shear modulus $\mu$.
   ```

   And:

   ```tex
   tension, bending, and saturation. The corresponding energy functional ...
   ```

   to something like:

   ```tex
   tension, shear, and bending. The corresponding energy functional ...
   ```

2. If this file reuses Eq.~\eqref{eq:energy-density} conceptually, make sure the **text** describes it without mentioning saturation.

---

### 4.2 Poisson equation: explicitly mark as an assumption

In the section where you write the effective Poisson equation (search around the `4\pi G` / Poisson relation, e.g. label `eq:poisson-from-brane`):

1. You already say something like “at the level of this paper we assume that … obeys a Poisson-type relation…”.

2. Make the “assumption” status visually explicit. For example, change the paragraph header to:

   ```tex
   \paragraph{Assumption: sourcing the potential by brane energy.}
   ```

   and in the text before/after the Poisson equation, add a sentence like:

   ```tex
   This Poisson-like relation is an explicit assumption at the
   current stage; a derivation of the exact proportionality from
   the microscopic brane Lagrangian is left for future work.
   ```

This makes it unambiguous for a referee that this is not claimed as a derivation.

---

### 4.3 Electromagnetism: clearly label the Maxwell sector as phenomenological

In the subsection that introduces the four-potential (A_\mu) and the Maxwell wave equation (around `\eqref{eq:maxwell-wave}`):

1. You already have wording like “introduced phenomenologically”.

2. Strengthen this by inserting a short lead-in phrase, e.g. change:

   ```tex
   We now promote the scalar potential $\Phi$ to a four-potential
   $A_\mu$ and postulate that it obeys a Lorentz-covariant wave equation ...
   ```

   to:

   ```tex
   \textbf{Phenomenological step.}
   We now promote the scalar potential $\Phi$ to a four-potential
   $A_\mu$ and postulate that it obeys a Lorentz-covariant wave equation ...
   ```

3. Keep the existing sentences that say (J_\mu) and the detailed form of the current are *introduced phenomenologically* and are not yet derived from the brane microphysics.

---

### 4.4 Remove the analog-gravity examples

1. In the subsection `Analogy: emergent relativity in physical media` (label `\label{analogy-emergent-relativity-in-physical-media}`), you have a paragraph/`itemize` starting with something like:

   ```tex
   \paragraph{Analog gravity models}, such as:
   \begin{itemize}
     \item Bose–Einstein condensates (BECs) ...
     \item Acoustic metrics in fluids ...
     \item Electromagnetic waves in metamaterials ...
   \end{itemize}
   ```

2. Delete this entire `\paragraph{Analog gravity models}` and its `itemize` environment so that the section no longer lists concrete analog-gravity examples.

3. You can keep a short sentence in the main text like:

   ```tex
   Similar emergent-relativity behaviour appears in various
   condensed-matter and optical systems, but here we treat the
   elastic brane as ontologically fundamental rather than as an
   effective medium.
   ```

   (optional, keep or shorten as you like).

---

## 5. `introduction.tex`

### 5.1 Remove analog-gravity *examples* from the intro

1. In the paragraph that cites Sakharov and mentions analog-gravity in condensed-matter systems (search for “condensed-matter systems such as” or similar), you likely have a sentence like:

   ```tex
   ... emergent relativity and gravity can arise in condensed-matter
   systems such as Bose–Einstein condensates, superfluid helium,
   and optical analogs ...
   ```

2. Replace that part with a version that **does not enumerate examples**, e.g.:

   ```tex
   ... emergent relativity and gravity can arise in a variety of
   condensed-matter and optical systems.
   ```

   So the reader still gets the idea that emergent-relativity media exist, but without a long example list.

---

## 6. `discussion.tex`

### 6.1 Explicitly acknowledge entanglement / Bell as an open question

You already mention entanglement and Bell-inequality violations in this file; we’ll make the open-problem status crystal clear.

1. **Bullet on “Particle quantum statistics and entanglement”**

   Find the bullet item that starts with something like:

   ```tex
   \item \textbf{Particle quantum statistics and entanglement}:
   ```

   Replace the entire text of that `\item` with something like:

   ```tex
   \item \textbf{Particle quantum statistics and entanglement}:
     While the present model reproduces wave-like propagation and
     classical uncertainty relations, it does not yet provide a
     concrete mechanism for quantum entanglement or for the
     experimentally observed Bell-inequality violations. At this
     stage we regard the reproduction of quantum statistics and
     nonlocal correlations within the brane framework as an open
     question, and we do not claim to have a satisfying solution.
   ```

2. **Paragraph where “entanglement correlations and Bell inequality violations” are mentioned**

   In the earlier discussion (search for that exact phrase), immediately after the first sentence that mentions entanglement and Bell, add:

   ```tex
   At present we do not have a satisfying explanation of these
   phenomena within the brane model; they remain an open problem.
   ```

   This line can be adjusted to avoid redundancy if it is very close to the bullet above.

---

## 7. Optional global notation clean-ups (for later)

These are not strictly required for your new requests, but they implement some of the earlier critique:

1. **Mass density notation**

   * Search across all `.tex` files for instances where `\rho` is clearly used as “mass density of the brane” and change them to `\rho_m` for consistency (leave things like `\rho_\mathrm{eff}`, `\rho_{\text{free}}` etc. alone).

2. **Remove dead references to saturation parameters**

   * After the edits above, do a project-wide search for:

     * `W_\text{sat}`
     * `\epsilon_\text{cr}`
     * “saturation”
   * Any remaining occurrences that still refer to a **material** saturation potential should either be deleted or rephrased so that the paper consistently says: **only geometric nonlinearity is used in the present work**.

---

If you’d like, next step I can draft specific replacement text for any paragraph you’re unsure about (e.g. the Sakharov/analog paragraph in the intro) so Claude code can just paste it in.
