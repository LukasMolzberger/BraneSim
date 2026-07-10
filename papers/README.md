# Papers

A five-paper series developing a deterministic 4D brane-lattice substrate hypothesis.
Each paper is self-contained but imports interface equations from earlier papers via `\externaldocument`.

## Structure

| Directory | Paper | Topic |
|-----------|-------|-------|
| `core/` | Paper I | Ontological foundation, lattice action, continuum limit, linearized wave structure, Bell commitment |
| `lorentz_gravity/` | Paper II | Emergent Lorentz symmetry (dual-observer argument) and the gravity channel |
| `gauge_color/` | Paper III | Berry connection as EM gauge potential; Wilczek–Zee connection as color gauge field; kinematic color confinement |
| `matter_mass/` | Paper IV | Soliton eigenproblem (VSH ansatz), Derrick stability, emergent rest mass as self-confined wave loop |
| `bell/` | Paper V | Bell constraint uniqueness argument; retrocausal worldtube interpretation; matter/antimatter as opposite-chirality worldtubes |
| `field_strength/` | Paper VII | Drill-down of Paper III: Faraday tensor `F_μν` and QCD field strength `G^a_μν` as continuum limits of plaquette holonomies on the 8-link spring stencil; equivalence to Berry/Wilczek–Zee curvature |

## Building

Each paper compiles independently from its own directory:

```bash
cd papers/core          # or lorentz_gravity, gauge_color, matter_mass, bell
pdflatex paper.tex
bibtex paper
pdflatex paper.tex
pdflatex paper.tex
```

Cross-paper label resolution (`\externaldocument`) requires the referenced paper's `.aux` file to exist.
Build order: Paper I → Paper III → Paper IV (Paper II and V can build in any order after Paper I).

## Shared files

`Definitions/` — MDPI journal class files (`mdpi.cls`, `mdpi.bst`, etc.), shared by all papers via a symlink in each paper directory.
Each paper's `references.bib` contains only the entries actually cited in that paper.