# Substrate Bridge — Problem Statement

The discrete 4D brane action must be shown to reduce, in the long-wavelength limit, to an isotropic hyperelastic continuum action on the spacelike slice.

## Specific open sub-problems

**(a) Cauchy relation and isotropy.** Show that the Cauchy relation on the 6-neighbor stencil does NOT automatically produce an isotropic acoustic tensor — the shell weights need retuning. Specifically: the axial-only stencil gives C_1111 − C_1122 − 2C_1212 = (2α−1)k_s/a ≠ 0, a leading-order (zeroth in ka) cubic anisotropy. Cubic isotropy requires 2w_I = w_II + (16/9)w_III on shell weights; the current 1/|δ|² choice does not lie on this curve.

**(b) 4D ontology and 3D slice integration.** Show that the 4D-in-4D ontology and the 3D-slice working theory are rigorously integrated. The saddle character of the Lorentzian action S = Σ_l Δt(T−V) (unbounded below) raises well-posedness questions that are currently deferred: the foundational solver must root-find ∇S = 0, never gradient-descend S; and the two-time (past+future) BVP is not unconditionally well-posed (non-uniqueness at resonant time-extents NΔt).

**(c) Timelike-link parameter verification.** Verify that the timelike-link parameter r_t = α·β·dt gives the prestressed canonical vacuum. The temporal link is a central-force spring with rest length r_t; at r_t = 0 it reduces exactly (not approximately) to the Verlet kinetic term. The open piece is the α_t = α consistency (OPEN_PROBLEMS A4a).