Read and thoroughly internalize the file PROJECT_PRINCIPLES.md located at the root of this repository before proceeding with any physics-related code changes or discussions.

**CRITICAL ARCHITECTURAL CONSTRAINTS:**

1. **Substrate-Only Evolution**: The simulation evolves ONLY the microscopic brane degrees of freedom (node positions R_p(t) and velocities). Forces come exclusively from F_p = -∂U/∂R_p where U = U_str + U_bend.

2. **No Back-Reaction**: ALL emergent fields (induced metric g_ij, electromagnetic potential Φ_EM, electric field E, gravitational potential Φ_grav) are PURELY DIAGNOSTIC measurements. They are computed post-hoc from the brane configuration and NEVER fed back into the equations of motion as forces, constraints, or damping terms.

3. **Pure Geometric Coupling**: The default mode is PURE_GEOMETRY. Amplitude-lateral coupling arises automatically from the 4D Euclidean distance |R_q - R_p| in the edge stiffness calculation. No separate coupling field is needed or allowed.

4. **No Artificial Cutoffs**: No hard amplitude clamps, no strain clamps, no piecewise cutoffs. All nonlinearity enters through the smooth elastic energy functional. Thresholds must emerge from dynamics, not be imposed.

5. **Ontological Foundation**: The 3D brane embedded in 4D is ontically real - this is THE fundamental entity. Time t is an external evolution parameter. The 4th coordinate X⁴ represents amplitude deformation. Everything else (gravity, quantum behavior, relativity, EM, particles) is emergent.

**Key Reminders:**
- Particles = topologically stable solitons (e.g., toroidal electron)
- Gravity = lateral contraction from amplitude deformations via geometric coupling
- Relativity = effective symmetry from isotropic wave propagation
- Charge = signed time-averaged amplitude X̄⁴ with chirality
- Mass = localized energy in motion (E/c²)

**Before making any changes to physics implementation, verify:**
- ✓ Are you only modifying substrate dynamics?
- ✓ Are emergent fields staying diagnostic-only?
- ✓ Are you avoiding artificial cutoffs or external forces?
- ✓ Is the geometric coupling mechanism preserved?

Refer to PROJECT_PRINCIPLES.md for complete details on all aspects of the theoretical framework.