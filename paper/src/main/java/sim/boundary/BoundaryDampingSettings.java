package sim.boundary;

/**
 * Configuration for boundary damping, which absorbs outgoing waves
 * to prevent reflections at the edges of the simulation domain.
 *
 * This implements an absorbing boundary layer similar to perfectly matched
 * layers (PML) in wave simulations, where damping increases smoothly
 * towards the boundaries.
 */
public final class BoundaryDampingSettings {

    /**
     * Profile for how damping increases from the interior towards the boundary.
     */
    public enum DampingProfile {
        /** Linear ramp: γ ∝ s */
        LINEAR,
        /** Quadratic ramp: γ ∝ s² */
        QUADRATIC,
        /** Cubic ramp: γ ∝ s³ (smoother transition) */
        CUBIC
    }

    /** Whether boundary damping is enabled */
    public final boolean enabled;

    /** Thickness of the damping layer in cells */
    public final int thicknessCells;

    /** Maximum damping coefficient γ_max at the boundary [1/s] */
    public final double gammaMax;

    /** Profile shape for damping ramp */
    public final DampingProfile profile;

    /** Whether to apply damping to amplitude (w) velocity */
    public final boolean applyAmplitude;

    /** Whether to apply damping to lateral (x,y,z) velocities */
    public final boolean applyLateral;

    /**
     * Creates boundary damping settings.
     *
     * @param enabled whether damping is active
     * @param thicknessCells number of boundary layers with damping
     * @param gammaMax maximum damping coefficient [1/s]
     * @param profile shape of damping ramp
     * @param applyAmplitude whether to damp w velocity
     * @param applyLateral whether to damp x,y,z velocities
     */
    public BoundaryDampingSettings(boolean enabled,
                                   int thicknessCells,
                                   double gammaMax,
                                   DampingProfile profile,
                                   boolean applyAmplitude,
                                   boolean applyLateral) {
        this.enabled = enabled;
        this.thicknessCells = thicknessCells;
        this.gammaMax = gammaMax;
        this.profile = profile;
        this.applyAmplitude = applyAmplitude;
        this.applyLateral = applyLateral;
    }

    /**
     * @return settings with damping disabled
     */
    public static BoundaryDampingSettings disabled() {
        return new BoundaryDampingSettings(false, 0, 0.0,
                DampingProfile.LINEAR, true, false);
    }
}