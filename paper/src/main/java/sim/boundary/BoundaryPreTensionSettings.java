package sim.boundary;

/**
 * Configuration for boundary pre-tension, which creates initial stress
 * in the boundary region by displacing boundary points outward.
 *
 * Pre-tension helps stabilize the boundary by ensuring springs connecting
 * boundary points to interior points start in a stretched state, providing
 * restoring forces that resist disturbances at the edges.
 */
public final class BoundaryPreTensionSettings {

    /**
     * Profile for how displacement increases from the interior towards the boundary.
     */
    public enum Profile {
        /** Linear ramp: displacement ∝ s */
        LINEAR,
        /** Quadratic ramp: displacement ∝ s² */
        QUADRATIC,
        /** Cubic ramp: displacement ∝ s³ (smoother transition) */
        CUBIC
    }

    /** Whether pre-tension is enabled */
    public final boolean enabled;

    /** Thickness of the pre-tension layer in cells */
    public final int thicknessCells;

    /**
     * Maximum geometric offset at the boundary (in meters).
     * Boundary points are displaced outward by up to this amount,
     * creating stretched springs. Typically a small fraction of
     * the stretched length (e.g., 0.1 * stretchedLength).
     */
    public final double maxOffset;

    /** Profile shape for displacement ramp */
    public final Profile profile;

    /**
     * Creates boundary pre-tension settings.
     *
     * @param enabled whether pre-tension is active
     * @param thicknessCells number of boundary layers under pre-tension
     * @param maxOffset maximum outward displacement at boundary [m]
     * @param profile shape of displacement ramp
     */
    public BoundaryPreTensionSettings(boolean enabled,
                                      int thicknessCells,
                                      double maxOffset,
                                      Profile profile) {
        this.enabled = enabled;
        this.thicknessCells = thicknessCells;
        this.maxOffset = maxOffset;
        this.profile = profile;
    }

    /**
     * @return settings with pre-tension disabled
     */
    public static BoundaryPreTensionSettings disabled() {
        return new BoundaryPreTensionSettings(false, 0, 0.0, Profile.LINEAR);
    }
}