package sim;

/**
 * Defines the dimensionality of a brane simulation.
 * This determines which spatial dimensions are active for wave propagation
 * and boundary effects.
 */
public enum BraneDimensionality {
    /**
     * One-dimensional: Only x varies (ny = nz = 1).
     * Wave propagation along a single axis.
     */
    ONE_D,

    /**
     * Two-dimensional: x and y vary (nz = 1).
     * Wave propagation in a plane.
     */
    TWO_D,

    /**
     * Three-dimensional: x, y, and z all vary.
     * Full 3D wave propagation.
     */
    THREE_D;

    /**
     * @return true if this is a 1D simulation
     */
    public boolean is1D() {
        return this == ONE_D;
    }

    /**
     * @return true if this is a 2D simulation
     */
    public boolean is2D() {
        return this == TWO_D;
    }

    /**
     * @return true if this is a 3D simulation
     */
    public boolean is3D() {
        return this == THREE_D;
    }
}