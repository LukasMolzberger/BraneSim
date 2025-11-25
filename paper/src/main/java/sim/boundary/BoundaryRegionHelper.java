package sim.boundary;

import sim.BraneConfig;
import sim.BraneDimensionality;

/**
 * Utility for determining how close a point is to the boundary of the brane,
 * taking into account the dimensionality of the simulation.
 *
 * This helper is used by both pre-tension and damping systems to identify
 * boundary regions in a dimension-aware manner.
 */
public final class BoundaryRegionHelper {

    private BoundaryRegionHelper() {
        // Utility class - no instantiation
    }

    /**
     * Returns the minimum distance (in cells) from the given point to any *active*
     * boundary plane. Active planes depend on the brane's dimensionality:
     * - 1D: only x boundaries are active
     * - 2D: x and y boundaries are active
     * - 3D: x, y, and z boundaries are active
     *
     * @param ix x-index of the point
     * @param iy y-index of the point
     * @param iz z-index of the point
     * @param cfg brane configuration containing dimensions and dimensionality
     * @return minimum distance to any boundary, or Integer.MAX_VALUE if outside boundary layer
     */
    public static int distanceToAnyBoundary(int ix, int iy, int iz, BraneConfig cfg) {
        final int nx = cfg.getN();
        final int ny = cfg.getM();
        final int nz = cfg.getK();
        final BraneDimensionality dim = cfg.getDimensionality();

        int d = Integer.MAX_VALUE;

        // X dimension is always active
        int dx = Math.min(ix, nx - 1 - ix);
        d = Math.min(d, dx);

        // Y dimension is active only in 2D/3D and when ny > 1
        if (!dim.is1D() && ny > 1) {
            int dy = Math.min(iy, ny - 1 - iy);
            d = Math.min(d, dy);
        }

        // Z dimension is active only in 3D and when nz > 1
        if (dim.is3D() && nz > 1) {
            int dz = Math.min(iz, nz - 1 - iz);
            d = Math.min(d, dz);
        }

        return d;
    }
}