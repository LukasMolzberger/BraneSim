package sim.boundary;

import sim.BraneConfig;

/**
 * Precomputed pre-tension factor field for boundary stress initialization.
 *
 * The factor field is zero in the interior and increases smoothly towards
 * the boundaries according to the configured profile. This factor (in [0,1])
 * is multiplied by maxOffset to determine the actual displacement.
 *
 * The displacement is applied during brane initialization to create
 * pre-stressed boundary regions.
 */
public final class BoundaryPreTensionField {

    private final double[][][] factor;

    /**
     * Constructs and precomputes the pre-tension factor field for the given
     * brane configuration.
     *
     * @param cfg brane configuration including pre-tension settings
     */
    public BoundaryPreTensionField(BraneConfig cfg) {
        int nx = cfg.getN();
        int ny = cfg.getM();
        int nz = cfg.getK();

        BoundaryPreTensionSettings s = cfg.getBoundaryPreTension();
        factor = new double[nx][ny][nz];

        if (!s.enabled || s.thicknessCells <= 0 || s.maxOffset == 0.0) {
            // All zeros - no pre-tension
            return;
        }

        int T = s.thicknessCells;

        for (int ix = 0; ix < nx; ix++) {
            for (int iy = 0; iy < ny; iy++) {
                for (int iz = 0; iz < nz; iz++) {
                    int d = BoundaryRegionHelper.distanceToAnyBoundary(ix, iy, iz, cfg);

                    if (d >= T) {
                        // Outside pre-tension layer
                        factor[ix][iy][iz] = 0.0;
                        continue;
                    }

                    // Normalized distance: 0 at inner edge of layer, 1 at boundary
                    double sNorm = 1.0 - (double) d / (double) T;

                    // Apply profile shape
                    double ramp;
                    switch (s.profile) {
                        case QUADRATIC:
                            ramp = sNorm * sNorm;
                            break;
                        case CUBIC:
                            ramp = sNorm * sNorm * sNorm;
                            break;
                        case LINEAR:
                        default:
                            ramp = sNorm;
                            break;
                    }

                    factor[ix][iy][iz] = ramp; // in [0, 1]
                }
            }
        }
    }

    /**
     * Gets the pre-tension factor at the specified grid point.
     *
     * @param ix x-index
     * @param iy y-index
     * @param iz z-index
     * @return pre-tension factor in [0, 1]
     */
    public double getFactor(int ix, int iy, int iz) {
        return factor[ix][iy][iz];
    }
}