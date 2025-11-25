package sim.boundary;

import sim.BraneConfig;

/**
 * Precomputed damping coefficient field γ(x,y,z) for boundary absorption.
 *
 * The damping field is zero in the interior and increases smoothly towards
 * the boundaries according to the configured profile. This is computed once
 * and reused throughout the simulation for efficiency.
 *
 * The damping is applied to velocities as: v(t+dt) = v(t) * (1 - γ * dt)
 */
public final class BoundaryDampingField {

    private final double[][][] gamma;

    /**
     * Constructs and precomputes the damping field for the given brane configuration.
     *
     * @param cfg brane configuration including damping settings
     */
    public BoundaryDampingField(BraneConfig cfg) {
        int nx = cfg.getN();
        int ny = cfg.getM();
        int nz = cfg.getK();

        BoundaryDampingSettings s = cfg.getBoundaryDamping();
        gamma = new double[nx][ny][nz];

        if (!s.enabled || s.thicknessCells <= 0 || s.gammaMax <= 0.0) {
            // All zeros - no damping
            return;
        }

        int T = s.thicknessCells;
        double gammaMax = s.gammaMax;

        for (int ix = 0; ix < nx; ix++) {
            for (int iy = 0; iy < ny; iy++) {
                for (int iz = 0; iz < nz; iz++) {
                    int d = BoundaryRegionHelper.distanceToAnyBoundary(ix, iy, iz, cfg);

                    if (d >= T) {
                        // Outside damping layer
                        gamma[ix][iy][iz] = 0.0;
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

                    gamma[ix][iy][iz] = gammaMax * ramp;
                }
            }
        }
    }

    /**
     * Gets the damping coefficient at the specified grid point.
     *
     * @param ix x-index
     * @param iy y-index
     * @param iz z-index
     * @return damping coefficient γ [1/s]
     */
    public double getGamma(int ix, int iy, int iz) {
        return gamma[ix][iy][iz];
    }
}