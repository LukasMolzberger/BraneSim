import numpy as np

from branesim.core.conventions import ActionParams, LatticeParams
from branesim.core.lattice import SpacelikeLattice
from branesim.diagnostics.alpha_separability import projection_operators
from branesim.diagnostics.run_measurements import device_color_channels
from branesim.initialization.vortex_worldtube import (
    CARRIER_RE_COMPONENTS,
    CARRIER_RE_WEIGHTS,
    VortexParams,
    inject_vortex_worldtube,
    measure_winding_closure,
    vacuum_offsets,
)


def _small_trace_seed():
    n_slices = 8
    lattice = SpacelikeLattice(
        LatticeParams(
            grid_shape=(11, 11, 11),
            spacing=1.0,
            periodic_axes=(True, True, True),
        )
    )
    params = ActionParams(
        k_s=1.0,
        alpha=0.7,
        rho=1.0,
        dt=0.25,
        n_slices=n_slices,
        m_ambient=4,
        r_t=0.7 * 0.25,
        beta=1.0,
    )
    vp = VortexParams(A0=0.2, r0=3.0, w=1.2, l=1, m=1, n_t=1)
    world, meta = inject_vortex_worldtube(lattice, params, vp, n_slices)
    ref = lattice.reference_positions(4)
    return world, ref, lattice, params, meta


def test_vortex_seed_re_component_is_lateral_trace():
    world, ref, lattice, params, meta = _small_trace_seed()
    off_re, _off_im = vacuum_offsets(params.n_slices, params.r_t)
    trace_weights = np.asarray(CARRIER_RE_WEIGHTS)
    P_U1, P_SU3 = projection_operators()

    disp_lat0 = world[0, :, :3] - ref[:, :3] - off_re[0] * trace_weights
    d_u1 = disp_lat0 @ P_U1.T
    d_su3 = disp_lat0 @ P_SU3.T

    assert tuple(meta["carrier_re_components"]) == CARRIER_RE_COMPONENTS
    assert np.sum(d_u1**2) > 0.0
    assert np.sum(d_su3**2) < 1e-24

    winding = measure_winding_closure(world, lattice, slice_index=0, r_t=params.r_t)
    assert abs(winding["winding_through_z_normal"] - 1.0) < 0.1


def test_color_channel_device_accepts_pure_trace_seed(tmp_path):
    world, ref, lattice, params, _meta = _small_trace_seed()

    result = device_color_channels(world, ref, lattice, params, tmp_path)

    assert result["u1_fraction_mean"] > 0.999999
    assert result["su3_fraction_mean"] < 1e-12