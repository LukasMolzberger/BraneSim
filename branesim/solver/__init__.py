"""Solver module: the 4D block BVP solver (JFNK root-find of ‖R‖=0).

``bvp.solve_block`` is the foundational mode (Dirichlet JFNK + chiral Cauchy
fast-path); ``breather`` is the time-periodic eigen-solver; ``worldvolume``
holds the shared WorldVolume data container.  The forward-Verlet IVP march was
removed 2026-06-06 (it raised for the prestressed r_t>0 regime and only served
as a regression baseline; the block solver now warm-starts from the seed).
"""
