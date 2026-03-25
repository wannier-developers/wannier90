# Tutorials Using the pwscf Interface

The `pwscf` plane-wave, density-functional theory code,
which is available as part of the
`quantum-espresso` distribution
(<http://www.quantum-espresso.org>), is fully interfaced to
`wannier90` via the `pw2wannier90` post-processing code that is also
available as part of `quantum-espresso`.

Note that both the `pwscf` executable `pw.x` *and* `pw2wannier90.x` can
be run in parallel, which for large calculations can reduce the
computation time very significantly. This requires compiling the code in
its parallel version, using the MPI libraries. Refer to the
`quantum-espresso` package for the documentation on how to
do so.

Moreover we remind here that both the `wannier90` executable and
`postw90.x` can be run in parallel. In this case any number of
processors can be used, independently of the number used for `pw.x` and
`pw2wannier90.x`.
