# Generic Band interpolation

By setting $\verb#geninterp#=\verb#TRUE#$, `postw90` will calculate the
band energies (and possibly the band derivatives, if also
`geninterp_alsofirstder` is set to `TRUE`) on a generic list of $k$
points provided by the user.

The list of parameters of the Generic Band Interpolation module are
summarized in
Table [\[parameter_keywords_geninterp\]](#parameter_keywords_geninterp){reference-type="ref"
reference="parameter_keywords_geninterp"}. The list of input $k$ points
for which the band have to be calculated is read from the file named
`seedname_geninterp.kpt`. The format of this file is described below.

### Files

#### `seedname_geninterp.kpt`

INPUT. Read by `postw90` if `geninterp` is `true`.

The first line is a comment (its maximum allowed length is 500
characters).

The second line must contain `crystal` (or `frac`) if the $k$-point
coordinates are given in crystallographic units, i.e., in fractional
units with respect to the primitive reciprocal lattice vectors.
Otherwise, it must contain `cart` (or `abs`) if instead the $k-$point
coordinates are given in absolute coordinates (in units of 1/Å) along
the $k_x$, $k_y$ and $k_z$ axes.

*Note on units*: In the case of absolute coordinates, if $a_{lat}$ is
the lattice constant expressed in angstrom, and you want to represent
for instance the point $X=\frac {2\pi}{a_{lat}} [0.5, 0, 0]$, then you
have to input for its $x$ coordinate $k_x = 0.5 * 2 * \pi / a_{lat}$. As
a practical example, if $a_{lat}=4$Å, then $k_x = 0.78539816339745$ in
absolute coordinates in units of 1/Å.

The third line must contain the number $n$ of following $k$ points.

The following $n$ lines must contain the list of $k$ points in the
format

    kpointidx k1 k2 k3

where `kpointidx` is an integer identifying the given $k$ point, and
`k1`, `k2` and `k3` are the three coordinates of the $k$ points in the
chosen units.

#### `seedname_geninterp.dat` or ` seedname_geninterp_NNNNN.dat` {#sec:seedname.geninterp.dat}

OUTPUT. This file/these files contain the interpolated band energies
(and also the band velocities if the input flag `geninterp_alsofirstder`
is `true`).

If the flag `geninterp_single_file` is `true`, then a single file
`seedname_geninterp.dat` is written by the code at the end of the
calculation. If instead one sets `geninterp_single_file` to `false`,
each process writes its own output file, named
`seedname_geninterp_00000.dat`, ` seedname_geninterp_00001.dat`, ...

This flag is useful when one wants to parallelize the calculation on
many nodes, and it should be used especially for systems with a small
number of Wannier functions, when one wants to compute the bands on a
large number of $k$ points (if the flag `geninterp_single_file` is
`true`, instead, all the I/O is made by the root node, which is a
significant bottleneck).

**Important!** The files are not deleted before the start of a
calculation, but only the relevant files are overwritten. Therefore, if
one first performs a calculation and then a second one with a smaller
number of processors, care is needed to avoid to mix the results of the
older calculations with those of the new one. In case of doubt, either
check the date stamp in the first line of the
` seedname_geninterp_*.dat` files, or simply delete the
` seedname_geninterp_*.dat` files before starting the new calculation.

To join the files, on can simply use the following command:

    cat seedname_geninterp_*.dat > seedname_geninterp.dat

or, if one wants to remove the comment lines:

    rm seedname_geninterp.dat
    for i in seedname_geninterp_*.dat ; do grep -v \# "$i" >> \
    seedname_geninterp.dat ; done

The first few lines of each files are comments (starting with #),
containing a datestamp, the comment line as it is read from the input
file, and a header. The following lines contain the band energies (and
derivatives) for each band and $k$ point (the energy index runs faster
than the $k$-point index). For each of these lines, the first four
columns contain the $k$-point index as provided in the input, and the
$k$ coordinates (always in absolute coordinates, in units of 1/Å). The
fifth column contains the band energy.

If `geninterp_alsofirstder` is `true`, three further columns are
printed, containing the three first derivatives of the bands along the
$k_x$, $k_y$ and $k_z$ directions (in units of eV$\cdot$Å).

The $k$ point coordinates are in units of 1/Å, the band energy is in eV.
