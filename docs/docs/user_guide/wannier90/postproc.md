# `wannier90` as a post-processing tool

This is a description of how to use `wannier90` as a post-processing
tool.

The code must be run twice. On the first pass either the logical keyword
`postproc_setup` must be set to `.true.` in the input file
`seedname.win` or the code must be run with the command line option
`-pp`. Running the code then generates the file `seedname.nnkp` which
provides the information required to construct the
$M_{mn}^{(\mathbf{k,b})}$ overlaps (Ref. [@marzari-prb97], Eq. (25)) and
$A_{mn}^{(\mathbf{k})}$ (Ref. [@marzari-prb97], Eq. (62);
Ref. [@souza-prb01], Eq. (22)).

Once the overlaps and projection have been computed and written to files
`seedname.mmn` and `seedname.amn`, respectively, set `postproc_setup` to
`.false.` and run the code. Output is written to the file
`seedname.wout`.

## `seedname.nnkp` file

OUTPUT, if `postproc_setup = .true.`

The file `seedname.nnkp` provides the information needed to determine
the required overlap elements $M_{mn}^{(\mathbf{k,b})}$ and projections
$A_{mn}^{(\mathbf{k})}$. It is written automatically when the code is
invoked with the `-pp` command-line option (or when
`postproc_setup=.true.` in `seedname.win`. There should be no need for
the user to edit this file.

Much of the information in `seedname.nnkp` is arranged in blocks
delimited by the strings `begin block_name` ... `end block_name`, as
described below.

### Keywords

The first line of the file is a user comment, e.g., the date and time:

`File written on 12Feb2006 at 15:13:12`

The only logical keyword is `calc_only_A`, eg,

`calc_only_A  :  F`

### `Real_lattice` block

```vi title="Input file"
begin real_lattice
 2.250000   0.000000   0.000000
 0.000000   2.250000   0.000000
 0.000000   0.000000   2.250000
end real_lattice
```

The real lattice vectors in units of Angstrom.

### `Recip_lattice` block

```vi title="Input file"
begin recip_lattice
 2.792527   0.000000   0.000000
 0.000000   2.792527   0.000000
 0.000000   0.000000   2.792527
end recip_lattice
```

The reciprocal lattice vectors in units of inverse Angstrom.

### `Kpoints` block

```vi title="Input file"
begin kpoints
  8
  0.00000   0.00000   0.00000
  0.00000   0.50000   0.00000
  .
  .
  .
  0.50000   0.50000   0.50000
end kpoints
```

The first line in the block is the total number of k-points `num_kpts`.
The subsequent `num_kpts` lines specify the k-points in crystallographic
co-ordinates relative to the reciprocal lattice vectors.

### `Projections` block

```vi title="Input file"
begin projections
   n_proj
   centre   l  mr  r
     z-axis   x-axis   zona
   centre   l  mr  r
     z-axis   x-axis   zona
   .
   .
end projections
```

Notes:

`n_proj`: integer; the number of projection centres, equal to the number
of MLWF `num_wann`.

`centre`: three real numbers; projection function centre in
crystallographic co-ordinates relative to the direct lattice vectors.

`l  mr  r`: three integers; $l$ and $m_\mathrm{r}$ specify the angular
part $\Theta_{lm_{\mathrm{r}}}(\theta,\varphi)$, and $\mathrm{r}$
specifies the radial part $R_{\mathrm{r}}(r)$ of the projection function
(see Tables [Angular functions](projections.md#angular-functions),
[Hybrids](projections.md#hybrids) and
[Radial functions](projections.md#radial-functions)).

`z-axis`: three real numbers; default is `0.0 0.0 1.0`; defines the axis
from which the polar angle $\theta$ in spherical polar coordinates is
measured.

`x-axis`: three real numbers; must be orthogonal to `z-axis`; default is
`1.0 0.0 0.0` or a vector perpendicular to `z-axis` if `z-axis` is
given; defines the axis from with the azimuthal angle $\varphi$ in
spherical polar coordinates is measured.

`zona`: real number; the value of $\frac{Z}{a}$ associated with the
radial part of the atomic orbital. Units are in reciprocal Angstrom.

### `spinor_projections` block

```vi title="Input file"
begin spinor_projections
   n_proj
   centre   l  mr  r
    z-axis   x-axis   zona
     spin spn_quant
   centre   l  mr  r
    z-axis   x-axis   zona
     spin spn_quant
   .
   .
end spinor_projections
```

Notes: Only one of projections and spinor_projections should be defined.
Variables are the same as the projections block with the addition of
`spin` and `spn_quant`.

`spin`: integer. '1' or '-1' to denote projection onto up or down
states.

`spn_quant`: three real numbers. Defines the spin quantisation axis in
Cartesian coordinates.

### `nnkpts` block

```vi title="Input file"
begin nnkpts
  10
  1   2   0  0  0
  .
  .
end nnkpts
```

First line: `nntot`, the number of nearest neighbours belonging to each
k-point of the Monkhorst-Pack mesh

Subsequent lines: `nntot`$\times$`num_kpts` lines, ie, `nntot` lines of
data for each k-point of the mesh.

Each line of consists of 5 integers. The first is the k-point number
`nkp`. The second to the fifth specify it's nearest neighbours
$\mathbf{k+b}$: the second integer points to the k-point that is the
periodic image of the $\mathbf{k+b}$ that we want; the last three
integers give the G-vector, in reciprocal lattice units, that brings the
k-point specified by the second integer (which is in the first BZ) to
the actual $\mathbf{k+b}$ that we need.

### `exclude_bands` block

```vi title="Input file"
begin exclude_bands 
  8 
  1 
  2 
  .
  .
end exclude_bands
```

To exclude bands (independent of k-point) from the calculation of the
overlap and projection matrices, for example to ignore shallow-core
states. The first line is the number of states to exclude, the following
lines give the states for be excluded.

### `auto_projections` block

```vi title="Input file"
begin auto_projections
   8
   0
end auto_projections
```

This block is only printed if `auto_projections=true` in the input. The
choice of an additional block has been made in order to maintain
back-compatibility with codes that interface with `wannier90`, e.g.
`pw2wannier90`. The first entry in the block (in the example above, `8`)
is the total number of target projections and it is equal to the number
of sought Wannier functions.

The second entry is a reserved flag with the value of zero. The
implementations of the interface codes MUST check for this value to be
zero and stop otherwise. In the future, one possible extension that we
plan is to combine the automatic generation of initial projections with
the selection of projections via a projections block. This will allow
the user to specify only a subset of initial projections in the
projections block and leave the interface code to automatically generate
the remaining ones. In that case the constraint on the second entry will
be lifted, so that it can take on the meaning of the number of
projections that need to be generated automatically.

The selected columns of the density matrix (SCDM)
method [@LinLin-ArXiv2017] is one way of generating the initial
$A_{mn}^{(\mathbf{k})}$ in an automatic way. This has been implemented
in the `pw2wannier90` interface code (you need v6.3 with the files
provided in the `pwscf` folder of Wannier90, or v6.4), see for instance
Tutorial [27](../../tutorials/tutorial_27.md) in the `wannier90` tutorial
that shows how to use it.

Moreover, also the automatic generation of initial projections with
spinor WFs is implemented in the `pw2wannier90` interface. See Tutorial
[31](../../tutorials/tutorial_31.md) in the `wannier90` tutorial that shows
how to use it.

Another automatic projection method is projectability-disentangled
Wannier function (PDWF) [@Qiao2023-pdwf], which uses pseudo-atomic
orbitals inside pseudopotentials as initial guesses.
<!-- TODO: Add a tutorial for PDWF. -->
<!-- See Tutorial [34](../../tutorials/tutorial_34.md) and 35. -->

### An example of projections

As a concrete example: one wishes to have a set of four sp$^3$
projection orbitals on, say, a carbon atom at (0.5,0.5,0.5) in
fractional co-ordinates relative to the direct lattice vectors. In this
case `seedname.win` will contain the following lines:

```vi title="Input file"
begin projections
 C:l=-1
end projections
```

and `seedname.nnkp`, generated on the first pass of `wannier90` (with
`postproc_setup=T`), will contain:

```vi title="Input file"
begin projections
   4
   0.50000    0.50000    0.50000    -1  1  1
     0.000  0.000  1.000   1.000  0.000  0.000   2.00
   0.50000    0.50000    0.50000    -1  2  1
     0.000  0.000  1.000   1.000  0.000  0.000   2.00
   0.50000    0.50000    0.50000    -1  3  1
     0.000  0.000  1.000   1.000  0.000  0.000   2.00
   0.50000    0.50000    0.50000    -1  4  1
     0.000  0.000  1.000   1.000  0.000  0.000   2.00
end projections
```

where the first line tells us that in total four projections are
specified, and the subsequent lines provide the projection centre, the
angular and radial parts of the orbital (see Section
[Orbital Definitions](projections.md#orbital-definitions) for definitions),
the $z$ and $x$ axes,
and the diffusivity and cut-off radius for the projection orbital.

`pwscf`, or any other *ab initio* electronic structure code,
then reads `seedname.nnkp` file, calculates the projections and writes
them to `seedname.amn`.

## `seedname.mmn` file

INPUT.

The file `seedname.mmn` contains the overlaps $M_{mn}^{(\mathbf{k,b})}$.

First line: a user comment, e.g., the date and time

Second line: 3 integers: `num_bands`, `num_kpts`, `nntot`

Then: `num_kpts * nntot` blocks of data:

First line of each block: 5 integers. The first specifies the
$\mathbf{k}$ (i.e., gives the ordinal corresponding to its position in
the list of k-points in `seedname.win`). The 2nd to 5th integers specify
$\mathbf{k+b}$. The 2nd integer, in particular, points to the k-point on
the list that is a periodic image of $\mathbf{k+b}$, and in particular
is the image that is actually mentioned in the list. The last three
integers specify the $\mathbf{G}$ vector, in reciprocal lattice units,
that brings the k-point specified by the second integer, and that thus
lives inside the first BZ zone, to the actual $\mathbf{k+b}$ that we
need.

Subsequent `num_bands * num_bands` lines of each block:
two real numbers per line. These are the real and imaginary parts,
respectively, of the actual scalar product $M_{mn}^{(\mathbf{k,b})}$ for
$m,n \in [1,\texttt{num_bands}]$. The order of these elements is such that
the first index $m$ is fastest.

## `seedname.amn` file

INPUT.

The file `seedname.amn` contains the projection $A_{mn}^{(\mathbf{k})}$.

First line: a user comment, e.g., the date and time

Second line: 3 integers: `num_bands`, `num_kpts`, `num_wann`

Subsequently
$\texttt{num_bands} \times \texttt{num_wann} \times \texttt{num_kpts}$ lines: 3
integers and 2 real numbers on each line. The first two integers are the
band index $m$ and the projection index $n$, respectively. The third
integer specifies the $\mathbf{k}$ by giving the ordinal corresponding
to its position in the list of $k$-points in `seedname.win`. The real
numbers are the real and imaginary parts, respectively, of the actual
$A_{mn}^{(\mathbf{k})}$.

## `seedname.dmn` file

INPUT.

The file `seedname.dmn` contains the data needed to construct
symmetry-adapted Wannier functions [@sakuma-prb13]. Required if
`site_symmetry = .true.`

First line: a user comment, e.g., the date and time

Second line: 4 integers: `num_bands`, `nsymmetry`, `nkptirr`,
`num_kpts`.

- `nsymmetry`: the number of symmetry operations
- `nkptirr`: the number of irreducible k-points

Blank line

`num_kpts` integers: Mapping between full k- and irreducible k-points.
Each k-point is related to some k-point in the irreducible BZ. The
information of this mapping is written. Each entry corresponds to a
k-point in the full BZ, in the order in which they appear in the k-point
list in `seedname.win` file. The (integer) value of each entry is the
k-point index in the IBZ to which the k-point maps. The number of unique
values is equal to the number of k-points in the IBZ. The data is
written 10 values per line.

Blank line

`nkptirr` integers: List of irreducible k-points. Each entry corresponds
to a k-point of the IBZ. The (integer) value of each entry is the
k-point index corresponding to the k-point list in `seedname.win` file.
The values should be between 1 and `num_kpts`. The data is written 10
values per line.

Blank line

`nkptirr` blocks of `nsymmetry` integer data (each block separated by a
blank line): List of k-points obtained by acting the symmetry operations
on the irreducible k-points. The data is written 10 values per line.

Blank line

$\texttt{nsymmetry} \times \texttt{nkptirr}$ blocks of data:
The information of $D$ matrix in Eq. (15) of Ref. [@sakuma-prb13]. Each
block contains $\texttt{num_wann} \times \texttt{num_wann}$ lines and is
separated by a blank line. The data are stored in
`d_matrix_wann(m,n,isym,ikirr)` with
$\texttt{m}, \texttt{n} \in [1,\texttt{num_wann}]$,
$\texttt{isym} \in [1,\texttt{nsymmetry}]$, and
$\texttt{ikirr} \in [1,\texttt{nkptirr}]$. The order of the elements is such
that left indices run faster than right indices (`m`: fastest, `ikirr`:
slowest).

Blank line

$\texttt{nsymmetry} \times \texttt{nkptirr}$ blocks of data:\
The information of $\tilde d$ matrix in Eq. (17) of
Ref. [@sakuma-prb13]. Each block contains
$\texttt{num_bands} \times \texttt{num_bands}$ lines and is separated by a
blank line. The data are stored in `d_matrix_band(m,n,isym,ikirr)` with
$\texttt{m}, \texttt{n} \in [1,\texttt{num_bands}]$,
$\texttt{isym} \in [1,\texttt{nsymmetry}]$, and
$\texttt{ikirr} \in [1,\texttt{nkptirr}]$. The order of the elements is such
that left indices run faster than right indices (`m`: fastest, `ikirr`:
slowest).

## `seedname.eig` file

INPUT.

Required if any of `disentanglement`, `plot_bands`, `plot_fermi_surface`
or `write_hr` are `.true.`

The file `seedname.eig` contains the Kohn-Sham eigenvalues
$\varepsilon_{n\mathbf{k}}$ (in eV) at each point in the Monkhorst-Pack
mesh.

Each line consist of two integers and a real number. The first integer
is the band index, the second integer gives the ordinal corresponding to
the $k$-point in the list of $k$-points in `seedname.win`, and the real
number is the eigenvalue.

E.g.,

```vi title="Input file"
            1           1  -6.43858831271328
            2           1   19.3977795287297
            3           1   19.3977795287297
            4           1   19.3977795287298
```

## Interface with pwscf

Interfaces between `wannier90` and many ab-initio codes such as
`pwscf`, `abinit` (<http://www.abinit.org>), `siesta`
(<http://www.icmab.es/siesta/>), `fleur`, `VASP` and `Wien2k`
(<http://www.wien2k.at>) are available. Here we describe the seamless
interface between `wannier90` and `pwscf`, a plane-wave DFT
code that comes as part of the Quantum ESPRESSO package (see
<http://www.quantum-espresso.org>). You will need to download and
compile `pwscf` (i.e., the `pw.x` code) and the
post-processing interface `pw2wannier90.x`. Please refer to the
documentation that comes with the Quantum ESPRESSO distribution for
instructions.

1. Run 'scf'/'nscf' calculation(s) with `pw`

2. Run `wannier90` with `postproc_setup` = `.true.` to generate
    `seedname.nnkp`

3. Run `pw2wannier90`. First it reads an input file, e.g.,
    `seedname.pw2wan`, which defines `prefix` and `outdir` for the
    underlying 'scf' calculation, as well as the name of the file
    `seedname.nnkp`, and does a consistency check between the direct and
    reciprocal lattice vectors read from `seedname.nnkp` and those
    defined in the files specified by `prefix`. `pw2wannier90` generates
    `seedname.mmn`, `seedname.amn` and `seedname.eig`. `seedname.dmn`
    and `seedname.sym` files are additionally created when
    `write_dmn = .true.` (see below).

4. Run `wannier90` with `postproc_setup` = `.false.` to disentangle
    bands (if required), localise MLWF, and use MLWF for plotting,
    bandstructures, Fermi surfaces etc.

Examples of how the interface with `pwscf` works are given
in the `wannier90` Tutorial.

### `seedname.pw2wan`

A number of keywords may be specified in the `pw2wannier90` input file:

- `outdir` -- Location to write output files. Default is `` `./' ``

- `prefix` -- Prefix for the `pwscf` calculation. Default
    is `` ` ' ``

- `seedname` -- Seedname for the `wannier90` calculation. Default is
    `` `wannier' ``

- `spin_component` -- Spin component. Takes values `` `up' ``,
    `` `down' `` or `` `none' `` (default).

- `wan_mode` -- Either `` `standalone' `` (default) or `` `library' ``

- `write_unk` -- Set to `.true.` to write the periodic part of the
    Bloch functions for plotting in `wannier90`. Default is `.false.`

- `reduce_unk` -- Set to `.true.` to reduce file-size (and resolution)
    of Bloch functions by a factor of 8. Default is `.false.` (only
    relevant if `write_unk=.true.`)

    !!! note
        Note that there is a small bug with this feature in v3.2 (and
        subsequent patches) of `quantum-espresso`. Please use a later
        version (if available) or the CVS version of `pw2wannier90.f90`,
        which has been fixed.

- `wvfn_formatted` -- Set to `.true.` to write formatted
    wavefunctions. Default is `.false.` (only relevant if
    `write_unk=.true.`)

- `write_amn` -- Set to `.false.` if $A_{mn}^{(\mathbf{k})}$ not
    required. Default is `.true.`

- `write_mmn` -- Set to `.false.` if $M_{mn}^{(\mathbf{k,b})}$ not
    required. Default is `.true.`

- `write_spn` -- Set to `.true.` to write out the matrix elements of
    $S$ between Bloch states (non-collinear spin calculation only).
    Default is `.false.`

- `spn_formatted` -- Set to `.true.` to write spn data as a formatted
    file. Default is `.false.` (only relevant if `write_spn=.true.`)

- `write_uHu` -- Set to `.true.` to write out the matrix elements

    $$
    \langle u_{n{\bf k}+{\bf b}_1}\vert
    H_{\bf k}\vert u_{m{\bf k}+{\bf b}_2}\rangle.
    $$

    Default is `.false.`

- `uHu_formatted` -- Set to `.true.` to write uHu data as a formatted
    file. Default is `.false.` (only relevant if `write_uHu=.true.`)

- `write_uIu` -- Set to `.true.` to write out the matrix elements of

    $$
    \langle  u_{n{\bf k}+{\bf b}_1}\vert
    u_{m{\bf k}+{\bf b}_2}\rangle.
    $$

    Default is `.false.`

- `uIu_formatted` -- Set to `.true.` to write uIu data as a formatted
    file. Default is `.false.` (only relevant if `write_uIu=.true.`)

- `write_unkg` -- Set to `.true.` to write the first few Fourier
    components of the periodic parts of the Bloch functions.

- `write_dmn` -- Set to `.true.` to construct symmetry-adapted Wannier
    functions. Default is `.false.`

- `read_sym` -- Set to `.true.` to customize symmetry operations to be
    used in symmetry-adapted mode. When `read_sym = .true.`, an
    additional input `seedname.sym` is required. Default is `.false.`
    (only relevant if `write_dmn=.true.`).

- `irr_bz` -- Set to `.true.` to use the irreducible BZ calculation
    (see [Section "Irreducible BZ calculation"](#irreducible-bz-calculation)).

- `atom_proj` -- Set to `.true.` to use pseudo-atomic orbitals for
    computing `amn`. Default is `.false.`.

- `atom_proj_exclude` -- A list of integers specifying the indices of
    pseudo-atomic projectors to be excluded from computing `amn`. Used
    only when `atom_proj = .true.`. No default.

- `atom_proj_ext` -- Set to `.true.` to use external pseudo-atomic
    orbitals for computing `amn`, will read data files from directory
    `atom_proj_dir`. Used only when `atom_proj = .true.`. Default is
    `.false.`.

- `atom_proj_dir` -- A string specifying the directory for external
    pseudo-atomic projectors. Used only when `atom_proj = .true.` and
    `atom_proj_ext = .true.`. No default.

For examples of use, refer to the `wannier90` Tutorial.

### `seedname.sym`

If `read_sym = .true.`, then this additional input file is required for
`pw2wannier90.x`
if `read_sym = .false.`, then this file is written by `pw2wannier90.x`
(only for reference -- it is not used in subsequent calculations)

The file `seedname.sym` contains the information of symmetry operations
used to create symmetry-adapted Wannier functions. If
`read_sym = .false.` (default), `pw2wannier90.x` uses the full symmetry
recognized by `pw.x`. If `read_sym = .true.`, you can specify symmetry
operations to be used in symmetry-adapted mode.

First line: an integer: `nsymmetry` (number of symmetry operations)

Second line: blank

Then: `nsymmetry` blocks of data. Each block (separated by a blank line)
consists of four lines. The order of the data in each block is as
follows:

```vi title="Input file"
      R(1,1)   R(2,1)   R(3,1)
      R(1,2)   R(2,2)   R(3,2)
      R(1,3)   R(2,3)   R(3,3)
       t(1)     t(2)     t(3)   
```

Here, $R$ is the rotational part of symmetry operations ($3\times3$
matrix), and $\bf t$ is the fractional translation in the unit of
"`alat`" (refer the definition of "`alat`" to the manual of
`pwscf`). Both data are given in Cartesian coordinates. The
symmetry operations act on a point $\bf r$ as ${\bf r} R - {\bf t}$.

## Irreducible BZ calculation

This section explains how to construct Wannier functions using wavefunctions
calculated only in the irreducible BZ (IBZ). Using this option, users only
have to calculate wavefunctions, overlap matrices ($M_{mn}$) and projection
matrices ($A_{mn}$) in the IBZ, which will reduce the computational cost.
Currently, this feature is implemented only in the interface code with
Quantum ESPRESSO and requires an additional python package,
[symWannier](https://github.com/wannier-utils-dev/symWannier). See also
the [website of symWannier](https://github.com/wannier-utils-dev/symWannier)
and Ref. [@Koretsune2023].

!!! note
    The `symWannier` package also has the ability to construct
    symmetry-adapted Wannier functions only from the `.mmn,`.amn` and `.eig`
    files in the IBZ and some symmetry information.

The procedure is as follows.

1. Run `scf`/`nscf` calculation(s) with Quantum ESPRESSO `pw.x`. In the `nscf`
    calculation, users only have to calculate k-points in the IBZ.

    E.g.:

    ```vi title="Input file"
    K_POINTS {automatic}
    4 4 4 0 0 0
    ```

    Users also have the option to skip the `nscf` calculation, if the grid
    of the `scf` step is already sufficient.

2. Run Wannier90 with `postproc_setup = .true.` (or equivalently with
    the `-pp` command line flag) to generate the file `seedname.nnkp`.

3. Run `pw2wannier90` with `irr_bz = .true.` in the input file.
    `pw2wannier90` then generates the files `seedname.immn`, `seedname.iamn`,
    `seedname.ieig` and `seedname.isym`.

4. Run the python code `write_full_data.py`, part of the `symWannier` package.
    You can first install the package from
    [GitHub](https://github.com/wannier-utils-dev/symWannier) and then
    run the script as

    ```bash title="Terminal"
    python write_full_data.py <seedname>
    ```

    (replace `<seedname>` with the seedname of your calculation).

    This python script generates the $M_{mn}$ and $A_{mn}$ matrices, and Kohn-Sham
    eigenvalues in the full BZ (`seedname.mmn`, `seedname.amn`, and
    `seedname.eig`) from those in the IBZ (`seedname.immn`,
    `seedname.iamn`, and `seedname.ieig`) by using the symmetry
    information (`seedname.isym`). The relations between these quantites
    in the IBZ and full BZ are explained in Ref. [@Koretsune2023].

5. Run `wannier90` (with `postproc_setup = .false.` or, equivalently,
    without the `-pp` command line option).

## `seedname.immn` file

INPUT

The format is the same as `seedname.mmn`, except that the number of
kpoints is the number of k-points in the IBZ. See also
[the section on the `seedname.mmn` file](#seednamemmn-file).

## `seedname.iamn` file

INPUT

The format is the same as `seedname.amn`, except that the number of
kpoints is the number of k-points in the IBZ. See also
See also [the section on the  `seedname.amn` file](#seednameamn-file).

## `seedname.ieig` file

INPUT

The format is the same as `seedname.eig`, except that the number of
kpoints is the number of k-points in the IBZ. See also
See also [the section on the  `seedname.eig` file](#seednameeig-file).

## `seedname.isym` file

INPUT

The file `seedname.isym` contains the data needed to construct Wannier
functions from imformation computed only in the IBZ.

First line: a user comment, e.g., the date and time

Second line: 2 integers: `nsymmetry`, `spinors`, where:

- `nsymmetry`: the number of symmetry operations
- `spinors`: `1` for spinor case or `0` for non-spinor case

Then: `nsymmetry` blocks of data. Each block describes a symmetry
operation, $\hat{g}$, and consists of 7 lines (11 lines for the spinor
case). The order of the data in each block is as follows:

Non-spinor case:

```vi title="Input file"
a comment
R(1,1)   R(2,1)   R(3,1)
R(1,2)   R(2,2)   R(3,2)
R(1,3)   R(2,3)   R(3,3)
 t(1)     t(2)     t(3)
 T
 invs
```

Spinor case:

```vi title="Input file"
a comment
R(1,1)   R(2,1)   R(3,1)
R(1,2)   R(2,2)   R(3,2)
R(1,3)   R(2,3)   R(3,3)
 t(1)     t(2)     t(3)
 T
u(1,1).real  u(1,1).imag
u(2,1).real  u(2,1).imag
u(1,2).real  u(1,2).imag
u(2,2).real  u(2,2).imag
 invs
```

Here, $R$ is the rotational part of symmetry operations, ($3 \times 3$
integer matrix), and $\bf t$ is the fractional translation in
crystal coordinates. The symmetry operations act on a point $\bf r$ in
crystal coordinates as $\hat{g} {\bf r} = {\bf r} R - {\bf t}$.
$T=1$ ($T=0$) indicates that this symmetry operation includes (does not
include) time-reversal operation, respectively.
$u$ is the SU(2) rotation matrix for
the spinor part.
`invs` represents the symmetry operation
corresponding to $\hat{g}^{-1}$.

2 blank lines

Then: An integer: `nk_ibz`, the number of k-points in the IBZ.

Then: `nk_ibz` lines, with 3 real numbers on each line corresponding
to a list of k-points in the IBZ in crystal coordinates.

2 blank lines

Then: 2 integers: `num_bands`, `nblocks`

Then: `nblocks` blocks of data for $h_{mn}({\bf k})$.
`nblocks` is a number of symmetry operations, $\hat{h}({\bf k})$, where
$\hat{h} \bf k = \bf k$ and $\bf k \in$ IBZ.

First line of each block: 3 integers. These specify the k-point in the
IBZ (a number from 1 to `nk_ibz`), the symmetry operation,
$\hat{h}({\bf k})$, (a number from 1 to `nsymmetry`) and a number
of non-zero elements of $\hat{h}({\bf k})$.

Subsequent lines of each block: 2 integers and 2 real numbers per line.
These specify band indeces, $m$ and $n$, and real and imaginary parts of
$h_{mn}({\bf k})$. Only non-zero components of $h_{mn}({\bf k})$ are
written.

2 blank lines

Then: An integer: `num_wann`: The number of Wannier functions.

Then: `nsymmetry` blocks of data.

First line of each block: 2 integers. These specify the symmetry
operation, $\hat{g}$, and a number of non-zero elements of the rotation
matrix, $D_{mn}(\hat{g})$.

Subsequent lines of each block: 2 integers and 2 real numbers per line.
These specify indeces of Wannier functions, $m$ and $n$, and real and
imaginary parts of $D_{mn}(\hat{g})$.
