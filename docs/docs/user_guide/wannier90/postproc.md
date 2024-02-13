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

### `seedname.nnkp` file

OUTPUT, if $\verb#postproc_setup#=\verb#.true.#$

The file `seedname.nnkp` provides the information needed to determine
the required overlap elements $M_{mn}^{(\mathbf{k,b})}$ and projections
$A_{mn}^{(\mathbf{k})}$. It is written automatically when the code is
invoked with the `-pp` command-line option (or when
`postproc_setup=.true.` in `seedname.win`. There should be no need for
the user to edit this file.

Much of the information in `seedname.nnkp` is arranged in blocks
delimited by the strings `begin block_name` ... `end block_name`, as
described below.

#### Keywords

The first line of the file is a user comment, e.g., the date and time:

`File written on 12Feb2006 at 15:13:12`

The only logical keyword is `calc_only_A`, eg,

`calc_only_A  :  F`

#### `Real_lattice` block

    begin real_lattice
     2.250000   0.000000   0.000000
     0.000000   2.250000   0.000000
     0.000000   0.000000   2.250000
    end real_lattice

The real lattice vectors in units of Angstrom.

#### `Recip_lattice` block

    begin recip_lattice
     2.792527   0.000000   0.000000
     0.000000   2.792527   0.000000
     0.000000   0.000000   2.792527
    end recip_lattice

The reciprocal lattice vectors in units of inverse Angstrom.

#### `Kpoints` block

    begin kpoints
      8
      0.00000   0.00000   0.00000
      0.00000   0.50000   0.00000
      .
      .
      .
      0.50000   0.50000   0.50000
    end kpoints

The first line in the block is the total number of k-points `num_kpts`.
The subsequent `num_kpts` lines specify the k-points in crystallographic
co-ordinates relative to the reciprocal lattice vectors.

#### `Projections` block

    begin projections
       n_proj
       centre   l  mr  r   
         z-axis   x-axis   zona
       centre   l  mr  r   
         z-axis   x-axis   zona
       .
       .
    end projections

Notes:

`n_proj`: integer; the number of projection centres, equal to the number
of MLWF `num_wann`.

`centre`: three real numbers; projection function centre in
crystallographic co-ordinates relative to the direct lattice vectors.

`l  mr  r`: three integers; $l$ and $m_\mathrm{r}$ specify the angular
part $\Theta_{lm_{\mathrm{r}}}(\theta,\varphi)$, and $\mathrm{r}$
specifies the radial part $R_{\mathrm{r}}(r)$ of the projection function
(see Tables [3.1](#tab:angular){reference-type="ref"
reference="tab:angular"}, [3.2](#tab:hybrids){reference-type="ref"
reference="tab:hybrids"} and [3.3](#tab:radial){reference-type="ref"
reference="tab:radial"}).

`z-axis`: three real numbers; default is `0.0 0.0 1.0`; defines the axis
from which the polar angle $\theta$ in spherical polar coordinates is
measured.

`x-axis`: three real numbers; must be orthogonal to `z-axis`; default is
`1.0 0.0 0.0` or a vector perpendicular to `z-axis` if `z-axis` is
given; defines the axis from with the azimuthal angle $\varphi$ in
spherical polar coordinates is measured.

`zona`: real number; the value of $\frac{Z}{a}$ associated with the
radial part of the atomic orbital. Units are in reciprocal Angstrom.

#### `spinor_projections` block

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

Notes: Only one of projections and spinor_projections should be defined.
Variables are the same as the projections block with the addition of
`spin` and `spn_quant`.

`spin`: integer. '1' or '-1' to denote projection onto up or down
states.

`spn_quant`: three real numbers. Defines the spin quantisation axis in
Cartesian coordinates.

#### `nnkpts` block

    begin nnkpts
      10
      1   2   0  0  0
      .
      .
    end nnkpts

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

#### `exclude_bands` block

    begin exclude_bands 
      8 
      1 
      2 
      .
      .
    end exclude_bands

To exclude bands (independent of k-point) from the calculation of the
overlap and projection matrices, for example to ignore shallow-core
states. The first line is the number of states to exclude, the following
lines give the states for be excluded.

#### []{#sec:auto-projections-block label="sec:auto-projections-block"}`auto_projections` block

    begin auto_projections
       8
       0
    end auto_projections

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
Example 27 in the `wannier90` tutorial that shows how to use it.

Moreover, also the automatic generation of initial projections with
spinor WFs is implemented in the `pw2wannier90` interface. See Example
31 in the `wannier90` tutorial that shows how to use it.

Another automatic projection method is projectability-disentangled
Wannier function (PDWF) [@Qiao2023-pdwf], which uses pseudo-atomic
orbitals inside pseudopotentials as initial guesses. See Example 34 and
35.

#### An example of projections {#sec:proj_example}

As a concrete example: one wishes to have a set of four sp$^3$
projection orbitals on, say, a carbon atom at (0.5,0.5,0.5) in
fractional co-ordinates relative to the direct lattice vectors. In this
case `seedname.win` will contain the following lines:

    begin projections
     C:l=-1
    end projections

and `seedname.nnkp`, generated on the first pass of `wannier90` (with
`postproc_setup=T`), will contain:

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

where the first line tells us that in total four projections are
specified, and the subsquent lines provide the projection centre, the
angular and radial parts of the orbital (see
Section [3.4](#sec:orbital-defs){reference-type="ref"
reference="sec:orbital-defs"} for definitions), the $z$ and $x$ axes,
and the diffusivity and cut-off radius for the projection orbital.

[pwscf]{.smallcaps}, or any other *ab initio* electronic structure code,
then reads `seedname.nnkp` file, calculates the projections and writes
them to `seedname.amn`.

### `seedname.mmn` file

INPUT.

The file `seedname.mmn` contains the overlaps $M_{mn}^{(\mathbf{k,b})}$.

First line: a user comment, e.g., the date and time

Second line: 3 integers: `num_bands`, `num_kpts`, `nntot`

Then: $\verb#num_kpts#\times\verb#nntot#$ blocks of data:

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

Subsequent $\verb#num_bands#\times\verb#num_bands#$ lines of each block:
two real numbers per line. These are the real and imaginary parts,
respectively, of the actual scalar product $M_{mn}^{(\mathbf{k,b})}$ for
$m,n \in [1,\verb#num_bands#]$. The order of these elements is such that
the first index $m$ is fastest.

### `seedname.amn` file

INPUT.

The file `seedname.amn` contains the projection $A_{mn}^{(\mathbf{k})}$.

First line: a user comment, e.g., the date and time

Second line: 3 integers: `num_bands`, `num_kpts`, `num_wann`

Subsequently
$\verb#num_bands#\times\verb#num_wann#\times\verb#num_kpts#$ lines: 3
integers and 2 real numbers on each line. The first two integers are the
band index $m$ and the projection index $n$, respectively. The third
integer specifies the $\mathbf{k}$ by giving the ordinal corresponding
to its position in the list of $k$-points in `seedname.win`. The real
numbers are the real and imaginary parts, respectively, of the actual
$A_{mn}^{(\mathbf{k})}$.

### `seedname.dmn` file

INPUT.

The file `seedname.dmn` contains the data needed to construct
symmetry-adapted Wannier functions [@sakuma-prb13]. Required if
`site_symmetry = .true.`

First line: a user comment, e.g., the date and time

Second line: 4 integers: `num_bands`, `nsymmetry`, `nkptirr`,
`num_kpts`.\
`nsymmetry`: the number of symmetry operations\
`nkptirr`: the number of irreducible k-points

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

$\verb#nsymmetry# \times \verb#nkptirr#$ blocks of data:\
The information of $D$ matrix in Eq. (15) of Ref. [@sakuma-prb13]. Each
block contains $\verb#num_wann# \times \verb#num_wann#$ lines and is
separated by a blank line. The data are stored in
`d_matrix_wann(m,n,isym,ikirr)` with
$\verb#m#, \verb#n# \in [1,\verb#num_wann#]$,
$\verb#isym# \in [1,\verb#nsymmetry#]$, and
$\verb#ikirr# \in [1,\verb#nkptirr#]$. The order of the elements is such
that left indices run faster than right indices (`m`: fastest, `ikirr`:
slowest).

Blank line

$\verb#nsymmetry# \times \verb#nkptirr#$ blocks of data:\
The information of $\tilde d$ matrix in Eq. (17) of
Ref. [@sakuma-prb13]. Each block contains
$\verb#num_bands# \times \verb#num_bands#$ lines and is separated by a
blank line. The data are stored in `d_matrix_band(m,n,isym,ikirr)` with
$\verb#m#, \verb#n# \in [1,\verb#num_bands#]$,
$\verb#isym# \in [1,\verb#nsymmetry#]$, and
$\verb#ikirr# \in [1,\verb#nkptirr#]$. The order of the elements is such
that left indices run faster than right indices (`m`: fastest, `ikirr`:
slowest).

### `seedname.eig` file

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

               1           1  -6.43858831271328
               2           1   19.3977795287297
               3           1   19.3977795287297
               4           1   19.3977795287298

### Interface with pwscf

Interfaces between `wannier90` and many ab-initio codes such as
[pwscf]{.smallcaps}, abinit (<http://www.abinit.org>), siesta
(<http://www.icmab.es/siesta/>), fleur, VASP and Wien2k
(<http://www.wien2k.at>) are available. Here we describe the seamless
interface between `wannier90` and [pwscf]{.smallcaps}, a plane-wave DFT
code that comes as part of the Quantum ESPRESSO package (see
<http://www.quantum-espresso.org>). You will need to download and
compile [pwscf]{.smallcaps} (i.e., the `pw.x` code) and the
post-processing interface `pw2wannier90.x`. Please refer to the
documentation that comes with the Quantum ESPRESSO distribution for
instructions.

1.  Run 'scf'/'nscf' calculation(s) with `pw`

2.  Run `wannier90` with `postproc_setup` = `.true.` to generate
    `seedname.nnkp`

3.  Run `pw2wannier90`. First it reads an input file, e.g.,
    `seedname.pw2wan`, which defines `prefix` and `outdir` for the
    underlying 'scf' calculation, as well as the name of the file
    `seedname.nnkp`, and does a consistency check between the direct and
    reciprocal lattice vectors read from `seedname.nnkp` and those
    defined in the files specified by `prefix`. `pw2wannier90` generates
    `seedname.mmn`, `seedname.amn` and `seedname.eig`. `seedname.dmn`
    and `seedname.sym` files are additionally created when
    `write_dmn = .true.` (see below).

4.  Run `wannier90` with `postproc_setup` = `.false.` to disentangle
    bands (if required), localise MLWF, and use MLWF for plotting,
    bandstructures, Fermi surfaces etc.

Examples of how the interface with [pwscf]{.smallcaps} works are given
in the `wannier90` Tutorial.

#### `seedname.pw2wan`

A number of keywords may be specified in the `pw2wannier90` input file:

-   `outdir` -- Location to write output files. Default is `` `./' ``

-   `prefix` -- Prefix for the [pwscf]{.smallcaps} calculation. Default
    is `` ` ' ``

-   `seedname` -- Seedname for the `wannier90` calculation. Default is
    `` `wannier' ``

-   `spin_component` -- Spin component. Takes values `` `up' ``,
    `` `down' `` or `` `none' `` (default).

-   `wan_mode` -- Either `` `standalone' `` (default) or `` `library' ``

-   `write_unk` -- Set to `.true.` to write the periodic part of the
    Bloch functions for plotting in `wannier90`. Default is `.false.`

-   `reduce_unk` -- Set to `.true.` to reduce file-size (and resolution)
    of Bloch functions by a factor of 8. Default is `.false.` (only
    relevant if `write_unk=.true.`)

    !!! note

        Note that there is a small bug with this feature in v3.2 (and
        subsequent patches) of ` quantum-espresso`. Please use a later
        version (if available) or the CVS version of `pw2wannier90.f90`,
        which has been fixed.

-   `wvfn_formatted` -- Set to `.true.` to write formatted
    wavefunctions. Default is `.false.` (only relevant if
    `write_unk=.true.`)

-   `write_amn` -- Set to `.false.` if $A_{mn}^{(\mathbf{k})}$ not
    required. Default is `.true.`

-   `write_mmn` -- Set to `.false.` if $M_{mn}^{(\mathbf{k,b})}$ not
    required. Default is `.true.`

-   `write_spn` -- Set to `.true.` to write out the matrix elements of
    $S$ between Bloch states (non-collinear spin calculation only).
    Default is `.false.`

-   `spn_formatted` -- Set to `.true.` to write spn data as a formatted
    file. Default is `.false.` (only relevant if `write_spn=.true.`)

-   `write_uHu` -- Set to `.true.` to write out the matrix elements
    $$\langle u_{n{\bf k}+{\bf b}_1}\vert
    H_{\bf k}\vert u_{m{\bf k}+{\bf b}_2}\rangle.$$ Default is `.false.`

-   `uHu_formatted` -- Set to `.true.` to write uHu data as a formatted
    file. Default is `.false.` (only relevant if `write_uHu=.true.`)

-   `write_uIu` -- Set to `.true.` to write out the matrix elements of
    $$\langle  u_{n{\bf k}+{\bf b}_1}\vert
    u_{m{\bf k}+{\bf b}_2}\rangle.$$ Default is `.false.`

-   `uIu_formatted` -- Set to `.true.` to write uIu data as a formatted
    file. Default is `.false.` (only relevant if `write_uIu=.true.`)

-   `write_unkg` -- Set to `.true.` to write the first few Fourier
    components of the periodic parts of the Bloch functions.

-   `write_dmn` -- Set to `.true.` to construct symmetry-adapted Wannier
    functions. Default is `.false.`

-   `read_sym` -- Set to `.true.` to customize symmetry operations to be
    used in symmetry-adapted mode. When `read_sym = .true.`, an
    additional input `seedname.sym` is required. Default is `.false.`
    (only relevant if `write_dmn=.true.`).

-   `atom_proj` -- Set to `.true.` to use pseudo-atomic orbitals for
    computing `amn`. Default is `.false.`.

-   `atom_proj_exclude` -- A list of integers specifying the indices of
    pseudo-atomic projectors to be excluded from computing `amn`. Used
    only when `atom_proj = .true.`. No default.

-   `atom_proj_ext` -- Set to `.true.` to use external pseudo-atomic
    orbitals for computing `amn`, will read data files from directory
    `atom_proj_dir`. Used only when `atom_proj = .true.`. Default is
    `.false.`.

-   `atom_proj_dir` -- A string specifying the directory for external
    pseudo-atomic projectors. Used only when `atom_proj = .true.` and
    `atom_proj_ext = .true.`. No default.

For examples of use, refer to the `wannier90` Tutorial.

#### `seedname.sym`

If `read_sym = .true.`, then this additional input file is required for
`pw2wannier90.x`\
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

      R(1,1)   R(2,1)   R(3,1)
      R(1,2)   R(2,2)   R(3,2)
      R(1,3)   R(2,3)   R(3,3)
       t(1)     t(2)     t(3)   

Here, $R$ is the rotational part of symmetry operations ($3\times3$
matrix), and $\bf t$ is the fractional translation in the unit of
"`alat`" (refer the definition of "`alat`" to the manual of
[pwscf]{.smallcaps}). Both data are given in Cartesian coordinates. The
symmetry operations act on a point $\bf r$ as ${\bf r} R - {\bf t}$.
