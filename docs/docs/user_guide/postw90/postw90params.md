# Parameters

## Introduction

The `wannier90.x` code described in
Part [wannier90.x](../wannier90/methodology.md) calculates the
maximally-localized Wannier functions.

The `postw90.x` executable contains instead a series of modules that
take the Wannier functions calculated by `wannier90.x` and use them to
calculate different properties. This executable is parallel (by means of
MPI libraries), so it can be run on multiple CPUs. The information on
the calculated Wannier functions is read from the checkpoint
`seedname.chk` file. Note that this is written in an unformatted
machine-dependent format. If you need to use this file on a different
machine, or you want to use a version of `postw90.x` compiled with a
different compiler, refer to
Sec. [w90chk2chk](../appendices/utilities.md#w90chk2chkx)
 in the Appendices for a description of how
to export/import this file.

## Usage

`postw90.x` can be run in parallel using MPI libraries to reduce the
computation time.

For serial execution use: `postw90.x [seedname]`

- `seedname`: If a seedname string is given the code will read its
    input from a file `seedname.win`. The default value is `wannier`.
    One can also equivalently provide the string `seedname.win` instead
    of `seedname`.

For parallel execution use: `mpirun -np NUMPROCS postw90.x [seedname]`

- `NUMPROCS`: substitute with the number of processors that you want
    to use.

Note that the `mpirun` command and command-line flags may be different
in your MPI implementation: read your MPI manual or ask your computer
administrator.

Note also that this requires that the `postw90.x` executable has been
compiled in its parallel version (follow the instructions in the file
`README.install` in the main directory of the wannier90 distribution)
and that the MPI libraries and binaries are installed and correctly
configured on your machine.

## `seedname.win` File

The `postw90.x` uses the same `seedname.win` input file of
`wannier90.x`. The input keywords of `postw90.x` must thus be added to
this file, using the same syntax described in
Sec. [wannier90.x](../wannier90/parameters.md#seednamewin-file).

Note that `wannier90.x` checks if the syntax of the input file is
correct, but then ignores the value of the flags that refer only to
modules of `postw90.x`, so one can safely run `wannier90.x` on a file
that contains also `postw90.x` flags.

Similarly, `postw90.x` ignores flags that refer only to `wannier90.x`
(as number of iterations, restart flags, ...). However, some parts of
the input file must be there, as for instance the number of Wannier
functions, etc.

The easiest thing to do is therefore to simply *add* the `postw90` input
keywords to the `seedname.win` file that was used to obtain the Wannier
functions.

## List of available modules
<!-- TODO: add link to tutorials -->

The currently available modules in `postw90.x` are:

- `dos`: Calculation of the density of states (DOS), projected density
    of states (PDOS), net spin etc.

- `kpath`: Calculation of $k$-space quantities such as energy bands,
    Berry curvature and Berry curvature-like term of spin Hall
    conductivity along a piecewise linear path in the BZ (see tutorials
    [17](../../tutorials/tutorial_17.md),
    [18](../../tutorials/tutorial_18.md)
    and [29](../../tutorials/tutorial_29.md) of the tutorial).

- `kslice`: Calculation of $k$-space quantities on a planar slice of
    the BZ (see tutorials [17](../../tutorials/tutorial_17.md),
    [18](../../tutorials/tutorial_18.md)
    and [29](../../tutorials/tutorial_29.md) of the tutorial).

- `berry`: Calculation of properties related to the BZ integral of the
    Berry curvature, Berry connection and Berry curvature-like term of
    spin Hall conductivity, including anomalous Hall conductivity,
    orbital magnetisation, optical conductivity, nonlinear shift current
    and spin Hall conductivity (see
    Chap. [Berry](berry.md) and tutorials
    [18](../../tutorials/tutorial_18.md),
    [19](../../tutorials/tutorial_19.md),
    [25](../../tutorials/tutorial_25.md),
    [29](../../tutorials/tutorial_29.md)
    and [30](../../tutorials/tutorial_30.md) of the
    tutorial). It also includes an option to compute $k\cdot p$
    expansion coefficients (see
    Sec. [kdotp](berry.md#berry_taskkdotp-k-p-coefficients) and tutorial
    [33](../../tutorials/tutorial_33.md) of the tutorial).

- `gyrotropic`: Calculation of gyrotropic properties, including
    natural and current0induced optical rotation, and the
    current-induced magnetization (see
    Chap. [Gyrotropic](gyrotropic.md) and tutorial
    [24](../../tutorials/tutorial_24.md) of the tutorial).

- `BoltzWann`: Calculation of electronic transport properties for bulk
    materials using the semiclassical Boltzmann transport equation (see
    Chap. [BoltzWann](boltzwann.md) and tutorial
    [16](../../tutorials/tutorial_16.md) of the tutorial).

- `geninterp` (Generic Band Interpolation): Calculation band energies
    (and band derivatives) on a generic list of $k$ points (see
    Chap. [Geninterp](geninterp.md)).

## Keyword List

On the next pages the list of available  input keywords is reported. In
particular,
Table [Global Parameters of `postw90`](#global-parameters-of-postw90) reports
keywords that affect the
generic behavior of all modules of . Often, these are "global" variables
that can be overridden by module-specific keywords (as for instance the
`kmesh` flag). The subsequent tables describe the input parameters for
each specific module.
A description of the behaviour of the global flags is described
Sec. [Global variables](#global-variables); the description of the flags
specific to the modules can be found in the following sections.

### Global Parameters of `postw90`

<!-- markdownlint-disable MD013 -->
{{ read_csv('docs/parameters/postw90-global-parameters.csv', colalign=('left', 'center', 'left')) }}
<!-- markdownlint-enable MD013 -->

`seedname.win` file keywords controlling the general behaviour of
  the modules in `postw90`. Argument types are represented by, I for a integer, R
  for a real number, P for a physical value, L for a logical value and S
  for a text string.

- The keyword `spin_moment` does not affect the behavior of the modules
  in `postw90`, and does not really belong to any of them. It is listed here for
  lack of a better place.

### `berry` Parameters

<!-- markdownlint-disable MD013 -->
{{ read_csv('docs/parameters/postw90-berry-parameters.csv', colalign=('left', 'center', 'left')) }}
<!-- markdownlint-enable MD013 -->

`seedname.win` file keywords controlling the `berry` module.
  Argument types are represented by, I for a integer, R for a real
  number, P for a physical value, L for a logical value and S for a text
  string.

### `dos` Parameters

<!-- markdownlint-disable MD013 -->
{{ read_csv('docs/parameters/postw90-dos-parameters.csv', colalign=('left', 'center', 'left')) }}
<!-- markdownlint-enable MD013 -->

`seedname.win` file keywords controlling the `dos` module.
Argument types are represented by, I for a integer, R for a real number,
P for a physical value, L for a logical value and S for a text string.

### `kpath` Parameters

<!-- markdownlint-disable MD013 -->
{{ read_csv('docs/parameters/postw90-kpath-parameters.csv', colalign=('left', 'center', 'left')) }}
<!-- markdownlint-enable MD013 -->

`seedname.win` file keywords controlling the `kpath` module.
  Argument types are represented by, I for a integer, R for a real
  number, P for a physical value, L for a logical value and S for a text
  string.

### `kslice` Parameters

<!-- markdownlint-disable MD013 -->
{{ read_csv('docs/parameters/postw90-kslice-parameters.csv', colalign=('left', 'center', 'left')) }}
<!-- markdownlint-enable MD013 -->

`seedname.win` file keywords controlling the `kslice` module. Argument types
are represented by, I for a integer, R for a real number, P for a physical
value, L for a logical value and S for a text string.

### `gyrotropic` Parameters

<!-- markdownlint-disable MD013 -->
{{ read_csv('docs/parameters/postw90-gyrotropic-parameters.csv', colalign=('left', 'center', 'left')) }}
<!-- markdownlint-enable MD013 -->

`seedname.win` file keywords controlling the `gyrotropic` module.
  Argument types are represented by, I for a integer, R for a real
  number, P for a physical value, L for a logical value and S for a text
  string.

### `BoltzWann` Parameters

<!-- markdownlint-disable MD013 -->
{{ read_csv('docs/parameters/postw90-boltzwann-parameters.csv', colalign=('left', 'center', 'left')) }}
<!-- markdownlint-enable MD013 -->

`seedname.win` file keywords controlling the `BoltzWann` module (calculation of
  the Boltzmann transport coefficients in the Wannier basis). Argument
  types are represented by, I for a integer, R for a real number, P for
  a physical value, L for a logical value and S for a text string.

### `geninterp` Parameters

<!-- markdownlint-disable MD013 -->
{{ read_csv('docs/parameters/postw90-geninterp-parameters.csv', colalign=('left', 'center', 'left')) }}
<!-- markdownlint-enable MD013 -->

`seedname.win` file keywords controlling the Generic Band Interpolation
(`geninterp`) module. Argument types are represented by,   I for a integer,
R for a real number, P for a physical value, L for a logical value and S for
a text string.

## Global variables

### `integer :: kmesh(:)`

Dimensions of the interpolation grid used in `postw90.x`.

*Not to be confused with the `mp_grid` input flag, which instead
specifies the Monkhorst--Pack grid used in the ab-initio calculation!*

If three integers $l$ $m$ $n$ are given, the reciprocal-space cell
subtended by the three primitive translations is sampled on a uniform
$l\times m\times n$ grid (including $\Gamma$). If only one integer $m$
is given, an $m\times m\times m$ grid is used.

If you use a module which needs a k-mesh, either `kmesh_spacing` or
`kmesh` must be defined.

### `real(kind=dp) :: kmesh_spacing`

An alternative way of specifying the interpolation grid. This flag
defines the minimum distance for neighboring $k$ points along each of
the three directions in $k$ space.

The units are Å$^{-1}$.

If you use a module which needs a k-mesh, either `kmesh_spacing` or
`kmesh` must be defined.

### `logical :: adpt_smr`

Determines whether to use an adaptive scheme for broadening the DOS and
similar quantities defined on the energy axis. If `true`, the values for
the smearing widths are controlled by the flag `adpt_smr_fac`.

The default value is `true`.

### `real(kind=dp) :: adpt_smr_fac`

The width $\eta_{n{\bf k}}$ of the broadened delta function used to
determine the contribution to the spectral property (DOS, \...) from
band $n$ at point ${\bf k}$ is calculated as
$\eta_{n{\bf k}}=\alpha\vert \nabla_{\bf k}
\varepsilon_{n{\bf k}}\vert \Delta k,$ where $\varepsilon_{n{\bf k}}$
is the energy eigenvalue and the dimensionless factor $\alpha$ is given
by `adpt_smr_fac`. $\Delta k$ is taken to be the largest of the mesh
spacings along the three reciprocal lattice vectors ${\bf b_1}$, ${\bf
  b_2}$, and ${\bf b_3}$. If the calculated value of $\eta_{n{\bf
    k}}$ exceeds `adpt_smr_max`, the latter value is used.

The default value is $\sqrt{2}$.

### `real(kind=dp) :: adpt_smr_max`

See description given immediately above.

The units are eV. The default value is 1.0.

### `character(len=120) :: smr_type`

Defines the analytical form used for the broadened delta function in the
computation of the DOS and similar quantities defined on the energy
axis.

- `gauss`: Gaussian smearing

- `m-pN`: derivative of the $N$-th order Methfessel-Paxton function
    ($N\geq 0$). Example: `m-p2` for the second-order Methfessel-Paxton
    function. If only `m-p` is provided, the first-order function is
    used, i.e., it is equivalent to `m-p1`.

- `m-v` or `cold`: derivative of the Marzari--Vanderbilt cold-smearing
    function

- `f-d`: derivative of the Fermi-Dirac distribution function

The default value is `gauss`.

### `logical :: smr_fixed_en_width`

Energy width for the smearing function for the DOS. Used only if
`adpt_smr` is `false`.

The units are eV. The default value is 0 eV. Note that if the width is
smaller than twice the energy step (e.g. `dos_energy_step` for the `dos`
module), the DOS will be unsmeared (thus the default is to have an
unsmeared properties when `adpt_smr` is set to `false`.).

### `integer :: num_elec_per_state`

Number of electrons per state. It can only take the values one or two.

The default value is 1 if `spinors=true`, 2 otherwise.

### `real(kind=dp) :: scissors_shift`

Scissors shift applied to the conduction bands.

!!! note
    This variable is deprecated and will be removed in future
    versions of the code. This applies the scissors shift only to the
    Hamiltonian, but also other matrices might need to be updated if a
    scissors shift is applied. If you are using BoltzWann, consider using
    `boltz_bandshift` instead. If you are calculating spin Hall
    conductivity, consider using `shc_bandshift` instead.

The units are eV. The default value is 0 eV (i.e., no scissors shift
applied).

### `integer :: num_valence_bands`

Number of valence bands of the system. Used in different modules and for
the scissors shift.

No default value.

### `logical :: spin_decomp`

If `true`, extra columns are added to some output files (such as
`seedname-dos.dat` for the `dos` module, and analogously for the `berry`
and `BoltzWann` modules).

For the `dos` and `BoltzWann` modules, two further columns are
generated, which contain the decomposition of the required property
(e.g., total or orbital-projected DOS) of a spinor calculation into
up-spin and down-spin parts (relative to the quantization axis defined
by the input variables `spin_axis_polar` and `spin_axis_azimuth`). For
the `berry` module with `berry_task = kubo`, three extra columns are
added to `seedname-jdos.dat`, containing the decomposition of the JDOS
into up $\rightarrow$ up, down $\rightarrow$ down, and spin-flip
transitions. In the same way, six extra columns are added to the data
files `seedname-kubo*.dat` where the complex optical conductivity is
stored.

The file `seedname.spn` must be present at input. Furthermore, if this
variable is set to `true` it requires `num_elec_per_state = 1`.

The default value is `false`.

### `real(kind=dp) :: spin_axis_polar`

Polar angle of the spin quantization axis.

The units are degrees. The default value is 0.

### `real(kind=dp) :: spin_axis_azimuth`

Azimuthal angle of the spin quantization axis.

The units are degrees. The default value is 0.

### `logical :: spin_moment`

Determines whether to evaluate the spin moment.

The default value is `false`.

### `logical :: uHu_formatted`

If `uHu_formatted`=`true`, then the uHu matrix elements will be read
from disk as formatted (ie ASCII) files; otherwise they will be read as
unformatted files.

The default value of this parameter is `false`.

### `logical :: spn_formatted`

If `spn_formatted`=`true`, then the spin matrix elements will be read
from disk as formatted (ie ASCII) files; otherwise they will be read as
unformatted files. Unformatted is generally preferable as the files will
take less disk space and I/O is significantly faster. However such files
will not be transferable between all machine architectures and formatted
files should be used if transferability is required (i.e., for test
cases).

The default value is `false`.

### `character(len=20) :: berry_curv_unit`

Unit in which the Berry curvature is specified at input (in
`berry_curv_adpt_kmesh_thresh`) or written to file (when
`kpath_task=curv` or `kpath_task=shc` or `kslice_task=curv` or
`kslice_task=shc`).

- `ang2`: Angstrom$^2$

- `bohr2`: Bohr$^2$ (atomic units)

The default value is `ang2`.

### `logical :: sc_use_eta_corr`

If `sc_use_eta_corr=true`, apply the finite-eta correction of the shift
current (Eq. (19) of Ref. [@Lihm_shift_eta_2021]). Without the
correction, the convention of the Bloch sum (determined by
`sc_phase_conv`) can affect the computed shift-current spectra. See Ref.
[@Lihm_shift_eta_2021] for details.

The default value is `true`.

## DOS

Note that the behavior of the `dos` module is also influenced by the
value of some global flags (listed in
Table [Global Parameters of postw90](#global-parameters-of-postw90)),
as `spin_decomp`,
`spin_axis_polar`, `spin_axis_azimuth`, `scissors_shift`, etc. Some of
the global flags can be possibly overridden by local flags of the DOS
module, listed below, which have the same name of the global flag but
are prefixed by `dos_`.

### `logical :: dos`

Determines whether to enter the DOS routines.

The default value is `false`.

### `character(len=20) :: dos_task`

The quantity to compute when `dos=true`

The valid options for this parameter are:

- `dos_plot` Density of states. An output data file `seedname-dos.dat`
    is created, containing the energy values in eV in the first column,
    and the total DOS per unit cell and unit energy range (in eV$^{-1}$)
    in the second. Two additional columns are present if
    `spin_decomp=true`

The default value is `dos_plot`.

### `real(kind=dp) :: dos_energy_min`

Lower limit of the energy range for computing the DOS. Units are eV.

The default value is the minimum value of the energy eigenvalues stored
in `seedname.eig`, minus 0.6667.

### `real(kind=dp) :: dos_energy_max`

Upper limit of the energy range for computing the DOS. Units are eV.

If an inner energy window was specified, the default value is the upper
bound of the innter energy window, plus 0.6667. Otherwise it is the
maximum value of the energy eigenvalues stored in `seedname.eig`, plus
0.6667.

### `real(kind=dp) :: dos_energy_step`

Energy step for the grid of energies used to plot the dos. Units are eV.

The default value is 0.01 eV.

### `integer :: dos_project(:)`

If present `postw90` computes, instead of the total DOS, the partial DOS
projected onto the WFs listed. The WFs are numbered according to the
file `seedname.wout`.

For example, to project onto WFs 2, 6, 7, 8, and 12:

`dos_project : 2, 6-8, 12`

The DOS projected onto a set ${\cal S}$ of orbitals is calculated as

$$
\begin{equation}
\begin{aligned}
\rho_{\cal S}(E)&=\frac{1}{N_k}\sum_{\bf k}\sum_n
\langle \psi_{n\bf k}^{({\rm H})}\vert \hat{P}_{\bf k}({\cal S})\vert
\psi_{n\bf k}^{({\rm H})}\rangle\delta(\varepsilon_{n\bf k}-E)\\
\hat{P}_{\bf k}({\cal S})&=\sum_{m\in{\cal S}}
\vert \psi_{n\bf k}^{({\rm W})}\rangle\langle \psi_{n\bf k}^{({\rm W})}\vert,
\end{aligned}
\end{equation}
$$

where $N_k$ is the number of mesh points used to sample the BZ, and the
superscript (H) and (W) refer to *Hamiltonian gauge* and *Wannier
gauge* [@wang-prb06].

### `integer :: dos_kmesh(:)`

Overrides the `kmesh` global variable (see
Sec. [Global variables](#global-variables)).

### `real(kind=dp) :: dos_kmesh_spacing`

Overrides the `kmesh_spacing` global variable (see
Sec. [Global variables](#global-variables)).

### `logical :: dos_adpt_smr`

Overrides the `adpt_smr` global variable (see
Sec. [Global variables](#global-variables)).

### `real(kind=dp) :: dos_adpt_smr_fac`

Overrides the `adpt_smr_fac` global variable (see
Sec. [Global variables](#global-variables)).

### `real(kind=dp) :: dos_adpt_smr_max`

Overrides the `adpt_smr_max` global variable (see
Sec. [Global variables](#global-variables)).

### `logical :: dos_smr_fixed_en_width`

Overrides the `smr_fixed_en_width` global variable (see
Sec. [Global variables](#global-variables)).

Note that if the width is smaller than twice the energy step
`dos_energy_step`, the DOS will be unsmeared (thus the default is to
have an unsmeared DOS).

### `character(len=20) :: dos_smr_type`

Overrides the `smr_type` global variable (see
Sec. [Global variables](#global-variables)).

## kpath

### `logical :: kpath`

Determines whether to enter the kpath routines.

The default value is `false`.

### `character(len=20) :: kpath_task`

The quantities to plot when `kpath=true`

The valid options for this parameter are:

- `bands` Energy bands, in eV. The following files are created:

    - `seedname-bands.dat` (data file)

    - `seedname-bands.gnu` (`gnuplot` script)

    - `seedname-bands.py` (`python` script)

    - `seedname-path.kpt` (list of $k$-points along the path, written
            in the `pwscf` format)

- `curv` Minus the Berry curvature given by
    [Berry Eq. \[anomalous Hall Conductivity\]](berry.md#mjx-eqn:eq:ahc) of
    Ch. [Berry](berry.md), in units of `berry_curv_unit`. The following
    files are created:

    - `seedname-curv.dat` (data file)

    - `seedname-curv_{x,y,z}.gnu` (`gnuplot` scripts)

    - `seedname-curv_{x,y,z}.py` (`python` scripts)

- `morb` The integrand of the $k$-space orbital magnetization formula
    \[[Berry Eq. \[orbital magnetization\]](berry.md#mjx-eqn:eq:morb) of
    Ch. [Berry](berry.md)\] in eV$\cdot$Å$^2$. Four output files are
    created:

    - `seedname-morb.dat` (data file)

    - `seedname-morb_{x,y,z}.gnu` (`gnuplot` scripts)

    - `seedname-morb_{x,y,z}.py` (`python` scripts)

- `shc` The band-projected Berry curvature-like term of spin Hall
    conductivity given by
    [Berry Eq. \[spin Hall conductivity\]](berry.md#mjx-eqn:eq:kubo_shc_berry) of
    Ch. [Berry](berry.md), in units of `berry_curv_unit`. The following
    files are created:

    - `seedname-shc.dat` (data file)

    - `seedname-shc.gnu` (`gnuplot` scripts)

    - `seedname-shc.py` (`python` scripts)

- Any combination of the above. The following combinations are of
    special interest

    - `kpath_task = bands+curv`

    - `kpath_task = bands+morb`

    - `kpath_task = bands+shc`

    They generate the following files:

    - `seedname-bands.dat` (data file)

    - `seedname-{curv,morb,shc}.dat` (data file)

    - `seedname-bands+{curv,morb}_{x,y,z}.py` or
            `seedname-bands+shc.py` (`python` scripts)

    Two-panel figures are produced, with the energy bands within $\pm
    0.65$ eV of the Fermi level in the top panel, and the Berry
    curvature (or $k$-space orbital magnetization, or $k$-resolved Berry
    curvature-like term of spin Hall conductivity) in the bottom panel.

The default value is `bands`.

### `integer :: kpath_num_points`

If `kpath`=`true`, then the number of points along the first
section of the bandstructure plot given by `kpoint_path`. Other sections
will have the same density of $k$-points.

The default value is 100.

### `character(len=20) :: kpath_bands_colour`

When `kpath_task=bands`, colour code the energy bands according to the
specified quantity.

The valid options for this parameter are:

- `spin` Spin projection (in units of $\hbar/2$) along the
    quantization axis defined by the variables `spin_axis_polar` and
    `spin_axis_azimuth`, for a spinor calculation

- `shc` Band-projected Berry curvature-like term of spin Hall
    conductivity (in units of
    `berry_curv_unit`) defined by the variables `shc_alpha`, `shc_beta`
    and `shc_gamma`, for a spinor calculation

- `none` no colour coding

The default value is `none`.

## kslice

### `logical :: kslice`

Determines whether to enter the kslice routines.

The default value is `false`.

### `character(len=20) :: kslice_task`

The quantity to plot when `kslice=true`

The valid options for this parameter are:

- `fermi_lines` Lines of intersection between constant-energy surfaces
    and the slice. The energy level is specified by the keyword
    `fermi_energy`. Output files:

    - `seedname-kslice-fermi-spn.dat` (data file when
            `kslice_fermi_lines_colour = spin`)

    - `seedname-bnd_n.dat` (`gnuplot` data files when
            `kslice_fermi_lines_colour = none`)

    - `seedname-kslice-coord.dat` (`python` data files when
            `kslice_fermi_lines_colour = none`)

    - `seedname-kslice-bands.dat` (`python` data file when
            `kslice_fermi_lines_colour = none`)

    - `seedname-kslice-fermi_lines.gnu` (`gnuplot` script)

    - `seedname-kslice-fermi_lines.py` (`python` script)

- `curv`\[+`fermi_lines`\] Heatmap of the Berry curvature of the
    occupied states \[together with the constant-energy contours\]. The
    unit of Berry curvature is `berry_curv_unit`.

    Output files:

    - `seedname-kslice-coord.dat` (data files)

    - `seedname-kslice-curv.dat` (data file)

    - `[seedname-kslice-bands.dat]` (data file)

    - `seedname-kslice-curv_{x,y,z}[+fermi_lines].py` (`python`
            scripts)

- `morb`\[+`fermi_lines`\] Heatmap of the $k$-space orbital
    magnetization in eV$\cdot$Å$^2$ \[together with the constant-energy
    contours\]. Output files:

    - `seedname-kslice-coord.dat` (data files)

    - `seedname-kslice-morb.dat` (data file)

    - `[seedname-kslice-bands.dat]` (data file)

    - `seedname-kslice-morb_{x,y,z}[+fermi_lines].py` (`python`
            scripts)

- `shc`\[+`fermi_lines`\] Heatmap of the Berry curvature-like term of
    the occupied states \[together with the constant-energy contours\].
    The unit of Berry curvature-like term is `berry_curv_unit`.

    Output files:

    - `seedname-kslice-coord.dat` (data files)

    - `seedname-kslice-shc.dat` (data file)

    - `[seedname-kslice-bands.dat]` (data file)

    - `seedname-kslice-shc[+fermi_lines].py` (`python` scripts)

The default value is `fermi_lines`.

!!! note
    When `kslice_fermi_lines_colour = none` the `gnuplot` scripts draw
    the $k$-slices with a square shape, even when `kslice_b1` and
    `kslice_b2` below are not at right angles, or do not have equal lengths.
    (The `python` scripts draw the slices with the correct parallelogram
    shape.)

### `real(kind=dp) :: kslice_corner(3)`

Reduced coordinates of the lower-left corner of the slice in k-space.

The default value is $(0.0,0.0,0.0)$

### `real(kind=dp) :: kslice_b1(3)`

Reduced coordinates of the first reciprocal-space vector defining the
slice.

The default value is $(1.0,0.0,0.0)$.

### `real(kind=dp) :: kslice_b2(3)`

Reduced coordinates of the second reciprocal-space vector defining the
slice.

The default value is $(0.0,1.0,0.0)$.

### `integer :: kslice_2dkmesh(2)`

Dimensions of the $k$-point grid covering the slice. If two integers $m$
$n$ are given, the slice is sampled on a uniform $m\times n$ grid. If
only one integer $m$ is given, an $m\times m$ grid is used.

The default value for `kslice_kmesh` is 50.

### `character(len=20) :: kslice_fermi_lines_colour`

When `kslice_task=fermi_lines` (but not when combined with `curv` or
`morb`), colour code the Fermi lines according to the specified
quantity.

The valid options for this parameter are:

- `spin` Spin projection (in units of $\hbar/2$) along the
    quantization axis defined by the variables `spin_axis_polar` and
    `spin_axis_azimuth`, for a spinor calculation

- `none` no colour coding

The default value is `none`.

## berry

### `logical :: berry`

Determines whether to enter the berry routines.

The default value is `false`.

### `character(len=120) :: berry_task`

The quantity to compute when `berry=true`

The valid options for this parameter are:

- `kubo` Complex optical conductivity and joint density of states.

    Output files:

    - `seedname-kubo-S_{xx,yy,zz,xy,xz,yz}.dat` (data files). First
        column: optical frequency $\hbar\omega$ in eV. Second and third
        columns: real and imaginary parts of the symmetric conductivity
        $\sigma^{\rm
            S}_{\alpha\beta}(\hbar\omega)=\sigma^{\rm
            S}_{\beta\alpha}(\hbar\omega)$ in S/cm. Six additional
        columns are present if `spin_decomp = true`.

    - `seedname-kubo-A_{yz,zx,xy}.dat` (data files). First column:
        optical frequency $\hbar\omega$ in eV. Second and third columns:
        real and imaginary parts of the antisymmetric conductivity
        $\sigma^{\rm A}_{\alpha\beta}(\hbar\omega)=-\sigma^{\rm
            A}_{\beta\alpha}(\hbar\omega)$ in S/cm. Six additional
        columns are present if `spin_decomp = true`.

    - `seedname-jdos.dat` (data file). First column: energy difference
        $\hbar\omega$ in eV between conduction ($c$) and valence ($v$)
        states with the same crystal momentum ${\bf
            k}$. Second column: joint density of states
        $\rho_{cv}(\hbar\omega)$ (number of states per unit cell per
        unit energy range, in eV$^{-1}$). Three additional columns are
        present if `spin_decomp = true`.

- `ahc` Anomalous Hall conductivity, in S/cm. The three independent
    components $\sigma_x=\sigma_{yz}$, $\sigma_y=\sigma_{zx}$, and
    $\sigma_z=\sigma_{xy}$ are computed.

    Output files:

    - `seedname-ahc-fermiscan.dat` (data file). The first column
        contains the Fermi level $\varepsilon_F$ in eV, and the
        following three column the values of
        $\sigma_{x,y,z}(\varepsilon_F)$. This file is written if a range
        of Fermi energies is specified via `fermi_energy_min` and
        `fermi_energy_max`. If a single Fermi energy is given, the AHC
        is printed in `seedname.wpout` only.

- `morb` Orbital magnetisation, in bohr magnetons per cell.

    Output files:

    - `seedname-morb-fermiscan.dat` (data file). The first column
        contains the Fermi level $\varepsilon_F$ in eV, and the
        following three column the values of $M^{\rm
            orb}_{x,y,z}(\varepsilon_F)$. This file is written if a
        range of Fermi energies is specified via `fermi_energy_min` and
        `fermi_energy_max`. If a single Fermi energy is given, ${\bf
            M}^{\rm orb}$ is printed in `seedname.wpout` only.

- `sc` Nonlinear shift current.

    Output files:

    - `seedname-sc_{xxx,xxy,xxz,...}.dat` (data files). The shift
        current is described by a $3\times3\times3$ tensor
        $\sigma^{abc}$. The program outputs a set of 18 files, and the 9
        remaining components can be obtained by taking into account that
        $\sigma^{abc}$ is symmetric under $b\leftrightarrow c$ index
        exchange. First column: optical frequency $\hbar\omega$ in eV.
        Second column: nonlinear shift current
        $\sigma^{abc}(\hbar\omega)$ in A/V$^{2}$.

- `shc` Spin Hall conductivity (SHC), in $(\hbar/e)$S/cm.

    Output files:

    - `seedname-shc-fermiscan.dat` (data file). The first column is
        the number of entries in the list, the second column contains
        the Fermi level $\varepsilon_F$ in eV, and the last column
        contains the values of
        $\sigma_{\alpha\beta}^{\text{spin}\gamma}(\varepsilon_F)$. This
        file is written if a range of Fermi energies is specified via
        `fermi_energy_min` and `fermi_energy_max`. If a single Fermi
        energy is given, the file will contain SHC at this specific
        energy.

    - `seedname-shc-freqscan.dat` (data file). The first column is the
        number of the entry in the list, the second column contains the
        frequency $\hbar\omega$ in eV, and the following two columns
        contain the values of the real part
        $\Re[\sigma_{\alpha\beta}^{\text{spin}\gamma}(\omega)]$ and
        imaginary part
        $\Im[\sigma_{\alpha\beta}^{\text{spin}\gamma}(\omega)]$ of ac
        SHC. This file is written if a range of frequencies is specified
        via `kubo_freq_min` and `kubo_freq_max`.

- `kdotp` $k\cdot p$ expansion coefficients.

    Output files:

    - `seedname-kdotp_{0,1,2}.dat` (data files); respectively, the
        zeroth, first and second order $k\cdot p$ expansion
        coefficients, in units of eV, eV$\cdot$Å, and eV$\cdot$Å$^{2}$.

There is no default value.

### `integer :: berry_kmesh(:)`

Overrides the `kmesh` global variable (see
Sec. [Global variables](#global-variables)).

### `real(kind=dp) :: berry_kmesh_spacing`

Overrides the `kmesh_spacing` global variable (see
Sec. [Global variables](#global-variables)).

### `integer :: berry_curv_adpt_kmesh`

If a positive integer $n$ is given and `berry_task=ahc`\[or
`berry_task=shc`\], an $n\times n\times n$ mesh is placed around points
on the uniform mesh (defined by either `berry_kmesh` or
`berry_kmesh_spacing`) where the magnitude of the $k$-space Berry
curvature\[$k$-space Berry curvature-like term of SHC\] exceeds the
threshold value specified in `berry_curv_adpt_kmesh_thresh`. This can
be used to densify the BZ integration mesh around spikes in the Berry
curvature\[Berry curvature-like term of SHC\].

The default value is 1.

### `real(kind=dp) :: berry_curv_adpt_kmesh_thresh`

Magnitude of the Berry curvature\[Berry curvature-like term of SHC\] (in
units of `berry_curv_unit`) that triggers adaptive mesh refinement when
`berry_task=ahc`\[`berry_task=shc`\].

The default value is 100.0.

### `real(kind=dp) :: kubo_freq_min`

Lower limit of the frequency range for computing the optical
conductivity, JDOS and ac SHC. Units are eV.

The default value 0.0.

### `real(kind=dp) :: kubo_freq_max`

Upper limit of the frequency range for computing the optical
conductivity, JDOS and ac SHC. Units are eV.

If an inner energy window was specified, the default value is
`dis_froz_max`-`fermi_energy`+0.6667. Otherwise it is the difference
between the maximum and the minimum energy eigenvalue stored in
`seedname.eig`, plus 0.6667.

### `real(kind=dp) :: kubo_freq_step`

Difference between consecutive values of the optical frequency between
`kubo_freq_min` and `kubo_freq_max`. Units are eV.

The default value is 0.01.

### `real(kind=dp) :: kubo_eigval_max`

Maximum energy eigenvalue of the eigenstates to be included in the
evaluation of the optical conductivity, JDOS and ac SHC. Units are eV.

If an inner energy window was specified, the default value is the upper
bound of the inner energy window plus 0.6667. Otherwise it is the
maximum energy eigenvalue stored in `seedname.eig` plus 0.6667.

### `logical :: kubo_adpt_smr`

Overrides the `adpt_smr` global variable (see
Sec. [Global variables](#global-variables)).

### `real(kind=dp) :: kubo_adpt_smr_fac`

Overrides the `adpt_smr_fac` global variable (see
Sec. [Global variables](#global-variables)).

### `real(kind=dp) :: kubo_adpt_smr_max`

Overrides the `adpt_smr_max` global variable (see
Sec. [Global variables](#global-variables)).

### `logical :: kubo_smr_fixed_en_width`

Overrides the `smr_fixed_en_width` global variable (see
Sec. [Global variables](#global-variables)).

### `character(len=120) :: kubo_smr_type`

Overrides the `smr_type` global variable (see
Sec. [Global variables](#global-variables)).

### `logical :: shc_freq_scan`

Determines whether to calculate the frequency scan (i.e. the ac SHC) or
the Fermi energy scan (i.e. the dc SHC) of the spin Hall conductivity.

The default value is `false`, which means dc SHC is calculated.

### `character(len=120) :: shc_method`

If it is `qiao/ryoo`, the ac or dc SHC is calculated using Junfeng
Qiao's/Jihoon Ryoo's method. To calculate the Kubo formula, the spin
current matrix elements are required, and the two methods differ in the
degree of approximation. For details, see the section
[shc](berry.md#sec:shc).

### `integer :: shc_alpha`

The $\alpha$ index of spin Hall conductivity
$\sigma_{\alpha\beta}^{\text{spin}\gamma}$, i.e. the direction of spin
current. Possible values are `1`, `2` and `3`, representing the `x`, `y`
and `z` directions respectively.

The default value is `1`.

### `integer :: shc_beta`

The $\beta$ index of spin Hall conductivity
$\sigma_{\alpha\beta}^{\text{spin}\gamma}$, i.e. the direction of
applied electric field. Possible values are `1`, `2` and `3`,
representing the `x`, `y` and `z` directions respectively.

The default value is `2`.

### `integer :: shc_gamma`

The $\gamma$ index of spin Hall conductivity
$\sigma_{\alpha\beta}^{\text{spin}\gamma}$, i.e. the spin direction of
spin current. Possible values are `1`, `2` and `3`, representing the
`x`, `y` and `z` directions respectively.

The default value is `3`.

If all the `shc_alpha`, `shc_beta` and `shc_gamma` are set as default
values, the $\sigma_{xy}^{\text{spin}z}$ is computed.

### `logical :: shc_bandshift`

Shift all conduction bands by a given amount (defined by
`shc_bandshift_energyshift`).

!!! note
    this flag slightly differs from the global `scissors_shift` flag:
    with `shc_bandshift`, an exact rigid shift is applied *after*
    interpolation; `scissors_shift` applies instead the shift *before*
    interpolation. As a consequence, results may slightly differ (and this
    is why we provide both possibilities). Note also that with
    `scissors_shift` you have to provide the number of valence bands
    `num_valence_bands`, while with `shc_bandshift` you should provide the
    first band to shift `shc_bandshift_firstband` = `num_valence_bands`$+1$.

The default value is `false`.

### `integer :: shc_bandshift_firstband`

Index of the first conduction band to shift.

That means that all bands with index
$i\ge {\tt shc\_bandshift\_firstband}$ will be shifted by
`shc_bandshift_energyshift`, if `shc_bandshift` is `true`.

The units are eV. No default value; if `shc_bandshift` is `true`, this
flag must be provided.

### `real(kind=dp) :: shc_bandshift_energyshift`

Energy shift of the conduction bands.

The units are eV. No default value; if `shc_bandshift` is `true`, this
flag must be provided.

### `real(kind=dp) :: sc_eta`

The width $\eta$ used to broaden energy differences in denominators of
the form

$$
\begin{equation}
\frac{1}{\varepsilon_{n\bf{k}}-\varepsilon_{m\bf{k}}}\rightarrow
\text{Re}\frac{1}{\varepsilon_{n\bf{k}}-\varepsilon_{m\bf{k}}+i\eta}.
\end{equation}
$$

The above is needed in shift-current calculations in order to avoid
numerical problems caused by near-degeneracies in the sum over virtual
states.

The units are eV. The default value is 0.4.

### `integer :: sc_phase_conv`

Convention for the expansion of the Bloch states in shift-current
calculations. It can only take the values one or two. We follow the
convention of Ref. [@pythtb]:

- 1: Include Wannier centre
    ${\bm \tau}_{n}=\langle w_{n{\bf 0}}|{\bf r}| w_{n{\bf 0}} \rangle$
    in the phase factor (so-called tight-binding convention):

    $$
    \begin{equation}
    |u_{n\bf{k}}\rangle = \sum_{\bf{R}} e^{-i{\bf k}({\bf r}-{\bf R}
        -{\bm \tau}_{n})}| w_{n\bf{R}} \rangle
    \end{equation}
    $$

- 2: Do not include Wannier centre in the phase factor (usual
    `Wannier90` convention):

    $$
    \begin{equation}
    |u_{n\bf{k}}\rangle = \sum_{\bf{R}} e^{-i\bf{k}(\bf{r}-\bf{R})}|
        w_{n\bf{R}} \rangle
    \end{equation}
    $$

If `sc_use_eta_corr=true`, the convention does not affect the full
shift-current matrix element, but it does affect the weights of the
internal components that compose it (see Ref.
[@ibanez-azpiroz_ab_2018]). If `sc_use_eta_corr=false`, the convention
can affect the full shift-current matrix element (see Ref.
[@Lihm_shift_eta_2021]).

The default value is 1.

### `real(kind=dp) :: sc_w_thr`

Parameter $\alpha_{t}$ for speeding up the frequency integration in
shift-current calculations. It settles the frequency threshold
$\omega_{t}=\alpha_{t}\eta_{n{\bf k}}$ (a factor times the broadening)
beyond which the delta functions are taken as zero.

The default value is 5.0.

### `real(kind=dp) :: kdotp_kpoint(3)`

Defines the $k$ point around which the $k\cdot p$ expansion is
performed.

The default value is 0.0 0.0 0.0 ($\Gamma$).

### `integer :: kdotp_num_bands`

Number of bands forming the $k\cdot p$ basis set.

No default value.

### `integer :: kdotp_bands(kdotp_num_bands)`

Band indexes of bands belonging to $k\cdot p$ basis. Number of entries
must be equal to the integer defined in `kdotp_num_bands`. The band
labelling follows that of "wannierised" bands.

No default value.

## Gyrotropic

### `logical :: gyrotropic`

Determines whether to enter the gyrotropic routines.

The default value is `false`.

### `character(len=120) :: gyrotropic_task`

The quantity to compute when `gyrotropic=true`

May contain one or more of the following valid options (note that each
option starts with a '-'):

- `-D0` The Berry-curvature dipole tensor (dimensionless)

    Output file: `seedname-gyrotropic-D.dat` ( see Sec.
    [output data format](#output-data-format) for file format description)

- `-Dw` The finite-frequency Berry-curvature dipole tensor
    (dimensionless)

    Output file: `seedname-gyrotropic-tildeD.dat` ( see Sec.
    [output data format](#output-data-format) for file format description)

- `-C` The ohmic conductivity tensor (Ampere/cm)

    Output file: `seedname-gyrotropic-C.dat` ( see Sec.
    [output data format](#output-data-format) for file format description)

- `-K` The orbital contribution to the kME tensor (Ampere)

    Output file: `seedname-gyrotropic-K_orb.dat` ( see Sec.
    [output data format](#output-data-format) for file format description)

    - `-spin` : if this task is present, compute also the spin
        contribution.

        Output file: `seedname-gyrotropic-K_spin.dat`

- `-NOA` The orbital contribution to the NOA (Å)

    Output file: `seedname-gyrotropic-NOA_orb.dat` ( see Sec.
    [output data format](#output-data-format) for file format description)

    - `-spin` : if this task is present, compute also the spin
        contribution.

        Output file: `seedname-gyrotropic-NOA_spin.dat`

- `-dos` the density of states

    Output file:
    `seedname-gyrotropic-DOS.dat`. First column - energy (eV), second
    column - DOS ($1/(\mathrm{eV}\times\mathring{A}^3)$)

There is no default value.

### output data format

The calculated tensors are written as functions of Fermi level $E_F$
(first column) and frequency $\omega$ (second column). If the tensor
does not denend on $\omega$, the second column is filled by zeros. Data
is grouped in blocks of the same $\omega$ separated by two blank lines.
In case of natural optical activity the columns 3 to 11 contain the
independent components of $\gamma_{abc}$ (antisymmetric in $ab$): $yzx$,
$zxy$ ,$xyz$, $yzy$, $yzz$, $zxz$, $xyy$, $yzz$ and $zxx$. For tensors
$C_{ab}$, $D_{ab}$, $\widetilde D_{ab}$, $K_{ab}$ the symmetric and
antisymmetric components are writted. Thus, the columns 3 to 11 are
marked as $xx$, $yy$, $zz$, $xy$, $xz$, $yz$, $x$, $y$, $z$, wich
correspond ,e.g., for $D_{ab}$ to $D_{xx}$, $D_{yy}$, $D_{zz}$,
$(D_{xy}+D_{yx})/2$, $(D_{xz}+D_{zx})/2$, $(D_{yz}+D_{zy})/2$,
$(D_{yz}-D_{zy})/2$, $(D_{zx}-D_{xz})/2$, $(D_{xy}-D_{yx})/2$

### `integer :: gyrotropic_kmesh(:)`

Overrides the `kmesh` global variable (see
Sec. [Global variables](#global-variables)).

### `real(kind=dp) :: gyrotropic_kmesh_spacing`

Overrides the `kmesh_spacing` global variable (see
Sec. [Global variables](#global-variables)).

### `real(kind=dp) :: gyrotropic_freq_min`

Lower limit of the frequency range for computing the optical activity.

Units are eV. The default value 0.0.

### `real(kind=dp) :: gyrotropic_freq_max`

Upper limit of the frequency range for computing the optical activity.
Units are eV.

If an inner energy window was specified, the default value is
`dis_froz_max`-`fermi_energy`+0.6667. Otherwise it is the difference
between the maximum and the minimum energy eigenvalue stored in
`seedname.eig`, plus 0.6667.

### `real(kind=dp) :: gyrotropic_freq_step`

Difference between consecutive values of the optical frequency between
`gyrotropic_freq_min` and `gyrotropic_freq_max`.

Units are eV. The default value is 0.01.

### `real(kind=dp) :: gyrotropic_eigval_max`

Maximum energy eigenvalue of the eigenstates to be included in the
evaluation of the Natural optical activity. Units are eV.

If an inner energy window was specified, the default value is the upper
bound of the inner energy window plus 0.6667. Otherwise it is the
maximum energy eigenvalue stored in `seedname.eig` plus 0.6667.

### `logical :: gyrotropic_smr_fixed_en_width`

Overrides the `smr_fixed_en_width` global variable (see
Sec. [Global variables](#global-variables)).

### `character(len=120) :: gyrotropic_smr_type`

Overrides the `smr_type` global variable (see
Sec. [Global variables](#global-variables)).

### `character(len=120) :: gyrotropic_degen_thresh`

The threshould to eliminate degenerate bands from the calculation in
order to avoid divergences.

Units are eV. The dfault value is 0.

### `character(len=120) :: gyrotropic_box_center`

\- three real numbers. Optionally the integration may be restricted to a
parallelogram, centered at `gyrotropic_box_center` and defined by
vectors `gyrotropic_box_b{1,2,3}`

In reduced coordinates. Default value is 0.5 0.5 0.5

### `character(len=120) :: gyrotropic_box_b1`

\- three real numbers. In reduced coordinates. Default value is 1.0 0.0
0.0

### `character(len=120) :: gyrotropic_box_b2`

\- three real numbers. In reduced coordinates. Default value is 0.0 1.0
0.0

### `character(len=120) :: gyrotropic_box_b3`

\- three real numbers. In reduced coordinates. Default value is 0.0 0.0
1.0

## BoltzWann

### `logical :: boltzwann`

Determines whether to enter the  routines.

The default value is `false`.

### `integer :: boltz_kmesh(:)`

It determines the interpolation $k$ mesh used to calculate the TDF (from
which the transport coefficient are calculated). If
`boltz_calc_also_dos` is `true`, the same $k$ mesh is used also for the
DOS. Overrides the `kmesh` global variable (see
Sec. [Global variables](#global-variables)).

### `real(kind=dp) :: boltz_kmesh_spacing`

Overrides the `kmesh_spacing` global variable (see
Sec. [Global variables](#global-variables)).

### `character(len=4) :: boltz_2d_dir` {#sec:boltz2ddir}

For two-dimensional systems, the direction along which the system is
non-periodic. It can assume the following values: `x` for a 2D system on
the $yz$ plane, `y` for a 2D system on the $xz$ plane, `z` for a 2D
system on the $xy$ plane, or `no` for a 3D system with periodicity along
all threee directions.

This value is used when calculating the Seebeck coefficient, where the
electrical conductivity tensor needs to be inverted. If the value is
different from zero, only the relevant $2\times 2$ sub-block of the
electrical conductivity is inverted.

The default value is `no`.

### `real(kind=dp) :: boltz_relax_time`

The relaxation time to be used for the calculation of the TDF and the
transport coefficients.

The units are fs. The default value is 10 fs.

### `real(kind=dp) :: boltz_mu_min`

Minimum value for the chemical potential $\mu$ for which we want to
calculate the transport coefficients.

The units are eV. No default value.

### `real(kind=dp) :: boltz_mu_max`

Maximum value for the chemical potential $\mu$ for which we want to
calculate the transport coefficients.

The units are eV. No default value.

### `real(kind=dp) :: boltz_mu_step`

Energy step for the grid of chemical potentials $\mu$ for which we want
to calculate the transport coefficients.

The units are eV. No default value.

### `real(kind=dp) :: boltz_temp_min`

Minimum value for the temperature $T$ for which we want to calculate the
transport coefficients.

The units are K. No default value.

### `real(kind=dp) :: boltz_temp_max`

Maximum value for the temperature $T$ for which we want to calculate the
transport coefficients.

The units are K. No default value.

### `real(kind=dp) :: boltz_temp_step`

Energy step for the grid of temperatures $T$ for which we want to
calculate the transport coefficients.

The units are K. No default value.

### `real(kind=dp) :: boltz_tdf_energy_step`

Energy step for the grid of energies for the TDF.

The units are eV. The default value is 0.001 eV.

### `character(len=120) :: boltz_tdf_smr_type`

The type of smearing function to be used for the TDF. The available
strings are the same of the global `smr_type` input flag.

The default value is the one given via the `smr_type` input flag (if
defined).

### `real(kind=dp) :: boltz_tdf_smr_fixed_en_width`

Energy width for the smearing function. Note that for the TDF, a
standard (non-adaptive) smearing scheme is used.

The units are eV. The default value is 0 eV. Note that if the width is
smaller than twice the energy step `boltz_tdf_energy_step`, the TDF will
be unsmeared (thus the default is to have an unsmeared TDF).

### `logical :: boltz_calc_also_dos`

Whether to calculate also the DOS while calculating the TDF.

If one needs also the DOS, it is faster to calculate the DOS using this
flag instead of using independently the routines of the `dos` module,
since in this way the interpolation on the $k$ points will be performed
only once.

The default value is `false`.

### `real(kind=dp) :: boltz_dos_energy_min`

The minimum value for the energy grid for the calculation of the DOS.

The units are eV. The default value is `minval(eigval)-0.6667`, where
`minval(eigval)` i s the minimum eigenvalue returned by the ab-initio
code on the ab-initio $q$ me sh.

### `real(kind=dp) :: boltz_dos_energy_max`

The maximum value for the energy grid for the calculation of the DOS.

The units are eV. The default value is `maxval(eigval)+0.6667`, where
`maxval(eigval)` i s the maximum eigenvalue returned by the ab-initio
code on the ab-initio $q$ me sh.

### `real(kind=dp) :: boltz_dos_energy_step`

Energy step for the grid of energies for the DOS.

The units are eV. The default value is 0.001 eV.

### `character(len=120) :: boltz_dos_smr_type`

Overrides the `smr_type` global variable (see
Sec. [Global variables](#global-variables)).

### `logical :: boltz_dos_adpt_smr`

Overrides the `adpt_smr` global variable (see
Sec. [Global variables](#global-variables)).

### `real(kind=dp) :: boltz_dos_adpt_smr_fac`

Overrides the `adpt_smr_fac` global variable (see
Sec. [Global variables](#global-variables)).

### `real(kind=dp) :: boltz_dos_adpt_smr_max`

Overrides the `adpt_smr_max` global variable (see
Sec. [Global variables](#global-variables)).

### `logical :: boltz_dos_smr_fixed_en_width`

Overrides the `smr_fixed_en_width` global variable (see
Sec. [Global variables](#global-variables)).

### `logical :: boltz_bandshift`

Shift all conduction bands by a given amount (defined by
`boltz_bandshift_energyshift`).

!!! note
    this flag slightly differs from the global `scissors_shift` flag:
    with `boltz_bandshift`, an exact rigid shift is applied *after*
    interpolation; `scissors_shift` applies instead the shift *before*
    interpolation. As a consequence, results may slightly differ (and this
    is why we provide both possibilities). Note also that with
    `scissors_shift` you have to provide the number of valence bands
    `num_valence_bands`, while with `boltz_bandshift` you should provide the
    first band to shift `boltz_bandshift_firstband` =
    `num_valence_bands`$+1$.

The default value is `false`.

### `integer :: boltz_bandshift_firstband`

Index of the first conduction band to shift.

That means that all bands with index
$i\ge {\tt boltz\_bandshift\_firstband}$ will be shifted by
`boltz_bandshift_energyshift`, if `boltz_bandshift` is `true`.

The units are eV. No default value; if `boltz_bandshift` is `true`, this
flag must be provided.

### `real(kind=dp) :: boltz_bandshift_energyshift`

Energy shift of the conduction bands.

The units are eV. No default value; if `boltz_bandshift` is `true`, this
flag must be provided.

## Generic Band Interpolation

### `logical :: geninterp`

Determines whether to enter the Generic Band Interpolation routines.

The default value is `false`.

### `logical :: geninterp_alsofirstder`

Whether to calculate also the first derivatives of the bands at the
given $k$ points.

The default value is `false`.

### `logical :: geninterp_single_file`

Whether to write a single `seedname_geninterp.dat` file (all I/O is done
by the root node); or instead multiple files (one for each node) with
names `seedname_geninterp_NNNNN.dat`, where `NNNNN` is the node number.
See also the discussion in
Sec. [`seedname_geninterp.dat` or `seedname_geninterp_NNNNN.dat`](geninterp.md#sec:seedname.geninterp.dat).

The default value is `true`.
