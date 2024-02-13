# Parameters

## Usage

`wannier90.x` can be run in parallel using MPI libraries to reduce the
computation time.

For serial execution use: `wannier90.x [-pp] [seedname]`

- `seedname`: If a seedname string is given the code will read its
    input from a file `seedname.win`. The default value is `wannier`.
    One can also equivalently provide the string `seedname.win` instead
    of `seedname`.

- `-pp`: This optional flag tells the code to generate a list of the
    required overlaps and then exit. This information is written to the
    file `seedname.nnkp`.

For parallel execution use:
`mpirun -np NUMPROCS wannier90.x [-pp] [seedname]`

- `NUMPROCS`: substitute with the number of processors that you want
    to use.

Note that the `mpirun` command and command-line flags may be different
in your MPI implementation: read your MPI manual or ask your computer
administrator.

Note also that this requires that the `wannier90.x` executable has been
compiled in its parallel version (follow the instructions in the file
`README.install` in the main directory of the wannier90 distribution)
and that the MPI libraries and binaries are installed and correctly
configured on your machine.

## `seedname.win` File

The `wannier90` input file `seedname.win` has a flexible free-form
structure.

The ordering of the keywords is not significant. Case is ignored (so
`num_bands` is the same as `Num_Bands`). Characters after !, or \# are
treated as comments. Most keywords have a default value that is used
unless the keyword is given in `seedname.win`. Keywords can be set in
any of the following ways

```
num_wann 4
num_wann = 4
num_wann : 4
```

A logical keyword can be set to `true` using any of the following
strings: `T`, `true`, `.true.`.

For further examples see Section [10.1](#winfile){reference-type="ref"
reference="winfile"} and the the `wannier90` Tutorial.

## Keyword List

  ------------------- ------ ----------------------------------------------------------------------------------
        Keyword        Type  Description
                             
   System Parameters         
       num_wann         I    Number of WF
       num_bands        I    Number of bands passed to the code
    unit_cell_cart      P    Unit cell vectors in Cartesian coordinates
     atoms_cart \*      P    Positions of atoms in Cartesian coordinates
     atoms_frac \*      R    Positions of atoms in fractional coordinates with respect to the lattice vectors
        mp_grid         I    Dimensions of the Monkhorst-Pack grid of k-points
        kpoints         R    List of k-points in the Monkhorst-Pack grid
      gamma_only        L    Wavefunctions from underlying ab initio calculation are manifestly real
        spinors         L    WF are spinors
      shell_list        I    Which shells to use in finite difference formula
     search_shells      I    The number of shells to search when determining finite difference formula
     skip_B1_tests      L    Check the condition B1 of Ref. [@marzari-prb97]
        nnkpts          I    Explicit list of nearest-neighbour k-points.
       kmesh_tol        R    The tolerance to control if two kpoint belong to the same shell
  ------------------- ------ ----------------------------------------------------------------------------------

  : `seedname.win` file keywords defining the system. Argument types are
  represented by, I for a integer, R for a real number, P for a physical
  value, L for a logical value and S for a text string.\
  \* atoms_cart and atoms_frac may not both be defined in the same input
  file.
:::
:::

::: center
::: {#parameter_keywords2}
  ---------------------- ------ -------------------------------------------------------------------------------------------------
         Keyword          Type  Description
                                
       Job Control              
      postproc_setup       L    To output the `seedname.nnkp` file
      exclude_bands        I    List of bands to exclude from the calculation
    select_projections     I    List of projections to use in Wannierisation
     auto_projections      L    To automatically generate initial projections
         restart           S    Restart from checkpoint file
          iprint           I    Output verbosity level
       length_unit         S    System of units to output lengths
      wvfn_formatted       L    Read the wavefunctions from a (un)formatted file
           spin            S    Which spin channel to read
        devel_flag         S    Flag for development use
       timing_level        I    Determines amount of timing information written to output
       optimisation        I    Optimisation level
   translate_home_cell     L    To translate final Wannier centres to home unit cell when writing xyz file
        write_xyz          L    To write atomic positions and final centres in xyz file format
      write_vdw_data       L    To write data for futher processing by w90vdw utility
      write_hr_diag        L    To write the diagonal elements of the Hamiltonian in the Wannier basis to seedname.wout (in eV)
  ---------------------- ------ -------------------------------------------------------------------------------------------------

  : `seedname.win` file keywords defining job control. Argument types
  are represented by, I for a integer, R for a real number, P for a
  physical value, L for a logical value and S for a text string.
  translate_home_cell only relevant if write_xyz is `.true.`
:::
:::

::: center
::: {#parameter_keywords4}
  ---------------------------- ------ ---------------------------------------------------------------------------------------
            Keyword             Type  Description
                                      
   Disentanglement Parameters         
          dis_win_min            P    Bottom of the outer energy window
          dis_win_max            P    Top of the outer energy window
          dis_froz_min           P    Bottom of the inner (frozen) energy window
          dis_froz_max           P    Top of the inner (frozen) energy window
         dis_froz_proj           L    To activate projectability disentanglement
          dis_proj_min           P    Lower threshold for projectability disentanglement
          dis_proj_max           P    Upper threshold for projectability disentanglement
          dis_num_iter           I    Number of iterations for the minimisation of $\Omega_{\mathrm{I}}$
         dis_mix_ratio           R    Mixing ratio during the minimisation of $\Omega_{\mathrm{I}}$
          dis_conv_tol           R    The convergence tolerance for finding $\Omega_{\mathrm{I}}$
        dis_conv_window          I    The number of iterations over which convergence of $\Omega_{\mathrm{I}}$ is assessed.
        dis_spheres_num          I    Number of spheres in k-space where disentaglement is performed
     dis_spheres_first_wann      I    Index of the first band to be considered a Wannier function
          dis_spheres            R    List of centres and radii, for disentanglement only in spheres
  ---------------------------- ------ ---------------------------------------------------------------------------------------

  : `seedname.win` file keywords controlling the disentanglement.
  Argument types are represented by, I for a integer, R for a real
  number, P for a physical value, L for a logical value and S for a text
  string.
:::
:::

::: center
::: {#parameter_keywords5}
  ----------------------- ------ ------------------------------------------------------------------------------------------------------------------
          Keyword          Type  Description
                                 
   Wannierise Parameters         
         num_iter           I    Number of iterations for the minimisation of $\Omega$
       num_cg_steps         I    During the minimisation of $\Omega$ the number of Conjugate Gradient steps before resetting to Steepest Descents
        conv_window         I    The number of iterations over which convergence of $\Omega$ is assessed
         conv_tol           P    The convergence tolerance for finding $\Omega$
          precond           L    Use preconditioning
      conv_noise_amp        R    The amplitude of random noise applied towards end of minimisation procedure
      conv_noise_num        I    The number of times random noise is applied
      num_dump_cycles       I    Control frequency of check-pointing
     num_print_cycles       I    Control frequency of printing
        write_r2mn          L    Write matrix elements of $r^2$ between WF to file
      guiding_centres       L    Use guiding centres
     num_guide_cycles       I    Frequency of guiding centres
     num_no_guide_iter      I    The number of iterations after which guiding centres are used
       trial_step \*        R    The trial step length for the parabolic line search during the minimisation of $\Omega$
       fixed_step \*        R    The fixed step length to take during the minimisation of $\Omega$, instead of doing a parabolic line search
   use_bloch_phases \*\*    L    To use phases for initial projections
    site_symmetry\*\*\*     L    To construct symmetry-adapted Wannier functions
   symmetrize_eps\*\*\*     R    The convergence tolerance used in the symmetry-adapted mode
         slwf_num           I    The number of objective WFs for selective localization
      slwf_constrain        L    Whether to constrain the centres of the objective WFs
        slwf_lambda         R    Value of the Lagrange multiplier for constraining the objective WFs
       slwf_centres         P    The centres to which the objective WFs are to be constrained
  ----------------------- ------ ------------------------------------------------------------------------------------------------------------------

  : `seedname.win` file keywords controlling the wannierisation.
  Argument types are represented by, I for a integer, R for a real
  number, P for a physical value, L for a logical value and S for a text
  string. \* fixed_step and trial_step may not both be defined in the
  same input file. \*\*Cannot be used in conjunction with
  disentanglement. \*\*\*Cannot be used in conjunction with the inner
  (frozen) energy window.
:::
:::

::: {#parameter_keywords6}
  ------------------------------------------ ------ -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                   Keyword                    Type  Description
                                                    
               Plot Parameters                      
                 wannier_plot                  L    Plot the WF
              wannier_plot_list                I    List of WF to plot
            wannier_plot_supercell             I    Size of the supercell for plotting the WF
             wannier_plot_format               S    File format in which to plot the WF
              wannier_plot_mode                S    Mode in which to plot the WF, molecule or crystal
             wannier_plot_radius               R    Cut-off radius of WF\*
              wannier_plot_scale               R    Scaling parameter for cube files
           wannier_plot_spinor_mode            S    Quantity to plot for spinor WF
          wannier_plot_spinor_phase            L    Include the "phase" when plotting spinor WF
                  bands_plot                   L    Plot interpolated band structure
                 kpoint_path                   P    K-point path for the interpolated band structure
               bands_num_points                I    Number of points along the first section of the k-point path
              bands_plot_format                S    File format in which to plot the interpolated bands
              bands_plot_project               I    WF to project the band structure onto
               bands_plot_mode                 S    Slater-Koster type interpolation or Hamiltonian cut-off
                bands_plot_dim                 I    Dimension of the system
              fermi_surface_plot               L    Plot the Fermi surface
           fermi_surface_num_points            I    Number of points in the Fermi surface plot
                 fermi_energy                  P    The Fermi energy
               fermi_energy_min                P    Lower limit of the Fermi energy range
               fermi_energy_max                P    Upper limit of the Fermi energy range
              fermi_energy_step                R    Step for increasing the Fermi energy in the specified range
          fermi_surface_plot_format            S    File format for the Fermi surface plot
        [hr_plot]{style="color: red"}          L    [This parameter is not used anymore. Use write_hr instead.]{style="color: red"}
       [write_hr]{style="color: blue"}         L    [Write the Hamiltonian in the WF basis]{style="color: blue"}
      [write_rmn ]{style="color: blue"}        L    [Write the position operator in the WF basis]{style="color: blue"}
      [write_bvec ]{style="color: blue"}       L    [Write to file the matrix elements of the bvectors and their weights]{style="color: blue"}
       [write_tb ]{style="color: blue"}        L    [Write lattice vectors, Hamiltonian, and position operator in WF basis]{style="color: blue"}
                  hr_cutoff                    P    Cut-off for the absolute value of the Hamiltonian
                 dist_cutoff                   P    Cut-off for the distance between WF
               dist_cutoff_mode                S    Dimension in which the distance between WF is calculated
           translation_centre_frac             R    Centre of the unit cell to which final WF are translated
   [use_ws_distance ]{style="color: blue"}     L    [Improve interpolation using minimum distance between WFs, see Chap. [\[chap:interpolation\]](#chap:interpolation){reference-type="ref" reference="chap:interpolation"}]{style="color: blue"}
   [ws_distance_tol ]{style="color: blue"}     R    [Absolute tolerance for the distance to equivalent positions.]{style="color: blue"}
    [ws_search_size ]{style="color: blue"}     I    [Maximum extension in each direction of the super-cell of the Born-von Karmann cell to search for points inside the Wigner-Seitz cell]{style="color: blue"}
   [write_u_matrices ]{style="color: blue"}    L    [Write $\mathbf{U}^{(\mathbf{k})}$ and $\mathbf{U}^{\mathrm{dis}(\mathbf{k})}$ matrices to files]{style="color: blue"}
                                                    
  ------------------------------------------ ------ -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  : `seedname.win` file keywords controlling the plotting. Argument
  types are represented by, I for a integer, R for a real number, P for
  a physical value, L for a logical value and S for a text string. \*
  Only applies when wannier_plot_format is `cube`.
:::

::: center
::: {#parameter_keywords7}
  -------------------------- ------ -----------------------------------------------------------
           Keyword            Type  Description
                                    
     Transport Parameters           
          transport            L    Calculate quantum conductance and density of states
        transport_mode         S    Bulk or left-lead_conductor_right-lead calculation
         tran_win_min          P    Bottom of the energy window for transport calculation
         tran_win_max          P    Top of the energy window for transport calculation
       tran_energy_step        R    Sampling interval of the energy values
         fermi_energy          R    The Fermi energy
         tran_num_bb           I    Size of a bulk Hamiltonian
         tran_num_ll           I    Size of a left-lead Hamiltonian
         tran_num_rr           I    Size of a right-lead Hamiltonian
         tran_num_cc           I    Size of a conductor Hamiltonian
         tran_num_lc           I    Number of columns in a left-lead_conductor Hamiltonian
         tran_num_cr           I    Number of rows in a conductor_right-lead Hamiltonian
       tran_num_cell_ll        I    Number of unit cells in PL of left lead
       tran_num_cell_rr        I    Number of unit cells in PL of right lead
        tran_num_bandc         I    Half-bandwidth+1 of a band-diagonal conductor Hamiltonian
        tran_write_ht          L    Write the Hamiltonian for transport calculation
         tran_read_ht          L    Read the Hamiltonian for transport calculation
      tran_use_same_lead       L    Left and right leads are the same
     tran_group_threshold      R    Distance that determines the grouping of WFs
          hr_cutoff            P    Cut-off for the absolute value of the Hamiltonian
         dist_cutoff           P    Cut-off for the distance between WF
       dist_cutoff_mode        S    Dimension in which the distance between WF is calculated
         one_dim_axis          S    Extended direction for a one-dimensional system
   translation_centre_frac     R    Centre of the unit cell to which final WF are translated
  -------------------------- ------ -----------------------------------------------------------

  : `seedname.win` file keywords controlling transport. Argument types
  are represented by, I for a integer, R for a real number, P for a
  physical value, L for a logical value and S for a text string.
:::
:::

## System

### `integer :: num_wann`

Number of WF to be found.

No default.

### `integer :: num_bands`

Total number of bands passed to the code in the `seedname.mmn` file.

Default `num_bands`=`num_wann`

### Cell Lattice Vectors

The cell lattice vectors should be specified in Cartesian coordinates.

`begin unit_cell_cart`\
`[units]` $$\begin{array}{ccc}
A_{1x} & A_{1y} & A_{1z} \\
A_{2x} & A_{2y} & A_{2z} \\
A_{3x} & A_{3y} & A_{3z}
\end{array}$$ `end unit_cell_cart`

Here $A_{1x}$ is the $x$-component of the first lattice vector
$\mathbf{A}_1$, $A_{2y}$ is the $y$-component of the second lattice
vector $\mathbf{A}_2$, etc.

`[units]` specifies the units in which the lattice vectors are defined:
either `Bohr` or `Ang`.

The default value is `Ang`.

### Ionic Positions

The ionic positions may be specified in fractional coordinates relative
to the lattice vectors of the unit cell, or in absolute Cartesian
coordinates. Only one of `atoms_cart` and `atoms_frac` may be given in
the input file.

#### Cartesian coordinates

`begin atoms_cart`\
`[units]` $$\begin{array}{cccc}
P  & R^{P}_{x} & R^{P}_{y} & R^{P}_{z} \\
Q  & R^{Q}_{x} & R^{Q}_{y} & R^{Q}_{z} \\
\vdots
\end{array}$$ `end atoms_cart`

The first entry on a line is the atomic symbol. The next three entries
are the atom's position $\mathbf{R}=(R_x , R_y, R_z)$ in Cartesian
coordinates. The first line of the block, `[units]`, specifies the units
in which the coordinates are given and can be either `bohr` or `ang`. If
not present, the default is `ang`.

#### Fractional coordinates

`begin atoms_frac` $$\begin{array}{cccc}
P  & F^{P}_{1} & F^{P}_{2} & F^{P}_{3} \\
Q  & F^{Q}_{1} & F^{Q}_{2} & F^{Q}_{3} \\
\vdots
\end{array}$$ `end atoms_frac`

The first entry on a line is the atomic symbol. The next three entries
are the atom's position in fractional coordinates $\mathbf{F} = F_1
\mathbf{A}_{1} + F_2 \mathbf{A}_{2} + F_3 \mathbf{A}_{3}$ relative to
the cell lattice vectors $\mathbf{A}_i$, $i\in [1,3]$.

### `integer, dimension :: mp_grid(3)`

Dimensions of the regular (Monkhorst-Pack) k-point mesh. For example,
for a $2\times2\times2$ grid:

`mp_grid : 2  2  2`

No default.

### K-points

Each line gives the coordinate $\mathbf{K}=K_1 \mathbf{B}_{1} + K_2
\mathbf{B}_{2} + K_3 \mathbf{B}_3$ of a k-point in relative
(crystallographic) units, i.e., in fractional units with respect to the
primitive reciprocal lattice vectors $\mathbf{B}_{i}$, $i \in [1,3]$.
The position of each k-point in this list assigns its numbering; the
first k-point is k-point 1, the second is k-point 2, and so on.

`begin kpoints`\
$$\begin{array}{ccc}
 K^{1}_{1} & K^{1}_{2} & K^{1}_{3} \\
 K^{2}_{1} & K^{2}_{2} & K^{2}_{3} \\
\vdots
\end{array}$$ `end kpoints`

There is no default.

**Note**: There is an utility provided with `wannier90`, called
`kmesh.pl`, which helps to generate the explicit list of $k$ points
required by `wannier90`. See Sec. [1.1](#sec:kmesh){reference-type="ref"
reference="sec:kmesh"}.

### `logical :: gamma_only`

If `gamma_only=true`, then `wannier90` uses a branch of algorithms for
disentanglement and localisation that exploit the fact that the Bloch
eigenstates obtained from the underlying ab initio calculation are
manifestly real. This can be the case when only the $\Gamma$-point is
used to sample the Brillouin zone. The localisation procedure that is
used in the $\Gamma$-only branch is based on the method of
Ref. [@gygi-cpc03].

The default value is `false`.

### `logical :: spinors`

If `spinors=true`, then `wannier90` assumes that the WF correspond to
singularly occupied spinor states and `num_elec_per_state=1`.

The default value is `false`.

### Shells

The MV scheme requires a finite difference expression for
$\nabla_{\bf k}$ defined on a uniform Monkhorst-Pack mesh of k-points.
The vectors $\{{\bf b}\}$ connect each mesh-point ${\bf k}$ to its
nearest neighbours. $N_{\mathrm{sh}}$ shells of neighbours are included
in the finite-difference formula, with $M_s$ vectors in the
$s^{\mathrm{th}}$ shell. For $\nabla_{{\bf k}}$ to be correct to linear
order, we require that the following equation is satisfied (Eq. B1 of
Ref. [@marzari-prb97]): $$\label{eq:B1}
\sum_{s}^{N_{\mathrm{sh}}} w_s \sum_i^{M_{\mathrm{s}}}
b_{\alpha}^{i,s} b_{\beta}^{i,s} = \delta_{\alpha\beta}\:,$$ where
${\bf b}^{i,s}$, $i\in[1,M_s]$, is the $i^{\mathrm{th}}$ vector
belonging to the $s^{\mathrm{th}}$ shell with associated weight $w_s$,
and $\alpha$ and $\beta$ run over the three Cartesian indices.

### `integer :: shell_list(:)`

`shell_list` is vector listing the shells to include in the finite
difference expression. If this keyword is absent, the shells are chosen
automatically.

### `integer :: search_shells`

Specifies the number of shells of neighbours over which to search in
attempting to determine an automatic solution to the B1 condition
Eq. [\[eq:B1\]](#eq:B1){reference-type="ref" reference="eq:B1"}. Larger
values than the default may be required in special cases e.g. for very
long thin unit cells.

The default value is 36.

### `logical :: skip_B1_tests`

If set to `.true.`, does not check the B1 condition
Eq. [\[eq:B1\]](#eq:B1){reference-type="ref" reference="eq:B1"}. This
should *only* be used if one knows why the B1 condition should not be
verified. A typical use of this flag is in conjunction with the Z2PACK
code: <http://www.physics.rutgers.edu/z2pack/>.

The default value is `.false.`.

### `integer, dimension(:, 5) :: nnkpts`

Specifies the nearest-neighbour k-points which are written to the
`.nnkp` file. This can be used to explicitly specify which overlap
matrices should be calculated.

    begin nnkpts
    1   2   0  0  0
    .
    .
    end nnkpts

Each nearest neighbour $\mathbf{k + b}$ is given by a line of 5
integers. The first specifies the k-point number `nkp` of $\mathbf{k}$.
The second is the k-point number of the neighbour. The final three
integers specify the reciprocal lattice vector which brings the k-point
specified by the second integer to $\mathbf{k + b}$.

This format is the same as in the `.nnkp` file, except that the number
of neighbours per k-point is not specified. However, the number of
neighbours still needs to be a multiple of the number of k-points.

This input parameter can be used only if `postproc_setup = .true.`, and
is not intended to be used with a full Wannier90 run. It can be used
also if the k-points do not describe a regular mesh.

### `real(kind=dp) :: kmesh_tol`

Two kpoints belong to the same shell if the distance between them is
less than `kmesh_tol`. Units are Ang.

The default value is 0.000001 Ang.

## Projection

The projections block defines a set of localised functions used to
generate an initial guess for the unitary transformations. This data
will be written in the `seedname.nnkp` file to be used by a
first-principles code.

`begin projections`\
.\
.\
`end projections`

If `guiding_centres`=`true`, then the projection centres are used as the
guiding centres in the Wannierisation routine.

For details see Section [3.1](#sec:proj){reference-type="ref"
reference="sec:proj"}.

## Job Control

### `logical :: postproc_setup`

If `postproc_setup`=`true`, then the wannier code will write
`seedname.nnkp` file and exit. If `wannier90` is called with the option
`-pp`, then `postproc_setup` is set to `true`, over-riding its value in
the `seedname.win` file.

The default value is `false`.

### `integer :: iprint`

This indicates the level of verbosity of the output from 0 ("low"), the
bare minimum, to 3 ("high"), which corresponds to full debugging output.

The default value is 1.

### `integer :: optimisation`

This indicates the level of optimisation used in the code. This is a
trade between speed and memory. A positive number indicates fastest
execution time at the cost of more memory. Zero or negative numbers
indicates a smaller memory footprint - at increased execution time.

At the moment the only values that have an effect are `optimisation<=0`
(low memory) and `optimisation>0` (fast)

The default value is 3.

### `character(len=20) :: length_unit`

The length unit to be used for writing quantities in the output file
`seedname.wout`.

The valid options for this parameter are:

-   `Ang` (default)

-   `Bohr`

### `character(len=50) :: devel_flag`

Not a regular keyword. Its purpose is to allow a developer to pass a
string into the code to be used inside a new routine as it is developed.

No default.

### `integer :: exclude_bands(:)`

A k-point independent list of states to excluded from the calculation of
the overlap matrices; for example to select only valence states, or
ignore semi-core states. This keyword is passed to the first-principles
code via the `seedname.nnkp` file. For example, to exclude bands 2, 6,
7, 8 and 12:

`exclude_bands : 2, 6-8, 12`

### `integer :: select_projections(:)`

A list of projections to be included in the wannierisation procedure. In
the case that `num_proj` is greater than `num_wann`, this keyword allows
a subset of the projections in the projection matrices to be used. For
example, to select the projections given by the indices 2, 6, 7, 8 and
12:

`select_projections : 2, 6-8, 12`

### `logical :: auto_projections`

If `.true.` and no projections block is defined, then `wannier90` writes
an additional block in the `.nnkp` file during the pre-processing step,
to instruct the interface code to automatically generate the
$A_{mn}^{(\mathbf{k})}$.

For additional information on the behavior and on the added block, see
Sec. [\[sec:auto-projections-block\]](#sec:auto-projections-block){reference-type="ref"
reference="sec:auto-projections-block"}.

**Note:** the interface code (e.g. `pw2wannier90.x`) must have at least
one implementation of a method to automatically generate initial
projections in order for this option to be usable.

The default value of this parameter is $\verb#false#$.

### `character(len=20) :: restart`

If `restart` is present the code will attempt to restart the calculation
from the `seedname.chk ` file. The value of the parameter determines the
position of the restart

The valid options for this parameter are:

-   `default`. Restart from the point at which the check file
    `seedname.chk` was written

-   `wannierise`. Restart from the beginning of the wannierise routine

-   `plot`. Go directly to the plotting phase

-   `transport`. Go directly to the transport routines

### `character(len=20) :: wvfn_formatted`

If `wvfn_formatted`=`true`, then the wavefunctions will be read from
disk as formatted (ie ASCII) files; otherwise they will be read as
unformatted files. Unformatted is generally preferable as the files will
take less disk space and I/O is significantly faster. However such files
will not be transferable between all machine architectures and formatted
files should be used if transferability is required (i.e., for test
cases).

The default value of this parameter is $\verb#false#$.

### `character(len=20) :: spin`

For bands from a spin polarised calculation `spin` determines which set
of bands to read in, either `up` or `down`.

The default value of this parameter is `up`.

### `integer :: timing_level`

Determines the amount of timing information regarding the calculation
that will be written to the output file. A value of 1 produces the least
information.

The default value is 1.

### `logical :: translate_home_cell`

Determines whether to translate the final Wannier centres to the home
unit cell at the end of the calculation. Mainly useful for molecular
systems in which the molecule resides entirely within the home unit cell
and user wants to write an xyz file (`write_xyz=.true.`) for the WF
centres to compare with the structure.

The default value is `false`.

### `logical :: write_xyz`

Determines whether to write the atomic positions and final Wannier
centres to an xyz file, `seedname_centres.xyz`, for subsequent
visualisation.

The default value is `false`.

### `logical :: write_vdw_data`

Determines whether to write `seedname.vdw` for subsequent
post-processing by the `w90vdw` utility (in the `utility/w90vdw/`
directory of the distribution) for calculating van der Waals energies.
Brillouin zone sampling must be at the Gamma-point only.

The default value is `false`.

## Disentanglement

These keywords control the disentanglement routine of
Ref. [@souza-prb01], i.e., the iterative minimisation of
$\Omega_{\mathrm{I}}$. This routine will be activated if
`num_wann`$\:<\:$`num_bands`.

### `real(kind=dp) :: dis_win_min`

The lower bound of the outer energy window for the disentanglement
procedure. Units are eV.

The default is the lowest eigenvalue in the system.

### `real(kind=dp) :: dis_win_max`

The upper bound of the outer energy window for the disentanglement
procedure. Units are eV.

The default is the highest eigenvalue in the given states (i.e., all
states are included in the disentanglement procedure).

### `real(kind=dp) :: dis_froz_min`

The lower bound of the inner energy window for the disentanglement
procedure. Units are eV.

If `dis_froz_max` is given, then the default for `dis_froz_min` is
`dis_win_min`.

### `real(kind=dp) :: dis_froz_max`

The upper bound of the inner (frozen) energy window for the
disentanglement procedure. If `dis_froz_max` is not specified, then
there are no frozen states. Units are eV.

No default.

### `logical :: dis_froz_proj`

To activate projectability disentanglement procedure, which selectively
discard/disentangle/freeze state $\vert n \mathbf{k}\rangle$ based on
its projectability onto some localized atomic orbitals.

Note: this requires the `amn` file is properly normalized, i.e.,
projectability computed from $A A^\dagger$ must be smaller than or equal
to 1. The pseudo-atomic projection satisfies such requirement, see
[3.6](#sec:proj_pdwf){reference-type="ref" reference="sec:proj_pdwf"}.

Additionally, one can combine projectability disentanglement with energy
disentanglement, i.e., enable both `dis_proj_min/max` and
`dis_froz_min/max` simultaneously in the `win` file. These settings will
freeze the union of inner energy window and high-projectability states,
and exclude the union of states outside outer energy window and having
low projectability.

### `real(kind=dp) :: dis_proj_min`

The lower bound for the projectability disentanglement procedure.

For states with projectabilities smaller than `dis_proj_min`, they will
be discarded in the disentanglement procedure, i.e., similar to the case
of outside of the outer energy window.

For states with projectabilities larger than or equal to `dis_proj_min`,
they will be included in the disentanglement procedure, i.e., similar to
the case of inside the outer energy window.

No unit.

The default value is 0.95.

### `real(kind=dp) :: dis_proj_max`

The upper bound for the projectability disentanglement procedure. For
states with projectability larger than or equal to `dis_proj_max`, they
will be freezed in the disentanglement procedure, i.e., similar to the
case of inside the inner energy window.

No unit.

The default value is 0.01.

### `integer :: dis_num_iter`

In the disentanglement procedure, the number of iterations used to
extract the most connected subspace.

The default value is 200.

### `real(kind=dp) :: dis_mix_ratio`

In the disentanglement procedure, the mixing parameter to use for
convergence (see pages 4-5 of Ref. [@souza-prb01]). A value of 0.5 is a
'safe' choice. Using 1.0 (i.e., no mixing) often gives faster
convergence, but may cause the minimisation of $\Omega_{\mathrm{I}}$ to
be unstable in some cases.

Restriction: $0.0<\:$`dis_mix_ratio`$\:\leq 1.0$

The default value is 0.5

### `real(kind=dp) :: dis_conv_tol`

In the disentanglement procedure, the minimisation of
$\Omega_{\mathrm{I}}$ is said to be converged if the fractional change
in the gauge-invariant spread between successive iterations is less than
`dis_conv_tol` for `dis_conv_window` iterations. Units are Å$^2$.

The default value is 1.0E-10

### `integer :: dis_conv_window`

In the disentanglement procedure, the minimisation is said to be
converged if the fractional change in the spread between successive
iterations is less than `dis_conv_tol` for `dis_conv_window` iterations.

The default value of this parameter is 3.

### `integer :: dis_spheres_num`

Number of spheres in reciprocal space where the k-dependent
disentanglement is performed. No disentanglement is performed for those
k-points that are not included in any of the spheres.

The default is 0, which means disentangle at every k-point in the full
BZ (the standard mode in Wannier90).

### `integer :: dis_spheres_first_wann`

Index of the first band that has to be considered as a Wannier function.
Used only if `dis_spheres_num` is greater than zero. At k-points where
disentanglement is not performed the bands from `dis_spheres_first_wann`
to `dis_spheres_first_wann+num_wann` are used to wannierise. The bands
excluded using `exclude_bands` should not be counted.

The default is 1, the band at the lowest energy.

### dis_spheres

Each line gives the coordinate $\mathbf{K}=K_1 \mathbf{B}_{1} + K_2
\mathbf{B}_{2} + K_3 \mathbf{B}_3$ of a k-point representing the center
of one of the spheres used for k-dependent disentanglement. The same
crystallographic units as for `kpoints` are used here. Each k-point
coordinate $\mathbf{K}^i$ must the followed by the respectice sphere
radius $r_{i}$ in inverse angstrom (on the same line).

The number of lines must be equal to `dis_spheres_num`.

`begin dis_spheres` $$\begin{array}{cccc}
 K^{1}_{1} & K^{1}_{2} & K^{1}_{3} & r_{1} \\
 K^{2}_{1} & K^{2}_{2} & K^{2}_{3} & r_{2} \\
\vdots
\end{array}$$ `end dis_spheres`

There is no default.

## Wannierise

Iterative minimisation of $\widetilde{\Omega}$, the non-gauge-invariant
part of the spread functional.

### `integer :: num_iter`

Total number of iterations in the minimisation procedure. Set
`num_iter=0` if you wish to generate projected WFs rather than
maximally-localized WFs (see Example 8 in the Tutorial).

The default value is 100

### `integer :: num_cg_steps`

Number of conjugate gradient steps to take before resetting to steepest
descents.

The default value is 5

### `integer :: conv_window`

If `conv_window`$\:>1$, then the minimisation is said to be converged if
the change in $\Omega$ over ` conv_window` successive iterations is less
than ` conv_tol`. Otherwise, the minimisation proceeds for num_iter
iterations (default).

The default value is -1

### `real(kind=dp) :: conv_tol`

If `conv_window`$\:>1$, then this is the convergence tolerance on
$\Omega$, otherwise not used. Units are Å$^2$.

The default value is 1.0E-10

### `logical :: precond`

Whether or not to use preconditioning to speed up the minimization of
the spreads. This is based on the same idea as the classical
Tetter-Payne-Allan preconditionning for DFT and dampens the
high-frequency oscillations of the gradient due to contributions from
large real lattice vectors. It is useful when the optimization is slow,
especially on fine grids. When `optimisation<3`, this uses a slower
algorithm to save memory.

The default value is `false`.

### `real(kind=dp) :: conv_noise_amp`

If `conv_noise_amp`$\:>0$, once convergence (as defined above) is
achieved, some random noise $f$ is added to the search direction, and
the minimisation is continued until convergence is achieved once more.
If the same value of $\Omega$ as before is arrived at, then the
calculation is considered to be converged. If not, then random noise is
added again and the procedure repeated up to a maximum of
` conv_noise_num` times. `conv_noise_amp` is the amplitude of the random
noise $f$ that is added to the search direction:
$0 < |f| <\:$`conv_noise_amp`. This functionality requires
` conv_window`$\:>1$. If `conv_window` is not specified, it is set to
the value 5 by default.

If `conv_noise_amp`$\:\leq 0$, then no noise is added (default).

The default value is -1.0

### `integer :: conv_noise_num`

If `conv_noise_amp`$\:>0$, then this is the number of times in the
minimisation that random noise is added.

The default value is 3

### `integer :: num_dump_cycles`

Write sufficient information to do a restart every `num_dump_cycles`
iterations.

The default is 100

### `integer :: num_print_cycles`

Write data to the master output file `seedname.wout` every
`num_print_cycles` iterations.

The default is 1

### `logical :: write_r2mn`

If $\verb#write_r2mn#=\verb#true#$, then the matrix elements
$\langle m|r^2|n\rangle$ (where $m$ and $n$ refer to WF) are written to
file `seedname.r2mn` at the end of the Wannierisation procedure.

The default value of this parameter is `false`.

### `logical :: guiding_centres`

Use guiding centres during the minimisation, in order to avoid local
minima.

`wannier90` uses a logarithm definition of the spread functional. As we
are taking the log of a complex argument there is a possibility that the
algorithm might make inconsistent choices for the branch cut. This
manifests itself as complex WF with a large spread. By using guiding
centres the code will attempt to make a consistent choice of branch cut.
Experience shows that with `guiding_centres` set to true this problem is
avoided and doing so does not cause any problems. For this reason we
recommend setting `guiding_centres` to true where possible (it is only
not possible if an explicit projection block is not defined).

The default value is `false`.

### `integer :: num_guide_cycles`

If `guiding_centres` is set to true, then the guiding centres are used
only every `num_guide_cycles`.

The default value is 1.

### `integer :: num_no_guide_iter`

If `guiding_centres` is set to true, then the guiding centres are used
only after `num_no_guide_iter` minimisation iterations have been
completed.

The default value is 0.

### `real(kind=dp) :: trial_step`

The value of the trial step for the parabolic fit in the line search
minimisation used in the minimisation of the spread function. Cannot be
used in conjunction with `fixed_step` (see below). If the minimisation
procedure doesn't converge, try decreasing the value of `trial_step` to
give a more accurate line search.

The default value is 2.0

### `real(kind=dp) :: fixed_step`

If this is given a value in the input file, then a fixed step of length
`fixed_step` (instead of a parabolic line search) is used at each
iteration of the spread function minimisation. Cannot be used in
conjunction with `trial_step`. This can be useful in cases in which
minimisation with a line search fails to converge.

There is no default value.

### `logical :: use_bloch_phases`

Determines whether to use the Bloch functions as the initial guess for
the projections. Can only be used if `disentanglement = false`.

The default value is `false`.

### `logical :: site_symmetry`

Construct symmetry-adapted Wannier functions. For the detail of the
theoretical background, see Ref. [@sakuma-prb13]. Cannot be used in
conjunction with the inner (frozen) energy window.

The default value is `false`.

### `real(kind=dp) :: symmetrize_eps`

Convergence threshold to check whether the symmetry condition (Eq. (19)
in Ref. [@sakuma-prb13]) on the unitary matrix
$\mathbf{U}^{(\mathbf{k})}$ is satisfied or not. See also Eq. (29) in
Ref. [@sakuma-prb13]. Used when `site_symmetry = .true`.

The default value is 1.0E-3.

### `integer :: slwf_num`

The number of objective Wannier functions for selective localisation in
the selectively localised Wannier function (SLWF) method of
Ref. [@Marianetti]. These functions are obtained by minimising the
spread functional only with respect to the degrees of freedom of a
subset of `slwf_num` $<$ `num_wann` functions. At convergence, the
objective WFs will have a minimum cumulative spread, whereas the
remaining `num_wann` $-$ `slwf_num` functions are left unoptimised. The
initial guesses for the objective WFs are given by the first `slwf_num`
orbitals in the `projections` block. If `slwf_num = num_wann` no
selective minimisation is performed. In this case, `wannier90` will
simply generate a set of `num_wann` MLWFs.

The default is `num_wann`.

### `logical :: slwf_constrain`

If `slwf_constrain=true`, then the centres of the objective Wannier
functions are constrained to either the centres of the first `slwf_num`
orbitals in the `projections` block or to new positions specified in the
`slwf_centres` block (see
Sec. [2.8.22](#sec:centre_constraints){reference-type="ref"
reference="sec:centre_constraints"}). In this case, a modified spread
functional, $\Omega_c$, with the addition of a constraint term, as
described in Ref. [@Marianetti].

The default is `false`

### `real(kind=dp) :: slwf_lambda`

The value of the Lagrange multiplier $\lambda$ for the constraint term
in term added to modify the spread functional:
$\lambda \sum_{n=1}^{J'} \left(\overline{\mathbf{r}}_n - \mathbf{r}_{0n}\right)^2$,
where $J'$ is `slwf_num`, and $\overline{\mathbf{r}}_{n}$ and
$\mathbf{r}_{0n}$ are the centre and target centre, respectively, for
the $n^{\text{th}}$ objective WF.

The default is `0.0`.

### Constraints on centres {#sec:centre_constraints}

If `slwf_constrain=true`, then by default the centres to which the
`slwf_num` objective Wannier function centres are constrained are given
by the first `slwf_num` rows of the `projections` block.

Optionally, the `slwf_centres` block may be used to define alternative
target centres for some or all of the `slwf_num` objective Wannier
functions.

The block below shows an example of how to set the constraints:

`begin slwf_centres`\
`   2  0.0   0.0  0.0`\
`   4  0.25  0.0  0.0`\
`end slwf_centres`

-   The first line sets the constraint for the centre of objective WF
    number 2 (as defined by the order of WFs in the `projections` block)
    to (0.0,0.0,0.0) in fractional co-ordinates.

-   The second line sets the constraint for the centre of objective WF
    number 4 (as defined by the order of WFs in the `projections` block)
    to (0.25,0.0,0.0) in fractional co-ordinates.

-   The target centres of all other objective Wannier functions remain
    as the centres given in the corresponding rows of the `projections`
    block.

## Post-Processing

Capabilities:

-   Plot the WF

-   Plot the interpolated band structure

-   Plot the Fermi surface

-   Output the Hamiltonian in the WF basis

-   Transport calculation (quantum conductance and density of states)

### `logical :: wannier_plot`

If $\verb#wannier_plot#=\verb#true#$, then the code will write out the
Wannier functions in a format specified by `wannier_plot_format`

The default value of this parameter is `false`.

### `integer :: wannier_plot_list(:)`

A list of WF to plot. The WF numbered as per the `seedname.wout` file
after the minimisation of the spread.

The default behaviour is to plot all WF. For example, to plot WF 4, 5, 6
and 10:

`wannier_plot_list : 4-6, 10`

### `integer :: wannier_plot_supercell`

The code generates the WFs on a grid corresponding to a
'super-unit-cell'. If `wannier_plot_supercell` is provided as a single
integer, then the size of the super-unit-cell is
`wannier_plot_supercell` times the size of the unit cell along all three
linear dimensions (the 'home' unit cell is kept approximately in the
middle); otherwise, if three integers are provided, the size of the
super-unit-cell is `wannier_plot_supercell(i)` times the size of the
unit cell along the $i-$th linear dimension.

The default value is 2.

### `character(len=20) :: wannier_plot_format`

WF can be plotted in either XCrySDen (xsf) format or Gaussian cube
format. The valid options for this parameter are:

-   `xcrysden` (default)

-   `cube`

If `wannier_plot_format=xsf`: the code outputs the WF on the entire
super-unit-cell specified by `wannier_plot_supercell`.

If `wannier_plot_format=cube`: the code outputs the WF on a grid that is
smaller than the super-unit-cell specified by `wannier_plot_supercell`.
This grid is determined by `wannier_plot_mode`, `wannier_plot_radius`
and `wannier_plot_scale`, described in detail below.

The code is able to output Gaussian cube files for systems with
non-orthogonal lattice vectors. Many visualisation programs (including
XCrySDen), however, are only able to handle cube files for systems with
*orthogonal* lattice vectors. One visualisation program that is capable
of dealing with non-orthogonal lattice vectors is VESTA
(<http://jp-minerals.org/vesta/en/>).

!!! note

    It's worth noting that another visualisation program, VMD
    (<http://www.ks.uiuc.edu/Research/vmd>), is able to deal with
    certain special cases of non-orthogonal lattice vectors; see
    <http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/cubeplugin.html>
    for details.

### `character(len=20) :: wannier_plot_mode`

Choose the mode in which to plot the WF, either as a molecule or as a
crystal.

The valid options for this parameter are:

-   `crystal` (default)

-   `molecule`

If `wannier_plot_format=cube`:

-   if `wannier_plot_mode = molecule`, then wherever the WF centre sits
    in the supercell, the origin of the cube is shifted (for the purpose
    of plotting only, ie, nothing is done to the U matrices etc) to
    coincide with the centre of mass of the atomic positions specified
    by the user in the `.win` input file. These atomic positions are
    also written to the cube file, so when it is visualised, the WF
    appears superimposed on the molecular structure.

-   if `wannier_plot_mode = crystal`, then the WF is not shifted, but
    instead the code searches for atoms that are within a radius of
    `wannier_plot_scale` $\times$ `wannier_plot_radius` of the WF centre
    and writes the coordinates of these atoms to the cube file. In this
    way, when the cube file is visualised, the WF appears superimposed
    on the nearest atoms to the WF centre.

-   `crystal` mode can be used for molecules, and `molecule` mode can be
    used for crystals.

### `real(kind=dp) :: wannier_plot_radius`

If `wannier_plot_format=cube`, then ` wannier_plot_radius` is the radius
of the sphere that must fit inside the parallelepiped in which the WF is
plotted. `wannier_plot_radius` must be greater than 0. Units are Å.

The default value is 3.5.

### `real(kind=dp) :: wannier_plot_scale`

If `wannier_plot_format=cube` and `wannier_plot_mode=crystal`, then the
code searches for atoms that are within a radius of `wannier_plot_scale`
$\times$ `wannier_plot_radius` of the WF centre and writes the
coordinates of these atoms to the cube file. In this way, when the cube
file is visualised, the WF appears superimposed on the nearest atoms to
the WF centre. `wannier_plot_scale` must be greater than 0. This
parameter is dimensionless.

The default value is 1.0.

### `character(len=20) :: wannier_plot_spinor_mode`

If $\verb#spinors#=\verb#true#$ then this parameter controls the
quantity to plot. For a spinor WF with components $[\phi,\psi]$ the
quatity plotted is

-   `total` (default). $\sqrt{[|\phi|^2+|\psi|^2}$

-   `up`. $|\phi|\times sign(Re\{\phi\})$ if
    $\verb#wannier_plot_spinor_phase#=\verb#true#$, otherwise $|\phi|$

-   `down`. $|\psi|\times sign(Re\{\psi\})$ if
    $\verb#wannier_plot_spinor_phase#=\verb#true#$, otherwise $|\psi|$

Note: making a visual representation of a spinor WF is not as
straightforward as for a scalar WF. While a scalar WF is typically a
real valued function, a spinor WF is a complex, two component spinor.
`wannier90` is able to plot several different quantities derived from a
spinor WF which should give you a good idea of the nature of the WF.

### `logical :: wannier_plot_spinor_phase`

If $\verb#wannier_plot_spinor_phase#=\verb#true#$ phase information will
be taken into account when plotting a spinor WF.

### `logical :: bands_plot`

If $\verb#bands_plot#=\verb#true#$, then the code will calculate the
band structure, through Wannier interpolation, along the path in k-space
defined by `bands_kpath` using `bands_num_points` along the first
section of the path and write out an output file in a format specified
by `bands_plot_format`.

The default value is `false`.

### kpoint_path

Defines the path in k-space along which to calculate the bandstructure.
Each line gives the start and end point (with labels) for a section of
the path. Values are in fractional coordinates with respect to the
primitive reciprocal lattice vectors.

`begin kpoint_path` $$\begin{array}{cccccccc}
G & 0.0 & 0.0 & 0.0 & L & 0.0 & 0.0 & 1.0 \\
L & 0.0 & 0.0 & 1.0 & N & 0.0 & 1.0 & 1.0 \\
\vdots
\end{array}$$ `end kpoint_path`

There is no default

### `integer :: bands_num_points`

If $\verb#bands_plot#=\verb#true#$, then the number of points along the
first section of the bandstructure plot given by `kpoint_path`. Other
sections will have the same density of k-points.

The default value for `bands_num_points` is 100.

### `character(len=20) :: bands_plot_format`

Format in which to plot the interpolated band structure. The valid
options for this parameter are:

-   `gnuplot` (default)

-   `xmgrace`

Note: it is possible to request output in both formats eg
$\verb#bands_format#=\verb#gnuplot xmgrace#$

### `integer :: bands_plot_project(:)`

If present `wannier90` will compute the contribution of this set of WF
to the states at each point of the interpolated band structure. The WF
are numbered according to the seedname.wout file. The result is written
in the `seedname_band.dat` file, and a corresponding gnuplot script to
`seedname_band_proj.dat` .

For example, to project on to WFs 2, 6, 7, 8 and 12:

`bands_plot_project : 2, 6-8, 12`

### `character(len=20) :: bands_plot_mode`

To interpolate the band structure along the k-point path, either use the
Slater-Koster interpolation scheme or truncate the Hamiltonian matrix in
the WF basis. Truncation criteria are provided by `hr_cutoff` and
`dist_cutoff`.

The valid options for this parameter are:

-   `s-k` (default)

-   `cut`

### `integer :: bands_plot_dim`

Dimension of the system. If $\verb#bands_plot_dim#<\:$`<!-- -->`{=html}3
and $\verb#bands_plot_mode#=\verb#cut#$, lattice vector
$\mathbf{R}=N_1 \mathbf{A}_{1} + N_2 \mathbf{A}_{2} + N_3 \mathbf{A}_3$,
where $N_i=0$ if $\mathbf{A}_i$ is parallel to any of the confined
directions specified by `one_dim_axis`, are exclusively used in the band
structure interpolation.

The valid options for this parameter are:

-   3 (default)

-   2

-   1

### `logical :: fermi_surface_plot`

If $\verb#fermi_surface_plot#=\verb#true#$, then the code will
calculate, through Wannier interpolation, the eigenvalues on a regular
grid with `fermi_surface_num_points` in each direction. The code will
write a file in bxsf format which can be read by XCrySDen in order to
plot the Fermi surface.

The default value is `false`.

### `integer :: fermi_surface_num_points`

If $\verb#fermi_surface_plot#=\verb#true#$, then the number of divisions
in the regular k-point grid used to calculate the Fermi surface.

The default value for `fermi_surface_num_points` is 50.

### `real(kind=dp) :: fermi_energy`

The Fermi energy in eV. This parameter is written into the bxsf file. If
`fermi_energy` is specified, ` fermi_energy_min`, `fermi_energy_max`,
and ` fermi_energy_step` should not be specified, and vice-versa.

The default value is 0.0

### `real(kind=dp) :: fermi_energy_min`

Instead of specifyfing a single Fermi energy, it is possible to scan the
Fermi level over a range of values, and recompute certain quantities for
each $\varepsilon_F$. This is the minimum value in the range (in
eV).

!!! note

    Scanning the Fermi level is currently supported only by the
    `postw90` module `berry`, for `berry_task=ahc,morb`. For all other
    functionalities that require a knowledge of $\varepsilon_F$, use
    `fermi_energy` instead.

There is no default value.

### `real(kind=dp) :: fermi_energy_max`

The maximum value in the range of Fermi energies. Units are eV.

The default value is `fermi_energy_min`+1.0.

### `real(kind=dp) :: fermi_energy_step`

Difference between consecutive values of the Fermi energy when scanning
from `fermi_energy_min` to ` fermi_energy_max`. Units are eV.

The default value is 0.01.

### `character(len=20) :: fermi_surface_plot_format`

Format in which to plot the Fermi surface. The valid options for this
parameter are:

-   `xcrysden` (default)

### `logical :: write_hr`

If $\verb#write_hr#=\verb#true#$, then the Hamiltonian matrix in the WF
basis will be written to a file `seedname_hr.dat`.

The default value is `false`.

### `logical :: write_rmn`

If $\verb#write_rmn#=\verb#true#$, then the position operator in the WF
basis will be written to a file `seedname_r.dat`.

The default value is `false`.

### `logical :: write_bvec`

If $\verb#write_bvec#=\verb#true#$, then the the matrix elements of
bvector and their weights will be written to a file `seedname.bvec`.

The default value is `false`.

### `logical :: write_tb`

If $\verb#write_tb#=\verb#true#$, then the lattice vectors, together
with the Hamiltonian and position-operator matrices in the WF basis,
will be written to a file `seedname_tb.dat`, in units of Angstrom and
eV.

The default value is `false`.

### `logical :: transport`

If $\verb#transport#=\verb#true#$, then the code will calculate quantum
conductance and density of states of a one-dimensional system. The
results will be written to files `seedname_qc.dat` and
`seedname_dos.dat`, respectively. Since both quantities are a function
of energy, they will be evaluated from `tran_win_min` to `tran_win_max`
with an interval of `tran_energy_step`.

The default value of this parameter is `false`.

### `character(len=20) :: transport_mode`

If $\verb#transport_mode#=\verb#bulk#$, quantum conductance and density
of states are calculated for a perfectly-periodic one-dimensional
system. In this case, the transport part can either use the Hamiltonian
matrix in the WF basis generated by `wannier90` or a Hamiltonian matrix
provided by the external file `seedname_htB.dat`.

If $\verb#transport_mode#=\verb#lcr#$, quantum conductance and density
of states are calculated for a system where semi-infinite, left and
right leads are connected through a central conductor region. In this
case, the transport part will work independently from the
disentanglement and wannierise procedure. Details of the method is
described in Ref. [@nardelli-prb99].

If $\verb#tran_read_ht# = \verb#true#$ then the Hamiltonian matrices
must be provided by the five external files:
`seedname_htL.dat, seedname_htLC.dat, seedname_htC.dat, seedname_htCR.dat, seedname_htR.dat`.
If $\verb#tran_read_ht# = \verb#false#$ then the Hamiltonian matrices
are found automatically provided the supercell adheres to conditions
outlined in Section [7.3](#sec:2c2){reference-type="ref"
reference="sec:2c2"}.

The valid options for this parameter are:

-   `bulk` (default)

-   `lcr`

### `real(kind=dp) :: tran_win_min`

The lower bound of the energy window for the transport calculation.
Units are eV.

The default value is -3.0.

### `real(kind=dp) :: tran_win_max`

The upper bound of the energy window for the transport calculation.
Units are eV.

The default value is 3.0.

### `real(kind=dp) :: tran_energy_step`

Sampling interval of the energy values from `tran_win_min` to
`tran_win_max`. Units are eV.

The default value is 0.01.

### `real(kind=dp) :: fermi_energy`

The Fermi energy in eV. The energy axis of the quantum conductance and
density of states data will be shifted rigidly by this amount.

The default value is 0.0

### `integer :: tran_num_bb`

Size of a bulk Hamiltonian matrix. This number is equal to the number of
WFs in one principal layer.

A one-dimensional system can be viewed as an array of principal layers
which are defined in a way that localized basis functions inside a
certain principal layer only interact with those in the nearest neighbor
principal layer. In `wannier90` a principal layer will be an integer
multiple of a unit cell, and the size is determined by `hr_cutoff`
and/or `dist_cutoff`. The criterion is rather arbitrary when WFs are
adopted as a localized basis set, and it is up to a user's choice.

The default value is 0.

### `integer :: tran_num_ll`

Size of a left-lead Hamiltonian matrix. If
$\verb#transport_mode# = \verb#lcr#$ and
$\verb#tran_read_ht# = \verb#false#$ then `tran_num_ll` is the number of
Wannier functions in a principal layer.

The default value is 0.

### `integer :: tran_num_rr`

Size of a right-lead Hamiltonian matrix.

The default value is 0.

### `integer :: tran_num_cc`

Size of a conductor Hamiltonian matrix.

The default value is 0.

### `integer :: tran_num_lc`

Number of columns in a left-lead_conductor Hamiltonian matrix. Number of
rows must be equal to `tran_num_ll`.

The default value is 0.

### `integer :: tran_num_cr`

Number of rows in a conductor_right-lead Hamiltonian matrix. Number of
columns must be equal to `tran_num_rr`.

The default value is 0.

### `integer :: tran_num_cell_ll`

Number of unit cells in one principal layer of left lead. Used if
$\verb#transport_mode# = \verb#lcr#$ and
$\verb#tran_read_ht# = \verb#false#$.

The default value is 0.

### `integer :: tran_num_cell_rr`

Number of unit cells in one principal layer of right lead. Not used at
present.

The default value is 0.

### `integer :: tran_num_bandc`

Half-bandwidth+1 of a band-diagonal conductor Hamiltonian matrix.

The Hamiltonian matrix of a central conductor part, which is read from
`seedname_htC.dat`, will be diagonally dominant when `tran_num_cc` is
very large. `tran_num_bandc` is used to construct a compact matrix which
contains the non-zero band-diagonal part of a full conductor Hamiltonian
matrix. Setting this parameter is only meaningful when `tran_num_bandc`
is greater than `tran_num_lc` and `tran_num_cr`.

The default value is 0.

### `logical :: tran_write_ht`

If $\verb#tran_write_ht#=\verb#true#$, then the Hamiltonian matrix
formatted for the transport calculation will be written to a file
`seedname_htB.dat`.

The default value is `false`.

### `logical :: tran_read_ht`

If $\verb#tran_write_ht#=\verb#true#$, then the Hamiltonian matrix
formatted for the transport calculation will be read from a set of files
described in the parameter `transport_mode`. Set
$\verb#tran_write_ht#=\verb#false#$ to perform automated lcr
calculations (see Section [7.3](#sec:2c2){reference-type="ref"
reference="sec:2c2"}).

The default value is `false`.

### `logical :: tran_use_same_lead`

If $\verb#tran_use_same_lead#=\verb#true#$, then the left and the right
leads are the same. In this case, `seedname_htR.dat` is not required.

The default value is `true`.

### `real(kind=dp) :: tran_group_threshold`

Used to group and sort Wannier functions according to the positions of
their centres. Wannier functions in a group are within
`tran_group_threshold` from one another in `x,y` and `z` directions.
Units are Å

The default is 0.15

### `real(kind=dp) :: translation_centre_frac(3)`

Centre of the unit cell to which the final Wannier centres are
translated. Numbers are in fractional coordinates with respect to the
lattice vectors.

The default value is (0.0,0.0,0.0).

### `logical :: use_ws_distance`

Improves the interpolation of the k-space Hamiltonian, by applying a
translation to each WF by a basis vector of the super-lattice that
minimises the distance between their centres. The translation is
dependent on both WF and on the unit cell vector to which they belong,
i.e., translate function $W_j({\bf r}-{\bf R})$ inside the Wigner-Seitz
cell centred on WF $W_i({\bf r})$.

For a longer explanation, see
Chapter [\[chap:interpolation\]](#chap:interpolation){reference-type="ref"
reference="chap:interpolation"}.

If `false` the code puts all the WF in the home cell, only possible
choice until wannier90 v2.0.1.

The default value is `true` (default changed since v.3.0). Introduced in
v2.1.

### `real(kind=dp) :: ws_distance_tol`

Tolerance when determining whether two values
$\|\mathbf{d}_{ij\mathbf{R}} + \tilde{\mathbf{R}}_{nml} \|$ and
$\|\mathbf{d}_{ij\mathbf{R}} + \tilde{\mathbf{R}}_{n'm'l'} \|$ (as
defined in
chapter [\[chap:interpolation\]](#chap:interpolation){reference-type="ref"
reference="chap:interpolation"}) for the shortest distance between two
Wannier functions are equivalent. If the difference in distance (in
Angstrom) is less than `ws_distance_tol`, they are taken to be
equivalent.

The default value is $10^{-5}$.

### `:: ws_search_size`

Maximum absolute value for the integers $n,m,l$ that identify the
super-lattice vectors $\tilde{\mathbf{R}}_{nml}$ (see
chapter [\[chap:interpolation\]](#chap:interpolation){reference-type="ref"
reference="chap:interpolation"}) when searching for points inside the
Wigner-Seitz cell. If `ws_search_size` is provided as a single integer,
then the number of repetitions of the Born-von Karman cell is the same
along all three linear dimensions; otherwise, if three integers are
provided, the number of repetitions along the $i-$th linear dimension is
`ws_search_size(i)`. The variable is used both in `hamiltonian.F90` and
in `ws_distance.F90`. In the latter case, its value is incremented by
one in order to account for WFs whose centre wanders away from the
original reference unit cell.\
The default value is generally sufficient, but might need to be
increased in case of elongated cells.

The default value is 2.

### `logical :: write_u_matrices`

Write the $\mathbf{U}^{(\mathbf{k})}$ and
$\mathbf{U}^{\mathrm{dis}(\mathbf{k})}$ matrices obtained at the end of
wannierization to files `seedname_u.mat` and `seedname_u_dis.mat`,
respectively.

The default value is `false`.

### `real(kind=dp) :: hr_cutoff`

The absolute value of the smallest matrix element of the Hamiltonian in
the WF basis. If $h_{mn}(\mathbf{R})>\:$` hr_cutoff`, then the matrix
element $h_{mn}(\mathbf{R})$ is retained and used in the band structure
interpolation (when $\verb#bands_plot_mode#=\verb#cut#$) or in the
transport calculation. Otherwise it is deemed to be insignificant and is
discarded. Units are eV.

The default value is 0.0.

### `real(kind=dp) :: dist_cutoff`

The largest distance between two WFs for which the Hamiltonian matrix
element is retained and used in the band interpolation (when
$\verb#bands_plot_mode#=\verb#cut#$) or in the transport calculation.
Units are Å.

The default value is 1000.0.

### `character(len=20) :: dist_cutoff_mode`

Dimension in which the distance between two WFs is calculated. The
vector connecting two WFs may be projected to a line (`one_dim`) or a
plane (`two_dim`). The size of the projected vector is calculated, and
`dist_cutoff` is applied. When `one_dim` or `two_dim` is used,
`one_dim_axis` must be given to specify extended or confined direction.

The valid options for this parameter are:

-   `three_dim` (default)

-   `two_dim`

-   `one_dim`

### `character(len=20) :: one_dim_axis`

Extended direction for a one-dimensional system or confined direction
for a two-dimensional system. This direction must be parallel to one of
the Cartesian axes.

The valid options for this parameter are:

-   `x`

-   `y`

-   `z`

No default.
