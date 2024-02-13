
## Utilities

The `wannier90` code is shipped with a few utility programs that may be
useful in some occasions. In this chapter, we describe their use.

### `kmesh.pl`

The `wannier90` code requires the definition of a full Monkhorst--Pack
grid of $k$ points. In the input file the size of this mesh is given by
means of the `mp_grid` variable. E.g., setting

    mp_grid = 4 4 4

tells `wannier90` that we want to use a $4\times 4\times 4$ $k$ grid.

One has then to specify (inside the `kpoints` block in the the
`seedname.win` file) the list of $k$ points of the grid. Here, the
`kmesh.pl` Perl script becomes useful, being able to generate the
required list.

The script can be be found in the `utility` directory of the
`wannier90` distribution. To use it, simply type:

    ./kmesh.pl nx ny nz

where `nx`, `ny` and `nz` define the size of the Monkhorst--Pack grid
that we want to use (for instance, in the above example of the
$4\times 4\times 4$ $k$ grid, `nx`$=$`ny`$=$`nz`$=$`<!-- -->`{=html}4).

This produces on output the list of $k$ points in Quantum Espresso
format, where (apart from a header) the first three columns of each line
are the $k$ coordinates, and the fourth column is the weight of each $k$
point. This list can be used to create the input file for the ab-initio
`nscf` calculation.

If one wants instead to generate the list of the $k$ coordinates without
the weight (in order to copy and paste the output inside the
`seedname.win` file), one simply has to provide a fourth argument on the
command line. For instance, for a $4\times 4\times 4$ $k$ grid, use

    ./kmesh.pl 4 4 4 wannier

and then copy the output inside the in the `kpoints` block in the
`seedname.win` file.

We suggest to always use this utility to generate the $k$ grids. This
allows to provide the $k$ point coordinates with the accuracy required
by `wannier90`, and moreover it makes sure that the $k$ grid used in the
ab-initio code and in `wannier90` are the same.

### $\tt{w90chk2chk.x}$[]{#sec:w90chk2chk label="sec:w90chk2chk"}

During the calculation of the Wannier functions, `wannier90` produces a
`.chk` file that contains some information to restart the calculation.

This file is also required by the `postw90` code. In particular, the
`postw90` code requires at least the `.chk` file, the `.win` input file,
and (almost always) the `.eig` file. Specific modules may require
further files: see the documentation of each module.

However, the `.chk` file is written in a machine-dependent format. If
one wants to run `wannier90` on a machine, and then continue the
calculation with `postw90` on a different machine (or with
`postw90` compiled with a different compiler), the file has to be
converted first in a machine-independent "formatted" format on the first
machine, and then converted back on the second machine.

To this aim, use the `w90chk2chk.x` executable. Note that this
executable is not compiled by default: you can obtain it by executing

    make w90chk2chk

in the main `wannier90` directory.

A typical use is the following:

1.  Calculate the Wannier functions with `wannier90`

2.  At the end of the calculation you will find a `seedname.chk` file.
    Run (in the folder with this file) the command

        w90chk2chk.x -export seedname

    or equivalently

        w90chk2chk.x -u2f seedname

    (replacing `seedname` with the seedname of your calculation).

    This command reads the `seedname.chk` file and creates a formatted
    file `seedname.chk.fmt` that is safe to be transferred between
    different machines.

3.  Copy the `seedname.chk.fmt` file (together with the `seedname.win`
    and `seedname.eig` files) on the machine on which you want to run
    `postw90`.

4.  On this second machine (after having compiled `w90chk2chk.x`) run

        w90chk2chk.x -import seedname

    or equivalently

        w90chk2chk.x -f2u seedname

    This command reads the `seedname.chk.fmt` file and creates an
    unformatted file `seedname.chk` ready to be used by `postw90`.

5.  Run the `postw90` code.

### $\tt{PL\_assessment}$ {#sec:pl_assessment}

The function of this utility is to assess the length of a principal
layer (in the context of a Landauer-Buttiker quantum conductance
calculation) of a periodic system using a calculation on a single unit
cell with a dense k-point mesh.

The utility requires the real-space Hamiltonian in the MLWF basis,
`seedname_hr.dat`.

The `seedname_hr.dat` file should be copied to a directory containing
executable for the utility. Within that directory, run:

    \$> ./PL_assess.x  nk1 nk2 nk3 num_wann 

where:

`nk1` is the number of k-points in x-direction `nk2` is the number of
k-points in y-direction `nk3` is the number of k-points in z-direction
`num_wann` is the number of wannier functions of your system

e.g.,

    \$> ./PL_assess.x  1 1 20 16

Note that the current implementation only allows for a single k-point in
the direction transverse to the transport direction.

When prompted, enter the seedname.

The programme will return an output file `seedname_pl.dat`, containing
four columns

1.  Unit cell number, $R$

2.  Average 'on-site' matrix element between MLWFs in the home unit
    cell, and the unit cell $R$ lattice vectors away

3.  Standard devaition of the quantity in (2)

4.  Maximum absolute value in (2)

### $\tt{w90vdw}$ {#sec:w90vdw}

This utility provides an implementation of a method for calculating van
der Waals energies based on the idea of density decomposition via MLWFs.

For theoretical details, please see the following publication and
references therein:

Lampros Andrinopoulos, Nicholas D. M. Hine and Arash A. Mostofi,
"Calculating dispersion interactions using maximally localized Wannier
functions", *J. Chem. Phys.* **135**, 154105 (2011).

For further details of this program, please see the documentation in
`utility/w90vdw/doc/` and the related examples in
`utility/w90vdw/examples/`.

### $\tt{w90pov}$ {#sec:w90pov}

An utility to create Pov files (to render the Wannier functions using
the PovRay utility) is provided inside `utility/w90pov`.

Please refer to the documentation inside `utility/w90pov/doc` for more
information.

### $\tt{k\_mapper.py}$ {#sec:k_mapper}

The `wannier90` code requires the definition of a full Monkhorst--Pack
grid of $\mathbf{k}$-vectors, which can be obtained by means of the
`kmesh.pl` utility. In order to perform a GW calculation with the Yambo
code, you need to perform a nscf calculation on a grid in the
irreducible BZ. Moreover, you may need a finer grid to converge the GW
calculation than what you need to interpolate the band structure. The
`k_mapper.py` tools helps in finding the $\mathbf{k}$-vectors indexes of
a full grid needed for interpolation into the reduced grid needed for
the GW calculation with Yambo. Usage:

`path/k_mapper.py nx ny nz QE_nscf_output`

where `path` is the path of `utility` folder, `QE_nscf_output` is the
path of the QE nscf output file given to Yambo.

### $\tt{gw2wannier90.py}$

This utility allows to sort in energy the input data of `wannier90`
(e.g. overlap matrices and energy eigenvalues). `gw2wannier90.py` allows
to use `wannier90` at the $G_0W_0$ level, where quasi-particle
corrections can change the energy ordering of eigenvalues (Some
`wannier90` modules require states to be ordered in energy). Usage:

`path/gw2wannier90.py seedname options` where `path` is the path of
`utility` folder.

Available options are:

    mmn, amn, spn, unk, uhu, uiu,
    spn_formatted, unk_formatted, uhu_formatted, uiu_formatted,
    write_formatted

If no options are specified, all the files
(`mmn, amn, spn, UNK, uHu, uIu`) are considered.

Binary (unformatted Fortran) files are supported, though not
reccommended, since they are compiler-dependent. A safer choice is to
use (bigger) formatted files, with options:

`spn_formatted, uiu_formatted, uhu_formatted, unk_formatted`

In default, the output format is the same as the input format. To
generate formatted files with unformatted input, use option:
`write_formatted` []{#sec:w90aaa label="sec:w90aaa"}

### $\tt{w90spn2spn.x}$[]{#sec:w90spn2spn label="sec:w90spn2spn"}

The interface between ab-initio code and `wannier90` (e.g.
`pw2wannier90.x`) can produce a `.spn` file that is used by `postw90` to
calculate some spin related quantities.

The `.spn` file can be written in a machine-dependent or a
machine-independent format depending on the input parameter
`spn_formatted` (the default is `false` which means the `.spn` file is
machine-dependent) of the `pw2wannier90.x`. If a `.spn` file has been
generated on a machine with machine-dependent format, and then one wants
to continue the calculation with `postw90` on a different machine (or
with `postw90` compiled with a different compiler), the file has to be
converted first in a machine-independent "formatted" format on the first
machine.

To this aim, use the `w90spn2spn.x` executable. Note that this
executable is not compiled by default: you can obtain it by executing

    make w90spn2spn

in the main `wannier90` directory.

A typical use is the following:

1.  Calculate the `.spn` file, e.g. by `pw2wannier90.x`

2.  At the end of the calculation you will find a `seedname.spn` file.
    If the file is "unformatted", run (in the folder with this file) the
    command

        	w90spn2spn.x -export seedname
        	

    or equivalently

        	w90spn2spn.x -u2f seedname
        	

    (replacing `seedname` with the seedname of your calculation).

    This command reads the `seedname.spn` file and creates a formatted
    file `seedname.spn.fmt` that is safe to be transferred between
    different machines.

3.  Copy the `seedname.spn.fmt` file on the machine on which you want to
    run `postw90`.

4.  On this second machine (after having compiled `w90spn2spn.x`) run

        	w90spn2spn.x -import seedname
        	

    or equivalently

        	w90spn2spn.x -f2u seedname
        	

    This command reads the `seedname.spn.fmt` file and creates an
    unformatted file `seedname.spn` ready to be used by `postw90`.

5.  Run the `postw90` code.

Note if `spn_formatted` is set to `true` in both `pw2wannier90.x` and
`postw90` input files, then the `.spn` file will be written and read as
"formatted", so `w90spn2spn.x` is not needed. However, if an
"unformatted" `seedname.spn` has been created and you do not want to
rerun `pw2wannier90.x`, then `w90spn2spn.x` can be useful. Also, once a
"formatted" `seedname.spn` has been generated, the step 4 can be skipped
if `spn_formatted` is set to `true` in `postw90` input file
`seedname.win`.

### `write_pdwf_projectors.py`

A python script to extract projectors from a `UPF` file and write them
into a `pw2wannier90.x`-recognizable `dat` file, which can be used to
compute `amn` using pseudo-atomic orbital projection.

Usage:

`path/write_pdwf_projectors.py UPF_filename`

where `path` is the path of `utility` folder, `UPF_filename` is the path
of a `UPF` file.

The script serves as a reference for writing the `dat` file, you can
generate your own pseudo-atomic orbitals by some other codes and use
them to compute `amn`.
