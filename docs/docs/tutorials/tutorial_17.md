# 17: Iron &#151; Spin-orbit-coupled bands and Fermi-surface contours {#iron-spin-orbit-coupled-bands-and-fermi-surface-contours .unnumbered}

Note: It is recommended that you go through Tutorial 8 first (bcc Fe
without spin-orbit).

Note: This tutorial requires a recent version of the `pw2wannier90`
interface.

-   Outline: *Plot the spin-orbit-coupled bands of ferromagnetic bcc Fe.
    Plot the Fermi-surface contours on a plane in the Brillouin zone.*

-   Directory: `tutorials/tutorial17/` *Files can be downloaded from [here](https://github.com/wannier-developers/wannier90/tutorials/tutorial17)*

-   Input files

    -    `Fe.scf` *The `pwscf` input file for ground state
        calculation*

    -    `Fe.nscf` *The `pwscf` input file to obtain Bloch
        states on a uniform grid*

    -    `Fe.pw2wan` *The input file for `pw2wannier90`*

    -    `Fe.win` *The `wannier90` and `postw90` input file*

Note that `num_wann =18` in `Fe.win`, but only nine trial orbitals are
provided. The line

```vi title="Input file"
spinors = true
```

tells `wannier90` to use in step 3 below the specified trial orbitals on
both the up- and down-spin channels, effectively doubling their number.

1.  Run `pwscf` to obtain the ferromagnetic ground state of
    iron[^2]

    ```bash title="Terminal"
    pw.x < Fe.scf > scf.out
    ```

2.  Run `pwscf` to obtain the Bloch states on a uniform
    k-point grid

    ```bash title="Terminal"
    pw.x < Fe.nscf > nscf.out
    ```

3.  Run `wannier90` to generate a list of the required overlaps (written
    into the `Fe.nnkp` file)

    ```bash title="Terminal"
    wannier90.x -pp Fe
    ```

4.  Run `pw2wannier90` to compute:

    -   The overlaps $\langle u_{n{\bf k}}\vert u_{m{\bf
                k}+{\bf b}}\rangle$ between *spinor* Bloch states
        (written in the `Fe.mmn` file)

    -   The projections for the starting guess (written in the `Fe.amn`
        file)

    -   The spin matrix elements $\langle \psi_{n{\bf
                k}}\vert \sigma_i\vert \psi_{m{\bf k}}\rangle$,
        $i=x,y,z$ (written in the `Fe.spn` file)

    ```bash title="Terminal"
    pw2wannier90.x < Fe.pw2wan > pw2wan.out
    ```

5.  Run `wannier90` to compute the MLWFs.\

    ```bash title="Terminal"
    wannier90.x Fe
    ```

6.  Run `postw90` to compute the energy eigenvalues and spin expectation
    values.

    ```bash title="Terminal"
    postw90.x Fe # (1)! 
    mpirun -np 8 postw90.x Fe # (2)!
    ```

    1.    serial execution
    2.    example of parallel execution with 8 MPI processes

In this tutorial we use the module `kpath` to plot the energy bands
coloured by the expectation value of the spin along \[001\]:

```vi title="Input file"
kpath = true

kpath_task = bands

kpath_bands_colour = spin

kpath_num_points=500
```

To plot the bands using `gnuplot` (version 4.2 or higher) issue

```bash title="Terminal"
gnuplot
```

```gnuplot title="Gnuplot shell"
load 'Fe-bands.gnu'
```

or, using `python`,


```bash title="Terminal"
python Fe-bands.py
```

Next we plot the Fermi-surface contours on the (010) plane $k_y=0$,
using the `kslice` module. Set `kpath = false` and uncomment the
following instructions in `Fe.win`,


```vi title="Input file"
kslice = true

kslice_task = fermi_lines

fermi_energy = [insert your value here]

kslice_corner = 0.0 0.0 0.0

kslice_b1 = 0.5 -0.5 -0.5

kslice_b2 = 0.5 0.5 0.5

kslice_2dkmesh = 200 200
```

taking the Fermi level value from `scf.out`. The energy eigenvalues are
computed on a $200\times 200$ $k$-point grid covering the BZ slice. The
lines of intersection between the Fermi surface and the (010) plane can
be visualized with the `gnuplot` or ` python` scripts generated at
runtime,


```bash title="Terminal"
gnuplot
```

```gnuplot title="Gnuplot shell"
load 'Fe-kslice-fermi_lines.gnu'
```

or


```bash title="Terminal"
python Fe-kslice-fermi_lines.py
```

The Fermi lines can be colour-coded by the spin expectation value
$\langle S_z\rangle$ of the states on the Fermi surface. Add to
` Fe.win` the line


```vi title="Input file"
kslice_fermi_lines_colour = spin
```

and re-run `postw90`. The names of the `gnuplot` and ` python` scripts
generated at runtime are unchanged. (However, the plotting algorithm is
different in this case, and the lines are not as smooth as before. You
may want to increase `kslice_2dkmesh`.)

## Further ideas {#further-ideas-5 .unnumbered}

-   Redraw the Fermi surface contours on the (010) plane starting from a
    calculation without spin-orbit coupling, by adding to the input
    files `iron_{up,down}.win` in Tutorial 8 the lines

    ```vi title="Input file"
    kslice = true
    
    kslice_task = fermi_lines
    
    fermi_energy = \[insert your value here\]
    
    kslice_corner = 0.0 0.0 0.0
    
    kslice_b1 = 0.5 -0.5 -0.5
    
    kslice_b2 = 0.5 0.5 0.5
    
    kslice_2dkmesh = 200 200
    ```

    before running `postw90`,

    ```vi title="Input file"
    postw90.x iron_up
    
    postw90.x iron_dn
    ```

    The `python` scripts generated at runtime draw the up- and down-spin
    Fermi lines on separate figures. To draw them together, use the
    script `iron_updn-kslice-fermi_lines.py` provided with Tutorial 17
    (or merge the two generated scripts). Compare the Fermi lines with
    and without spin-orbit, and note the spin-orbit-induced avoided
    crossings.

-   In Tutorial [8](tutorial_8.md#iron-spin-polarized-wfs-dos-projeced-wfs-versus-mlwfs) we obtained MLWFs separately for the up- and down-spin
    channels of bcc Fe without spin-orbit. The Wannier-interpolated DOS
    was therefore automatically separated into minority and majority
    contributions. For a spinor calculation we can still spin-decompose
    the DOS, using

    ```vi title="Input file"
    dos = true
    
    spin_decomp = true
    
    dos_kmesh = 25 25 25
    ```

    The data file `Fe-dos.dat` created by `postw90` contains the up-spin
    and down-spin contributions in the third and fourth columns,

    ```bash title="Terminal"
    gnuplot
    ```

    ```gnuplot title="Gnuplot shell"
    plot 'Fe-dos.dat' u (-\$3):(\$1-12.6285) w
    l,'Fe-dos.dat' u (\$4):(\$1-12.6285) w l
    ```

    (You should replace 12.6285 with your value of the Fermi energy). An
    alternative approach is to project the DOS onto the up-spin and
    down-spin WFs separately. To find the DOS projected onto the up-spin
    (odd-numbered) WFs replace `spin_decomp = true` with

    ```vi title="Input file"
    dos_project = 1,3,5,7,9,11,13,15,17
    ```

    and re-run `postw90`. This approach has the advantage that it does
    not require the `Fe.spn` file.

[^2]: Please note the following counterintuitive feature in `pwscf`: in order to obtain a ground state with magnetization
along the positive z-axis, one should use a negative value for the variable `starting_magnetization`.
