# 29: Platinum &#151; Spin Hall conductivity {#platinum-spin-hall-conductivity .unnumbered}

-   Outline: *Calculate spin Hall conductivity (SHC) and plot Berry
    curvature-like term of fcc Pt considering spin-orbit coupling. To
    gain a better understanding of this tutorial, it is suggested to read
    Ref. [@qiao-prb2018] for a detailed description of the theory and
    Ch. 12.5 of the User Guide.*

-   Directory: `tutorials/tutorial29/` *Files can be downloaded [here](https://github.com/wannier-developers/wannier90/tutorials/tutorial29)*

-   Input files

    -    `Pt.scf` *The `pwscf` input file for ground state
        calculation*

    -    `Pt.nscf` *The `pwscf` input file to obtain Bloch
        states on a uniform grid*

    -    `Pt.pw2wan` *The input file for `pw2wannier90`*

    -    `Pt.win` *The `wannier90` and `postw90` input file*

&nbsp;

1.  Run `pwscf` to obtain the ground state of platinum
    
    ```bash title="Terminal"
    pw.x < Pt.scf > scf.out
    ```

2.  Run `pwscf` to obtain the Bloch states on a uniform
    $k$-point grid
 
    ```bash title="Terminal"
    pw.x < Pt.nscf > nscf.out
    ```

3.  Run `wannier90` to generate a list of the required overlaps (written
    into the `Pt.nnkp` file)
    
    ```bash title="Terminal"
    wannier90.x -pp Pt
    ```

4.  Run `pw2wannier90` to compute the overlaps between Bloch states and
    the projections for the starting guess (written in the `Pt.mmn` and
    `Pt.amn` files)

    ```bash title="Terminal"
    pw2wannier90.x < Pt.pw2wan > pw2wan.out
    ```

5.  Run `wannier90` to compute the MLWFs

    ```bash title="Terminal"
    wannier90.x Pt
    ```

6.  Run `postw90`

    ```bash title="Terminal"
    postw90.x Pt # (1)! 
    mpirun -np 8 postw90.x Pt # (2)! 
    ```

    1.     serial execution
    2.     example of parallel execution with 8 MPI processes

## Spin Hall conductivity {#spin-hall-conductivity .unnumbered}

The intrinsic spin Hall conductivity
$\sigma_{\alpha\beta}^{\text{spin}\gamma}$ is proportional to the BZ
integral of the Berry curvature-like term. To evaluate the SHC using a
$25\times
25\times 25$ $k$-point mesh, set the following lines in `Pt.win`,

```vi title="Input file"
berry = true

berry_task = shc

berry_kmesh = 25 25 25
```

When calculating SHC, adaptive smearing can be used by commenting the
following two lines,

```vi title="Input file"
#kubo_adpt_smr = false

#kubo_smr_fixed_en_width = 1
```

Then set the Fermi energy $\varepsilon_F$ to a specific value

```vi title="Input file"
fermi_energy = [insert your value here]
```

or invoke Fermi energy scan by setting

```vi title="Input file"
fermi_energy_min = [insert here your lower range]

fermi_energy_max = [insert here your upper range]

fermi_energy_step = [insert here your step]
```

and re-run `postw90`. The SHC is written in the output file
`Pt-shc-fermiscan.dat`. If only `fermi_energy` is set, the output file
will contain SHC at this specific energy; if a list of Fermi energies
are set, the output file will contain SHC calculated at each energy
point in the list: we call this the "Fermi energy scan" of SHC.

To plot the Fermi energy scan of SHC $\sigma_{xy}^{\text{spin}z}$ versus
$\varepsilon_F$, issue

```bash title="Terminal"
gnuplot
```

```gnuplot title="Gnuplot shell"
plot 'Pt-shc-fermiscan.dat' u 2:3 w lp
```

As a result of the strong and rapid variations of the Berry
curvature-like term across the BZ, the SHC converges rather slowly with
$k$-point sampling, and a $25\times 25\times 25$ kmesh does not yield a
well-converged value.

-   Increase the kmesh density by changing ` berry_kmesh`.

-   To accelerate the convergence, adaptively refine the kmesh around
    spikes in the Berry curvature-like term, by adding to ` Pt.win` the
    lines This adds a $5\times 5\times 5$ fine mesh around those points
    where
    $\vert{\Omega_{\alpha\beta}^{\text{spin}\gamma}}({\bm k})\vert$
    exceeds 100 `berry_curv_unit`. The percentage of points triggering
    adaptive refinement is reported in `Pt.wpout`.

Compare the converged SHC value with those obtained in
Refs. [@qiao-prb2018] and [@guo-prl2008].

Note some rough estimations of computation progress and time are
reported in `Pt.wpout` (see the SHC part of the Solution Booklet). These
may be helpful if the computation time is very long.

## Notes {#notes .unnumbered}

-   Since the Kubo formula of SHC involves unoccupied bands, we need to
    include some unoccupied bands and construct more MLWF. Thus the
    following parameters should be increased accordingly:

-   Normally we calculate the SHC $\sigma_{xy}^{\text{spin}z}$, i.e.
    $\alpha = x, \beta = y, \gamma = z$. To calculate other components,
    the following parameters can be set as `1, 2, 3` with `1, 2, 3`
    standing for `x, y, z` respectively.

## Berry curvature-like term plots {#berry-curvature-like-term-plots .unnumbered}

The band-projected Berry curvature-like term
$\Omega_{n,\alpha\beta}^{\text{spin} \gamma}({\bm k})$ is defined in
this [equation](../../user_guide/postw90/berry#mjx-eqn:eq:kubo_shc) of the User Guide. 
The following lines in `Pt.win` are used
to calculate the energy bands colored by the band-projected Berry
curvature-like term
$\Omega_{n,\alpha\beta}^{\text{spin} \gamma}({\bm k})$ (in Å$^2$), as
well as the $k$-resolved Berry curvature-like term
$\Omega_{\alpha\beta}^{\text{spin} \gamma}({\bm k})$ along high-symmetry
lines in $k$-space, i.e. the `kpath` plot. First comment the line
`berry = true` and then set

```vi title="Input file"
kpath = true

kpath_task = bands+shc

kpath_bands_colour = shc

kpath_num_points = 400

kubo_adpt_smr = false

kubo_smr_fixed_en_width = 1

fermi_energy = [insert your value here]

berry_curv_unit = ang2
```

After executing `postw90`, four files are generated: `Pt-bands.dat`,
`Pt-path.kpt`, `Pt-shc.dat` and `Pt-bands+shc.py`. Then plot the
band-projected Berry curvature-like term
$\Omega_{n,\alpha\beta}^{\text{spin}\gamma}({\bm k})$ using the script
generated at runtime,

```bash title="Terminal"
python Pt-bands+shc.py
```

and compare with Fig. 2 of Ref. [@qiao-prb2018]. Note a large fixed
smearing of 1 eV is used to recover the result in Ref. [@qiao-prb2018].
You can adjust the `kubo_smr_fixed_en_width` as you like to draw a
visually appealing figure. A `kpath` plot of 0.05 eV smearing is shown
in the Solution Booklet.

Besides, you can set `kpath_task = shc` to only draw $k$-resolved term
$\Omega_{\alpha\beta}^{\text{spin} \gamma}({\bm k})$ (the lower panel of
the figure), or set `kpath_task = bands` and `kpath_bands_colour = shc`
to only draw energy bands colored by the band-projected term
$\Omega_{n,\alpha\beta}^{\text{spin} \gamma}({\bm k})$ (the upper panel
of the figure).

Similar to that of AHC, we can get a heatmap plot of the $k$-resolved
Berry curvature-like term
$\Omega_{\alpha\beta}^{\text{spin}\gamma}({\bm k})$, i.e. the `kslice`
plot. To move forward, set `kpath = false` and uncomment the following
lines in `Pt.win`, Note the `kslice_b2` is actually
$(\frac{\sqrt{2}}{4},   \frac{3\sqrt{2}}{4},0.0)$ which leads to a
square slice in the BZ, making it easier to plot in the generated
`python` script. Re-run `postw90`, and issue

```bash title="Terminal"
python Pt-kslice-shc+fermi_lines.py
```

Compare the generated figure with Fig. 3 in Ref. [@qiao-prb2018], or the
Solution Booklet.

## Notes {#notes-1 .unnumbered}

-   Adaptive smearing depends on a uniform kmesh, so when running
    `kpath` and `kslice` plots adaptive smearing should not be used. A
    fixed smearing is needed to avoid near zero number in the
    denominator of the [Kubo formula](../../user_guide/postw90/berry#mjx-eqn:eq:kubo_shc) 
    in the User Guide. To
    add a fixed smearing of 0.05 eV, add the following keywords in the
    `Pt.win`,

## Input parameters for SHC {#input-parameters-for-shc .unnumbered}

Finally, we provide a complete list of input parameters that can be used
to control the SHC calculation, including the calculation of alternating
current (ac) SHC which will be introduced in the next tutorial.

-   general controls for SHC

-   kmesh

-   ac SHC

-   smearing

-   Fermi energy

-   kpath

-   kslice

Their meanings and usages can be found in Ch. 11.5 of the User Guide.


