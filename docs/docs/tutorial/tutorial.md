# 26: Gallium Arsenide -- Selective localization and constrained centres {#gallium-arsenide-selective-localization-and-constrained-centres .unnumbered}

-   Outline: *Application of the selectively localised Wannier function
    (SLWF) method to gallium arsenide (GaAs), following the example in
    Ref. [@Marianetti], which is essential reading for this tutorial
    example.*

-   Directory: `examples/example26/`

-   Input files:

    -   `GaAs.scf` *The `PWSCF` input file for ground state calculation*

    -   `GaAs.nscf` *The `PWSCF` input file to obtain Bloch states on a
        uniform grid*

    -   `GaAs.pw2wan` *The input file for* `pw2wannier90`

    -   `GaAs.win` *The* `wannier90` *and* `postw90` *input file*

    1.  Run `PWSCF` to obtain the ground state of Gallium Arsenide

        `pw.x < GaAs.scf > scf.out`

    2.  Run `PWSCF` to obtain the ground state of Gallium Arsenide

        `pw.x < GaAs.nscf > nscf.out`

    3.  Run `Wannier90` to generate a list of the required overlaps
        (written into the `GaAs.nnkp` file)

        `wannier90.x -pp GaAs`

    4.  Run `pw2wannier90` to compute:

        -   The overlaps $\langle u_{n\bf{k}}|u_{n\bf{k+b}}\rangle$
            between Bloch states (written in the `GaAs.mmn` file)

        -   The projections for the starting guess (written in the
            `GaAs.amn` file)

        `pw2wannier90.x < GaAs.pw2wan > pw2wan.out`

    5.  Inspect the `.win` file.

        -   Make sure you understand the new keywords corresponding to
            the selective localisation algorithm.

        -   Run `wannier90` to compute the SLWFs, in this case using one
            objective Wannier function.

        `wannier90.x GaAs`

    To constrain the centre of the SLWF you need to add
    `slwf_constrain = true` and\
    `slwf_lambda = 1` to the input file and uncomment the `slwf_centres`
    block. This will add a penalty functional to the total spread, which
    will try to constrain the centre of the SLWF to be on the As atom
    (as explained in Ref. [@Marianetti], particularly from Eq. 24 to
    Eq. 35).

    Look at the value of the penalty functional, is this what you would
    expect at convergence? Does the chosen value of the Lagrange
    multiplier `slwf_lambda` give a SLWF function centred on the As
    atom?

    Alternatively, you can modify the `slwf_centres` block to constrain
    the centre of the SLWF to be on the Ga atom. Do you need a different
    value of `slwf_lambda` in this case to converge? Take a look at the
    result in Vesta and explain what you see. Do these functions
    transform like the identity under the action of the $T_d$ group?

# 27: Silicon -- Selected columns of density matrix algorithm for automated MLWFs {#silicon-selected-columns-of-density-matrix-algorithm-for-automated-mlwfs .unnumbered}

Note: This example requires a recent version of the `pw2wannier90.x`
post-processing code of `Quantum ESPRESSO` (v6.4 or above).

-   Outline: *For bulk crystalline Silicon, generate the $A_{mn}$
    matrices via the selected columns of density matrix (SCDM) algorithm
    and the corresponding MLWFs for 1) Valence bands 2) Valence bands
    and 4 low-lying conduction bands 3) Conduction bands only. To better
    understand the input files and the results of these calculations, it
    is crucial that the Reader has familiarized with the concepts and
    methods explained in Ref. [@LinLin-ArXiv2017]. More info on the
    keywords related to the SCDM method may be found in the user_guide.*

-   Directory: `examples/example27/`

-   Input Files: `input_files`, and in the three subfolders
    `isolated, erfc` and `gaussian`. The `input_files` folder contains:

    -    `si.scf` *The [pwscf]{.smallcaps} input file for the ground
        state calculation*

    -    `si_4bands.nscf ` *The [pwscf]{.smallcaps} input file to obtain
        Bloch states on a uniform grid for 4 bands.*

    -    `si_12bands.nscf ` *The [pwscf]{.smallcaps} input file to
        obtain Bloch states on a uniform grid for 12 bands.*

-    Whereas the three subfolders `isolated, erfc` and `gaussian`
    contain the `si.win` `wannier90`  input files and `si.pw2wan`
    `pw2wannier90` input files each corresponding to one of the
    scenarios listed in the outline.

```{=html}
<!-- -->
```
-   Valence bands: In this case we will compute 4 localized WFs
    corresponding to the 4 valence bands of Silicon. These 4 bands
    constitute a manifold that is separated in energy from other bands.
    In this case the columns of the density matrix are already localized
    in real space and no extra parameter is required.

    1.  Copy the input files `si.scf` and `si_4bands.nscf` from the
        `input_files` directory into the `isolated` folder

    2.  Run [pwscf]{.smallcaps} to obtain the ground state charge of
        bulk Silicon.\
        `pw.x < si.scf > scf.out`

    3.  Run [pwscf]{.smallcaps} to obtain the Bloch states on a uniform
        k-point grid of 4x4x4 for 4 bands.\
        `pw.x < si_4bands.nscf > nscf.out`

    4.  Inspect the `si.win` input file and make sure that the
        `auto_projections` flag is set to `.true.`. Also, make sure that
        no projections block is present.

    5.  Run `wannier90` to generate a list of the required overlaps and
        also info on the SCDM method (written into the `si.nnkp` file).\
        `wannier90.x -pp si`

    6.  Inspect the `si.nnkp` file and make sure you find the
        `auto_projections` block and that no projections have been
        written in the `projections` block.

    7.  Inspect the `.pw2wan` input file. You will find two new
        keywords, i.e. `scdm_proj` and `scdm_entanglement`. The former,
        will instruct `pw2wannier90.x` to use the SCDM method when
        generating the $A_{mn}$ matrix. The latter, defines which
        formula to adopt for the function $f(\varepsilon_{n\mathbf{k}})$
        (see [@LinLin-ArXiv2017] and point below).

    8.  Run `pw2wannier90` to compute the overlap between Bloch states
        and the projections via the SCDM method (written in the `si.mmn`
        and `si.amn` respectively).\
        `pw2wannier90.x < si.pw2wan > pw2wan.out`

    9.  Run `wannier90` to compute the MLWFs.\
        `wannier90.x si`\
        At this point, you should have obtained 4 Wannier functions and
        the interpolated valence bands for Silicon. Inspect the output
        file `si.wout`. In particular, look at the geometric centres of
        each WF, do they lie at the centre of the Si-Si bond as for the
        MLWFs computed from user-defined initial $s$-like projections
        (see Example11)? Plot these WFs using Vesta. Do they show the
        $\sigma$ character one would expect from chemical arguments?

-   Valence bands + conduction bands: In this case we will compute 8
    localized WFs corresponding to the 4 valence bands and 4 low-lying
    conduction bands. Here, we don't have a separate manifold, since the
    conduction bands are entangled with other high-energy bands and the
    columns of the density matrix are not exponentially localized by
    construction. A modified density matrix is required in this
    case[@LinLin-ArXiv2017], and it is defined as:
    $$P(\mathbf{r},\mathbf{r}') = \sum_{n,\mathbf{k}} \psi_{n\mathbf{k}}(\mathbf{r})f(\varepsilon_{n,\mathbf{k}})\psi_{n\mathbf{k}}^\ast(\mathbf{r}'),$$
    where $\psi_{n\mathbf{k}}$ and $\varepsilon_{n,\mathbf{k}}$ are the
    energy eigestates and eigenvalues from the first-principle
    calculation respectively. The function
    $f(\varepsilon_{n,\mathbf{k}})$ contains two free parameters $\mu$
    and $\sigma$ and is defined as a complementary error function:
    $$f(\varepsilon_{n,\mathbf{k}}) = \frac{1}{2}\mathrm{erfc}\left(\frac{\varepsilon_{n,\mathbf{k}} - \mu}{\sigma}\right).$$

    1.  Copy the input files `si.scf` and `si_12bands.nscf` from the
        `input_files` folder into the `erfc` folder

    2.  Run [pwscf]{.smallcaps} to obtain the ground state charge of
        bulk Silicon.\
        `pw.x < si.scf > scf.out`

    3.  Run [pwscf]{.smallcaps} to obtain the Bloch states on a uniform
        k-point grid of 4x4x4 for 12 bands this time.\
        `pw.x < si_12bands.nscf > nscf.out`

    4.  Inspect the `si.win` input file and make sure that the
        `auto_projections` flag is set to `.true.`. Also, make sure that
        no projection block is present.

    5.  Run `wannier90` to generate a list of the required overlaps and
        also info on the SCDM method (written into the `si.nnkp` file).\
        `wannier90.x -pp si`

    6.  Inspect the `si.nnkp` file and make sure you find the
        `auto_projections` block and that no projections have been
        written in the `projections` block.

    7.  Inspect the `.pw2wan` input file. You will find other two new
        keywords, i.e. `scdm_mu` and `scdm_sigma`. These are the values
        in eV of $\mu$ and $\sigma$ in $f(\varepsilon_{n,\mathbf{k}})$,
        respectively.

    8.  Run `pw2wannier90` to compute the overlap between Bloch states
        and the projections via the SCDM method (written in the `si.mmn`
        and `si.amn` respectively).\
        `pw2wannier90.x < si.pw2wan > pw2wan.out`

    9.  Run `wannier90` to compute the MLWFs.\
        `wannier90.x si`\
        At this point, you should have obtained 8 localized Wannier
        functions and the interpolated valence and conduction bands for
        Silicon. Again, compare the results for the geometric centres
        and the individual spreads with the ones from Example11. Is the
        final value of total spread bigger or smaller than the one from
        Example11? Look at the WFs with Vesta. Can you explain what you
        see? Where do the major lobes of the $sp3$-like WFs point in
        this case?

-   Conduction bands only: In this case we will compute 4 localized WFs
    corresponding to the 4 low-lying conduction bands only. As for the
    previous point, we need to define a modified density
    matrix[@LinLin-ArXiv2017]. Since we are only interested in a subset
    of the conduction states, within a bounded energy region, a good
    choice for $f(\varepsilon_{n,\mathbf{k}})$ is:
    $$f(\varepsilon_{n,\mathbf{k}}) = \exp\left(-\frac{(\varepsilon_{n,\mathbf{k}} - \mu)^2}{\sigma^2}\right).$$

    1.  Copy the input files `si.scf` and `si_12bands.nscf` from the
        `input_files` directory into the `gaussian` folder

    2.  Run [pwscf]{.smallcaps} to obtain the ground state charge of
        bulk Silicon.\
        `pw.x < si.scf > scf.out`

    3.  Run [pwscf]{.smallcaps} to obtain the Bloch states on a uniform
        k-point grid of 4x4x4 for 12 bands this time.\
        `pw.x < si_12bands.nscf > nscf.out`

    4.  Inspect the `si.win` input file and make sure that the
        `auto_projections` flag is set to `.true.`. Also, make sure that
        no projections block is present.

    5.  Run `wannier90` to generate a list of the required overlaps and
        also info on the SCDM method (written into the `si.nnkp` file).\
        `wannier90.x -pp si`

    6.  Inspect the `si.nnkp` file and make sure you find the
        `auto_projections` block and that no projections have been
        written in the `projections` block.

    7.  Run `pw2wannier90` to compute the overlap between Bloch states,
        the projections for the starting guess via the SCDM method
        (written in the `si.mmn` and `si.amn` respectively).\
        `pw2wannier90.x < si.pw2wan > pw2wan.out`

    8.  Run `wannier90` to compute the MLWFs.\
        `wannier90.x si`\
        At this point, you should have obtained 4 localized Wannier
        functions and the interpolated conduction bands for Silicon.
        From chemical intuition, we would expect these functions to be
        similar to anti-bonding orbitals of molecules with tetrahedral
        symmetry. Plot the WFs and check if this is confirmed.

# 28: Diamond -- plotting of MLWFs using Gaussian cube format and VESTA {#diamond-plotting-of-mlwfs-using-gaussian-cube-format-and-vesta .unnumbered}

-   Outline: *Obtain MLWFs for the valence bands of diamond and output
    them in Gaussian cube format*

-   Directory: `examples/example28/` The input files for this examples
    are the same as the ones in example05

-   Input Files

    -    `diamond.scf` *The [pwscf]{.smallcaps} input file for ground
        state calculation*

    -    `diamond.nscf` *The [pwscf]{.smallcaps} input file to obtain
        Bloch states on a uniform grid*

    -    `diamond.pw2wan` *The input file for `pw2wannier90`*

    -    `diamond.win` *The `wannier90` input file*

1.  Run [pwscf]{.smallcaps} to obtain the ground state of diamond\
    `pw.x < diamond.scf > scf.out`

2.  Run [pwscf]{.smallcaps} to obtain the Bloch states on a uniform
    k-point grid\
    `pw.x < diamond.nscf > nscf.out`

3.  Run `wannier90` to generate a list of the required overlaps (written
    into the `diamond.nnkp` file).\
    `wannier90.x -pp diamond`

4.  Run `pw2wannier90` to compute the overlap between Bloch states and
    the projections for the starting guess (written in the `diamond.mmn`
    and `diamond.amn` files).\
    `pw2wannier90.x < diamond.pw2wan > pw2wan.out`

5.  When the lattice vectors are non-orthogonal, not all the
    visualisation programs are capable to plot volumetric data in the
    Gaussian cube format. One program that can read volumetric data for
    these systems is VESTA. To instruct `wannier90` to output the MLWFs
    data in Gaussian cube format you need to add the following lines to
    the `.win` file

        wannier_plot           = .true.
        wannier_plot_supercell = 3
        wannier_plot_format    = cube
        wannier_plot_mode      = crystal
        wannier_plot_radius    = 2.5
        wannier_plot_scale     = 1.0

    Run `wannier90` to compute the MLWFs and output them in the Gaussian
    cube file.\
    `wannier90.x diamond`

6.  Plot the first MLWF with VESTA `vesta diamond_00001.cube`

Extra: Instead of using `wannier_plot_mode = crystal` try to use the
molecule mode as `wannier_plot_mode = molecule` (see the user guide for
the definition of this keyword). Add the following line to the `.win`
file:

    restart = plot

and re-run `wannier90`. Use VESTA to plot the resulting MLWFs, do you
see any difference from the `crystal` mode case? Can you explain why?
Try to change the size of the supercell from 3 to 5, do you expect the
results to be different? (*Hint:* When using the Gaussian cube format
the code outputs the WF on a grid that is smaller than the super
unit-cell. The size of the grid is specified by `wannier_plot_scale` and
`wannier_plot_radius`.)

# 29: Platinum -- Spin Hall conductivity {#platinum-spin-hall-conductivity .unnumbered}

-   Outline: *Calculate spin Hall conductivity (SHC) and plot Berry
    curvature-like term of fcc Pt considering spin-orbit coupling. To
    gain a better understanding of this example, it is suggested to read
    Ref. [@qiao-prb2018] for a detailed description of the theory and
    Ch. 12.5 of the User Guide.*

-   Directory: `examples/example29/`

-   Input files

    -    `Pt.scf` *The [pwscf]{.smallcaps} input file for ground state
        calculation*

    -    `Pt.nscf` *The [pwscf]{.smallcaps} input file to obtain Bloch
        states on a uniform grid*

    -    `Pt.pw2wan` *The input file for `pw2wannier90`*

    -    `Pt.win` *The `wannier90` and `postw90` input file*

1.  Run [pwscf]{.smallcaps} to obtain the ground state of platinum\
    `pw.x < Pt.scf > scf.out`

2.  Run [pwscf]{.smallcaps} to obtain the Bloch states on a uniform
    $k$-point grid\
    `pw.x < Pt.nscf > nscf.out`

3.  Run `wannier90` to generate a list of the required overlaps (written
    into the `Pt.nnkp` file)\
    `wannier90.x -pp Pt`

4.  Run `pw2wannier90` to compute the overlaps between Bloch states and
    the projections for the starting guess (written in the `Pt.mmn` and
    `Pt.amn` files)\
    `pw2wannier90.x < Pt.pw2wan > pw2wan.out`

5.  Run `wannier90` to compute the MLWFs\
    `wannier90.x Pt`

6.  Run `postw90` \
    `postw90.x Pt` (serial execution)\
    `mpirun -np 8 postw90.x Pt` (example of parallel execution with 8
    MPI processes)

## Spin Hall conductivity {#spin-hall-conductivity .unnumbered}

The intrinsic spin Hall conductivity
$\sigma_{\alpha\beta}^{\text{spin}\gamma}$ is proportional to the BZ
integral of the Berry curvature-like term. To evaluate the SHC using a
$25\times
25\times 25$ $k$-point mesh, set the following lines in `Pt.win`,

` `

> berry = true
>
> berry_task = shc
>
> berry_kmesh = 25 25 25

When calculating SHC, adaptive smearing can be used by commenting the
following two lines,

` `

> #kubo_adpt_smr = false
>
> #kubo_smr_fixed_en_width = 1

Then set the Fermi energy $\varepsilon_F$ to a specific value

` `

> fermi_energy = \[insert your value here\]

or invoke Fermi energy scan by setting

` `

> fermi_energy_min = \[insert here your lower range\]
>
> fermi_energy_max = \[insert here your upper range\]
>
> fermi_energy_step = \[insert here your step\]

and re-run `postw90`. The SHC is written in the output file
`Pt-shc-fermiscan.dat`. If only `fermi_energy` is set, the output file
will contain SHC at this specific energy; if a list of Fermi energies
are set, the output file will contain SHC calculated at each energy
point in the list: we call this the "Fermi energy scan" of SHC.

To plot the Fermi energy scan of SHC $\sigma_{xy}^{\text{spin}z}$ versus
$\varepsilon_F$, issue

` `

> myshell\> gnuplot
>
> gnuplot\> plot 'Pt-shc-fermiscan.dat' u 2:3 w lp

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
Eq. (12.22) of the User Guide. The following lines in `Pt.win` are used
to calculate the energy bands colored by the band-projected Berry
curvature-like term
$\Omega_{n,\alpha\beta}^{\text{spin} \gamma}({\bm k})$ (in Å$^2$), as
well as the $k$-resolved Berry curvature-like term
$\Omega_{\alpha\beta}^{\text{spin} \gamma}({\bm k})$ along high-symmetry
lines in $k$-space, i.e. the `kpath` plot. First comment the line
`berry = true` and then set

` `

> kpath = true
>
> kpath_task = bands+shc
>
> kpath_bands_colour = shc
>
> kpath_num_points = 400
>
> kubo_adpt_smr = false
>
> kubo_smr_fixed_en_width = 1
>
> fermi_energy = \[insert your value here\]
>
> berry_curv_unit = ang2

After executing `postw90`, four files are generated: `Pt-bands.dat`,
`Pt-path.kpt`, `Pt-shc.dat` and `Pt-bands+shc.py`. Then plot the
band-projected Berry curvature-like term
$\Omega_{n,\alpha\beta}^{\text{spin}\gamma}({\bm k})$ using the script
generated at runtime,

` `

> myshell\> python Pt-bands+shc.py

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

` `

> myshell\> python Pt-kslice-shc+fermi_lines.py

Compare the generated figure with Fig. 3 in Ref. [@qiao-prb2018], or the
Solution Booklet.

## Notes {#notes-1 .unnumbered}

-   Adaptive smearing depends on a uniform kmesh, so when running
    `kpath` and `kslice` plots adaptive smearing should not be used. A
    fixed smearing is needed to avoid near zero number in the
    denominator of the Kubo formula, Eq. (12.22) in the User Guide. To
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

# 30: Gallium Arsenide -- Frequency-dependent spin Hall conductivity {#gallium-arsenide-frequency-dependent-spin-hall-conductivity .unnumbered}

-   Outline: *Calculate the alternating current (ac) spin Hall
    conductivity of gallium arsenide considering spin-orbit coupling. To
    gain a better understanding of this example, it is suggested to read
    Ref. [@qiao-prb2018] for a detailed description of the theory and
    Ch. 12.5 of the User Guide.*

-   Directory: `examples/example30/`

-   Input files

    -    `GaAs.scf` *The [pwscf]{.smallcaps} input file for ground state
        calculation*

    -    `GaAs.nscf` *The [pwscf]{.smallcaps} input file to obtain Bloch
        states on a uniform grid*

    -    `GaAs.pw2wan` *The input file for `pw2wannier90`*

    -    `GaAs.win` *The `wannier90` and `postw90` input file*

1.  Run [pwscf]{.smallcaps} to obtain the ground state of gallium
    arsenide\
    `pw.x < GaAs.scf > scf.out`

2.  Run [pwscf]{.smallcaps} to obtain the Bloch states on a uniform
    $k$-point grid\
    `pw.x < GaAs.nscf > nscf.out`

3.  Run `wannier90` to generate a list of the required overlaps (written
    into the `GaAs.nnkp` file)\
    `wannier90.x -pp GaAs`

4.  Run `pw2wannier90` to compute the overlaps between Bloch states and
    the projections for the starting guess (written in the `GaAs.mmn`
    and `GaAs.amn` files)\
    `pw2wannier90.x < GaAs.pw2wan > pw2wan.out`

5.  Run `wannier90` to compute the MLWFs\
    `wannier90.x GaAs`

6.  Run `postw90` \
    `postw90.x GaAs` (serial execution)\
    `mpirun -np 8 postw90.x GaAs` (example of parallel execution with 8
    MPI processes)

## ac spin Hall conductivity {#ac-spin-hall-conductivity .unnumbered}

The spin Hall conductivity is also dependent on the frequency $\omega$
in the Eq. (12.22) of the User Guide. The direct current (dc) SHC
calculated in the previous example corresponds to
$\sigma_{\alpha\beta}^{\text{spin}\gamma}$ in the limit
$\omega\rightarrow
0$ and it is a real number. At finite frequency
$\sigma_{\alpha\beta}^{\text{spin}\gamma}$ acquires an imaginary part.

To compute the ac spin Hall conductivity for $\hbar\omega$ up to 8 eV,
add the lines

` `

> shc_freq_scan = true
>
> kubo_freq_min = 0.0
>
> kubo_freq_max = 8.0
>
> kubo_freq_step = 0.01

and re-run `postw90`. The file `GaAs-shc-freqscan.dat` contains the
calculated ac SHC. Reasonably converged spectra can be obtained with a
$250\times 250\times 250$ $k$-point mesh. To plot the ac SHC, issue the
following commands

` `

> myshell\> gnuplot
>
> gnuplot\> plot 'GaAs-shc-freqscan.dat' u 2:3 w l title 'Re',
> 'GaAs-shc-freqscan.dat' u 2:4 w l title 'Im'

and then compare the result with Fig. 4 in Ref. [@qiao-prb2018] or the
Solution Booklet.

## Notes {#notes-2 .unnumbered}

-   When calculating ac SHC, adaptive smearing can be used by add the
    following keywords in the `GaAs.win`,

-   Adaptive kmesh refinement is not implemented for ac SHC calculation.

-   The first 10 semi-core states are excluded from the calculation by
    using the following keywords and in the case of GaAs disentanglement
    is not adopted so

-   Since the band gap is often under estimated by LDA/GGA calculations,
    a scissors shift is applied to recover the experimental band gap by
    using the following keywords or by

# 31: Platinum -- Selected columns of density matrix algorithm for spinor wavefunctions {#platinum-selected-columns-of-density-matrix-algorithm-for-spinor-wavefunctions .unnumbered}

Note: This example requires a recent version of the `pw2wannier90.x`
post-processing code of `Quantum ESPRESSO` (v6.3 or above).

-   Outline: *For bulk crystalline platinum with spin-orbit coupling,
    generate the $A_{mn}$ matrices via the selected columns of density
    matrix (SCDM) algorithm and the corresponding spinor-MLWFs. To
    better understand the input files and the results of these
    calculations, it is crucial that the Reader has familiarized with
    the concepts and methods explained in Ref. [@LinLin-ArXiv2017]. More
    info on the keywords related to the SCDM method may be found in the
    user_guide.*

    This example focuses on the use of the SCDM method for
    spin-noncollinear systems. For the overview of the use of SCDM
    method to spinless systems, please refer to example27.

-   Directory: `examples/example31/`

    The input files for this examples are similar to the ones in example
    29, except that a coarser k-point grid is used and that the keywords
    related to `postw90.x` are removed.

-   Input Files:

    -   `Pt.scf` *The [pwscf]{.smallcaps} input file for the ground
        state calculation*

    -   `Pt.nscf` *The [pwscf]{.smallcaps} input file to obtain Bloch
        states on a uniform grid*

    -   `Pt.pw2wan` *The input file for `pw2wannier90` with keywords
        related to the SCDM method*

    -   `Pt.win` *The `wannier90` input file*

We will compute 18 localized WFs. Since the band structure of platinum
is metallic, the low-lying bands are entangled with other high-energy
bands, and the columns of the density matrix are not exponentially
localized by construction. Thus, we use a modified density matrix
[@LinLin-ArXiv2017], with the function $f(\varepsilon_{n,\mathbf{k}})$
defined as a complementary error function. Refer to example 27 for the
definition of the modified density matrix and the functional form of
$f(\varepsilon_{n,\mathbf{k}})$.

1.  Run [pwscf]{.smallcaps} to obtain the ground state of platinum\
    `pw.x < Pt.scf > scf.out`

2.  Run [pwscf]{.smallcaps} to obtain the Bloch states on a uniform
    $7\times 7\times 7$ $k$-point grid\
    `pw.x < Pt.nscf > nscf.out`

3.  Inspect the `Pt.win` input file and make sure that the
    `auto_projections` flag is set to `.true.`. Also, make sure that no
    projection block is present.

4.  Run `wannier90` to generate a list of the required overlaps (written
    into the `Pt.nnkp` file)\
    `wannier90.x -pp Pt`

5.  Inspect the `Pt.nnkp` file and make sure you find the
    `auto_projections` block and that no projections have been written
    in the `projections` block.

6.  Inspect the `Pt.pw2wan` input file. You will find four SCDM-related
    keywords: `scdm_proj`, `scdm_entanglement`, `scdm_mu` and
    `scdm_sigma`. In particular, the keyword `scdm_proj` will instruct
    `pw2wannier90.x` to use the SCDM method when generating the $A_{mn}$
    matrix. The remaining three keywords defines the formula and
    parameters to define the function $f(\varepsilon_{n\mathbf{k}})$
    (see Ref. [@LinLin-ArXiv2017] and example 27).

7.  Run `pw2wannier90` to compute the overlap between Bloch states and
    the projections via the SCDM method (written in the `Pt.mmn` and
    `Pt.amn` respectively).\
    `pw2wannier90.x < Pt.pw2wan > pw2wan.out`

8.  Inspect the `pw2wan.out` output file. Compared to the spinless case,
    you will find the following two additional lines.

             Number of pivot points with spin up  :     9
             Number of pivot points with spin down:     9

    These lines give information on the pivots obtained by the QR
    decomposition with column pivoting (QRCP) in the SCDM algorithm.
    Each pivot determines a point in the real-space grid and a spin
    state. The basis of the spin state is determined by the basis used
    in the electronic structure code. In [pwscf]{.smallcaps}, the basis
    states are spin up and down states along the Cartesian $z$-axis.

9.  Run `wannier90` to compute the MLWFs\
    `wannier90.x Pt`\

# 32: Tungsten --- SCDM parameters from projectability {#tungsten-scdm-parameters-from-projectability .unnumbered}

-   Outline: *Compute the Wannier interpolated band structure of
    tungsten (W) using the SCDM method to calculate the initial guess
    (see Example 27 for more details). The free parameters in the SCDM
    method, i.e., $\mu$ and $\sigma$, are obtained by fitting a
    complementary error function to the projectabilities. The number of
    MLWFs is given by the number of pseudo-atomic orbitals (PAOs) in the
    pseudopotential, $21$ in this case. All the steps shown in this
    example have been automated in the AiiDA[@Pizzi_AiiDA] workflow that
    can be downloaded from the MaterialsCloud
    website[@MaterialsCloudArchiveEntry]*.

-   Directory: `examples/example31/`

-   Input files

    -    `W.scf` *The [pwscf]{.smallcaps} input file for ground state
        calculation*

    -    `W.nscf` *The [pwscf]{.smallcaps} input file to obtain Bloch
        states on a uniform grid*

    -    `W.pw2wan` *The input file for `pw2wannier90`*

    -    `W.proj` *The input file for `projwfc`*

    -    `generate_weights.sh` *The bash script to extract the
        projectabilities from the output of `projwfc`*

    -    `W.win` *The `wannier90` input file*

1.  Run [pwscf]{.smallcaps} to obtain the ground state of tungsten\
    `pw.x -in W.scf > scf.out`

2.  Run [pwscf]{.smallcaps} to obtain the Bloch states on a
    $10\times10\times10$ uniform $k$-point grid\
    `pw.x -in W.nscf > nscf.out`

3.  Run `wannier90` to generate a list of the required overlaps (written
    into the `W.nnkp` file)\
    `wannier90.x -pp W`

4.  Run `projwfc` to compute the projectabilities of the Bloch states
    onto the Bloch sums obtained from the PAOs in the pseudopotential\
    `projwfc.x -in W.proj > proj.out`

5.  Run `generate_weights` to extract the projectabilitites from
    `proj.out` in a format suitable to be read by Xmgrace or gnuplot\
    `./generate_weights.sh`

6.  Plot the projectabilities and fit the data with the complementary
    error function
    $$f(\epsilon;\mu,\sigma) = \frac{1}{2}\mathrm{erfc}(-\frac{\mu - \epsilon}{\sigma}).$$
    We are going Xmgrace to plot the projectabilities and perform the
    fitting. Open Xmgrace\
    `xmgrace `

    To Import the `p_vs_e.dat` file, click on `Data` from the top bar
    and then `Import -> ASCII...`. At this point a new window
    `Grace: Read sets` should pop up. Select `p_vs_e.dat` in the `Files`
    section, click `Ok` at the bottom and close the window. You should
    now be able to see a quite noisy function that is bounded between 1
    and 0. You can modify the appearence of the plot by clicking on
    `Plot` in the top bar and then `Set appearance...`. In the `Main`
    section of the pop-up window change the symbol type from `None` to
    `Circle`. Change the line type from straight to none, since the
    lines added by default by Xmgrace are not meaningful. For the
    fitting, go to
    `Data -> Transformations -> Non-linear curve fitting`. In this
    window, select the source from the `Set` box and in the `Formula`
    box insert the following\
    `y = 0.5 * erfc( ( x - A0 ) / A1 )`\
    Select 2 as number of parameters, give 40 as initial condition for
    `A0` and 7 for `A1`. Click `Apply`. A new window should pop up with
    the stats of the fitting. In particular you should find a
    `Correlation coefficient` of 0.96 and a value of $39.9756$ for `A0`
    and $6.6529$ for `A1`. These are the value of $\mu_{fit}$ and
    $\sigma_{fit}$ we are going to use for the SCDM method. In
    particular, $\mu_{SCDM} = \mu_{fit} - 3\sigma_{fit} = 20.0169$ eV
    and $\sigma_{SCDM} = \sigma_{fit} = 6.6529$ eV. The motivation for
    this specific choice of $\mu_{fit}$ and $\sigma_{fit}$ may be found
    in Ref. [@Vitale2019automated], where the authors also show
    validation of this approach on a dataset of 200 materials. You
    should now see the fitting function, as well as the
    projectabilities, in the graph (see Fig.
    [10](#fig:W_fit){reference-type="ref" reference="fig:W_fit"}-(a)).\

7.  Open `W.pw2wan` and append the following lines

    scdm_entanglement = erfc

    scdm_mu = 20.0169

    scdm_proj = .true.

    scdm_sigma = 6.6529

    /

8.  Run `pw2wannier90` to compute the overlaps between Bloch states and
    the projections for the starting guess (written in the `W.mmn` and
    `W.amn` files)\
    `pw2wannier90.x -in W.pw2wan > pw2wan.out`

9.  Run `wannier90` to obtain the interpolated bandstructure (see Fig.
    [10](#fig:W_fit){reference-type="ref" reference="fig:W_fit"}-(b)).\
    `wannier90.x W`

    Please cite Ref. [@Vitale2019automated] in any publication employing
    the procedure outlined in this example to obtain $\mu$ and $\sigma$.

<figure id="fig:W_fit">

<figcaption> a) Each blue dot represents the projectability as defined
in Eq. (22) of Ref. <span class="citation"
data-cites="Vitale2019automated"></span> of the state <span
class="math inline">|<em>n</em><strong>k</strong>⟩</span> as a function
of the corresponding energy <span
class="math inline"><em>ϵ</em><sub><em>n</em><strong>k</strong></sub></span>
for tungsten. The yellow line shows the fitted complementary error
function. The vertical red line represents the value of <span
class="math inline"><em>σ</em><sub><em>f</em><em>i</em><em>t</em></sub></span>
while the vertical green line represents the optimal value of <span
class="math inline"><em>μ</em><sub><em>S</em><em>C</em><em>D</em><em>M</em></sub></span>,
i.e. <span
class="math inline"><em>μ</em><sub><em>S</em><em>C</em><em>D</em><em>M</em></sub> = <em>μ</em><sub><em>f</em><em>i</em><em>t</em></sub> − 3<em>σ</em><sub><em>f</em><em>i</em><em>t</em></sub></span>.
b) Band structure of tungsten on the <span
class="math inline"><em>Γ</em></span>-H-N-<span
class="math inline"><em>Γ</em></span> path from DFT calculations (solid
black) and Wannier interpolation using the SCDM method to construct the
initial guess (red dots).</figcaption>
</figure>

# 33: Monolayer BC$_2$N -- $k\cdot p$ expansion coefficients {#monolayer-bc_2n-kcdot-p-expansion-coefficients .unnumbered}

-   Outline: *Calculate $k\cdot p$ expansion coefficients monolayer
    BC$_2$N using quasi-degenerate (Löwdin) perturbation theory. In
    preparation for this example it may be useful to read Ref.
    [@ibanez-azpiroz-ArXiv2019]*

-   Directory: `examples/example33/`

-   Input files:

    -   `bc2n.scf` *The `PWSCF` input file for ground state calculation*

    -   `bc2n.nscf` *The `PWSCF` input file to obtain Bloch states on a
        uniform grid*

    -   `bc2n.pw2wan` *The input file for* `pw2wannier90`

    -   `bc2n.win` *The* `wannier90` *and* `postw90` *input file*

    1.  Run `PWSCF` to obtain the ground state of Gallium Arsenide

        `pw.x < bc2n.scf > scf.out`

    2.  Run `PWSCF` to obtain the ground state of Gallium Arsenide

        `pw.x < bc2n.nscf > nscf.out`

    3.  Run `Wannier90` to generate a list of the required overlaps
        (written into the `GaAs.nnkp` file)

        `wannier90.x -pp bc2n`

    4.  Run `pw2wannier90` to compute:

        -   The overlaps $\langle u_{n\bf{k}}|u_{n\bf{k+b}}\rangle$
            between spinor Bloch states (written in the `bc2n.mmn` file)

        -   The projections for the starting guess (written in the
            `bc2n.amn` file)

        `pw2wannier90.x < bc2n.pw2wan > pw2wan.out`

    5.  Run `wannier90` to compute MLWFs

        `wannier90.x bc2n`

    6.  Run `postw90` to compute expansion coefficients

        `postw90.x bc2n`

    ## Expansion coefficients {#expansion-coefficients .unnumbered}

    For computing $k\cdot p$ expansion coefficients as given by
    quasi-degenerate (Löwdin) perturbation theory, set

        berry = true
        berry_task = kdotp

    Select the k-point around which the expansion coefficients will be
    computed, *e.g.*, the S point

        kdotp_kpoint  =  0.5000 0.0000 0.5000

    Set number of bands that should be taken into account for the
    $k\cdot p$ expansion, as well as their band indexes within the
    Wannier basis

        kdotp_num_bands = 2
        kdotp_bands =  2,3

    Since no k-space integral is needed, set

        berry_kmesh = 1 1 1

    Although not used, we also need to input the value of the Fermi
    level in eV

        fermi_energy = [insert your value here]

    On output, the program generates three files, namely
    `SEED-kdotp_0.dat`, `SEED-kdotp_1.dat` and `SEED-kdotp_2.dat`, which
    correspond to the zeroth, first and second order expansion
    coefficients, respectively. The dimension of the matrix contained in
    each file is $3^{l}\times N^{2}$, where $N$ is the number of bands
    set by `kdotp_num_bands`, and $l$ is the order of the expansion term
    (currently $l=0,1$ or $2$).

    These coefficients can be used, among other things, to compute the
    energy dispersion of the bands of interest around the chosen
    k-point. The $k\cdot p$ band dispersion can be computed and plotted
    along $k_x$ (from S to X) using python and the file `kdotp_plot.py`
    provided in the example folder

        python kdotp_plot.py

    For comparison, the exact band structure calculated using Wannier90
    (file `bc2n_band.dat`, generated automatically) is also plotted
    along (see Fig. [11](#fig:bc2n-bnd){reference-type="ref"
    reference="fig:bc2n-bnd"}).

    <figure id="fig:bc2n-bnd">
    <div class="center">
    <embed src="kdotp_bands_SX.pdf" style="width:12cm" />
    </div>
    <figcaption>Band dispersion of monolayer BC<span
    class="math inline"><sub>2</sub></span>N around <span
    class="math inline"><em>S</em></span> point. Exact results (solid dots)
    are compared to first-order (blue) and second-order (red) <span
    class="math inline"><em>k</em> ⋅ <em>p</em></span> model results for
    valence and conduction bands.</figcaption>
    </figure>

[^1]: Once `XCrySDen` starts, click on `Tools` $\rightarrow$ `Data Grid`
    in order to specify an isosurface value to plot.

[^2]: Please note the following counterintuitive feature in `pwscf`: in
    order to obtain a ground state with magnetization along the
    *positive* $z$-axis, one should use a *negative* value for the
    variable `starting_magnetization`.

[^3]: The calculation of the AHC using `berry_task = kubo` involves a
    truncation of the sum over empty states in the Kubo-Greenwood
    formula: see description of the keyword `kubo_eigval_max` in the
    User Guide. As discussed around Eq. (11.17) of the User Guide, no
    truncation is done with `berry_task = ahc`.
