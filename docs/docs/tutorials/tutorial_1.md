# 1: Gallium Arsenide &#151; MLWFs for the valence bands {#gallium-arsenide-mlwfs-for-the-valence-bands .unnumbered}

-   Outline: *Obtain and plot MLWFs for the four valence bands of GaAs.*

-   Generation details: *From `pwscf`, using norm-conserving
    pseudopotentials and a <br> 
    2$\times$2$\times$2 k-point grid. Starting guess: four bond-centred Gaussians.*

-   Directory: `tutorials/tutorial01/` *Files can be
    downloaded from [here](https://github.com/wannier-developers/wannier90/tutorials/tutorial01)*

-   Input Files

    -    `gaas.win` *The master input file*

    -    `gaas.mmn` *The overlap matrices
        $\mathbf{M}^{(\mathbf{k},\mathbf{b})}$*

    -    `gaas.amn` *Projection $\mathbf{A}^{(\mathbf{k})}$ of the Bloch
        states onto a set of trial localised orbitals*

    -    `UNK00001.1` *The Bloch states in the real space unit cell. For
        plotting only.*

1.  Run `wannier90` to minimise the MLWFs spread

    ```bash title="Terminal"
    wannier90.x  gaas
    ```

    Inspect the output file `gaas.wout`. The total spread converges to
    its minimum value after just a few iterations. Note that the
    geometric centre of each MLWF lies along a Ga-As bond, slightly
    closer to As than Ga. Note also that the memory requirement for the
    minimisation of the spread is very low as the MLWFs are defined at
    each k-point by just the 4$\times$4 unitary
    matrices $\mathbf{U}^{(\mathbf{k})}$.

2.  Plot the MLWFs by adding the following keywords to the input file
    `gaas.win`

    ```vi title="Input file"
    wannier_plot = true
    ```

    and re-running `wannier90`. To visualise the MLWFs we must represent
    them explicitly on a real space grid (see the [User guide page](../../user_guide/wannier90/methodology#methodology)). As a
    consequence, plotting the MLWFs is slower and uses more memory than
    the minimisation of the spread. The four files that are created
    (`gaas_00001.xsf`, etc.) can be viewed using `XCrySDen`[^1],
    e.g.,

    ```bash title="Terminal"
    xcrysden --xsf gaas_00001.xsf
    ```

    For large systems, plotting the MLWFs may be time consuming and
    require a lot of memory. Use the keyword `wannier_plot_list` to plot
    a subset of the MLWFs. E.g., to plot the 1st and 3rd MLWFs use

    ```vi title="Input file"
    wannier_plot_list = 1 3
    ```

    The MLWFs are plotted in a supercell of the unit cell. The size of
    this supercell is set through the keyword ` wannier_plot_supercell`.
    The default value is 2 (corresponding to a supercell with eight
    times the unit cell volume). We recommend not using values great
    than 3 as the memory and computational cost scales cubically with
    supercell size.

    Plot the 3rd MLWFs in a supercell of size 3. Choose a low value for
    the isosurface (say 0.5). Can you explain what you see?

    *Hint:* For a finite k-point mesh, the MLWFs are in fact periodic
    and the period is related to the spacing of the k-point mesh. For
    mesh with $n$ divisions in the $i^{\mathrm{th}}$ direction in the
    Brillouin zone, the MLWFs "live" in a supercell $n$ times the unit
    cell.
    [^1]: http://www.xcrysden.org/
