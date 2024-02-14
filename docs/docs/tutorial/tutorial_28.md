# 28: Diamond &#151; plotting of MLWFs using Gaussian cube format and VESTA {#diamond-plotting-of-mlwfs-using-gaussian-cube-format-and-vesta .unnumbered}

-   Outline: *Obtain MLWFs for the valence bands of diamond and output
    them in Gaussian cube format*

-   Directory: `tutorial/tutorial28/` *The input files for this examples
    are the same as the ones in example05 and can be downloaded [here]{}*

-   Input Files

    \-    `diamond.scf` *The `pwscf` input file for ground
        state calculation*

    \-    `diamond.nscf` *The `pwscf` input file to obtain
        Bloch states on a uniform grid*

    \-    `diamond.pw2wan` *The input file for `pw2wannier90`*

    \-    `diamond.win` *The `wannier90` input file*

1.  Run `pwscf` to obtain the ground state of diamond

    ```bash title="Terminal"
    pw.x < diamond.scf > scf.out
    ```

2.  Run `pwscf` to obtain the Bloch states on a uniform
    k-point grid

    ```bash title="Terminal"
    pw.x < diamond.nscf > nscf.out
    ```

3.  Run `wannier90` to generate a list of the required overlaps (written
    into the `diamond.nnkp` file).

    ```bash title="Terminal"
    wannier90.x -pp diamond
    ```

4.  Run `pw2wannier90` to compute the overlap between Bloch states and
    the projections for the starting guess (written in the `diamond.mmn`
    and `diamond.amn` files).

    ```bash title="Terminal"
    pw2wannier90.x < diamond.pw2wan > pw2wan.out
    ```

5.  When the lattice vectors are non-orthogonal, not all the
    visualisation programs are capable to plot volumetric data in the
    Gaussian cube format. One program that can read volumetric data for
    these systems is VESTA. To instruct `wannier90` to output the MLWFs
    data in Gaussian cube format you need to add the following lines to
    the `.win` file

    ```vi title="Input file"
    wannier_plot           = .true.
    wannier_plot_supercell = 3
    wannier_plot_format    = cube
    wannier_plot_mode      = crystal
    wannier_plot_radius    = 2.5
    wannier_plot_scale     = 1.0
    ```

    Run `wannier90` to compute the MLWFs and output them in the Gaussian
    cube file.

    ```bash title="Terminal"
    wannier90.x diamond
    ```

6.  Plot the first MLWF with VESTA `vesta diamond_00001.cube`

Extra: Instead of using `wannier_plot_mode = crystal` try to use the
molecule mode as `wannier_plot_mode = molecule` (see the user guide for
the definition of this keyword). Add the following line to the `.win`
file:

```vi title="Input file"
restart = plot
```

and re-run `wannier90`. Use VESTA to plot the resulting MLWFs, do you
see any difference from the `crystal` mode case? Can you explain why?
Try to change the size of the supercell from 3 to 5, do you expect the
results to be different? (*Hint:* When using the Gaussian cube format
the code outputs the WF on a grid that is smaller than the super
unit-cell. The size of the grid is specified by `wannier_plot_scale` and
`wannier_plot_radius`.)


