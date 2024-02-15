# 6: Copper &#151; Fermi surface {#copper-fermi-surface .unnumbered}

-   Outline: *Obtain MLWFs to describe the states around the Fermi-level
    in copper*

-   Directory: `tutorials/tutorial06/` *Files can be downloaded from [here](https://github.com/wannier-developers/wannier90/tutorials/tutorial06)*

-   Input Files

    -    `copper.scf` *The `pwscf` input file for ground
        state calculation*

    -    `copper.nscf` *The `pwscf` input file to obtain
        Bloch states on a uniform grid*

    -    `copper.pw2wan` *Input file for `pw2wannier90`*

    -    `copper.win` *The `wannier90` input file*

1.  Run `pwscf` to obtain the ground state of copper

    ```bash title="Terminal"
    pw.x < copper.scf > scf.out
    ```

2.  Run `pwscf` to obtain the Bloch states on a uniform
    k-point grid

    ```bash title="Terminal"
    pw.x < copper.nscf > nscf.out
    ```

3.  Run `wannier90` to generate a list of the required overlaps (written
    into the `copper.nnkp` file).

    ```bash title="Terminal"
    wannier90.x -pp copper
    ```

4.  Run `pw2wannier90` to compute the overlap between Bloch states and
    the projections for the starting guess (written in the `copper.mmn`
    and `copper.amn` files).

    ```bash title="Terminal"
    pw2wannier90.x < copper.pw2wan > pw2wan.out
    ```

5.  Run `wannier90` to compute the MLWFs.

    ```bash title="Terminal"
    wannier90.x copper
    ```

Inspect the output file `copper.wout`.

1.  Use Wannier interpolation to obtain the Fermi surface of copper.
    Rather than re-running the whole calculation we can use the unitary
    transformations obtained in the first calculation and restart from
    the plotting routine. Add the following keywords to the
    ` copper.win` file:

    ```vi title="Input file"
    restart = plot
    
    fermi_energy = [insert your value here]
    
    fermi_surface_plot = true
    ```

    and re-run `wannier90`. The value of the Fermi energy can be
    obtained from the initial first principles calculation.
    `wannier90` calculates the band energies, through Wannier
    interpolation, on a dense mesh of k-points in the Brillouin zone.
    The density of this grid is controlled by the keyword
    `fermi_surface_num_points`. The default value is 50 (i.e., 50$^3$
    points). The Fermi surface file `copper.bxsf` can be viewed using
    `XCrySDen`, e.g.,
    
    ```bash title="Terminal"
    xcrysden --bxsf copper.bxsf
    ```

2.  Plot the interpolated bandstructure. A suitable path in k-space is

    ```vi title="Input file"
    begin kpoint_path
    G 0.00 0.00 0.00 X 0.50 0.50 0.00
    X 0.50 0.50 0.00 W 0.50 0.75 0.25
    W 0.50 0.75 0.25 L 0.00 0.50 0.00
    L 0.00 0.50 0.00 G 0.00 0.00 0.00
    G 0.00 0.00 0.00 K 0.00 0.50 -0.50
    end kpoint_path
    ```

## Further ideas {#further-ideas .unnumbered}

-   Compare the Wannier interpolated bandstructure with the full
    `pwscf` bandstructure. Obtain MLWFs using a denser
    k-point grid. To plot the bandstructure you can use the
    `pwscf` tool `bands.x` or the small FORTRAN program
    available at <http://www.tcm.phy.cam.ac.uk/~jry20/bands.html>.

-   Investigate the effects of the outer and inner energy windows on the
    interpolated bands.

-   Instead of extracting a subspace of seven states, we could extract a
    nine dimensional space (i.e., with $s$, $p$ and $d$ character).
    Examine this case and compare the interpolated bandstructures.
