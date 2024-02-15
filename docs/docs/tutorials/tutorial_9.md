# 9: Cubic BaTiO$_3$ {#cubic-batio_3 .unnumbered}

-   Outline: *Obtain MLWFs for a perovskite*

-   Directory: `tutorials/tutorial09/` *Files can be downloaded from [here](https://github.com/wannier-developers/wannier90/tutorials/tutorial09)*

-   Input Files

    -    `batio3.scf` *The `pwscf` input file for ground
        state calculation*

    -    `batio3.nscf` *The `pwscf` input file to obtain
        Bloch states on a uniform grid*

    -    `batio3.pw2wan` *Input file for `pw2wannier90`*

    -    `batio3.win` *The `wannier90` input file*

To start with, we are going to obtain MLWFs for the oxygen 2p states.
From the bandstructure [@marzari-aip98], these form an isolated group
of bands. We use the `wannier90` keyword `exclude_bands` to remove all
but the 2p bands from the calculation of the overlap and projection
matrices (we don't have to do this, but it saves time).

1.  Run `pwscf` to obtain the ground state of BaTiO$_3$

    ```bash title="Terminal"
    pw.x < BaTiO3.scf > scf.out
    ```

2.  Run `pwscf` to obtain the Bloch states on a uniform
    k-point grid

    ```bash title="Terminal"
    pw.x < BaTiO3.nscf > nscf.out
    ```

3.  Run `wannier90` to generate a list of the required overlaps (written
    into the `BaTiO3.nnkp` file).

    ```bash title="Terminal"
    wannier90.x -pp BaTiO3
    ```

4.  Run `pw2wannier90` to compute the overlap between Bloch states and
    the projections for the starting guess (written in the `BaTiO3.mmn`
    and `BaTiO3.amn` files).

    ```bash title="Terminal"
    pw2wannier90.x < BaTiO3.pw2wan > pw2wan.out
    ```

5.  Run `wannier90` to compute the MLWFs.

    ```bash title="Terminal"
    wannier90.x BaTiO3
    ```

Inspect the output file `BaTiO3.wout`.

Plot the second MLWF, as described in Section 1, by adding the following
keywords to the input file `BaTiO3.win`


```vi title="Input file"
wannier_plot = true
restart = plot
wannier_plot_list = 2
wannier_plot_supercell = 3
```

and re-running `wannier90`. Visualise it using `XCrySDen`,

```bash title="Terminal"
xcrysden `--`xsf BaTiO3_00002.xsf
```

We can now simulate the ferroelectric phase by displacing the Ti atom.
Change its position to

```vi title="Input file"
Ti 0.505 0.5 0.5
```

and regenerate the MLWFs (i.e., compute the ground-state charge density
and Bloch states using `pwscf`, etc.) and look at the change
in the second MLWF.

## Further ideas {#further-ideas-1 .unnumbered}

-   Look at MLWFs for other groups of bands. What happens if you form
    MLWFs for the whole valence manifold?

-   Following Ref. [@marzari-aip98], compute the Born effective
    charges from the change in Wannier centres under an atomic
    displacement.


