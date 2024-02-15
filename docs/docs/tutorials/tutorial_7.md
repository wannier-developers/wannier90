# 7: Silane (SiH$_4$) &#151; Molecular MLWFs using $\Gamma$-point sampling {#silane-sih_4-molecular-mlwfs-using-gamma-point-sampling .unnumbered}

-   Outline: *Obtain MLWFs for the occupied states of molecular silane.
    $\Gamma$-point sampling*

-   Directory: `tutorials/tutorial07/` *Files can be downloaded from [here](https://github.com/wannier-developers/wannier90/tutorials/tutorial07)*

-   Input Files

    -    `silane.scf` *The `pwscf` input file for ground
        state calculation*

    -    `silane.nscf` *The `pwscf` input file to obtain
        Bloch states on a uniform grid*

    -    `silane.pw2wan` *Input file for `pw2wannier90`*

    -    `silane.win` *The `wannier90` input file*

1.  Run `pwscf` to obtain the ground state of silane

    ```bash title="Terminal"
    pw.x < silane.scf > scf.out
    ```

2.  Run `pwscf` to obtain the Bloch states on a uniform
    k-point grid

    ```bash title="Terminal"
    pw.x < silane.nscf > nscf.out
    ```

3.  Run `wannier90` to generate a list of the required overlaps (written
    into the `silane.nnkp` file).

    ```bash title="Terminal"
    wannier90.x -pp silane
    ```

4.  Run `pw2wannier90` to compute the overlap between Bloch states and
    the projections for the starting guess (written in the `silane.mmn`
    and `silane.amn` files).

    ```bash title="Terminal"
    pw2wannier90.x < silane.pw2wan > pw2wan.out
    ```

5.  Run `wannier90` to compute the MLWFs.

    ```bash title="Terminal"
    wannier90.x silane
    ```

