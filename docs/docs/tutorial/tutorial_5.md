# 5: Diamond &#151; MLWFs for the valence bands {#diamond-mlwfs-for-the-valence-bands .unnumbered}

-   Outline: *Obtain MLWFs for the valence bands of diamond*

-   Directory: `examples/example5/`

-   Input Files

    \-    `diamond.scf` *The `pwscf` input file for ground
        state calculation*

    \-    `diamond.nscf` *The `pwscf` input file to obtain
        Bloch states on a uniform grid*

    \-    `diamond.pw2wan` *The input file for `pw2wannier90`*

    \-    `diamond.win` *The `wannier90` input file*

1.  Run `pwscf` to obtain the ground state of diamond

    ```bash
    pw.x < diamond.scf > scf.out
    ```

2.  Run `pwscf` to obtain the Bloch states on a uniform
    k-point grid\

    ```bash
    pw.x < diamond.nscf > nscf.out
    ```

3.  Run `wannier90` to generate a list of the required overlaps (written
    into the `diamond.nnkp` file).
    
    ```bash
    wannier90.x -pp diamond
    ```

4.  Run `pw2wannier90` to compute the overlap between Bloch states and
    the projections for the starting guess (written in the `diamond.mmn`
    and `diamond.amn` files).

    ```bash
    pw2wannier90.x < diamond.pw2wan > pw2wan.out
    ```

5.  Run `wannier90` to compute the MLWFs.

    ```bash
    wannier90.x diamond
    ```


