# 21: Gallium Arsenide &#151; Symmetry-adapted Wannier functions {#gallium-arsenide-symmetry-adapted-wannier-functions .unnumbered}

Note: This example requires a recent version of the `pw2wannier90`
interface.

-   Outline: *Obtain symmetry-adapted Wannier functions out of four
    valence bands of GaAs. For the theoretical background of the
    symmetry-adapted Wannier functions, see R. Sakuma, Phys. Rev. B
    **87**, 235109 (2013).*

-   Directory: `tutorial/tutorial21/atom_centered_As_sp/` *Files can be downloaded from [here]{}*

-   Input Files

    \-    `GaAs.scf` *The `pwscf` input file for ground state
        calculation*

    \-    `GaAs.nscf` *The `pwscf` input file to obtain Bloch
        states on a uniform grid*

    \-    `GaAs.pw2wan` *The input file for `pw2wannier90`*

    \-    `GaAs.win` *The `wannier90` input file*

1.  Run `pwscf` to obtain the ground state of GaAs

    ```bash title="Terminal"
    pw.x < GaAs.scf > scf.out
    ```

2.  Run `pwscf` to obtain the Bloch states on a uniform
    k-point grid

    ```bash title="Terminal"
    pw.x < GaAs.nscf > nscf.out
    ```

3.  Run `wannier90` to generate a list of the required overlaps (written
    into the `GaAs.nnkp` file).

    ```bash title="Terminal"
    wannier90.x -pp GaAs
    ```

4.  Run `pw2wannier90` to compute the overlap between Bloch states, the
    projections for the starting guess, and the symmetry information
    needed for symmetry-adapted mode (written in the `GaAs.mmn`,
    `GaAs.amn`, and `GaAs.dmn` files, respectively).

    ```bash title="Terminal"
    pw2wannier90.x < GaAs.pw2wan > pw2wan.out
    ```

5.  Run `wannier90` to compute the MLWFs.

    ```bash title="Terminal"
    wannier90.x GaAs
    ```

Each directory creates different kind of symmetry-adapted Wannier
function. See more detail in `examples/example21/README`.


