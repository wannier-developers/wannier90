# 37: Bcc Iron &#151; Translationally-invariant Wannier Interpolation

- Outline: *Perform Wannier interpolation of orbital magnetization
    using translationally-invariant formulas.*

- Directory: [`tutorial/tutorial37/`](https://github.com/wannier-developers/wannier90/tree/develop/tutorials/tutorial37)

- Input Files

    - `Fe.scf` *The `pwscf` input file for ground
        state calculation*

    - `Fe.nscf` *The `pwscf` input file to obtain
        Bloch states on a uniform grid*

    - `Fe.pw2wan` *Input file for `pw2wannier90`*

    - `Fe.win` *The `wannier90` input file*

1. For both directories, run `pwscf` to obtain the ground state of silicon

    ```bash title="Terminal"
    pw.x < Fe.scf > Fe.out
    ```

2. Run `pwscf` to obtain the Bloch states on a uniform
    k-point grid.

    ```bash title="Terminal"
    pw.x < Fe.nscf > nscf.out
    ```

3. Run `wannier90` to generate a list of the required overlaps (written
    into the `Fe.nnkp` file).

    ```bash title="Terminal"
    wannier90.x -pp Fe
    ```

4. Run `pw2wannier90` to compute the overlap between Bloch states and
    the projections for the starting guess (written in the `Fe.mmn`,
    `Fe.amn`, and `Fe.uHu` files).

    ```bash title="Terminal"
    pw2wannier90.x < Fe.pw2wan > Fe.out
    ```

5. Run `wannier90` to compute the MLWFs.

    ```bash title="Terminal"
    wannier90.x Fe
    ```

6. Run `postw90` to perform Wannier interpolation of orbital magnetization
   with and without translationally-invariant formulas. You can use the option
   `transl_inv_full`.

    ```bash title="Terminal"
    postw90.x Fe
    ```

7. Do step 6 varying the fermi energy and obtain a graph for the relation
   between orbital magnetization and the fermi energy. Compare the results of
   `without_translation` and `with_translation` directories for both
   `transl_inv_full = F` and `transl_inv_full = T`.

## Further ideas

- Increase the grid and find out that the `transl_inv_full = T` case
  converges faster.
