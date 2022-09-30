# 10: Graphite

- Outline: *Obtain MLWFs for graphite (AB, Bernal)*

- Directory: `tutorials/tutorial10/` *Files can be downloaded from
    [here](https://github.com/wannier-developers/wannier90/tree/develop/tutorials/tutorial10)*

- Input Files

    - `graphite.scf` *The `pwscf` input file for ground
        state calculation*

    - `graphite.nscf` *The `pwscf` input file to obtain
        Bloch states on a uniform grid*

    - `graphite.pw2wan` *Input file for `pw2wannier90`*

    - `graphite.win` *The `wannier90` input file*

1. Run `pwscf` to obtain the ground state of graphite

    ```bash title="Terminal"
    pw.x < graphite.scf > scf.out
    ```

2. Run `pwscf` to obtain the Bloch states on a uniform
    k-point grid

    ```bash title="Terminal"
    pw.x < graphite.nscf > nscf.out
    ```

3. Run `wannier90` to generate a list of the required overlaps (written
    into the `graphite.nnkp` file).

    ```bash title="Terminal"
    wannier90.x -pp graphite
    ```

4. Run `pw2wannier90` to compute the overlap between Bloch states and
    the projections for the starting guess (written in the
    `graphite.mmn` and `graphite.amn` files).

    ```bash title="Terminal"
    pw2wannier90.x < graphite.pw2wan > pw2wan.out
    ```

5. Run `wannier90` to compute the MLWFs.

    ```bash title="Terminal"
    wannier90.x graphite
    ```
