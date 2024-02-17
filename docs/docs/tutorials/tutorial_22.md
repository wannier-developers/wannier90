# 22: Copper &#151; Symmetry-adapted Wannier functions

Note: This tutorial requires a recent version of the `pw2wannier90`
interface.

- Outline: *Obtain symmetry-adapted Wannier functions for Cu. By
    symmetry-adapted mode, for example, we can make atomic centered
    $s$-like Wannier function, which is not possible in the usual
    procedure to create maximally localized Wannier functions. For the
    theoretical background of the symmetry-adapted Wannier functions,
    see R. Sakuma, Phys. Rev. B **87**, 235109 (2013).*

- Directory: `tutorials/tutorial22/s_at_0.00/` *Files can be downloaded from [here](https://github.com/wannier-developers/wannier90/tree/develop/tutorials/tutorial22)*
    \

- Input Files

    - `Cu.scf` *The `pwscf` input file for ground state
        calculation*

    - `Cu.nscf` *The `pwscf` input file to obtain Bloch
        states on a uniform grid*

    - `Cu.pw2wan` *The input file for `pw2wannier90`*

    - `Cu.sym` *Used only in `tutorials/tutorial22/s_at_0.25/`.
        `pw2wannier90` reads this file when `“read_sym = .true.”` in
        `Cu.pw2wan`. By default,
        `“read_sym = .false.” and ``Cu.sym`` is the output of ``pw2wannier90`,
        in which the symmetry operations employed in the calculation are
        written for reference.*

    - `Cu.win` *The `wannier90` input file*

1. Run `pwscf` to obtain the ground state of Cu

    ```bash title="Terminal"
    pw.x < Cu.scf > scf.out
    ```

2. Run `pwscf` to obtain the Bloch states on a uniform
    k-point grid

    ```bash title="Terminal"
    pw.x < Cu.nscf > nscf.out
    ```

3. Run `wannier90` to generate a list of the required overlaps (written
    into the `Cu.nnkp` file).

    ```bash title="Terminal"
    wannier90.x -pp Cu
    ```

4. Run `pw2wannier90` to compute the overlap between Bloch states, the
    projections for the starting guess, and the symmetry information
    needed for symmetry-adapted mode (written in the `Cu.mmn`, `Cu.amn`,
    and `Cu.dmn` files, respectively).

    ```bash title="Terminal"
    pw2wannier90.x < Cu.pw2wan > pw2wan.out
    ```

5. Run `wannier90` to compute the MLWFs.

    ```bash title="Terminal"
    wannier90.x Cu
    ```

Each directory creates $s$-like symmetry-adapted Wannier function
centered at different position on top of atomic centered $d$-like
Wannier functions. See more detail in `tutorials/tutorial22/README`.
