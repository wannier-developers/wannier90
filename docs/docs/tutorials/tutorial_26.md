# 26: Gallium Arsenide &#151; Selective localization and constrained centres

- Outline: *Application of the selectively localised Wannier function
    (SLWF) method to gallium arsenide (GaAs), following the example in
    Ref. [@Marianetti], which is essential reading for this tutorial
    tutorial.*

- Directory: `tutorials/tutorial26/` *Files can be downloaded from [here](https://github.com/wannier-developers/wannier90/tree/develop/tutorials/tutorial26)*

- Input files:

    - `GaAs.scf` *The `pwscf` input file for ground state calculation*

    - `GaAs.nscf` *The `pwscf` input file to obtain Bloch states on a
        uniform grid*

    - `GaAs.pw2wan` *The input file for* `pw2wannier90`

    - `GaAs.win` *The* `wannier90` *and* `postw90` *input file*

&nbsp;

1. Run `pwscf` to obtain the ground state of Gallium Arsenide

    ```bash title="Terminal"
    pw.x < GaAs.scf > scf.out
    ```

2. Run `pwscf` to obtain the ground state of Gallium Arsenide

    ```bash title="Terminal"
    pw.x < GaAs.nscf > nscf.out
    ```

3. Run `Wannier90` to generate a list of the required overlaps
    (written into the `GaAs.nnkp` file)

    ```bash title="Terminal"
    wannier90.x -pp GaAs
    ```

4. Run `pw2wannier90` to compute:

    - The overlaps $\langle u_{n\bf{k}}|u_{n\bf{k+b}}\rangle$
        between Bloch states (written in the `GaAs.mmn` file)

    - The projections for the starting guess (written in the
        `GaAs.amn` file)

    ```bash title="Terminal"
    pw2wannier90.x < GaAs.pw2wan > pw2wan.out
    ```

5. Inspect the `.win` file.

    - Make sure you understand the new keywords corresponding to
        the selective localisation algorithm.

    - Run `wannier90` to compute the SLWFs, in this case using one
        objective Wannier function.

    ```bash title="Terminal"
    wannier90.x GaAs
    ```

    To constrain the centre of the SLWF you need to add
    `slwf_constrain = true` and\
    `slwf_lambda = 1` to the input file and uncomment the `slwf_centres`
    block. This will add a penalty functional to the total spread, which
    will try to constrain the centre of the SLWF to be on the As atom
    (as explained in Ref. [@Marianetti], particularly from Eq. 24 to
    Eq. 35).

    Look at the value of the penalty functional, is this what you would
    expect at convergence? Does the chosen value of the Lagrange
    multiplier `slwf_lambda` give a SLWF function centred on the As
    atom?

    Alternatively, you can modify the `slwf_centres` block to constrain
    the centre of the SLWF to be on the Ga atom. Do you need a different
    value of `slwf_lambda` in this case to converge? Take a look at the
    result in Vesta and explain what you see. Do these functions
    transform like the identity under the action of the $T_d$ group?
