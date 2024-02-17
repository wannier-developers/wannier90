# 19: Iron &#151; Orbital magnetization

Note: This tutorial requires a recent version of the `pw2wannier90`
interface.

- Outline: *Calculate the orbital magnetization of ferromagnetic bcc
    Fe by Wannier interpolation.*

- Directory: `tutorials/tutorial19/` *Files can be downloaded from [here](https://github.com/wannier-developers/wannier90/tree/develop/tutorials/tutorial19)*

- Input files

    - `Fe.scf` *The `pwscf` input file for ground state
        calculation*

    - `Fe.nscf` *The `pwscf` input file to obtain Bloch
        states on a uniform grid*

    - `Fe.pw2wan` *The input file for `pw2wannier90`*

    - `Fe.win` *The `wannier90` and `postw90` input file*

The sequence of steps below is the same of Tutorials [17](tutorial_17.md) and [18](tutorial_18.md). If you
have already run one of those tutorials, you can reuse the output files
from steps 1 &#151; 3 and 5. Steps 4 and 6 should be carried out again using
the new input files `Fe.pw2wan` and `Fe.win`.

1. Run `pwscf` to obtain the ground state of iron

    ```bash title="Terminal"
    pw.x < Fe.scf > scf.out
    ```

2. Run `pwscf` to obtain the Bloch states on a uniform
    k-point grid

    ```bash title="Terminal"
    pw.x < Fe.nscf > nscf.out
    ```

3. Run `wannier90` to generate a list of the required overlaps (written
    into the `Fe.nnkp` file).

    ```bash title="Terminal"
    wannier90.x -pp Fe
    ```

4. Run `pw2wannier90` to compute:

    - The overlaps $\langle u_{n{\bf k}}\vert u_{m{\bf k}+{\bf
                  b}}\rangle$ (written in the `Fe.mmn` file)

    - The projections for the starting guess (written in the `Fe.amn`
        file)

    - The matrix elements $\langle u_{n{\bf k}+{\bf b}_1}\vert
              H_{\bf k}\vert u_{m{\bf k}+{\bf b}_2}\rangle$ (written in
        the `Fe.uHu` file)

    ```bash title="Terminal"
    pw2wannier90.x < Fe.pw2wan > pw2wan.out
    ```

5. Run `wannier90` to compute the MLWFs.

    ```bash title="Terminal"
    wannier90.x Fe
    ```

6. Run `postw90` to compute the orbital magnetization.

    ```bash title="Terminal"
    postw90.x Fe # (1)! 
    mpirun -np 8 postw90.x Fe  # (2)!
    ```

    1. serial execution
    2. example of parallel execution with 8 MPI processes

The orbital magnetization is computed as the BZ integral of the quantity
${\bf M}^{\rm orb}({\bf k})$ defined in this [equation](../user_guide/postw90/berry.md#mjx-eqn:eq:morb)
of the User Guide.
The relevant lines in `Fe.win` are

```vi title="Input file"
berry = true

berry_task = morb

berry_kmesh = 25 25 25

fermi_energy = [insert your value here]
```

After running `postw90`, compare the value of the orbital magnetization
reported in `Fe.wpout` with the spin magnetization in `scf.out`. Set
`iprint = 2` to report the decomposition of ${\bf M}^{\rm orb}$ into the
terms $J0$, $J1$, and $J2$ defined in Ref. [@lopez-prb12].

To plot $M_z^{\rm orb}({\bf k})$ along high-symmetry lines set
`berry = false` and uncomment in `Fe.win` the block of instructions
containing

```vi title="Input file"
kpath = true

kpath_task = bands+morb
```

After running `postw90`, issue

```bash title="Terminal"
python Fe-bands+morb_z.py
```

Compare with Fig. 2 of Ref. [@lopez-prb12], bearing in mind the factor
of $-1/2$ difference in the definition of ${\bf M}^{\rm
  orb}({\bf k})$ (see Ch. 11 in the User Guide).

To plot $M_z^{\rm orb}({\bf k})$ together with the Fermi contours on the
(010) BZ plane set `kpath = false`, uncomment in `Fe.win` the block of
instructions containing

```vi title="Input file"
kslice = true

kslice_task = morb+fermi_lines
```

re-run `postw90`, and issue

```bash title="Terminal"
python Fe-kslice-morb_z+fermi_lines.py
```

$M_z^{\rm orb}({\bf k})$ is much more evenly distributed in $k$-space
than the Berry curvature (see Tutorial [18](tutorial_18.md)). As a result, the integrated
orbital magnetization converges more rapidly with the BZ sampling.
