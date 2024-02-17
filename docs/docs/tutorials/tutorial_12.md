# 12: Benzene &#151; Valence and low-lying conduction states

## Valence States

- Outline: *Obtain MLWFs for the valence states of benzene*

- Directory: `tutorials/tutorial12/` *Files can be downloaded from [here](https://github.com/wannier-developers/wannier90/tree/develop/tutorials/tutorial12)*

- Input Files

    - `benzene.scf` *The `pwscf` input file for ground
        state calculation*

    - `benzene.pw2wan` *Input file for `pw2wannier90`*

    - `benzene.win` *The `wannier90` input file*

1. Run `pwscf` to obtain the ground state of benzene

    ```bash title="Terminal"
    pw.x < benzene.scf > scf.out
    ```

2. Run `wannier90` to generate a list of the required overlaps (written
    into the `benzene.nnkp` file).

    ```bash title="Terminal"
    wannier90.x -pp benzene
    ```

3. Run `pw2wannier90` to compute the overlap between Bloch states and
    the projections for the starting guess (written in the `benzene.mmn`
    and `benzene.amn` files).

    ```bash title="Terminal"
    pw2wannier90.x < benzene.pw2wan > pw2wan.out
    ```

4. Run `wannier90` to compute the MLWFs.

    ```bash title="Terminal"
    wannier90.x benzene
    ```

Inspect the output file `benzene.wout`. The total spread converges to
its minimum value after just a few iterations.

Plot the MLWFs by adding the following keywords to the input file
`benzene.win`

```vi title="Input file"
restart = plot
wannier_plot = true
wannier_plot_format = cube
wannier_plot_list = 2-4
```

and re-running `wannier90`. Visualise them using, e.g., `XCrySDen`.

## Valence + Conduction States

- Outline: *Obtain MLWFs for the valence and low-lying conduction
    states of benzene.*

- Input Files

    - `benzene.scf` *The `pwscf` input file for ground
        state calculation*

    - `benzene.nscf` *The `pwscf` input file to obtain
        Bloch states for the conduction states*

    - `benzene.pw2wan` *Input file for `pw2wannier90`*

    - `benzene.win` *The `wannier90` input file*

In order to form localised WF we use the disentanglement procedure. The
position of the inner energy window is set to lie in the energy gap; the
outer energy window is set to 4.0 eV. Modify the input file
appropriately.

1. Run `pwscf` and `wannier90`.\
    Inspect the output file `benzene.wout`. The minimisation of the
    spread occurs in a two-step procedure. First, we minimise
    $\Omega_{\rm
      I}$. Then, we minimise $\Omega_{\rm O}+\Omega_{{\rm OD}}$.

2. Plot the MLWFs by adding the following commands to the input file
    `benzene.win`

    ```vi title="Input file"
    restart = plot
    wannier_plot = true
    wannier_plot_format = cube
    wannier_plot_list = 1,7,13
    ```

    and re-running `wannier90`. Visualise them using, e.g., `XCrySDen`.
