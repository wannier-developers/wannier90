# 11: Silicon &#151; Valence and low-lying conduction states {#silicon-valence-and-low-lying-conduction-states .unnumbered}

## Valence States {#valence-states .unnumbered}

-   Outline: *Obtain MLWFs for the valence bands of silicon.*

-   Directory: `tutorials/tutorial11/` *Files can be downloaded from [here](https://github.com/wannier-developers/wannier90/tutorials/tutorial11)*

-   Input Files

    -    `silicon.scf` *The `pwscf` input file for ground
        state calculation*

    -    `silicon.nscf` *The `pwscf` input file to obtain
        Bloch states on a uniform grid*

    -    `silicon.pw2wan` *Input file for `pw2wannier90`*

    -    `silicon.win` *The `wannier90` input file*

1.  Run `pwscf` to obtain the ground state of silicon

    ```bash title="Terminal"
    pw.x < silicon.scf > scf.out
    ```

2.  Run `pwscf` to obtain the Bloch states on a uniform
    k-point grid. Note that we request the lower 4 (valence) bands

    ```bash title="Terminal"
    pw.x < silicon.nscf > nscf.out
    ```

3.  Run `wannier90` to generate a list of the required overlaps (written
    into the `silicon.nnkp` file).

    ```bash title="Terminal"
    wannier90.x -pp silicon
    ```

4.  Run `pw2wannier90` to compute the overlap between Bloch states and
    the projections for the starting guess (written in the `silicon.mmn`
    and `silicon.amn` files).

    ```bash title="Terminal"
    pw2wannier90.x < silicon.pw2wan > pw2wan.out
    ```

5.  Run `wannier90` to compute the MLWFs.

    ```bash title="Terminal"
    wannier90.x silicon
    ```

Inspect the output file `silicon.wout`. The total spread converges to
its minimum value after just a few iterations. Note that the geometric
centre of each MLWF lies at the centre of the Si-Si bond. Note also that
the memory requirement for the minimisation of the spread is very low as
the MLWFs are defined by just the 4$\times$4 unitary
matrices $\mathbf{U}^{(\mathbf{k})}$.

Plot the MLWFs by adding the following keywords to the input file
` silicon.win`

```vi title="Input file"
wannier_plot = true
```

and re-running `wannier90`. Visualise them using `XCrySDen`, e.g.,

```bash title="Terminal"
xcrysden --xsf silicon_00001.xsf
```

## Valence + Conduction States {#valence-conduction-states .unnumbered}

-   Outline: *Obtain MLWFs for the valence and low-lying conduction-band
    states of Si. Plot the interpolated bandstructure. Apply a scissors
    correction to the conduction bands.*

-   Input Files

    -    `silicon.scf` *The `pwscf` input file for ground
        state calculation*

    -    `silicon.nscf` *The `pwscf` input file to obtain
        Bloch states on a uniform grid*

    -    `silicon.pw2wan` *Input file for `pw2wannier90`*

    -    `silicon.win` *The `wannier90` input file*

The valence and lower conduction states can be represented by MLWFs with
$sp^3$-like symmetry. The lower conduction states are not separated by
an energy gap from the higher states. In order to form localised WF we
use the disentanglement procedure introduced in Ref. [@souza-prb01]. The
position of the inner and outer energy windows are shown in
[this plot](../tutorial_3#fig:si.bnd){reference-type="ref" reference="fig:si.bnd"}.

1.  Modify the input file and run `pwscf` and `wannier90`.\
    Inspect the output file `silicon.wout`. The minimisation of the
    spread occurs in a two-step procedure. First, we minimise
    $\Omega_{\rm
      I}$ -- this is the extraction of the optimal subspace in the
    disentanglement procedure. Then, we minimise $\Omega_{\rm
      O}+\Omega_{{\rm OD}}$.

2.  Plot the bandstructure by adding the following commands to the input
    file `silicon.win`

    ```vi title="Input file"
    restart = plot
    
    bands_plot = true
    ```

    and re-running `wannier90`. The files `silicon_band.dat` and
    ` silicon_band.gnu` are created. To plot the bandstructure using
    gnuplot

    ```bash title="Terminal"
    gnuplot
    ```

    ```gnuplot title="Gnuplot shell"
    load 'silicon_band.gnu'
    ```

    The k-point path for the bandstructure interpolation is set in the
    ` kpoint_path` block. Try plotting along different paths.

## Further ideas {#further-ideas-2 .unnumbered}

-   Compare the Wannier-interpolated bandstructure with the full
    `pwscf` bandstructure. Recompute the MLWFs using a finer
    $k$-point grid (e.g.,
    6$\times$6$\times$6 or
    8$\times$8$\times$8) and note how
    the accuracy of the interpolation increases [@yates-prb07].

-   Compute four MLWFs spanning the low-lying conduction states (see
    Ref. [@souza-prb01]).
