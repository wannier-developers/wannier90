# 3: Silicon &#151; Disentangled MLWFs {#silicon-disentangled-mlwfs .unnumbered}

-   Outline: *Obtain disentangled MLWFs for the valence and low-lying
    conduction states of Si. Plot the interpolated bandstructure*

-   Generation Details: *From `pwscf`, using norm-conserving
    pseudopotentials and a <br>
    4$\times$4$\times$4 k-point grid.
    Starting guess: atom-centred sp$^3$ hybrid orbitals*

-   Directory: `tutorials/tutorial03/` *Files can be downloaded from [here](https://github.com/wannier-developers/wannier90/tutorials/tutorial03)*

-   Input Files

    -    `silicon.win` *The master input file*

    -    `silicon.mmn` *The overlap matrices
        $\mathbf{M}^{(\mathbf{k},\mathbf{b})}$*

    -    `silicon.amn` *Projection $\mathbf{A}^{(\mathbf{k})}$ of the
        Bloch states onto a set of trial localised orbitals*

    -    `silicon.eig` *The Bloch eigenvalues at each k-point*

The valence and lower conduction states can be represented by MLWFs with
$sp^3$-like symmetry. The lower conduction states are not separated from
the higher states by an energy gap. In order to form localised WF, we
use the disentanglement procedure introduced in Ref. [@souza-prb01]. The
position of the inner and outer energy windows are shown in the bandstructure
[plot](#fig:si.bnd).

1.  Run `wannier90`.

    ```bash title="Terminal"
    wannier90.x silicon
    ```

    Inspect the output file `silicon.wout`. The minimisation of the
    spread occurs in a two-step procedure [@souza-prb01]. First, we
    minimise $\Omega_{\rm I}$ -- this is the extraction of the optimal
    subspace in the disentanglement procedure. Then, we minimise
    $\Omega_{\rm D} +
    \Omega_{{\rm OD}}$.

2.  Plot the energy bands by adding the following commands to the input
    file `silicon.win`

    ```vi title="Input file"
    restart = plot
    
    bands_plot = true
    ```

    and re-running `wannier90`. The files `silicon_band.dat` and
    ` silicon_band.gnu` are created. To plot the bandstructure run
    gnuplot
    
    ```bash title="Terminal"
    gnuplot
    ```

    and within the gnuplot shell type
   
    ```gnuplot title="Gnuplot shell"
    load 'silicon_band.gnu'
    ```

    The k-point path for the bandstructure interpolation is set in the
    `kpoint_path` block. Try plotting along different paths.

<figure markdown="span" id="fig:si.bnd">
![Image title](./si.webp){ width="500" }
<figcaption> Bandstructure of silicon showing the position of the outer
and inner energy windows.</figcaption>
</figure>
