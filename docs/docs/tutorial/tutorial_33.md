# 33: Monolayer BC$_2$N &#151; $k\cdot p$ expansion coefficients {#monolayer-bc_2n-kcdot-p-expansion-coefficients .unnumbered}

-   Outline: *Calculate $k\cdot p$ expansion coefficients monolayer
    BC$_2$N using quasi-degenerate (Löwdin) perturbation theory. In
    preparation for this example it may be useful to read Ref.
    [@ibanez-azpiroz-ArXiv2019]*

-   Directory: `tutorial/tutorial33/` *Files can be downloaded from [here](https://github.com/wannier-developers/wannier90/tutorials/tutorial33)*

-   Input files:

    \-   `bc2n.scf` *The `pwscf` input file for ground state calculation*

    \-   `bc2n.nscf` *The `pwscf` input file to obtain Bloch states on a
        uniform grid*

    \-   `bc2n.pw2wan` *The input file for* `pw2wannier90`

    \-   `bc2n.win` *The* `wannier90` *and* `postw90` *input file*

&nbsp;

1.  Run `pwscf` to obtain the ground state

    ```bash title="Terminal"
    pw.x < bc2n.scf > scf.out
    ```

2.  Run `pwscf` to obtain the ground state
    ```bash title="Terminal"
    pw.x < bc2n.nscf > nscf.out
    ```
3.  Run `Wannier90` to generate a list of the required overlaps
    (written into the `bc2n.nnkp` file)
    ```bash title="Terminal"
    wannier90.x -pp bc2n
    ```
4.  Run `pw2wannier90` to compute:
    -   The overlaps $\langle u_{n\bf{k}}|u_{n\bf{k+b}}\rangle$
        between spinor Bloch states (written in the `bc2n.mmn`file)
    -   The projections for the starting guess (written in the
        `bc2n.amn` file)
        
    ```bash title="Terminal"
    pw2wannier90.x < bc2n.pw2wan > pw2wan.out
    ```

5.  Run `wannier90` to compute MLWFs
    ```bash title="Terminal"
    wannier90.x bc2n
    ```

6.  Run `postw90` to compute expansion coefficients
    ```bash title="Terminal"
    postw90.x bc2n
    ```

## Expansion coefficients {#expansion-coefficients .unnumbered}

For computing $k\cdot p$ expansion coefficients as given by
quasi-degenerate (Löwdin) perturbation theory, set

```vi title="Input file"
berry = true
berry_task = kdotp
```

Select the k-point around which the expansion coefficients will be
computed, *e.g.*, the S point

```vi title="Input file"
kdotp_kpoint  =  0.5000 0.0000 0.5000
```

Set number of bands that should be taken into account for the
$k\cdot p$ expansion, as well as their band indexes within the
Wannier basis

```vi title="Input file"
kdotp_num_bands = 2
kdotp_bands =  2,3
```

Since no k-space integral is needed, set

```vi title="Input file"
berry_kmesh = 1 1 1
```

Although not used, we also need to input the value of the Fermi
level in eV

```vi title="Input file"
fermi_energy = [insert your value here]
```

On output, the program generates three files, namely
`SEED-kdotp_0.dat`, `SEED-kdotp_1.dat` and `SEED-kdotp_2.dat`, which
correspond to the zeroth, first and second order expansion coefficients, respectively. The dimension of the matrix contained in each file is $3^{l}\times N^{2}$, where $N$ is the number of bands set by `kdotp_num_bands`, and $l$ is the order of the expansion term (currently $l=0,1$ or $2$).

These coefficients can be used, among other things, to compute the
energy dispersion of the bands of interest around the chosen
k-point. The $k\cdot p$ band dispersion can be computed and plotted
along $k_x$ (from S to X) using python and the file `kdotp_plot.py`
provided in the example folder

```bash title="Terminal"
    python kdotp_plot.py
```

For comparison, the exact band structure calculated usingWannier90 (file `bc2n_band.dat`, generated automatically) is also plotted along (see the band dispersion [plot](#fig:bc2n-bnd)).

<figure markdown="span" id="fig:bc2n-bnd">
![Image title](./kdotp_bands_SX.webp){ width="500"}
<figcaption>Band dispersion of monolayer BC<sub>2</sub>N around <em>S</em> point. Exact results (solid dots) are compared to first-order (blue) and second-order (red) <em>k</em>⋅ <em>p</em> model results for valence and conduction bands.</figcaption>
</figure>

&nbsp;

[^1]: Once `XCrySDen` starts, click on `Tools` $\rightarrow$ `Data Grid`
    in order to specify an isosurface value to plot.

[^2]: Please note the following counterintuitive feature in `pwscf`: in
    order to obtain a ground state with magnetization along the
    *positive* $z$-axis, one should use a *negative* value for the
    variable `starting_magnetization`.

[^3]: The calculation of the AHC using `berry_task = kubo` involves a
    truncation of the sum over empty states in the Kubo-Greenwood
    formula: see description of the keyword [`kubo_eigval_max`]() in the
    User Guide. As discussed around [the formula for anomalous Hall conductivity]() of the User Guide, no truncation is done with `berry_task = ahc`.