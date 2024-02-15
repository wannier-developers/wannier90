# 27: Silicon &#151; Selected columns of density matrix algorithm for automated MLWFs {#silicon-selected-columns-of-density-matrix-algorithm-for-automated-mlwfs .unnumbered}

Note: This tutorial requires a recent version of the `pw2wannier90.x`
post-processing code of `Quantum ESPRESSO` (v6.4 or above).

-   Outline: *For bulk crystalline Silicon, generate the $A_{mn}$
    matrices via the selected columns of density matrix (SCDM) algorithm
    and the corresponding MLWFs for 1) Valence bands 2) Valence bands
    and 4 low-lying conduction bands 3) Conduction bands only. To better
    understand the input files and the results of these calculations, it
    is crucial that the Reader has familiarized with the concepts and
    methods explained in Ref. [@LinLin-ArXiv2017]. More info on the
    keywords related to the SCDM method may be found in the user_guide.*

-   Directory: `tutorials/tutorial27/` *Files can be downloaded from [here](https://github.com/wannier-developers/wannier90/tutorials/tutorial27)*

-   Input Files: `input_files`, and in the three subfolders
    `isolated, erfc` and `gaussian`. The `input_files` folder contains:

    -    `si.scf` *The `pwscf` input file for the ground
        state calculation*

    -    `si_4bands.nscf ` *The `pwscf` input file to obtain
        Bloch states on a uniform grid for 4 bands.*

    -    `si_12bands.nscf ` *The `pwscf` input file to
        obtain Bloch states on a uniform grid for 12 bands.*

-    Whereas the three subfolders `isolated, erfc` and `gaussian`
    contain the `si.win` `wannier90`  input files and `si.pw2wan`
    `pw2wannier90` input files each corresponding to one of the
    scenarios listed in the outline.


###Valence bands 
In this case we will compute 4 localized WFs
corresponding to the 4 valence bands of Silicon. These 4 bands
constitute a manifold that is separated in energy from other bands.
In this case the columns of the density matrix are already localized
in real space and no extra parameter is required.

1.  Copy the input files `si.scf` and `si_4bands.nscf` from the
    `input_files` directory into the `isolated` folder

2.  Run `pwscf` to obtain the ground state charge of
    bulk Silicon.

    ```bash title="Terminal"
    pw.x < si.scf > scf.out
    ```

3.  Run `pwscf` to obtain the Bloch states on a uniform
    k-point grid of 4x4x4 for 4 bands.

    ```bash title="Terminal"
    pw.x < si_4bands.nscf > nscf.out
    ```

4.  Inspect the `si.win` input file and make sure that the
    `auto_projections` flag is set to `.true.`. Also, make sure that
    no projections block is present.

5.  Run `wannier90` to generate a list of the required overlaps and
    also info on the SCDM method (written into the `si.nnkp` file).

    ```bash title="Terminal"
    wannier90.x -pp si
    ```

6.  Inspect the `si.nnkp` file and make sure you find the
    `auto_projections` block and that no projections have been
    written in the `projections` block.

7.  Inspect the `.pw2wan` input file. You will find two new
    keywords, i.e. `scdm_proj` and `scdm_entanglement`. The former,
    will instruct `pw2wannier90.x` to use the SCDM method when
    generating the $A_{mn}$ matrix. The latter, defines which
    formula to adopt for the function $f(\varepsilon_{n\mathbf{k}})$
    (see [@LinLin-ArXiv2017] and point below).

8.  Run `pw2wannier90` to compute the overlap between Bloch states
    and the projections via the SCDM method (written in the `si.mmn`
    and `si.amn` respectively).

    ```bash title="Terminal"
    pw2wannier90.x < si.pw2wan > pw2wan.out
    ```

9.  Run `wannier90` to compute the MLWFs.

    ```bash title="Terminal"
    wannier90.x si
    ```

At this point, you should have obtained 4 Wannier functions and
the interpolated valence bands for Silicon. Inspect the output
file `si.wout`. In particular, look at the geometric centres of
each WF, do they lie at the centre of the Si-Si bond as for the
MLWFs computed from user-defined initial $s$-like projections
(see Tutorial [11](../tutorial_11#silicon-valence-and-low-lying-conduction-states))? Plot these WFs using Vesta. Do they show the
$\sigma$ character one would expect from chemical arguments?

###Valence bands + conduction bands 
In this case we will compute 8
localized WFs corresponding to the 4 valence bands and 4 low-lying
conduction bands. Here, we don't have a separate manifold, since the
conduction bands are entangled with other high-energy bands and the
columns of the density matrix are not exponentially localized by
construction. A modified density matrix is required in this
case[@LinLin-ArXiv2017], and it is defined as:

$$
P(\mathbf{r},\mathbf{r}') = \sum_{n,\mathbf{k}} \psi_{n\mathbf{k}}(\mathbf{r})f(\varepsilon_{n,\mathbf{k}})\psi_{n\mathbf{k}}^\ast(\mathbf{r}'),
$$

where $\psi_{n\mathbf{k}}$ and $\varepsilon_{n,\mathbf{k}}$ are the
energy eigestates and eigenvalues from the first-principle
calculation respectively. The function
$f(\varepsilon_{n,\mathbf{k}})$ contains two free parameters $\mu$
and $\sigma$ and is defined as a complementary error function:

$$
f(\varepsilon_{n,\mathbf{k}}) = \frac{1}{2}\mathrm{erfc}\left(\frac{\varepsilon_{n,\mathbf{k}} - \mu}{\sigma}\right).
$$

1.  Copy the input files `si.scf` and `si_12bands.nscf` from the
    `input_files` folder into the `erfc` folder

2.  Run `pwscf` to obtain the ground state charge of
    bulk Silicon.

    ```bash title="Terminal"
    pw.x < si.scf > scf.out
    ```

3.  Run `pwscf` to obtain the Bloch states on a uniform
    k-point grid of 4x4x4 for 12 bands this time.

    ```bash title="Terminal"
    pw.x < si_12bands.nscf > nscf.out
    ```

4.  Inspect the `si.win` input file and make sure that the
    `auto_projections` flag is set to `.true.`. Also, make sure that
    no projection block is present.

5.  Run `wannier90` to generate a list of the required overlaps and
    also info on the SCDM method (written into the `si.nnkp` file).

    ```bash title="Terminal"
    wannier90.x -pp si
    ```

6.  Inspect the `si.nnkp` file and make sure you find the
    `auto_projections` block and that no projections have been
    written in the `projections` block.

7.  Inspect the `.pw2wan` input file. You will find other two new
    keywords, i.e. `scdm_mu` and `scdm_sigma`. These are the values
    in eV of $\mu$ and $\sigma$ in $f(\varepsilon_{n,\mathbf{k}})$,
    respectively.

8.  Run `pw2wannier90` to compute the overlap between Bloch states
    and the projections via the SCDM method (written in the `si.mmn`
    and `si.amn` respectively).

    ```bash title="Terminal"
    pw2wannier90.x < si.pw2wan > pw2wan.out
    ```

9.  Run `wannier90` to compute the MLWFs.

    ```bash title="Terminal"
    wannier90.x si
    ```

At this point, you should have obtained 8 localized Wannier
functions and the interpolated valence and conduction bands for
Silicon. Again, compare the results for the geometric centres
and the individual spreads with the ones from Tutorial [11](../tutorial_11#silicon-valence-and-low-lying-conduction-states). Is the
final value of total spread bigger or smaller than the one from
Tutorial [11](../tutorial_11#silicon-valence-and-low-lying-conduction-states)? 
Look at the WFs with Vesta. Can you explain what you
see? Where do the major lobes of the $sp3$-like WFs point in
this case?

###Conduction bands only 
In this case we will compute 4 localized WFs
corresponding to the 4 low-lying conduction bands only. As for the
previous point, we need to define a modified density
matrix[@LinLin-ArXiv2017]. Since we are only interested in a subset
of the conduction states, within a bounded energy region, a good
choice for $f(\varepsilon_{n,\mathbf{k}})$ is:

$$
f(\varepsilon_{n,\mathbf{k}}) = \exp\left(-\frac{(\varepsilon_{n,\mathbf{k}} - \mu)^2}{\sigma^2}\right).
$$


1.  Copy the input files `si.scf` and `si_12bands.nscf` from the
    `input_files` directory into the `gaussian` folder

2.  Run `pwscf` to obtain the ground state charge of
    bulk Silicon.

    ```bash title="Terminal"
    pw.x < si.scf > scf.out
    ```

3.  Run `pwscf` to obtain the Bloch states on a uniform
    k-point grid of 4x4x4 for 12 bands this time.

    ```bash title="Terminal"
    pw.x < si_12bands.nscf > nscf.out
    ```

4.  Inspect the `si.win` input file and make sure that the
    `auto_projections` flag is set to `.true.`. Also, make sure that
    no projections block is present.

5.  Run `wannier90` to generate a list of the required overlaps and
    also info on the SCDM method (written into the `si.nnkp` file).

    ```bash title="Terminal"
    wannier90.x -pp si
    ```

6.  Inspect the `si.nnkp` file and make sure you find the
    `auto_projections` block and that no projections have been
    written in the `projections` block.

7.  Run `pw2wannier90` to compute the overlap between Bloch states,
    the projections for the starting guess via the SCDM method
    (written in the `si.mmn` and `si.amn` respectively).

    ```bash title="Terminal"
    pw2wannier90.x < si.pw2wan > pw2wan.out
    ``` 

8.  Run `wannier90` to compute the MLWFs.

    ```bash title="Terminal"
    wannier90.x si
    ```

At this point, you should have obtained 4 localized Wannier
functions and the interpolated conduction bands for Silicon.
From chemical intuition, we would expect these functions to be
similar to anti-bonding orbitals of molecules with tetrahedral
symmetry. Plot the WFs and check if this is confirmed.


