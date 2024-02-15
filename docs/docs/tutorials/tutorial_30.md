# 30: Gallium Arsenide &#151; Frequency-dependent spin Hall conductivity {#gallium-arsenide-frequency-dependent-spin-hall-conductivity .unnumbered}

-   Outline: *Calculate the alternating current (ac) spin Hall
    conductivity of gallium arsenide considering spin-orbit coupling. To
    gain a better understanding of this tutorial, it is suggested to read
    Ref. [@qiao-prb2018] for a detailed description of the theory and the
    [berry_task=shc: spin Hall conductivity](../../user_guide/postw90/berry#sec:shc) 
    chapter of the User Guide.*

-   Directory: `tutorials/tutorial30/`

-   Input files

    -    `GaAs.scf` *The `pwscf` input file for ground state
        calculation*

    -    `GaAs.nscf` *The `pwscf` input file to obtain Bloch
        states on a uniform grid*

    -    `GaAs.pw2wan` *The input file for `pw2wannier90`*

    -    `GaAs.win` *The `wannier90` and `postw90` input file*

1.  Run `pwscf` to obtain the ground state of gallium
    arsenide

    ```bash title="Terminal"
    pw.x < GaAs.scf > scf.out
    ```

2.  Run `pwscf` to obtain the Bloch states on a uniform
    $k$-point grid

    ```bash title="Terminal"
    pw.x < GaAs.nscf > nscf.out
    ```

3.  Run `wannier90` to generate a list of the required overlaps (written
    into the `GaAs.nnkp` file)

    ```bash title="Terminal"
    wannier90.x -pp GaAs
    ```

4.  Run `pw2wannier90` to compute the overlaps between Bloch states and
    the projections for the starting guess (written in the `GaAs.mmn`
    and `GaAs.amn` files)

    ```bash title="Terminal"
    pw2wannier90.x < GaAs.pw2wan > pw2wan.out
    ```

5.  Run `wannier90` to compute the MLWFs

    ```bash title="Terminal"
    wannier90.x GaAs
    ```

6.  Run `postw90`

    ```bash title="Terminal"
    postw90.x GaAs # (1)! 
    mpirun -np 8 postw90.x GaAs # (2)! 
    ```

    1.     serial execution
    2.     example of parallel execution with 8 MPI processes

## ac spin Hall conductivity {#ac-spin-hall-conductivity .unnumbered}

The spin Hall conductivity is also dependent on the frequency $\omega$
in this [equation](../../user_guide/postw90/berry#mjx-eqn:eq:kubo_shc) of the User Guide. 
The direct current (dc) SHC calculated in the previous tutorial corresponds to
$\sigma_{\alpha\beta}^{\text{spin}\gamma}$ in the limit
$\omega\rightarrow
0$ and it is a real number. At finite frequency
$\sigma_{\alpha\beta}^{\text{spin}\gamma}$ acquires an imaginary part.

To compute the ac spin Hall conductivity for $\hbar\omega$ up to 8 eV,
add the lines

```vi title="Input file"
shc_freq_scan = true

kubo_freq_min = 0.0

kubo_freq_max = 8.0

kubo_freq_step = 0.01
```

and re-run `postw90`. The file `GaAs-shc-freqscan.dat` contains the
calculated ac SHC. Reasonably converged spectra can be obtained with a
$250\times 250\times 250$ $k$-point mesh. To plot the ac SHC, issue the
following commands

```bash title="Terminal"
gnuplot
```

```gnuplot title="Gnuplot shell"
plot 'GaAs-shc-freqscan.dat' u 2:3 w l title 'Re',
'GaAs-shc-freqscan.dat' u 2:4 w l title 'Im'
```

and then compare the result with Fig. 4 in Ref. [@qiao-prb2018] or the
Solution Booklet.

## Notes {#notes-2 .unnumbered}

-   When calculating ac SHC, adaptive smearing can be used by add the
    following keywords in the `GaAs.win`,

-   Adaptive kmesh refinement is not implemented for ac SHC calculation.

-   The first 10 semi-core states are excluded from the calculation by
    using the following keywords and in the case of GaAs disentanglement
    is not adopted so

-   Since the band gap is often under estimated by LDA/GGA calculations,
    a scissors shift is applied to recover the experimental band gap by
    using the following keywords or by


