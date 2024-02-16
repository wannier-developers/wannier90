# 18: Iron &#151; Berry curvature, anomalous Hall conductivity and optical conductivity {#iron-berry-curvature-anomalous-hall-conductivity-and-optical-conductivity .unnumbered}

Note: This tutorial requires a recent version of the `pw2wannier90`
interface.

-   Outline: *Calculate the Berry curvature, anomalous Hall
    conductivity, and (magneto)optical conductivity of ferromagnetic bcc
    Fe with spin-orbit coupling. In preparation for this tutorial it may
    be useful to read Ref. [@yao-prl04] and Ch. 11 of the User Guide.*

-   Directory: `tutorials/tutorial18/` *Files can be downloaded from [here](https://github.com/wannier-developers/wannier90/tree/develop/tutorials/tutorial18)*

-   Input files

    -    `Fe.scf` *The `pwscf` input file for ground state
        calculation*

    -    `Fe.nscf` *The `pwscf` input file to obtain Bloch
        states on a uniform grid*

    -    `Fe.pw2wan` *The input file for `pw2wannier90`*

    -    `Fe.win` *The `wannier90` and `postw90` input file*

The sequence of steps below is the same of Tutorial [17](tutorial_17.md#iron-spin-orbit-coupled-bands-and-fermi-surface-contours). If you have
already run that example, you can reuse the output files from steps
1 &#151; 5, and only step 6 must be carried out again using the new input file
`Fe.win`.

1.  Run `pwscf` to obtain the ground state of iron

    ```bash title="Terminal"
    pw.x < Fe.scf > scf.out
    ```

2.  Run `pwscf` to obtain the Bloch states on a uniform
    k-point grid

    ```bash title="Terminal"
    pw.x < Fe.nscf > nscf.out
    ```

3.  Run `wannier90` to generate a list of the required overlaps (written
    into the `Fe.nnkp` file)

    ```bash title="Terminal"
    wannier90.x -pp Fe
    ```

4.  Run `pw2wannier90` to compute the overlaps between Bloch states and
    the projections for the starting guess (written in the `Si.mmn` and
    `Si.amn` files)

    ```bash title="Terminal"
    pw2wannier90.x < Fe.pw2wan > pw2wan.out
    ```

5.  Run `wannier90` to compute the MLWFs

    ```bash title="Terminal"
    wannier90.x Fe
    ```

6.  Run `postw90` 

    ```bash title="Terminal"
    postw90.x Fe # (1)! 
    mpirun -np 8 postw90.x Fe # (2)!
    ```

    1.    serial execution
    2.    example of parallel execution with 8 MPI processes

## Berry curvature plots {#berry-curvature-plots .unnumbered}

The Berry curvature $\Omega_{\alpha\beta}({\bf k})$ of the occupied
states is defined in this [equation](../user_guide/postw90/berry.md#mjx-eqn:eq:ahc)  of the User Guide. The following lines
in `Fe.win` are used to calculate the energy bands and the Berry
curvature (in bohr$^2$) along high-symmetry lines in $k$-space.

```vi title="Input file"
fermi_energy = [insert your value here]

berry_curv_unit = bohr2

kpath = true

kpath_task = bands+curv

kpath_bands_colour = spin

kpath_num_points = 1000
```

After executing `postw90`, plot the Berry curvature component
$\Omega_z({\bf k})=\Omega_{xy}({\bf k})$ along the magnetization
direction using the script generated at runtime,

```bash title="Terminal"
python Fe-bands+curv_z.py
```

and compare with Fig. 2 of Ref. [@yao-prl04].

In Tutorial [17](tutorial_17.md#iron-spin-orbit-coupled-bands-and-fermi-surface-contours) we plotted the Fermi lines on the (010) plane $k_y=0$. To
combine them with a heatmap plot of (minus) the Berry curvature set
`kpath = false`, uncomment the following lines in `Fe.win`, re-run
`postw90`, and issue

```bash title="Terminal"
python Fe-kslice-curv_z+fermi_lines.py
```

Compare with Fig. 3 in Ref. [@yao-prl04]. Note how the Berry curvature
"hot-spots" tend to occur near spin-orbit-induced avoided crossings (the
Fermi lines with and without spin-orbit were generated in Tutorial [17](tutorial_17.md#iron-spin-orbit-coupled-bands-and-fermi-surface-contours)).

## Anomalous Hall conductivity {#anomalous-hall-conductivity .unnumbered}

The intrinsic anomalous Hall conductivity (AHC) is proportional to the
BZ integral of the Berry curvature. In bcc Fe with the magnetization
along $\hat{\bf z}$, the only nonzero components are
$\sigma_{xy}=-\sigma_{yx}$. To evaluate the AHC using a $25\times
25\times 25$ $k$-point mesh, set `kslice = false`, uncomment the
following lines in `Fe.win`,

```vi title="Input file"
berry = true

berry_task = ahc

berry_kmesh = 25 25 25
```

and re-run `postw90`. The AHC is written in the output file `Fe.wpout`
in vector form. For bcc Fe with the magnetization along \[001\], only
the $z$-component $\sigma_{xy}$ is nonzero.

As a result of the strong and rapid variations of the Berry curvature
across the BZ, the AHC converges rather slowly with $k$-point sampling,
and a $25\times 25\times 25$ does not yield a well-converged value.

\-   Increase the BZ mesh density by changing ` berry_kmesh`.

\-   To accelerate the convergence, adaptively refine the mesh around
    spikes in the Berry curvature, by adding to ` Fe.win` the lines

This adds a $5\times 5\times 5$ fine mesh around those points where
$\vert{\bm \Omega}({\bf k})\vert$ exceeds 100 bohr$^2$. The percentage
of points triggering adaptive refinement is reported in `Fe.wpout`.

Compare the converged AHC value with those obtained in
Refs. [@wang-prb06] and [@yao-prl04].

The Wannier-interpolation formula for the Berry curvature comprises
three terms, denoted $D$-$D$, $D$-$\overline{A}$, and
$\overline{\Omega}$ in Ref. [@wang-prb06], and $J2$, $J1$, and $J0$ in
Ref. [@lopez-prb12]. To report in `Fe.wpout` the decomposition of the
total AHC into these three terms, set ` iprint` (verbosity level) to a
value larger than one in ` Fe.win`.

## Optical conductivity {#optical-conductivity .unnumbered}

The optical conductivity tensor of bcc Fe with magnetization along
$\hat{\bf z}$ has the form

$$
\bm{\sigma}=\bm{\sigma}^{\rm S}+\bm{\sigma}^{\rm A}=
\left(
\begin{array}{ccc}
\sigma_{xx} & 0 & 0\\
0 & \sigma_{xx} & 0\\
0 & 0 & \sigma_{zz}
\end{array}
\right)+
\left(
\begin{array}{ccc}
0 & \sigma_{xy} & 0 \\
-\sigma_{xy} & 0 & 0\\
0 & 0 & 0
\end{array}
\right),
$$

where "S" and "A" stand for the symmetric and antisymmetric
parts and $\sigma_{xx}=\sigma_{yy}\not=\sigma_{zz}$. The dc AHC
calculated earlier corresponds to $\sigma_{xy}$ in the limit
$\omega\rightarrow
0$. At finite frequency $\sigma_{xy}=-\sigma_{yx}$ acquires an imaginary
part which describes magnetic circular dichoism (MCD).

To compute the complex optical conductivity for $\hbar\omega$ up to
7 eV, replace

```vi title="Input file"
berry_task = ahc
```

with


```vi title="Input file"
berry_task = kubo
```

add the line

```vi title="Input file"
kubo_freq_max = 7.0
```

and re-run `postw90`. Reasonably converged spectra can be obtained with
a $125\times 125\times 125$ $k$-point mesh. Let us first plot the ac AHC
in S/cm, as in the lower panel of Fig. 5 in Ref. [@yao-prl04],

```bash title="Terminal"
gnuplot
```

```gnuplot title="Gnuplot shell"
plot 'Fe-kubo_A_xy.dat' u 1:2 w l
```

Comapare the $\omega\rightarrow 0$ limit with the result obtained
earlier by integrating the Berry curvature.[^3]

Next we plot the MCD spectrum. Following Ref. [@yao-prl04], we plot
${\rm Im}[\omega\sigma_{xy}(\hbar\omega)]$, in units of
$10^{29}$ sec$^{-2}$. The needed conversion factor is $9\times
10^{-18}\times e/\hbar\simeq 0.0137$ ($e$ and $\hbar$ in SI units),


```gnuplot title="Gnuplot shell"
set yrange[-5:15]
plot 'Fe-kubo_A_xy.dat' u 1:(\$1)\*(\$3)\*0.0137 w l
```

## Further ideas {#further-ideas-6 .unnumbered}

\-   Recompute the AHC and optical spectra of bcc Fe using projected $s$,
    $p$, and $d$-type Wannier functions instead of the hybridrized MLWFs
    (see Tutorial [8](tutorial_8.md#iron-spin-polarized-wfs-dos-projected-wfs-versus-mlwfs)), 
    and compare the results.

\-   A crude way to model the influence of heterovalent alloying on the
    AHC is to assume that its only effect is to donate or deplete
    electrons, i.e., to shift the Fermi level of the pure
    crystal [@yao-prb07]. Recalculate the AHC of bcc Fe for a range of Fermi energies within
    $\pm 0.5$ eV of the true Fermi level. This calculation can be
    streamlined by replacing in `Fe.win`

```vi title="Input file"
fermi_energy = [insert your value here]
```

with

```vi title="Input file"
fermi_energy_min = [insert here your value minus 0.5]

fermi_energy_max = [insert here your value plus 0.5]
```

Use a sufficiently dense BZ mesh with adaptive refinement. To plot
$\sigma_{xy}$ versus $\varepsilon_F$, issue

```bash title="Terminal"
gnuplot
```

```gnuplot title="Gnuplot shell"
plot 'Fe-ahc-fermiscan.dat' u 1:4 w lp
```
