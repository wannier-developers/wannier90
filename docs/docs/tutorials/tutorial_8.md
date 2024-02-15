# 8: Iron &#151; Spin-polarized WFs, DOS, projected WFs versus MLWFs {#iron-spin-polarized-wfs-dos-projected-wfs-versus-mlwfs .unnumbered}

-   Outline: *Generate both maximally-localized and projected Wannier
    functions for ferromagnetic bcc Fe. Calculate the total and
    orbital-projected density of states by Wannier interpolation.*

-   Directory: `tutorials/tutorial08/` *Files can be downloaded from [here](https://github.com/wannier-developers/wannier90/tutorials/tutorial08)*

-   Input Files

    -    `iron.scf` *The `pwscf` input file for the
        spin-polarized ground state calculation*

    -    `iron.nscf` *The `pwscf` input file to obtain Bloch
        states on a uniform grid*

    -    `iron_{up,down}.pw2wan` *Input files for `pw2wannier90`*

    -    `iron_{up,down}.win` *Input files for `wannier90` and
        ` postw90`*

-   Note that in a spin-polarized calculation the spin-up and spin-down
    MLWFs are computed separately. (The more general case of spinor WFs
    will be treated in Tutorial [17](../tutorial_17#iron-spin-orbit-coupled-bands-and-fermi-surface-contours).

1.  Run `pwscf` to obtain the ferromagnetic ground state of
    bcc Fe

    ```bash title="Terminal"
    pw.x < iron.scf > scf.out
    ```

2.  Run `pwscf` to obtain the Bloch states on a uniform
    k-point grid

    ```bash title="Terminal"
    pw.x < iron.nscf > nscf.out
    ```

3.  Run `wannier90` to generate a list of the required overlaps (written
    into the `.nnkp` files).

    ```bash title="Terminal"
    wannier90.x -pp iron_up
    wannier90.x -pp iron_dn
    ```

4.  Run `pw2wannier90` to compute the overlap between Bloch states and
    the projections for the starting guess (written in the `.mmn` and
    `.amn` files).

    ```bash title="Terminal"
    pw2wannier90.x < iron_up.pw2wan > pw2wan_up.out
    pw2wannier90.x < iron_dn.pw2wan > pw2wan_dn.out
    ```

5.  Run `wannier90` to compute the MLWFs.

    ```bash title="Terminal"
    wannier90.x iron_up
    wannier90.x iron_dn
    ```

### Density of states {#density-of-states .unnumbered}

To compute the DOS using a $25\times 25 \times 25$ $k$-point grid add to
the two `.win` files

```vi title="Input file"
dos = true
dos_kmesh = 25
```

run `postw90`,

```vi title="Input file"
postw90.x iron_up
postw90.x iron_dn
```

and plot the DOS with `gnuplot`,

```bash title="Terminal"
gnuplot
```

```gnuplot title="Gnuplot shell"
plot 'iron_up_dos.dat' u (-\$2):(\$1-12.6256) w
l,'iron_dn_dos.dat' u 2:(\$1-12.6256) w l
```

Energies are referred to the Fermi level (12.6256 eV, from ` scf.out`).
Note the exchange splitting between the up-spin and down-spin DOS. Check
the convergence by repeating the DOS calculations with more $k$-points.

### Projected versus maximally-localized Wannier functions {#projected-versus-maximally-localized-wannier-functions .unnumbered}

In the calculations above we chose $s$, $p$, and $d$-type trial orbitals
in the `.win` files,

```vi title="Input file"
Fe:s;p;d
```

Let us analyze the evolution of the WFs during the gauge-selection step.
Open one of the `.wout` files and search for "`Initial state`"; those
are the *projected* WFs. As expected they are atom-centred, with spreads
organized in three groups, 1+3+5: one $s$, three $p$, and five $d$. Now
look at the final state towards the end of the file. The Wannier spreads
have re-organized in two groups, 6+3; moreover, the six more diffuse WFs
are off-centred: the initial atomic-like orbitals hybridized with one
another, becoming more localized in the process. It is instructive to
visualize the final-state MLWFs using `XCrySDen`, following Tutorial 
[1](../tutorial_1#gallium-arsenide-mlwfs-for-the-valence-bands).
For more details, see Sec. IV.B of Ref. [@wang-prb06].

Let us plot the evolution of the spread functional $\Omega$,

```bash title="Terminal"
grep SPRD iron_up.wout > sprd_up

gnuplot
```

```gnuplot title="Gnuplot shell"
plot 'sprd_up' u 6 w l
```


The first plateau corresponds to atom-centred WFs of separate $s$, $p$,
and $d$ character, and the sharp drop signals the onset of the
hybridization. With hindsight, we can redo steps 4 and 5 more
efficiently using trial orbitals with the same character as the final
MLWFs,


```vi title="Input file"
Fe:sp3d2;dxy;dxz,dyz
```

With this choice the minimization converges much more rapidly.

Any reasonable set of localized WFs spanning the states of interest can
be used to compute physical quantities (they are "gauge-invariant"). Let
us recompute the DOS using, instead of MLWFs, the WFs obtained by
projecting onto $s$, $p$, and $d$-type trial orbitals, without further
iterative minimization of the spread functional. This can be done by
setting

```vi title="Input file"
num_iter = 0
```

But note that we still need to do disentanglement! Recalculate the DOS
to confirm that it is almost identical to the one obtained earlier using
the hybridized set of MLWFs. Visualize the projected WFs using
`XCrySDen`, to see that they retain the pure orbital character of the
individual trial orbitals.

### Orbital-projected DOS and exchange splitting {#orbital-projected-dos-and-exchange-splitting .unnumbered}

With projected WFs the total DOS can be separated into $s$, $p$ and $d$
contributions, in a similar way to the orbital decomposition of the
energy bands in Tutorial [4](../tutorial_4#copper-fermi-surface-orbital-character-of-energy-bands).

In order to obtain the partial DOS projected onto the $p$-type WFs, add
to the `.win` files

```vi title="Input file"
dos_project = 2,3,4
```

and re-run `postw90`. Plot the projected DOS for both up- and down-spin
bands. Repeat for the $s$ and $d$ projections.

Projected WFs can also be used to quantify more precisely the exchange
splitting between majority and minority states. Re-run `wannier90` after
setting `dos=false` and adding to the `.win` files

```vi title="Input file"
write_hr_diag = true
```

This instructs `wannier90` to print in the output file the on-site
energies $\langle {\bf 0}n\vert H\vert {\bf 0}n\rangle$. The difference
between corresponding values in `iron_up.wout` and in `iron_dn.wout`
gives the exchange splittings for the individual orbitals. Compare their
magnitudes with the splittings displayed by the orbital-projected DOS
plots. In agreement with the Stoner criterion, the largest exchange
splittings occur for the localized $d$-states, which contribute most of
the density of states at the Fermi level.

<figure markdown="span" id="fig:Fe-sprd">
![Image title](./Fe-spread.webp){ width="500" }
<figcaption markdown="span"> Fig.3: Evolution of the Wannier spread $\Omega$ of the minority (spin-up) bands of
bcc Fe during the iterative minimization of $\widetilde{\Omega}$, starting from s, p and
d-type trial orbitals.</figcaption>
</figure>
