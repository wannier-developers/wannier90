# Projections

## Specification of projections in `seedname.win` {#sec:proj}

Here we describe the projection functions used to construct the initial
guess $A_{mn}^{(\mathbf{k})}$ for the unitary transformations.

Each projection is associated with a site and an angular momentum state
defining the projection function. Optionally, one may define, for each
projection, the spatial orientation, the radial part, the diffusivity,
and the volume over which real-space overlaps $A_{mn}$ are calculated.

The code is able to

1.  project onto s,p,d and f angular momentum states, plus the hybrids
    sp, sp$^2$, sp$^3$, sp$^3$d, sp$^3$d$^2$.

2.  control the radial part of the projection functions to allow higher
    angular momentum states, e.g., both 3s and 4s in silicon.

The atomic orbitals of the hydrogen atom provide a good basis to use for
constructing the projection functions: analytical mathematical forms
exist in terms of the good quantum numbers $n$, $l$ and $m$; hybrid
orbitals (sp, sp$^{2}$, sp$^{3}$, sp$^{3}$d etc.) can be constructed by
simple linear combination $|\phi\rangle =
\sum_{nlm} C_{nlm}|nlm\rangle$ for some coefficients $C_{nlm}$.

The angular functions that use as a basis for the projections are not
the canonical spherical harmonics $Y_{lm}$ of the hydrogenic Schrödinger
equation but rather the *real* (in the sense of non-imaginary) states
$\Theta_{lm_{\mathrm{r}}}$, obtained by a unitary transformation. For
example, the canonical eigenstates associated with $l=1$, $m=\{-1,0,1\}$
are not the real p$_{x}$, p$_{y}$ and p$_{z}$ that we want. See
Section [3.4](#sec:orbital-defs){reference-type="ref"
reference="sec:orbital-defs"} for our mathematical conventions regarding
projection orbitals for different $n$, $l$ and $m_{\mathrm{r}}$.

We use the following format to specify projections in `<seedname>.win`:

`Begin Projections`\
`[units]`\
`site:ang_mtm:zaxis:xaxis:radial:zona`\
`    `⋮\
`End Projections`

Notes:

`units`:\
Optional. Either `Ang` or `Bohr` to specify whether the projection
centres specified in this block (if given in Cartesian co-ordinates) are
in units of Angstrom or Bohr, respectively. The default value is `Ang`.

`site`:\
`C`, `Al`, etc. applies to all atoms of that type\
`f=0,0.50,0` -- centre on (0.0,0.5,0.0) in **f**ractional coordinates
(crystallographic units) relative to the direct lattice vectors\
`c=0.0,0.805,0.0` -- centre on (0.0,0.805,0.0) in **C**artesian
coordinates in units specified by the optional string `units` in the
first line of the projections block (see above).

`ang_mtm`:\
Angular momentum states may be specified by `l` and `mr`, or by the
appropriate character string. See
Tables [3.1](#tab:angular){reference-type="ref" reference="tab:angular"}
and [3.2](#tab:hybrids){reference-type="ref" reference="tab:hybrids"}.
Examples:\
`l=2,mr=1 ` or ` dz2` -- a single projection with $l=2$,
$m_{\textrm{r}}=1$ (i.e., d$_{z^{2}}$)\
`l=2,mr=1,4 ` or ` dz2,dx2-y2` -- two functions: d$_{z^{2}}$ and
d$_{xz}$\
`l=-3 ` or ` sp3` -- four sp$^{3}$ hybrids\
Specific hybrid orbitals may be specified as follows:\
`l=-3,mr=1,3 ` or ` sp3-1,sp3-3` -- two specific sp$^{3}$ hybrids\
Multiple states may be specified by separating with '`;`', e.g.,\
`sp3;l=0 ` or ` l=-3;l=0` -- four sp$^{3}$ hybrids and one s orbital

`zaxis` (optional):\
`z=1,1,1` -- set the $z$-axis to be in the (1,1,1) direction. Default is
`z=0,0,1`

`xaxis` (optional):\
`x=1,1,1` -- set the $x$-axis to be in the (1,1,1) direction. Default is
`x=1,0,0`

`radial` (optional):\
`r=2` -- use a radial function with one node (ie second highest
pseudostate with that angular momentum). Default is `r=1`. Radial
functions associated with different values of `r` should be orthogonal
to each other.

`zona` (optional):\
`zona=2.0` -- the value of $\frac{Z}{a}$ for the radial part of the
atomic orbital (controls the diffusivity of the radial function). Units
always in reciprocal Angstrom. Default is `zona=1.0`.

**Examples**

1\. CuO, s,p and d on all Cu; sp$^3$ hybrids on O:

`Cu:l=0;l=1;l=2 `

`O:l=-3 ` or ` O:sp3`

2\. A single projection onto a p$_z$ orbital orientated in the (1,1,1)
direction:

`c=0,0,0:l=1,mr=1:z=1,1,1 ` or ` c=0,0,0:pz:z=1,1,1`

3\. Project onto s, p and d (with no radial nodes), and s and p (with
one radial node) in silicon:

`Si:l=0;l=1;l=2`

`Si:l=0;l=1:r=2`

## Spinor Projections

When `spinors=.true.` it is possible to select a set of localised
functions to project onto 'up' states and a set to project onto 'down'
states where, for complete flexibility, it is also possible to set the
local spin quantisation axis.

Note, however, that this feature requires a recent version of the
interface between the ab-initio code and Wannier90 (i.e., written after
the release of the 2.0 version, in October 2013) supporting spinor
projections.

`Begin Projections`\
`[units]`\
`site:ang_mtm:zaxis:xaxis:radial:zona(spin)[quant_dir]`\
`    `⋮\
`End Projections`

`spin` (optional):\
Choose projection onto 'up' or 'down' states\
`u` -- project onto 'up' states.\
`d` -- project onto 'down' states.\
Default is `u,d`

`quant_dir` (optional):\
`1,0,0` -- set the spin quantisation axis to be in the (1,0,0)
direction. Default is `0,0,1`

**Examples**

-   18 projections on an iron site

    `Fe:sp3d2;dxy;dxx;dyz`

-   same as above

    `Fe:sp3d2;dxy;dxx;dyz(u,d)`

-   same as above

    `Fe:sp3d2;dxy;dxz;dyz(u,d)[0,0,1]`

-   same as above but quantisation axis is now x

    `Fe:sp3d2;dxy;dxz;dyz(u,d)[1,0,0]`

-   now only 9 projections onto up states

    `Fe:sp3d2;dxy;dxz;dyz(u)`

-   9 projections onto up-states and 3 on down

    `Fe:sp3d2;dxy;dxz;dyz(u) `\
    `Fe:dxy;dxz;dyz(d)`

-   projections onto alternate spin states for two lattice sites (Cr1,
    Cr2)

    `Cr1:d(u)`\
    `Cr2:d(d)`

## Short-Cuts

### Random projections

It is possible to specify the projections, for example, as follows:

`Begin Projections`\
`random`\
`C:sp3`\
`End Projections`

in which case `wannier90` uses four sp$^3$ orbitals centred on each C
atom and then chooses the appropriate number of randomly-centred s-type
Gaussian functions for the remaining projection functions. If the block
only consists of the string `random` and no specific projection centres
are given, then all of the projection centres are chosen randomly.

### Bloch phases

Setting `use_bloch_phases = true` in the input file absolves the user of
the need to specify explicit projections. In this case, the Bloch
wave-functions are used as the projection orbitals, namely
$A_{mn}^{(\mathbf{k})} =
\langle\psi_{m\mathbf{k}}|\psi_{n\mathbf{k}}\rangle = \delta_{mn}$.

## Orbital Definitions {#sec:orbital-defs}

The angular functions $\Theta_{lm_{\mathrm{r}}}(\theta,\varphi)$
associated with particular values of $l$ and $m_{\mathrm{r}}$ are given
in Tables [3.1](#tab:angular){reference-type="ref"
reference="tab:angular"} and [3.2](#tab:hybrids){reference-type="ref"
reference="tab:hybrids"}.

The radial functions $R_{\mathrm{r}}(r)$ associated with different
values of $r$ should be orthogonal. One choice would be to take the set
of solutions to the radial part of the hydrogenic Schrödinger equation
for $l=0$, i.e., the radial parts of the 1s, 2s, 3s... orbitals, which
are given in Table [3.3](#tab:radial){reference-type="ref"
reference="tab:radial"}.

::: center
::: {#tab:angular}
  ----- ------------------ -------------- ---------------------------------------------------------------------------------------------
                                          
   $l$   $m_{\mathrm{r}}$       Name                               $\Theta_{lm_{\mathrm{r}}}(\theta,\varphi)$
                                          
                                          
    0           1               `s`                                          $\frac{1}{\sqrt{4\pi}}$
                                          
                                          
    1           1               `pz`                                    $\sqrt{\frac{3}{4\pi}}\cos\theta$
                                          
    1           2               `px`                              $\sqrt{\frac{3}{4\pi}}\sin\theta\cos\varphi$
                                          
    1           3               `py`                              $\sqrt{\frac{3}{4\pi}}\sin\theta\sin\varphi$
                                          
                                          
    2           1              `dz2`                              $\sqrt{\frac{5}{16\pi}}(3\cos^{2}\theta -1)$
                                          
    2           2              `dxz`                         $\sqrt{\frac{15}{4\pi}}\sin\theta\cos\theta\cos\varphi$
                                          
    2           3              `dyz`                         $\sqrt{\frac{15}{4\pi}}\sin\theta\cos\theta\sin\varphi$
                                          
    2           4             `dx2-y2`                         $\sqrt{\frac{15}{16\pi}}\sin^{2}\theta\cos2\varphi$
                                          
    2           5              `dxy`                           $\sqrt{\frac{15}{16\pi}}\sin^{2}\theta\sin2\varphi$
                                          
                                          
    3           1              `fz3`                       $\frac{\sqrt{7}}{4\sqrt{\pi}}(5\cos^{3}\theta-3\cos\theta)$
                                          
    3           2              `fxz2`               $\frac{\sqrt{21}}{4\sqrt{2\pi}}(5\cos^{2}\theta-1)\sin\theta\cos\varphi$
                                          
    3           3              `fyz2`               $\frac{\sqrt{21}}{4\sqrt{2\pi}}(5\cos^{2}\theta-1)\sin\theta\sin\varphi$
                                          
    3           4           `fz(x2-y2)`               $\frac{\sqrt{105}}{4\sqrt{\pi}}\sin^{2}\theta\cos\theta\cos2\varphi$
                                          
    3           5              `fxyz`                 $\frac{\sqrt{105}}{4\sqrt{\pi}}\sin^{2}\theta\cos\theta\sin2\varphi$
                                          
    3           6           `fx(x2-3y2)`   $\frac{\sqrt{35}}{4\sqrt{2\pi}}\sin^{3}\theta(\cos^{2}\varphi-3\sin^{2}\varphi)\cos\varphi$
                                          
    3           7           `fy(3x2-y2)`   $\frac{\sqrt{35}}{4\sqrt{2\pi}}\sin^{3}\theta(3\cos^{2}\varphi-\sin^{2}\varphi)\sin\varphi$
                                          
  ----- ------------------ -------------- ---------------------------------------------------------------------------------------------

  : Angular functions $\Theta_{lm_{\mathrm{r}}}(\theta,\varphi)$
  associated with particular values of $l$ and $m_{\mathrm{r}}$ for
  $l\ge0$.
:::
:::

::: center
::: {#tab:hybrids}
  ----------------------- ------------------ ----------- -----------------------------------------------------------------------------
                                                         
            $l$            $m_{\mathrm{r}}$     Name                      $\Theta_{lm_{\mathrm{r}}}(\theta,\varphi)$
                                                         
                                                         
   $-$`<!-- -->`{=html}1          1            `sp-1`                  $\frac{1}{\sqrt{2}}$`s` $+\frac{1}{\sqrt{2}}$`px`
                                                         
   $-$`<!-- -->`{=html}1          2            `sp-2`                  $\frac{1}{\sqrt{2}}$`s` $-\frac{1}{\sqrt{2}}$`px`
                                                         
                                                         
   $-$`<!-- -->`{=html}2          1            `sp2-1`    $\frac{1}{\sqrt{3}}$`s` $-\frac{1}{\sqrt{6}}$`px` $+\frac{1}{\sqrt{2}}$`py`
                                                         
   $-$`<!-- -->`{=html}2          2            `sp2-2`    $\frac{1}{\sqrt{3}}$`s` $-\frac{1}{\sqrt{6}}$`px` $-\frac{1}{\sqrt{2}}$`py`
                                                         
   $-$`<!-- -->`{=html}2          3            `sp2-3`                 $\frac{1}{\sqrt{3}}$`s` $+\frac{2}{\sqrt{6}}$`px`
                                                         
                                                         
   $-$`<!-- -->`{=html}3          1            `sp3-1`                   $\frac{1}{2}$(`s` $+$ `px` $+$ `py` $+$ `pz`)
                                                         
   $-$`<!-- -->`{=html}3          2            `sp3-2`                   $\frac{1}{2}$(`s` $+$ `px` $-$ `py` $-$ `pz`)
                                                         
   $-$`<!-- -->`{=html}3          3            `sp3-3`                   $\frac{1}{2}$(`s` $-$ `px` $+$ `py` $-$ `pz`)
                                                         
   $-$`<!-- -->`{=html}3          4            `sp3-4`                   $\frac{1}{2}$(`s` $-$ `px` $-$ `py` $+$ `pz`)
                                                         
                                                         
   $-$`<!-- -->`{=html}4          1           `sp3d-1`    $\frac{1}{\sqrt{3}}$`s` $-\frac{1}{\sqrt{6}}$`px` $+\frac{1}{\sqrt{2}}$`py`
                                                         
   $-$`<!-- -->`{=html}4          2           `sp3d-2`    $\frac{1}{\sqrt{3}}$`s` $-\frac{1}{\sqrt{6}}$`px` $-\frac{1}{\sqrt{2}}$`py`
                                                         
   $-$`<!-- -->`{=html}4          3           `sp3d-3`                 $\frac{1}{\sqrt{3}}$`s` $+\frac{2}{\sqrt{6}}$`px`
                                                         
   $-$`<!-- -->`{=html}4          4           `sp3d-4`                $\frac{1}{\sqrt{2}}$`pz` $+\frac{1}{\sqrt{2}}$`dz2`
                                                         
   $-$`<!-- -->`{=html}4          5           `sp3d-5`               $-\frac{1}{\sqrt{2}}$`pz` $+\frac{1}{\sqrt{2}}$`dz2`
                                                         
                                                         
   $-$`<!-- -->`{=html}5          1           `sp3d2-1`             $\frac{1}{\sqrt{6}}\verb#s#-\frac{1}{\sqrt{2}}\verb#px#
                                                                   -\frac{1}{\sqrt{12}}\verb#dz2#+\frac{1}{2}\verb#dx2-y2#$
                                                         
   $-$`<!-- -->`{=html}5          2           `sp3d2-2`             $\frac{1}{\sqrt{6}}\verb#s#+\frac{1}{\sqrt{2}}\verb#px#
                                                                   -\frac{1}{\sqrt{12}}\verb#dz2#+\frac{1}{2}\verb#dx2-y2#$
                                                         
   $-$`<!-- -->`{=html}5          3           `sp3d2-3`             $\frac{1}{\sqrt{6}}\verb#s#-\frac{1}{\sqrt{2}}\verb#py#
                                                                   -\frac{1}{\sqrt{12}}\verb#dz2#-\frac{1}{2}\verb#dx2-y2#$
                                                         
   $-$`<!-- -->`{=html}5          4           `sp3d2-4`             $\frac{1}{\sqrt{6}}\verb#s#+\frac{1}{\sqrt{2}}\verb#py#
                                                                   -\frac{1}{\sqrt{12}}\verb#dz2#-\frac{1}{2}\verb#dx2-y2#$
                                                         
   $-$`<!-- -->`{=html}5          5           `sp3d2-5`             $\frac{1}{\sqrt{6}}\verb#s#-\frac{1}{\sqrt{2}}\verb#pz#
                                                                                +\frac{1}{\sqrt{3}}\verb#dz2#$
                                                         
   $-$`<!-- -->`{=html}5          6           `sp3d2-6`             $\frac{1}{\sqrt{6}}\verb#s#+\frac{1}{\sqrt{2}}\verb#pz#
                                                                                +\frac{1}{\sqrt{3}}\verb#dz2#$
                                                         
  ----------------------- ------------------ ----------- -----------------------------------------------------------------------------

  : Angular functions $\Theta_{lm_{\mathrm{r}}}(\theta,\varphi)$
  associated with particular values of $l$ and $m_{\mathrm{r}}$ for
  $l<0$, in terms of the orbitals defined in
  Table [3.1](#tab:angular){reference-type="ref"
  reference="tab:angular"}.
:::
:::

::: center
::: {#tab:radial}
  ---------- --------------------------------------------
             
     $r$                 $R_{\mathrm{r}}(r)$
             
             
      1            $2 \alpha^{3/2}\exp(-\alpha r)$
             
             
      2       $\frac{1}{2\sqrt{2}}\alpha^{3/2}(2-\alpha
                         r)\exp(-\alpha r/2)$
             
             
      3       $\sqrt{\frac{4}{27}}\alpha^{3/2}(1-2\alpha
              r/3+2\alpha^{2}r^{2}/27)\exp(-\alpha r/3)$
             
  ---------- --------------------------------------------

  :  One possible choice for the radial functions $R_{\mathrm{r}}(r)$
  associated with different values of $r$: the set of solutions to the
  radial part of the hydrogenic Schrödinger equation for $l=0$, i.e.,
  the radial parts of the 1s, 2s, 3s... orbitals, where
  $\alpha=Z/a={\tt zona}$.
:::
:::

## Projections via the SCDM-**k** method in pw2wannier90

For many systems, such as aperiodic systems, crystals with defects, or
novel materials with complex band structure, it may be extremely hard to
identify *a-priori* a good initial guess for the projection functions
used to generate the $A_{mn}^{(\mathbf{k})}$ matrices. In these cases,
one can use a different approach, known as the SCDM-**k**
method[@LinLin-ArXiv2017], based on a QR factorization with column
pivoting (QRCP) of the density matrix from the self-consistent field
calculation, which allows one to avoid the tedious step of specifying a
projection block altogether, hence to avoid . This method is robust in
generating well localised function with the correct spatial orientations
and in general in finding the global minimum of the spread functional
$\Omega$. Any electronic-structure code should in principle be able to
implement the SCDM-**k** method within their interface with Wannier90,
however at the moment (develop branch on the GitHub repository July
2019) only the Quantum ESPRESSO package has this capability implemented
in the `pw2wannier90` interface program. Moreover, the `pw2wannier90`
interface program supports also the SCDM-**k** method for
spin-noncollinear systems. The SCDM-**k** can operate in two modes:

1.  In isolation, i.e., without performing a subsequent Wannier90
    optimisation (not recommended). This can be achieved by setting
    `num_iter=0` and `dis_num_iter=0` in the `<seedname>.win` input
    file. The rationale behind this is that in general the projection
    functions obtained with the SCDM-**k** are already well localised
    with the correct spatial orientations. However, the spreads of the
    resulting functions are usually larger than the MLWFs ones.

2.  In combination with the Marzari-Vanderbilt (recommended option). In
    this case, the SCDM-**k** is only used to generate the initial
    $A_{mn}^{(\mathbf{k})}$ matrices as a replacement scheme for the
    projection block.

The following keywords need to be specified in the `pw2wannier90.x`
input file `<seedname>.pw2wan`: `scdm_proj` `scdm_entanglement`
`scdm_mu` `scdm_sigma`

## Projections via pseudo-atomic orbitals in pw2wannier90 {#sec:proj_pdwf}

When generating pseudopotentials, often the atomic wavefunctions of
isolated atom are pseudized and bundled together with the
pseudopotential files. These orbitals are often used for computing the
projectabilities, for instance, measuring orbital contributions to band
structures. Instead of manually specifying the initial projections in
the `projections` block, one can use these pseudo-atomic orbitals to
automate the initial projection process.

Currently (July 2023), this functionality is implemented in the
[quantum-espresso]{.smallcaps} interface, but in principle it can be
done in any other interface as well. In the following, we will use the
[quantum-espresso]{.smallcaps} interface as an example to illustrate the
whole procedure.

To activate pseudo-atomic orbital projection, one needs to set
`auto_projections = .true.` in the `win` file, and remove the
`projections` block.

Then in the `pw2wannier90` input file, one needs to add an additional
tag `atom_proj = .true.`. This will ask `pw2wannier90` to read the
pseudo-atomic orbitals from the pseudopotential files, and use them to
compute the `amn` file.

Some times, one may want to exclude semi-core states from
Wannierisation, for such cases, one can inspect the stdout of
`pw2wannier90`, which will print the orbitals used for computing `amn`,
e.g.,

` `

> -------------------------------------\
> \*\*\* Compute A with atomic projectors\
> -------------------------------------\
> Use atomic projectors from UPF\
> \
> (read from pseudopotential files):\
> state \# 1: atom 1 (C ), wfc 1 (l=0 m= 1)\
> state \# 2: atom 1 (C ), wfc 2 (l=1 m= 1)\
> state \# 3: atom 1 (C ), wfc 2 (l=1 m= 2)\
> state \# 4: atom 1 (C ), wfc 2 (l=1 m= 3)\
> state \# 5: atom 2 (C ), wfc 1 (l=0 m= 1)\
> state \# 6: atom 2 (C ), wfc 2 (l=1 m= 1)\
> state \# 7: atom 2 (C ), wfc 2 (l=1 m= 2)\
> state \# 8: atom 2 (C ), wfc 2 (l=1 m= 3)\

Here it shows that there are two carbon atoms, each with one $s$ and
three $p$ orbitals. If one wants to exclude specific orbital(s), there
is an additional input `atom_proj_exclude`, which accept a list of
integers, e.g.,

`atom_proj_exclude = 1 5`

which will exclude the two $s$ orbitals from computing `amn`.

#### Advanced usage

If the pseudopotential orbitals are not enough, one could also generate
a custom set of orbitals, and ask `pw2wannier90` to use them for
computing `amn`. This can be done by setting

`atom_proj_dir = ’./ext_proj’`

where the directory `ext_proj` contains the orbitals for all the atomic
species used in the calculation. For example, for a silicon calculation,
the directory `ext_proj` should contain a file named `Si.dat`. The
format of the file is:

1.  The first line contains two integers: the number of radial grid
    points ($n_g$) and the number of projectors ($n_p$), e.g.,

        1141 2

    which means the radial grid has $n_g = 1141$ points, and there are
    $n_p = 3$ projectors.

2.  The second line contains $n_p$ integers specifying the angular
    momentums of all the projectors, e.g.,

        0 1

    standing for the two projectors having $s$ and $p$ characters,
    respectively.

3.  The rest of the file contains $n_g$ rows of the radial wavefunctions
    of the projectors. There are $2+n_p$ columns: the first column is
    the $x$-grid, the second column is the $r$-grid in Bohr unit, and
    they are related by $r = \exp(x)$. The rest are $n_p$ columns of the
    radial wavefunctions of the projectors,

        -9.639057329615259 0.000065134426111 3.32211124436945e-05 1.86840239681223e-09
        -9.626557329615258 0.000065953716334 3.363898259696903e-05 1.915701228607072e-09
        -9.614057329615258 0.000066783311958 3.406210890972733e-05 1.964197436025957e-09
        ...

    Inside `pw2wannier90.x`, the radial wavefunction will be read and
    multiplied by spherical harmonics to form the actual projectors.

    For a practical example of extracting pseudo-atomic orbitals from
    UPF file and writing to a `pw2wannier90`-recognizable `.dat` file,
    see the script `utility/write_pdwf_projectors.py`.

    For an actual example of a `Si.dat` file for silicon, see the file
    `examples/example35/ext_proj/Si.dat`.
