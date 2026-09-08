---
author:
- Valerio Vitale
bibliography:
- biblio.bib
nocite:
- "[@PhysRevLett92]"
- "[@PhysRevB74]"
- "[@PhysRevLett92]"
- "[@PhysRevB74]"
- "[@PhysRevLett92]"
- "[@PhysRevB74]"
- "[@PhysRevLett92]"
- "[@PhysRevB85]"
- "[@PhysRevB74]"
- "[@PhysRevB85]"
- "[@PhysRevB85]"
- "[@PhysRevB85]"
- "[@Sakuma]"
- "[@Sakuma]"
- "[@qiao-prb2018]"
- "[@qiao-prb2018]"
- "[@guo-prl2008]"
- "[@qiao-prb2018]"
- "[@guo-prl2008]"
- "[@qiao-prb2018]"
- "[@guo-prl2008]"
- "[@qiao-prb2018]"
- "[@qiao-prb2018]"
title: "[wannier90]{.smallcaps} v3.1.0: Solution booklet"
---

# Preliminaries {#sec:preliminaries .unnumbered}

Welcome to [wannier90]{.smallcaps}! This is the solution booklet for the
examples in the [wannier90]{.smallcaps} v3.1.0 tutorial
<http://www.wannier.org/doc/tutorial.pdf>. Info on the installation
process and the theory of Maximally Localized Wannier Functions (MLWFs)
is not reported here as they can be found elsewhere[^1]. The solutions
in this booklet are for the v3.1.0 only! The following (open-source)
programs are required to reproduce the plots and figures in this
booklet:

-   `gnuplot` is used to plot bandstructures. It is available for many
    operating systems and is often installed by default on Unix/Linux
    distributions. In particular, we used gnuplot 4.6 patchlevel 6.\
    <http://www.gnuplot.info>

-   `Grace` is another plotting tool to visualise bandstructures.\
    <http://plasma-gate.weizmann.ac.il/Grace/>

-   `Vesta` is the default 3D visualisation program[@vesta] adopted in
    this booklet. It is used to visualise crystal structures, volumetric
    data (such as WFs and denisities). Download is available for several
    OS here: <http://jp-minerals.org/vesta/en/>.

-   `XCrySDen` is also used to visualise crystal structures and Fermi
    surfaces in particular. It is available for Unix/Linux, Windows
    (using cygwin), and OSX. To correctly display files from wannier90,
    version 1.4 or later must be used.\
    <http://www.xcrysden.org>

-   `VMD` may also be used to visualise crystal structures and
    3D-fields. It can also read a great variety of input formats and it
    comes with handy postprocessing tools.
    <http://www.ks.uiuc.edu/Research/vmd>

**Disclaimer:** All the band structure interpolations have been carried
out with `ws_distance = .false.`, which is the default value for the
version 2.1. However, in the new [wannier90]{.smallcaps} release,
corresponding to version 3.0, the default value of `ws_distance` has
been changed to `.true.`, as, to the best of our knowledge, the
Wigner-Seitz interpolation scheme never lowers the quality of the
interpolation and it is often superior to the default scheme.

# About this booklet {#about-this-booklet .unnumbered}

This solution manual consists of $24$ sections, each containing the
solutions, in the form of plots, tabs, and texts, to the corresponding
example in the [wannier90]{.smallcaps} v3.1.0 tutorial! For each
example, only the outline and key questions from the tutorial are
reported here. All of the [wannier90]{.smallcaps} input files have been
provided. From example 5 onwards, input files for the pwscf interface
([http://www.
quantum-espresso.org](http://www.
quantum-espresso.org){.uri}) to [wannier90]{.smallcaps} have also been
provided. You will need a recent working version of the
[quantum espresso]{.smallcaps} package (`v6.2` and above), to run these
examples. In particular, you will need `pw.x` and `pw2wannier90.x`, as
explained in the [wannier90]{.smallcaps} v3.1.0 tutorial. Please visit
<http://www.quantum-espresso.org> to download the package and follow the
instruction on the website for installation. Further details on how to
run the calculations for each example may be found in the corresponding
section of the [wannier90]{.smallcaps} v3.1.0 tutorial. There are
interfaces to a number of other electronic structure codes including:
[ABINIT]{.smallcaps} ([http://www.
abinit.org](http://www.
abinit.org){.uri}), [fleur]{.smallcaps} (<http://www.flapw.de>),
[OpenMX]{.smallcaps} (<http://www.openmx-square.org/>),
[GPAW]{.smallcaps} (<https://wiki.fysik.dtu.dk/gpaw/>),
[VASP]{.smallcaps} (<http://www.vasp.at>), and [Wien2K]{.smallcaps}
(<http://www.wien2k.at>).

All the tests were performed on an x86_64 octa-core Intel(R) Xeon(R) CPU
E5620 \@2.40GHz. Two packages `intel-suite/2015.3.187` and
`mkl/2015.3.187` were used for the compilation of
[wannier90]{.smallcaps}. We expect some of the numerical results to
depend on the architecture, the compiler distribution, e.g. `gcc` vs
`gfortran`, the version of the compiler and the libraries. However,
general trends should not be affected by these parameters and, in
principle, it should be safe to ignore the differences between different
set-ups. If you find that your results are significantly different from
the one given in this booklet, please report it to the
[wannier90]{.smallcaps} developers team by opening an issue on the
GitHub repository <https://github.com/wannier-developers/wannier90> or
by writing an email to the forum at `wannier@quantum-espresso.org` (we
strongly reccommend the first option). Moreover, if you know how to
solve the issue and you have a fix for it, you can open a pull-request
on the GitHub repo.

# Contact us {#contact-us .unnumbered}

[]{#sec:contacts label="sec:contacts"}

If you have any suggestions on how this solution manual may be improved
and for any other issue, open an issue on the official repository of the
[wannier90]{.smallcaps} code on GitHub
<https://github.com/wannier-developers/wannier90> or send an email on
the forum at `wannier@quantum-espresso.org` (we strongly reccommend the
former). For the forum note that you will need to be registered. Emails
from non-registered users will be deleted automatically. You can
register by following the links at <http://www.wannier.org/forum.html>.

# Gallium Arsenide --- MLWFs for the valence bands {#sec1:GaAs}

-   Outline: *Obtain and plot MLWFs for the four valence bands of GaAs.*

![Unit cell of GaAs crystal plotted with the [XCrySDen]{.smallcaps}
program.](figure/example01/GaAs.png){#fig1 width="0.25\\columnwidth"}

1.  *Inspect the output file `gaas.wout`.*

    Tab. [1](#tab1.1){reference-type="ref" reference="tab1.1"} shows the
    converged values (after 20 iterations) for a $2\times2\times2$
    $\mathbf{k}$ -point mesh of the total spread functional $\Omega$ and
    its three components, i.e. the gauge-invariant component
    $\Omega\ifmmode _{\mbox{\scriptsize{I}}} \else$ \_I $~\fi$, the
    off-diagonal component of the gauge-dependent part
    $\Omega\ifmmode _{\mbox{\scriptsize{OD}}} \else$ \_OD $~\fi$, and
    the diagonal component of the gauge-dependent part
    $\Omega\ifmmode _{\mbox{\scriptsize{D}}} \else$ \_D $~\fi$,
    respectively. These can be found at the end of the `gaas.wout` file,
    from the line starting with `Final State`. You will find the MLWF
    centres and their spreads together with the information on the
    spread functional components as reported below, and summarized in
    tab.[1](#tab1.1){reference-type="ref" reference="tab1.1"}.

    ::: tcolorbox
         Final State
          WF centre and spread    1  ( -0.866253,  1.973841,  1.973841 )     1.11672024
          WF centre and spread    2  ( -0.866253,  0.866253,  0.866253 )     1.11672024
          WF centre and spread    3  ( -1.973841,  1.973841,  0.866253 )     1.11672024
          WF centre and spread    4  ( -1.973841,  0.866253,  1.973841 )     1.11672024
          Sum of centres and spreads ( -5.680188,  5.680188,  5.680188 )     4.46688098

                 Spreads (Ang^2)       Omega I      =     3.956862958
                ================       Omega D      =     0.008030049
                                       Omega OD     =     0.501987969
            Final Spread (Ang^2)       Omega Total  =     4.466880976
          -----------------------------------------------------------------------------
    :::

    The geometric centre lies along the Ga-As bond, slightly closer to
    As than Ga. To see this, we introduce a measure $\beta$ defined as
    $$\beta \triangleq \frac{d\ifmmode _{\mbox{\scriptsize{Ga-MLWF}}} \else $ _{\mbox{\scriptsize{Ga-MLWF}}} $~\fi}{d\ifmmode _{\mbox{\scriptsize{Ga-As}}} \else $ _{\mbox{\scriptsize{Ga-As}}} $~\fi},$$

    where $d\ifmmode _{\mbox{\scriptsize{Ga-MLWF}}} \else$ \_Ga-MLWF
    $~\fi$ is the distance of the Ga atom placed in the origin and the
    MLWF centre (along the Ga-As bond), and
    $d\ifmmode _{\mbox{\scriptsize{Ga-As}}} \else$ \_Ga-As $~\fi$ is the
    Ga-As bond length, cf. Ref. [@marzarivanderbilt1997]. A value of 0.5
    corresponds to the MLWF centre being equidistant from the Ga atom
    and As atom. In our case, we find:
    $$\beta = \frac{d\ifmmode _{\mbox{\scriptsize{Ga-MLWF(2)}}} \else $ _{\mbox{\scriptsize{Ga-MLWF(2)}}} $~\fi}{d\ifmmode _{\mbox{\scriptsize{Ga-As}}} \else $ _{\mbox{\scriptsize{Ga-As}}} $~\fi} = \frac{0.866253\sqrt{3}}{1.4200\sqrt{3}} \approx 0.61 .$$
    Maximum RAM allocated for the wannierisation was 0.06Mb.

    ::: {#tab1.1}
      MP mesh             $\Omega$   $\Omega\ifmmode _{\mbox{\scriptsize{I}}} \else$ \_I $~\fi$   $\Omega\ifmmode _{\mbox{\scriptsize{OD}}} \else$ \_OD $~\fi$   $\Omega\ifmmode _{\mbox{\scriptsize{D}}} \else$ \_D $~\fi$   $\beta$
      ------------------- ---------- ------------------------------------------------------------ -------------------------------------------------------------- ------------------------------------------------------------ ---------
      $2\times2\times2$   4.467      3.957                                                        0.502                                                          0.008                                                        0.610

      : Converged values of the components of spread functional and
      their sum, given in Å$^2$. Here $\beta$ is the distance of the Ga
      atom placed in the origin and the MLWF centre (along the Ga-As
      bond) as a fraction of the Ga-As bond length $2.4595\si{~\AA}$.
    :::

    []{#tab1.1 label="tab1.1"}

2.  *Plot the MLWFs.*

    In Fig. [2](#fig1.2){reference-type="ref" reference="fig1.2"} are
    shown the four valence MLWFs as plotted by [XCrySDen]{.smallcaps},
    where we used the following parameters in the `Tools` $\mapsto$
    `Data Grid` section:

    ::: tcolorbox
    ` Degree of triCubic Spline = 3; Isovalue = 0.95; Render +/- isovalue = yes `
    :::

    <figure id="fig1.2">

    <figcaption>Four valence MLWFs for the Ga-As system plotted using the
    <span class="smallcaps">XCrySDen</span> visualisation
    program.</figcaption>
    </figure>

3.  *Plot the 3$rd$ MLWFs in a supercell of size 3. Choose a low value
    for the isosurface (say 0.5). Can you explain what you see?*

    With [XCrySDen]{.smallcaps} we can also plot the
    $3\ifmmode ^{\mbox{\scriptsize{rd}}} \else$ \^rd $~\fi$ MLWF to
    check its periodicity. The period in each direction is given by the
    spacing used to sample the first irreducible Brillouin zone, i.e.
    the $\mathbf{k}$ -point mesh. We used a $2\times2\times2$
    $\mathbf{k}$ -point mesh, hence we expect to find a periodic image
    of our MLWF in a supercell which is 2 times larger than the unit
    cell along each direction. This is shown in
    Fig. [3](#fig1.3){reference-type="ref" reference="fig1.3"}, where
    the $3\ifmmode ^{\mbox{\scriptsize{rd}}} \else$ \^rd $~\fi$ has been
    plotted using both the [XCrySDen]{.smallcaps} program and the
    [vesta]{.smallcaps} program.

<figure id="fig1.3">

<figcaption><span class="math inline">$3\ifmmode
^{\mbox{\scriptsize{rd}}} \else$</span> ^<span><span>rd</span></span>
<span class="math inline">$~\fi$</span> MLWF with a supercell value of 3
and for an isovalue of 0.5 using (a) <span
class="smallcaps">XCrySDen</span> and (b) <span
class="smallcaps">vesta</span>.</figcaption>
</figure>

# Lead --- Wannier-interpolated Fermi surface {#sec2:lead}

-   Outline: *Obtain MLWFs for the four lowest states in lead. Use
    Wannier interpolation to plot the Fermi surface.*

![Unit cell of lead crystal plotted with the [XCrySDen]{.smallcaps}
program.](figure/example02/lead.png){#fig2.0 width="0.25\\columnwidth"}

1.  *Inspect the output file `lead.wout`.*

    A summary of the wannierisation is given in
    tab.[2](#tab2.1){reference-type="ref" reference="tab2.1"}. At the
    end of the `.wout` file you should find the info on the final state
    of the minimization as

    ::: tcolorbox
         Final State
          WF centre and spread    1  (  0.397070,  0.397070,  0.397070 )     1.93781315
          WF centre and spread    2  (  0.397070, -0.397070, -0.397070 )     1.93781315
          WF centre and spread    3  ( -0.397070,  0.397070, -0.397070 )     1.93781315
          WF centre and spread    4  ( -0.397070, -0.397070,  0.397070 )     1.93781315
          Sum of centres and spreads (  0.000000, -0.000000, -0.000000 )     7.75125261

                 Spreads (Ang^2)       Omega I      =     6.039099038
                ================       Omega D      =     0.007065754
                                       Omega OD     =     1.705087819
            Final Spread (Ang^2)       Omega Total  =     7.751252611
         ------------------------------------------------------------------------------
    :::

    ::: {#tab2.1}
      MP mesh             $\Omega$   $\Omega\ifmmode _{\mbox{\scriptsize{I}}} \else$ \_I $~\fi$   $\Omega\ifmmode _{\mbox{\scriptsize{OD}}} \else$ \_OD $~\fi$   $\Omega\ifmmode _{\mbox{\scriptsize{D}}} \else$ \_D $~\fi$
      ------------------- ---------- ------------------------------------------------------------ -------------------------------------------------------------- ------------------------------------------------------------
      $4\times4\times4$   7.751      6.039                                                        0.007                                                          1.705

      : Converged values of the components of spread functional and
      their sum, given in Å$^2$.
    :::

    []{#tab2.1 label="tab2.1"}

2.  *Use Wannier interpolation to generate the Fermi surface of lead.*

    As can be seen from the bandstructure plot in the
    [wannier90]{.smallcaps} tutorial, that we report here, cf.
    Fig. [5](#fig2.1){reference-type="ref" reference="fig2.1"}, the four
    lower valence bands are separated in energy from the higher
    conduction states (there is however an indirect band gap). The Fermi
    level lies somewhere in the middle of the manifold making
    crystalline lead a metal. As a results, the states belonging to
    these manifold will have partial occupancy. We can see that 2 bands
    are entirely below and above the Fermi level, respectively. Hence,
    no Fermi surface can be plotted from these bands. On the other hand,
    the two central bands do cross the Fermi energy level and the
    corresponding Fermi surfaces are shown in
    Fig. [6](#fig2.2){reference-type="ref" reference="fig2.2"}.

    ![Bandstructure of lead showing the position of the Fermi level.
    Only the lowest four bands are included in the
    calculation.](figure/example02/lead.pdf){#fig2.1}

    <figure id="fig2.2">

    <figcaption>Fermi surfaces for band 2 and band 3 in lead. The value of
    the Fermi energy is 5.2676eV, and it was obtained from the first
    principle calculation, with a <span class="math inline">4 × 4 × 4</span>
    <span class="math inline"><strong>k</strong></span> -point mesh. To
    calculate the band energies and to plot Fermi surfaces, Wannier
    interpolation was employed on a dense mesh in the Brillouin zone
    consisting of <span class="math inline">50<sup>3</sup></span>
    points.</figcaption>
    </figure>

# Silicon --- Disentangled MLWFs {#sec3:silicon}

-   Outline: *Obtain disentangled MLWFs for the valence and low-lying
    conduction states of Si. Plot the interpolated bandstructure.*

![Unit cell of Silicon crystal plotted with the [XCrySDen]{.smallcaps}
program.](figure/example03/silicon.png){#fig3.0
width="0.25\\columnwidth"}

1.  *Inspect the output file `silicon.wout``.`*

    Starting from 4 $sp3$ orbitals on each Silicon atom we obtain two
    sets of 4 WFs, all with the same spread, that show the $sp3$
    character one would expect from symmetry considerations. A summary
    of the wannierisation is given in
    tab.[3](#tab3.1){reference-type="ref" reference="tab3.1"}. At the
    end of the `.wout` file you should find the info on the final state
    of the minimization, here we show an extract of the output file

    ::: tcolorbox
        	 Final State
          WF centre and spread    1  ( -0.460754, -0.460711, -0.460767 )     1.81241746
          WF centre and spread    2  ( -0.460743,  0.460722,  0.460718 )     1.81246400
          WF centre and spread    3  (  0.460703, -0.460761,  0.460685 )     1.81248774
          WF centre and spread    4  (  0.460704,  0.460724, -0.460764 )     1.81244947
          WF centre and spread    5  (  1.810128,  1.810112,  1.810113 )     1.81247628
          WF centre and spread    6  (  1.810097,  0.888662,  0.888617 )     1.81242110
          WF centre and spread    7  (  0.888640,  1.810140,  0.888660 )     1.81240614
          WF centre and spread    8  (  0.888643,  0.888652,  1.810090 )     1.81245230
          Sum of centres and spreads (  5.397417,  5.397539,  5.397353 )    14.49957450

                 Spreads (Ang^2)       Omega I      =    11.849193709
                ================       Omega D      =     0.105470243
                                       Omega OD     =     2.544910550
            Final Spread (Ang^2)       Omega Total  =    14.499574503
         ------------------------------------------------------------------------------
        	
    :::

    ::: {#tab3.1}
      MP mesh             $\Omega$   $\Omega\ifmmode _{\mbox{\scriptsize{I}}} \else$ \_I $~\fi$   $\Omega\ifmmode _{\mbox{\scriptsize{OD}}} \else$ \_OD $~\fi$   $\Omega\ifmmode _{\mbox{\scriptsize{D}}} \else$ \_D $~\fi$
      ------------------- ---------- ------------------------------------------------------------ -------------------------------------------------------------- ------------------------------------------------------------
      $4\times4\times4$   14.5       11.849                                                       2.545                                                          0.105

      : Converged values of the components of spread functional and
      their sum, given in Å$^2$.
    :::

2.  *Plot the energy bands.*

    As can be seen from DFT bandstructure plot in the
    [wannier90]{.smallcaps} tutorial, that we report here, cf.
    Fig. [8](#fig3.1){reference-type="ref" reference="fig3.1"}, the four
    lower valence bands are separated in energy from the higher
    conduction states (there is however an indirect band gap). The Fermi
    level lies inside the gap, making crystalline Silicon a
    semiconductor.

    The path in $\ifmmode  \mathbf{k}  \else$ $~\fi$-space given in the
    tutorial (L-$\Gamma$-X-K-$\Gamma$) and shown in
    Fig. [9](#fig3.2){reference-type="ref" reference="fig3.2"}-(a)-top
    is the following

    ` `

    > begin kpoint_path
    >
    > L 0.50000 0.50000 0.5000 G 0.00000 0.00000 0.0000
    >
    > G 0.00000 0.00000 0.0000 X 0.50000 0.00000 0.5000
    >
    > X 0.50000 -0.50000 0.0000 K 0.37500 -0.37500 0.0000
    >
    > K 0.37500 -0.37500 0.0000 G 0.00000 0.00000 0.0000
    >
    > end kpoint_path

    which gives the bandstructure shown in
    Fig. [9](#fig3.2){reference-type="ref"
    reference="fig3.2"}-(a)-bottom.

3.  *Try plotting along different paths.*

    Another path usually used for Silicon is W-$\Gamma$-X-W-L-$\Gamma$
    shown in Fig. [8](#fig3.1){reference-type="ref"
    reference="fig3.1"}-(b)-top and the corresponding bands are shown in
    Fig. [9](#fig3.2){reference-type="ref"
    reference="fig3.2"}-(b)-bottom. To obtain this path you need to
    replace the previous `kpoint_path` block with the following block

    ` `

    > begin kpoint_path
    >
    > W 0.25000 0.75000 0.50000 G 0.00000 0.00000 0.00000
    >
    > G 0.00000 0.00000 0.00000 X 0.50000 0.50000 0.00000
    >
    > X 0.50000 0.50000 0.00000 W -0.25000 0.25000 -0.25000
    >
    > W -0.25000 0.25000 -0.25000 L 0.00000 0.50000 0.00000
    >
    > L 0.00000 0.50000 0.00000 G 0.00000 0.00000 0.00000
    >
    > end kpoint_path

    ![Bandstructure of Silicon showing the position of the Fermi level
    and of the inner and outer windows. Both the 4 valence bands and the
    4 low-lying conduction bands are included in the
    calculation.](figure/example03/si.eps){#fig3.1
    width="0.7\\columnwidth"}

    <figure id="fig3.2">
    <p><br />
    </p>
    <figcaption>Bandstructure of Silicon showing the position of the Fermi
    level and of both the inner and outer windows. The 4 valence bands
    together with the 4 low-lying conduction bands are included in the
    calculation. Panel a) Interpolation with Wannier90 on the L-<span
    class="math inline"><em>Γ</em></span>-X-K-<span
    class="math inline"><em>Γ</em></span> path in <span
    class="math inline">$\ifmmode  \mathbf{k}  \else$</span> <span
    class="math inline">$~\fi$</span> space (red dots) and DFT reference
    bandstructure (solid black). Panel b) Interpolation with Wannier90 on
    the W-<span class="math inline"><em>Γ</em></span>-X-W-L-<span
    class="math inline"><em>Γ</em></span> path in <span
    class="math inline">$\ifmmode  \mathbf{k}  \else$</span> <span
    class="math inline">$~\fi$</span> space (red dots) and DFT reference
    bandstructure (solid black).</figcaption>
    </figure>

# Copper --- Fermi surface, orbital character of energy bands {#sec4:copper}

-   Outline: *Obtain MLWFs to describe the states around the Fermi-level
    in copper*

![Unit cell of Copper crystal plotted with the [XCrySDen]{.smallcaps}
program.](figure/example04/copper_crystal.png){#fig4.0
width="0.25\\columnwidth"}

1.  *Run [wannier90]{.smallcaps} to minimise the MLWFs spread. Inspect
    the output file `copper.wout.`*

    Starting from 5 $d$ orbitals centred on the Copper atom and 2 $s$
    orbitals in the interstitial regions of the FCC, we obtain the
    following spreads and centres after 200 iterations (extract from the
    `copper.wout`, a summary of the wannierisation is given in
    tab.[4](#tab4.1){reference-type="ref" reference="tab4.1"}.):

    ::: tcolorbox
         Final State
          WF centre and spread    1  ( -0.000000,  0.000000, -0.000000 )     0.40838932
          WF centre and spread    2  ( -0.000000, -0.000000, -0.000000 )     0.30784969
          WF centre and spread    3  ( -0.000000, -0.000000,  0.000000 )     0.30784979
          WF centre and spread    4  ( -0.000000, -0.000000,  0.000000 )     0.40838973
          WF centre and spread    5  (  0.000000, -0.000000, -0.000000 )     0.30784886
          WF centre and spread    6  ( -0.902512,  0.902512,  0.902512 )     1.14385632
          WF centre and spread    7  (  0.902512, -0.902512, -0.902512 )     1.14385635
          Sum of centres and spreads (  0.000000, -0.000000, -0.000000 )     4.02804006

                 Spreads (Ang^2)       Omega I      =     3.662691490
                ================       Omega D      =     0.001894482
                                       Omega OD     =     0.363454087
            Final Spread (Ang^2)       Omega Total  =     4.028040058
         ------------------------------------------------------------------------------
        	
    :::

    We can readily see that looking at the individual spreads we find
    two groups of MLWFs, a group of 5 $d$-like MLWFs centred on the
    Copper atom, whose spreads are $0.4084$Å$^2$and
    $0.3078$Å$^2$respectively, and a group of 2 $s$-like MLWFs centred
    on two opposite (with respect to the origin) interstitial points,
    whose spread is $1.1439$Å$^2$. The 3+2 $d$-like MLWFs are the basis
    of two representations of the $O_h$ group, with character $t2_g$ and
    $e_g$ respectively.

    ::: {#tab4.1}
      MP mesh             $\Omega$   $\Omega\ifmmode _{\mbox{\scriptsize{I}}} \else$ \_I $~\fi$   $\Omega\ifmmode _{\mbox{\scriptsize{OD}}} \else$ \_OD $~\fi$   $\Omega\ifmmode _{\mbox{\scriptsize{D}}} \else$ \_D $~\fi$
      ------------------- ---------- ------------------------------------------------------------ -------------------------------------------------------------- ------------------------------------------------------------
      $4\times4\times4$   4.028      3.66                                                         0.363                                                          0.002

      : Converged values of the components of spread functional and
      their sum, given in Å$^2$.
    :::

    []{#tab4.1 label="tab4.1"}

2.  *Plot the Fermi surface, it should look familiar! The Fermi energy
    is at 12.2103 eV.*

    As explained in example 2 of the tutorial, we need to add the
    following lines to the input file (`copper.win`):

    ` `

    > restart = plot\
    > fermi_energy = 12.2103\
    > fermi_surface_plot = true\
    > fermi_surface_num_points = 50

    and re-run the [wannier90]{.smallcaps} executable. The result will
    be a file named `copper.bxsf`, which contains volumetric data in a
    format suitable for `xcrysden`. There is only one band that crosses
    the Fermi level (12.2103eV), i.e. band 6, as shown in
    Fig. [11](#fig4.1){reference-type="ref" reference="fig4.1"}-(a). The
    corresponding Fermi surface is shown in
    Fig. [11](#fig4.1){reference-type="ref" reference="fig4.1"}-(b)

    <figure id="fig4.1">

    <figcaption>Fermi surfaces for band 6 in copper. The value of the Fermi
    energy is 12.2103eV, and it was obtained from DFT calculation, with a
    <span class="math inline">4 × 4 × 4</span> <span
    class="math inline"><strong>k</strong></span> -point mesh. To calculate
    the band energies and to plot Fermi surfaces, Wannier interpolation was
    employed on a dense mesh in the Brillouin zone consisting of <span
    class="math inline">50<sup>3</sup></span> points.</figcaption>
    </figure>

3.  *Plot the interpolated bandstructure.*

    Interpolated bandstructure, with path in k-space given in the
    tutorial, is shown in Fig. [12](#fig4.2){reference-type="ref"
    reference="fig4.2"}.

    ![Interpolated bandstructure of Copper (solid red) showing the
    position of the Fermi level (dashed red) and both inner and outer
    windows (dotted and dashed-dotted respectively). The reference DFT
    bandstructure (solid black) was obtained with Quantum ESPRESSO, see
    procedure in Example [6](#sec6:copper){reference-type="ref"
    reference="sec6:copper"}.](figure/example04/copper_bs_qe_w90.pdf){#fig4.2
    width="0.8\\columnwidth"}

4.  *Plot the contribution of the interstitial WF to the bandstructure.*

    The contribution of the 2 $s$-like MLWFs to the band structure is
    shown in Fig. [13](#fig4.3){reference-type="ref" reference="fig4.3"}

    ![Bandstructure of Copper showing the contribution from the 2
    $s$-like MLWFs in
    red.](figure/example04/copper_bs_projection.pdf){#fig4.3
    width="0.9\\columnwidth"}

5.  *Investigate the effect of the outer and inner energy window on the
    interpolated bands.*

    From now on, we will refer to the inner window energy level as
    $\varepsilon_{\mathrm{in}}$ and to the outer window energy level as
    $\varepsilon_{\mathrm{out}}$. The reference values are in this case
    $\varepsilon_{\mathrm{in}}^0 = 13$eV and
    $\varepsilon_{\mathrm{out}}^0 = 38$eV. Hereafter, we will use
    $\varepsilon_{\mathrm{min}}$ to refer to the minimum energy
    eigenvalue (for the chosen path). The value of
    $\varepsilon_{\mathrm{out}}^0$ must be chosen such that there are at
    least $N_{\mathrm{wannier}}$ (7 in this case) states inside the
    outer energy window for each $\mathbf{k}$ -point. This means that
    for a given path in the BZ there exists a lower bound to
    $\varepsilon_{\mathrm{out}}^0$. The actual value depends on the path
    in the BZ and on the choice of the zero for the pseudopotential. The
    result for several values of $\varepsilon_{\mathrm{in}}$ and
    $\varepsilon_{\mathrm{out}}$ are shown in
    Fig. [14](#fig4.4){reference-type="ref" reference="fig4.4"}.

    <figure id="fig4.4">
    <p><br />
    <br />
    </p>
    <figcaption>Interpolated bandstructure of Copper (solid red) with DFT
    reference (solid black) and different values of the inner window. Panel
    a) <span class="math inline">…</span> </figcaption>
    </figure>

# Diamond --- MLWFs for the valence bands {#sec5:diamond}

-   Outline: *Obtain MLWFs for the valence bands of diamond.*

![Unit cell of Diamond crystal plotted with the [XCrySDen]{.smallcaps}
program.](figure/example05/diamond.png){#fig5.0
width="0.25\\columnwidth"}

1.  *Run pwscf to obtain the ground state of diamond.*

    Convergence of the self-consistent field calculation in Quantum
    Espresso can be checked at the end of the `scf.out` file. At the
    very end of the file one should find the line confirming that the
    job has finished without crashing, e.g.

    ::: tcolorbox
        	=------------------------------------------------------------------------------=
        	   JOB DONE.
        	=------------------------------------------------------------------------------=
        	
    :::

    Depending on the output verbosity one may or may not find info about
    WALL times for the calls to the different routines. Just above this
    block, if present, one may find the info about the convergence of
    the SCF loop, such as the scf accuracy and the number of iterations
    to required to achieve it:

    ::: tcolorbox
            !    total energy              =     -22.58128615 Ry
                 Harris-Foulkes estimate   =     -22.58128615 Ry
                 estimated scf accuracy    <          1.0E-14 Ry


                 The total energy is the sum of the following terms:

                 one-electron contribution =      11.69117931 Ry
                 hartree contribution      =       1.57036314 Ry
                 xc contribution           =      -7.58421586 Ry
                 ewald contribution        =     -28.25861274 Ry

                 convergence has been achieved in   9 iterations
        	
    :::

2.  *Run pwscf to obtain the Bloch states on a uniform k-point grid.*

    Similarly for the non-scf calculation one can check that the
    calculation has been carried out without crashing by looking at the
    last three line of the `nscf.out` file. A useful information to
    check is the value of the highest eigenvalue (for insulators and
    semiconductors) or the value of the Fermi level (for metals). In the
    diamond we case, we find:

    ::: tcolorbox
        	highest occupied level (ev):    19.3978
        	
    :::

3.  *Run [wannier90]{.smallcaps} to compute the MLWFs.*

    The result of the wannierisation, after 20 iterations, may be found
    at the end of `diamond.wout` file:

    ::: tcolorbox
        	 Final State
          WF centre and spread    1  ( -0.000000,  0.000000, -0.000000 )     0.58022623
          WF centre and spread    2  ( -0.806995,  0.806995,  0.000000 )     0.58022623
          WF centre and spread    3  ( -0.000000,  0.806995,  0.806995 )     0.58022623
          WF centre and spread    4  ( -0.806995, -0.000000,  0.806995 )     0.58022623
          Sum of centres and spreads ( -1.613990,  1.613990,  1.613990 )     2.32090491

                 Spreads (Ang^2)       Omega I      =     1.954619859
                ================       Omega D      =     0.000000000
                                       Omega OD     =     0.366285054
            Final Spread (Ang^2)       Omega Total  =     2.320904912
         ------------------------------------------------------------------------------
        	
    :::

4.  *Plot the 4 MLWFs.*

    The resulting 4 $\sigma$-bonding MLWFs are shown in
    Fig. [16](#fig5.1){reference-type="ref" reference="fig5.1"}

<figure id="fig5.1">

<figcaption>4 MLWFs in diamond describing the valence bands plotted
using <span class="smallcaps">vesta</span>.</figcaption>
</figure>

# Copper --- Fermi surface {#sec6:copper}

-   Outline: *Obtain MLWFs to describe the states around the Fermi-level
    in copper.*

![Unit cell of Copper crystal plotted with the [XCrySDen]{.smallcaps}
program.](figure/example04/copper_crystal.png){#fig6.0
width="0.25\\columnwidth"}

After checking that the calculations have converged as shown in Example
[5](#sec5:diamond){reference-type="ref" reference="sec5:diamond"}, one
can proceed with other points in the example.

1.  *Use Wannier interpolation to obtain the Fermi surface of copper.*

    To obtain the value of the Fermi energy we can use the `grep`
    command (only for Linux/Unix systems) as:

    > \$ \> grep Fermi nscf.out

    The output should be:

    ` `

    > the Fermi energy is 12.9344 ev

    Alternatively, one can open the `nscf.out` file with the editor of
    choice and search for \"Fermi\" inside the file. We can then use
    this value in the `.win` file to compute the Fermi surface as done
    in Example [2](#sec2:lead){reference-type="ref"
    reference="sec2:lead"}. The interpolated Fermi surface is shown in
    Fig. [18](#fig6.1){reference-type="ref" reference="fig6.1"}-(a).

    <figure id="fig6.1">

    <figcaption>a) Fermi surface of Copper. b) Band structure of Copper
    along the <span class="math inline"><em>Γ</em></span>–X–W–L–<span
    class="math inline"><em>Γ</em></span>–K computed from a non-scf DFT
    calculation (solid black) and via Wannier interpolation (solid
    red).</figcaption>
    </figure>

2.  *Plot the interpolated bandstructure.*

    Bandstructure is shown in Fig. [18](#fig6.1){reference-type="ref"
    reference="fig6.1"}-(b). One way to obtain the DFT bandstructure on
    exactly the same path as the one in the `.win` input file, is given
    by the `bands.x` program available at
    <http://www.tcm.phy.cam.ac.uk/~jry20/bands.html>.

<figure id="fig6.2">

<figcaption>a) Bandstructure of Copper along the <span
class="math inline"><em>Γ</em></span>–X–W–L–<span
class="math inline"><em>Γ</em></span>–K from a non-scf DFT calculation
(solid black) and via Wannier interpolation using two different sets of
initial projections, namely <span class="math inline">2<em>s</em></span>
and 5<span class="math inline"><em>d</em></span> (<span
class="math inline"><em>N</em><sub><em>w</em></sub> = 7</span>) (solid
red) and 1<span class="math inline"><em>s</em></span> 3<span
class="math inline"><em>p</em></span> and 5<span
class="math inline"><em>d</em></span> (<span
class="math inline"><em>N</em><sub><em>w</em></sub> = 9</span>) (solid
blue). b) <span class="math inline"><em>p</em></span> character of bands
computed using <span><code>bands_plot_project = 2,3,4</code></span> in
the input file. A color scheme is used to measure the <span
class="math inline"><em>p</em></span> <span><em>character</em></span> of
the bands.</figcaption>
</figure>

::: tcolorbox
Here we summarize the main steps to produce the bandstructure with the
`bands.x` code:

-   Compilation:

    > eg. g95 -o bands.x bands.F90
    >
    > ifort -o bands.x bands.F90
    >
    > for NAG
    >
    > f95 -o bands.x bands.F90 -DNAG

-   Usage: First you need to generate an `copper.inp` file, with the
    following structure

            ! Input file for Copper
            !
            ! First the unit cell (in atomic units = Bohr)
            -3.411 0.000 3.411
             0.000 3.411 3.411
            -3.411 3.411 0.000

            !then the number of points along the 1st special path
            100

            ! then the special kpoints and their labels
            G 0.00  0.00  0.00    X 0.50  0.50  0.00
            X 0.50  0.50  0.00    W 0.50  0.75  0.25
            W 0.50  0.75  0.25    L 0.00  0.50  0.00
            L 0.00  0.50  0.00    G 0.00  0.00  0.00
            G 0.00  0.00  0.00    K 0.00  0.50 -0.50

    Then you need to generate the kpoint list by running the `bands.x`
    program with the `-pp` flag

    ` `

    > \$ \> ./bands.x -pp copper

    This will read data from `copper.inp` and write kpoints into
    `copper_band.kpt`.

    **WARNING**: if you already have a `copper_band.kpt` file from a
    previous Wannier90 calculation, running the above command will
    overwrite it.

    Now you need to calculate a non-scf or bands calculation with
    Quantum ESPRESSO on the k-points given in `copper_band.kpt`. To do
    so, copy the `copper.nscf` to `copper.bands` and modify it
    accordingly. Run a non-scf calculation

    ` `

    > \$ \> pw.x \< copper.bands \> copper.pwscf

    **WARNING**: the output file must terminate with `.pwscf` in order
    to be read by `bands.x`.

    Now extract the bands from the `copper.pwscf` file

    ` `

    > \$ \> bands.x copper

    The bands are written into `copper_band.dat`. WARNING: if you
    already have a `coppper_band.dat` file and a `copper_band.gnu` file
    from a previous Wannier90 calculation, running the above command
    will overwrite them.

    Plot with gnuplot

    ` `

    > \$ \> gnuplot --persist copper_band.gnu
:::

1.  *Compare the Wannier interpolated bandstructure with the full pwscf
    bandstructure. Obtain MLWFs using a denser k-point grid.*

2.  *Investigate the effects of the outer and inner energy windows on
    the interpolated bands.*

    The effect of different energy windows has already been discussed in
    the Example [4](#sec4:copper){reference-type="ref"
    reference="sec4:copper"}, so it won't be repeated here.

3.  *Instead of extracting a subspace of seven states, we could extract
    a nine dimensional space (i.e., with $s$, $p$ and $d$ character).
    Examine this case and compare the interpolated bandstructures.*

    Using $s,p$, and $d$ orbitals as initial guesses for the MLWFs of
    Copper, yields the bandstructure (solid blue) shown in
    Fig. [19](#fig6.2){reference-type="ref" reference="fig6.2"} (with
    reference values for the inner and outer windows). The bandstructure
    obtained starting from $2s$ and $5d$ orbitals is shown in red,
    whereas the DFT reference bandstructure, computed with the procedure
    described above, is in black. It is clear from
    Fig. [19](#fig6.2){reference-type="ref" reference="fig6.2"}(a), and
    Fig. [19](#fig6.2){reference-type="ref" reference="fig6.2"}(b) that
    the bands of interest have very little $p$ character, particularly
    the 5 flat bands, which are well very described by $d$ states.

# Silane (SiH 4 ) --- Molecular MLWFs using $\Gamma$-point sampling {#sec7:silane}

-   Outline: *Obtain MLWFs for the valence bands of silane.*

![Silane molecule in a periodic cell plotted with the
[XCrySDen]{.smallcaps} program.](figure/example07/silane.png){#fig7.0
width="0.25\\columnwidth"}

1.  Convergence of the self-consistent field calculation in Quantum
    Espresso can be checked at the end of the `scf.out` file. At the
    very end of the file one should find the line confirming that the
    job has finished without crashing, e.g.

    ::: tcolorbox
        	=------------------------------------------------------------------------------=
        	   JOB DONE.
        	=------------------------------------------------------------------------------=
        	
    :::

    Just above the block reporting the info about WALL times, if
    present, one may find the info about the convergence of the SCF
    loop, such as the scf accuracy and the number of iterations to
    required to achieve it:

    ::: tcolorbox
        !    total energy              =     -12.25602944 Ry
             Harris-Foulkes estimate   =     -12.25602944 Ry
             estimated scf accuracy    <          7.0E-11 Ry

                 The total energy is the sum of the following terms:

                 one-electron contribution =      11.69117931 Ry
                 hartree contribution      =       1.57036314 Ry
                 xc contribution           =      -7.58421586 Ry
                 ewald contribution        =     -28.25861274 Ry

                 convergence has been achieved in   9 iterations
        	
    :::

2.  Similarly for the non-scf calculation one can check that the
    calculation has been carried out without crashing by looking at the
    last three line of the `nscf.out` file. A useful information to
    check is the value of the highest eigenvalue (for insulators and
    semiconductors) or the value of the Fermi level (for metals). In the
    diamond we case, we find:

    ::: tcolorbox
             highest occupied level (ev):    -6.5316
        	
    :::

3.  The result of the wannierisation, after 20 iterations, may be found
    at the end of `silane.wout` file:

    ::: tcolorbox
         Final State
          WF centre and spread    1  (  0.762490,  0.762490,  0.762490 )     1.01124580
          WF centre and spread    2  (  0.762491, -0.762492, -0.762491 )     1.01124445
          WF centre and spread    3  ( -0.762491,  0.762490, -0.762491 )     1.01124473
          WF centre and spread    4  ( -0.762491, -0.762491,  0.762491 )     1.01124420
          Sum of centres and spreads ( -0.000002, -0.000003, -0.000001 )     4.04497917

                 Spreads (Ang^2)       Omega I      =     3.920639090
                ================       Omega D      =     0.000000000
                                       Omega OD     =     0.124340085
            Final Spread (Ang^2)       Omega Total  =     4.044979175
         ------------------------------------------------------------------------------
        	
    :::

# Iron --- Spin-polarized WFs, DOS, projected WFs versus MLWFs {#sec8:Iron}

-   Outline : *Generate both maximally-localized and projected Wannier
    functions for ferromagnetic bcc Fe. Calculate the total and
    orbital-projected density of states by Wannier interpolation.*

![Unit cell of Iron crystal plotted with the [XCrySDen]{.smallcaps}
program.](figure/example08/iron.png){#fig8.0 width="0.25\\columnwidth"}

1.  *Converged values for the total spread functional and its components
    for both spin channels are shown in
    Tab. [5](#tab8.1){reference-type="ref" reference="tab8.1"}.* The
    final state for spin-up MLWFs is

    ::: tcolorbox
        	 Final State
          WF centre and spread    1  (  0.709852,  0.000108,  0.000131 )     1.08935224
          WF centre and spread    2  (  0.000131,  0.000053, -0.709852 )     1.08935218
          WF centre and spread    3  ( -0.709852, -0.000108, -0.000131 )     1.08935221
          WF centre and spread    4  (  0.000108, -0.709852, -0.000053 )     1.08935218
          WF centre and spread    5  ( -0.000131, -0.000053,  0.709852 )     1.08935226
          WF centre and spread    6  (  0.000000,  0.000000,  0.000000 )     0.43234428
          WF centre and spread    7  ( -0.000000,  0.000000,  0.000000 )     0.43234429
          WF centre and spread    8  ( -0.000108,  0.709852,  0.000053 )     1.08935225
          WF centre and spread    9  (  0.000000,  0.000000, -0.000000 )     0.43234428
          Sum of centres and spreads (  0.000000, -0.000000, -0.000000 )     7.83314616

                 Spreads (Ang^2)       Omega I      =     5.948424630
                ================       Omega D      =     0.017027691
                                       Omega OD     =     1.867693841
            Final Spread (Ang^2)       Omega Total  =     7.833146162
         ------------------------------------------------------------------------------
        	
    :::

    and for spin-down MLWFs is

    ::: tcolorbox
        [sharp corners,boxrule=0.5pt]
        	 Final State
          WF centre and spread    1  ( -0.685467, -0.000123,  0.000259 )     1.10268580
          WF centre and spread    2  ( -0.000259, -0.000207, -0.685467 )     1.10268617
          WF centre and spread    3  (  0.685468,  0.000123, -0.000259 )     1.10268605
          WF centre and spread    4  ( -0.000123,  0.685467, -0.000207 )     1.10268595
          WF centre and spread    5  (  0.000259,  0.000207,  0.685467 )     1.10268552
          WF centre and spread    6  (  0.000000,  0.000000, -0.000000 )     0.41116646
          WF centre and spread    7  ( -0.000000,  0.000000, -0.000000 )     0.41116648
          WF centre and spread    8  (  0.000123, -0.685467,  0.000207 )     1.10268572
          WF centre and spread    9  (  0.000000,  0.000000,  0.000000 )     0.41116644
          Sum of centres and spreads (  0.000000, -0.000000,  0.000000 )     7.84961460

                 Spreads (Ang^2)       Omega I      =     5.946718376
                ================       Omega D      =     0.014524283
                                       Omega OD     =     1.888371944
            Final Spread (Ang^2)       Omega Total  =     7.849614603
         ------------------------------------------------------------------------------
        	
    :::

    As it is clear from the output file snippets above, the $s,p$ and
    $d$ orbitals hybridize to give rise to two groups of functions for
    both spin channels. A first group made of 6 MLWFs coming from the
    hybridisation of $sp^3$ and $d_{e_g}$ MLWFs, with a total spread of
    1.089(1.103)Å$^2$ for spin-up(down). A second group made of 3 MLWFs
    with a $d_{t_{2g}}$ character, with a total spread of
    0.432(0.4112)Å$^2$ for spin-up(down). Two sample MLWFs, one for each
    group, are shown in Fig. [22](#fig8.3){reference-type="ref"
    reference="fig8.3"}.

<figure id="fig8.3">

<figcaption>2 representative MLWFs from the wannierisation of 9 spin-up
bands of iron. a) A representative of the hybrid (<span
class="math inline"><em>s</em><em>p</em><sup>3</sup></span> and <span
class="math inline"><em>d</em><sub><em>e</em><sub><em>g</em></sub></sub></span>)
group of MLWFs. b) A representative of the <span
class="math inline"><em>d</em><sub><em>t</em><sub>2<em>g</em></sub></sub></span>
group of MLWFs.</figcaption>
</figure>

::: {#tab8.1}
  spin   $\Omega$   $\Omega\ifmmode _{\mbox{\scriptsize{I}}} \else$ \_I $~\fi$   $\Omega\ifmmode _{\mbox{\scriptsize{OD}}} \else$ \_OD $~\fi$   $\Omega\ifmmode _{\mbox{\scriptsize{D}}} \else$ \_D $~\fi$   $N_{\mathrm{iter}}$
  ------ ---------- ------------------------------------------------------------ -------------------------------------------------------------- ------------------------------------------------------------ ---------------------
  up     7.8331     5.9484                                                       1.8677                                                         0.0170                                                       400
  down   7.8496     5.9467                                                       1.8884                                                         0.0145                                                       400

  : Converged values of the components of spread functional and their
  sums for both spin chanels for ferromagnetic bcc Fe, given in Å$^2$.
:::

[]{#tab8.1 label="tab8.1"}

## Density of states {#density-of-states .unnumbered}

-   *run `postw90` and plot the DOS with `gnuplot`*

![Interpolated DOS of bcc iron on a $25\times25\times25$
$\mathbf{k}$-mesh. Up-spin channel (solid red). Down-spin channel (solid
blue).](figure/example08/DOS_iron_bcc.pdf){#fig8.1
width="0.7\\columnwidth"}

-   *Check the convergence by repeating the DOS calculations with more
    k-points.*

    Plots of the DOS calculated with different k-point mesh densities
    for the spin-down channel are shown in
    Fig. [24](#fig8.2){reference-type="ref" reference="fig8.2"}-(a). In
    Fig. [24](#fig8.2){reference-type="ref" reference="fig8.2"}-(b)-(c)
    and (d) we show the convergence of the DOS for the spin-down
    channel, spin-up channel and both spin channels respectively. The
    convergence is assessed by looking at the number of states $N$
    computed by integrating the DOS up to the Fermi level using the
    formula
    $$N_{\uparrow/\downarrow} = \int_{-\infty}^{\epsilon_F}\!\! \mathrm{d}\epsilon\,\, f\ifmmode _{\mbox{\scriptsize{MV}}} \else $ _{\mbox{\scriptsize{MV}}} $~\fi(\epsilon,\uparrow/\downarrow)\, g(\epsilon,\uparrow/\downarrow),$$
    where $f\ifmmode _{\mbox{\scriptsize{MV}}} \else$ \_MV
    $~\fi(\epsilon,\uparrow) = \int_{-\infty}^{\epsilon} \mathrm{d}\epsilon'\,\widetilde{\delta}(\epsilon')$
    is the Marzari-Vanderbilt occupation number function, with
    $$\widetilde{\delta}(x) = \frac{2}{\sqrt{\pi}}e^{-[x-(1/\sqrt{2})]^2}(2\,-\,\sqrt{2}x), \quad x=\frac{\mu-\epsilon}{\sigma},$$
    where $\epsilon_F$ is the Fermi energy ($12.6256$ eV) and $\sigma$
    is the smearing ($0.02$ eV). $g(\epsilon,\uparrow)$ is the DOS from
    [wannier90]{.smallcaps} interpolation.

    <figure id="fig8.2">
    <p><br />
    </p>
    <figcaption>Panel (a) interpolated DOS for the down-spin channel of bcc
    iron for different <span
    class="math inline"><strong>k</strong></span>-mesh sizes. Panel (b)
    corresponding integrated DOS. The integral of the DOS is used as a
    convergence criterion. <span
    class="math inline"><em>N</em><sub> ↑  + ↓</sub></span> has been scaled
    such as the final value is equal to the total number of
    electrons.</figcaption>
    </figure>

## Projected versus maximally-localized Wannier functions {#projected-versus-maximally-localized-wannier-functions .unnumbered}

-   *Open one of the `.wout` files and search for "Initial state"; those
    are the projected WFs.*

    For the spin-up channel one finds

    ::: tcolorbox
        	 Initial State
          WF centre and spread    1  ( -0.000000, -0.000000, -0.000000 )     2.25930561
          WF centre and spread    2  ( -0.000000,  0.000000, -0.000000 )     2.32454089
          WF centre and spread    3  (  0.000000,  0.000000, -0.000000 )     2.32428592
          WF centre and spread    4  (  0.000000, -0.000000, -0.000000 )     2.32428592
          WF centre and spread    5  ( -0.000000,  0.000000, -0.000000 )     0.54443303
          WF centre and spread    6  (  0.000000, -0.000000, -0.000000 )     0.51353680
          WF centre and spread    7  (  0.000000,  0.000000, -0.000000 )     0.51353680
          WF centre and spread    8  (  0.000000,  0.000000, -0.000000 )     0.54447716
          WF centre and spread    9  (  0.000000,  0.000000,  0.000000 )     0.51347734
          Sum of centres and spreads (  0.000000, -0.000000, -0.000000 )    11.86187946
        	
    :::

    It is clear from the spreads and the centres that these are the
    projected WFs. In particular, WF 1 is the $s$-projected WF. WF 2-4
    are the $p$-projected WFs and WF 5-9 are the $d$-projected WF, with
    $e_g$ (5,8) and $t2_g$ (6,7,9) charachter, respectively (see
    Fig. [25](#fig8.5){reference-type="ref" reference="fig8.5"}).

    <figure id="fig8.5">

    <figcaption>3 representative MLWFs from the wannierisation via
    projections of 9 spin-up bands of iron. a) MLWF from projection onto 1
    <span class="math inline"><em>s</em></span> orbital. b) A representative
    of the MLWFs from projection onto <span
    class="math inline"><em>p</em></span> orbitals. c) A representative of
    the MLWFs from projection onto <span
    class="math inline"><em>d</em></span> orbitals.</figcaption>
    </figure>

-   The Wannier spreads have re-organized in two groups, 6+3; moreover,
    the six more diffuse WFs are off-centred: the initial atomic-like
    orbitals hybridized with one another, becoming more localized in the
    process.

    ::: tcolorbox
        	 Final State
          WF centre and spread    1  ( -0.709852,  0.000191,  0.000015 )     1.08935227
          WF centre and spread    2  ( -0.000015, -0.000041, -0.709852 )     1.08935223
          WF centre and spread    3  (  0.709852, -0.000191, -0.000015 )     1.08935227
          WF centre and spread    4  ( -0.000191, -0.709852,  0.000041 )     1.08935226
          WF centre and spread    5  (  0.000015,  0.000041,  0.709852 )     1.08935227
          WF centre and spread    6  ( -0.000000, -0.000000,  0.000000 )     0.43234437
          WF centre and spread    7  (  0.000000,  0.000000,  0.000000 )     0.43234440
          WF centre and spread    8  (  0.000191,  0.709852, -0.000041 )     1.08935228
          WF centre and spread    9  ( -0.000000,  0.000000,  0.000000 )     0.43234438
          Sum of centres and spreads ( -0.000000, -0.000000, -0.000000 )     7.83314672
        	
    :::

-   *The first plateau corresponds to atom-centred WFs of separate s, p,
    and d character, and the sharp drop signals the onset of the
    hybridization. With hindsight, we can redo steps 4 and 5 more
    efficiently using trial orbitals with the same character as the
    final MLWFs,*

    Fe : sp3d2;dxy;dxz;dyz

    With this choice the minimization converges much more rapidly as can
    be seen in Fig. [26](#fig8.4){reference-type="ref"
    reference="fig8.4"}-(a).

-   *Let us recompute the DOS using, instead of MLWFs, the WFs obtained
    by projecting onto s, p, and d-type trial orbitals.*

    <figure id="fig8.4">

    <figcaption>a) Convergence of <span
    class="math inline"><em>Ω</em></span> for two different sets of initial
    projections: <span
    class="math inline"><em>s</em>; <em>p</em>; <em>d</em></span> (solid
    black) and <span
    class="math inline"><em>s</em><em>p</em><sub>3</sub><em>d</em><sub>2</sub>; <em>d</em><sub><em>x</em><em>y</em></sub>; <em>d</em><sub><em>x</em><em>z</em></sub>; <em>d</em><sub><em>y</em><em>z</em></sub></span>
    (solid red). b) DOS with MLWFs (solid black) and projected <span
    class="math inline"><em>s</em>; <em>p</em>; <em>d</em></span> Wannier
    functions (solid blue).</figcaption>
    </figure>

## Orbital--projected DOS and exchange splitting {#orbitalprojected-dos-and-exchange-splitting .unnumbered}

*In order to obtain the partial DOS projected onto the $p$--type WFs,
add to the `.win` files*

dos_project = 2,3,4

and re-run `postw90`.

-   *Plot the projected DOS for both up-- and down--spin bands. Repeat
    for the $s$ and $d$ projections.*

    Results are shown in figure below
    (Fig. [27](#fig8.7){reference-type="ref" reference="fig8.7"}).

    <figure id="fig8.7">

    <figcaption>Partial DOS projected onto a) 1 <span
    class="math inline"><em>s</em></span>-like WF, b) 3<span
    class="math inline"><em>p</em></span>-like WFs and c) 5<span
    class="math inline"><em>d</em></span>-like WFs.</figcaption>
    </figure>

-   *The difference between corresponding values of the on-site energies
    the on-site energies
    $\braket{\boldsymbol{0}n\vert H \vert \boldsymbol{0}n}$ in
    `iron_up.wout` and in `iron_dn.wout` gives the exchange splittings
    for the individual orbitals.*

    Results are shown in Tab. [6](#tab8.2){reference-type="ref"
    reference="tab8.2"}.

    ::: {#tab8.2}
      --- ----------- -------------------------------------------------------------------------- ------------------------------------------------------------------------- ----------
      n    character   $\braket{\boldsymbol{0}n\vert H \vert \boldsymbol{0}n}$ for $\downarrow$   $\braket{\boldsymbol{0}n \vert H \vert \boldsymbol{0}n}$ for $\uparrow$   $\Delta$
                                                        \[eV\]                                                                    \[eV\]                                     \[eV\]
      1       $s$                                     21.307132                                                                  22.074648                                  0.767516
      2       $p$                                     26.353088                                                                  26.817526                                  0.464438
      3       $p$                                     26.352956                                                                  26.817207                                  0.464251
      4       $p$                                     26.352956                                                                  26.817207                                  0.464251
      5       $d$                                     10.531720                                                                  13.206631                                  2.67491
      6       $d$                                     10.775917                                                                  12.808277                                  2.03236
      7       $d$                                     10.775917                                                                  12.808277                                  2.03236
      8       $d$                                     10.532108                                                                  13.207139                                  2.67503
      9       $d$                                     10.775177                                                                  12.807388                                  2.03221
      --- ----------- -------------------------------------------------------------------------- ------------------------------------------------------------------------- ----------

      : Exchange splittings for individual orbitals in eV.
    :::

    []{#tab8.2 label="tab8.2"}

-   Compare their magnitudes with the splittings displayed by the
    orbital-projected DOS plots

# Cubic BaTiO$_3$ {#sec9:BaTiO3}

-   Outline : *Obtain MLWFs for a perovskite.*

![Unit cell of cubic BaTiO$_3$ crystal plotted with the
[XCrySDen]{.smallcaps} program.](figure/example09/BaTiO3.png){#fig9.0
width="0.25\\columnwidth"}

-   *Compute the MLWFs.*

    Converged values for the total spread functional and its components
    are shown in Tab. [7](#tab9.1){reference-type="ref"
    reference="tab9.1"}.

::: {#tab9.1}
  $\Omega$   $\Omega\ifmmode _{\mbox{\scriptsize{I}}} \else$ \_I $~\fi$   $\Omega\ifmmode _{\mbox{\scriptsize{OD}}} \else$ \_OD $~\fi$   $\Omega\ifmmode _{\mbox{\scriptsize{D}}} \else$ \_D $~\fi$   $N_{\mathrm{iter}}$
  ---------- ------------------------------------------------------------ -------------------------------------------------------------- ------------------------------------------------------------ ---------------------
  12.7187    12.5662                                                      0.1525                                                         0.000                                                        50

  : Converged values of the components of spread functional and their
  sums for cubic BaTiO$_3$ in Å$^2$.
:::

[]{#tab9.1 label="tab9.1"}

-   *Plot the second MLWF.*

    The result is shown in Fig. [29](#fig9.1){reference-type="ref"
    reference="fig9.1"}-(a) and -(b).

    <figure id="fig9.1">

    <figcaption>Top-view (a) and side-view (b) of the second MLWF in
    BaTiO<span class="math inline"><sub>3</sub></span> </figcaption>
    </figure>

-   *We can now simulate the ferroelectric phase by displacing the Ti
    atom. Regenerate the MLWFs (i.e., compute the ground-state charge
    density and Bloch states using pwscf, etc.) and look at the change
    in the second MLWF.*

    The result is shown in Fig. [30](#fig9.2){reference-type="ref"
    reference="fig9.2"}-(a) and -(b).

    <figure id="fig9.2">

    <figcaption>Top-view (a) and side-view (b) of the second MLWF in
    BaTiO<span class="math inline"><sub>3</sub></span> with the Ti atom
    displaced.</figcaption>
    </figure>

## Further ideas {#further-ideas .unnumbered}

-   *Look at MLWFs for other groups of bands.*

    Plots of MLWFs for other group of bands are shown in
    Fig. [31](#fig9.3){reference-type="ref"
    reference="fig9.3"}-(a)-(b)-(c)-(d)-(e) and -(f).

    <figure id="fig9.3">
    <p><br />
    </p>
    <figcaption>MLWFs for other group of bands.</figcaption>
    </figure>

-   *What happens if you form MLWFs for the whole valence manifold?*

    Some representative MLWFs from the wannierisation of the whole
    valence bands are shown in Fig. [32](#fig9.4){reference-type="ref"
    reference="fig9.4"}.

    <figure id="fig9.4">
    <p><br />
    </p>
    <figcaption>MLWFs formed from the whole valence manifold, i.e. from 20
    bands.</figcaption>
    </figure>

# Graphite {#sec10:graphite}

-   Outline: *Obtain MLWFs for the graphite (AB, Bernal).*

![Unit cell of Graphite plotted with the [XCrySDen]{.smallcaps}
program.](figure/example10/graphite.png){#fig10.0
width="0.25\\columnwidth"}

-   *Compute the MLWFs.*

    Converged values for the total spread functional and its components
    are shown in Tab. [8](#tab10.1){reference-type="ref"
    reference="tab10.1"}. Three MLWFs, one $\sigma$ and two $p_z$ on
    different layers are shown in
    Fig. [34](#fig10.1){reference-type="ref" reference="fig10.1"}(a),(b)
    and (c) respectively.

<figure id="fig10.1">

<figcaption>MLWFs for graphite. (a) <span
class="math inline"><em>σ</em></span>–like MLWF centred on a C–C bond.
(b) <span
class="math inline"><em>p</em><sub><em>z</em></sub></span>–like MLWF on
the top layer. (c) <span
class="math inline"><em>p</em><sub><em>z</em></sub></span>–like MLWF on
the bottom layer.</figcaption>
</figure>

::: {#tab10.1}
  $\Omega$   $\Omega\ifmmode _{\mbox{\scriptsize{I}}} \else$ \_I $~\fi$   $\Omega\ifmmode _{\mbox{\scriptsize{OD}}} \else$ \_OD $~\fi$   $\Omega\ifmmode _{\mbox{\scriptsize{D}}} \else$ \_D $~\fi$   $N_{\mathrm{iter}}$
  ---------- ------------------------------------------------------------ -------------------------------------------------------------- ------------------------------------------------------------ ---------------------
  7.3809     5.7641                                                       1.5874                                                         0.0293                                                       100

  : Converged values of the components of spread functional and their
  sums for graphite (AB, Bernal) in Å$^2$.
:::

[]{#tab10.1 label="tab10.1"}

# Silicon --- valence and low lying conduction states {#sec11:silicon}

![Unit cell of Silicon crystal plotted with the [XCrySDen]{.smallcaps}
program.](figure/example11/silicon.png){#fig11.0
width="0.25\\columnwidth"}

## Valence States {#valence-states .unnumbered}

-   Outline: *Obtain MLWFs for the valence bands of silicon.*

```{=html}
<!-- -->
```
-   *Inspect the output file `silicon.wout`. The total spread converges
    to its minimum value after just a few iterations. Note that the
    geometric centre of each MLWF lies at the centre of the Si--Si bond.
    Note also that the memory requirement for the minimisation of the
    spread is very low as the MLWFs are defined by just the $4\times4$
    unitary matrices $U(\mathbf{k})$.*

    Below a snippet from the `silicon.wout` output file

    ::: tcolorbox
         Final State
          WF centre and spread    1  ( -0.674701,  0.674701, -0.674701 )     1.59185520
          WF centre and spread    2  ( -0.674701, -0.674701,  0.674701 )     1.59185520
          WF centre and spread    3  (  0.674701,  0.674701,  0.674701 )     1.59185520
          WF centre and spread    4  (  0.674701, -0.674701, -0.674701 )     1.59185520
          Sum of centres and spreads ( -0.000000,  0.000000,  0.000000 )     6.36742081

                 Spreads (Ang^2)       Omega I      =     5.801375426
                ================       Omega D      =     0.000000000
                                       Omega OD     =     0.566045385
            Final Spread (Ang^2)       Omega Total  =     6.367420811
         ------------------------------------------------------------------------------
    :::

    Memory estimates may be found in the `MEMORY ESTIMATE` section of
    the `silicon.wout` file.

    ::: tcolorbox
         *============================================================================*
         |                              MEMORY ESTIMATE                               |
         |         Maximum RAM allocated during each phase of the calculation         |
         *============================================================================*
         |                        Disentanglement            1.57 Mb                  |
         |                            Wannierise:            0.47 Mb                  |
    :::

    Converged values for the total spread functional and its components
    are shown in Tab. [9](#tab11.1){reference-type="ref"
    reference="tab11.1"}.

    ::: {#tab11.1}
      MP mesh             $\Omega$   $\Omega\ifmmode _{\mbox{\scriptsize{I}}} \else$ \_I $~\fi$   $\Omega\ifmmode _{\mbox{\scriptsize{OD}}} \else$ \_OD $~\fi$   $\Omega\ifmmode _{\mbox{\scriptsize{D}}} \else$ \_D $~\fi$
      ------------------- ---------- ------------------------------------------------------------ -------------------------------------------------------------- ------------------------------------------------------------
      $4\times4\times4$   6.3674     5.8014                                                       0.5660                                                         0.0000

      : Converged values of the components of spread functional and
      their sum, given in Å$^2$.
    :::

    []{#tab11.1 label="tab11.1"}

-   *Plot the MLWFs* The four MLWFs with $\sigma$ character describing
    the valence manifold of Si are shown in
    Fig. [36](#fig11.1){reference-type="ref" reference="fig11.1"}(a),(b)
    and (c) respectively.

    <figure id="fig11.1">

    <figcaption>Four MLWFs for the valence manifold of Si.</figcaption>
    </figure>

## Valence + Conduction States {#valence-conduction-states .unnumbered}

-   Outline: *Obtain MLWFs for the valence and low--lying
    conduction-band states of Si. Plot the interpolated bandstructure.
    Apply a scissors correction to the conduction bands.*

```{=html}
<!-- -->
```
-   *Inspect the output file `silicon.wout`. The minimisation of the
    spread occurs in a two-step procedure. First, we minimise
    $\Omega\ifmmode _{\mbox{\scriptsize{I}}} \else$ \_I $~\fi$ -- this
    is the extraction of the optimal subspace in the disentanglement
    procedure. Then, we minimise
    $\Omega\ifmmode _{\mbox{\scriptsize{D}}} \else$ \_D
    $~\fi + \Omega\ifmmode _{\mbox{\scriptsize{OD}}} \else$ \_OD
    $~\fi$.*

    Converged values for the total spread functional and its components
    are shown in Tab. [10](#tab11.2){reference-type="ref"
    reference="tab11.2"}. The two groups of four MLWFs with $sp3$
    character are shown in Fig. [37](#fig11.2){reference-type="ref"
    reference="fig11.2"}

    ::: tcolorbox
                           Extraction of optimally-connected subspace
                           ------------------------------------------
         +---------------------------------------------------------------------+<-- DIS
         |  Iter     Omega_I(i-1)      Omega_I(i)      Delta (frac.)    Time   |<-- DIS
         +---------------------------------------------------------------------+<-- DIS
               1      12.97640155      12.44630235       4.259E-02      0.00    <-- DIS
               .		.					.				.			 .

              79      12.33580893      12.33580893      -6.531E-11      0.23    <-- DIS
              80      12.33580893      12.33580893      -5.241E-11      0.23    <-- DIS

                     <<<      Delta < 1.000E-10  over  3 iterations     >>>
                     <<< Disentanglement convergence criteria satisfied >>>

                Final Omega_I    12.33580893 (Ang^2)

         +----------------------------------------------------------------------------+
        	
    :::

    ::: tcolorbox
        	 Final State
          WF centre and spread    1  (  1.807167,  1.807167,  1.807167 )     2.01695824
          WF centre and spread    2  (  1.807167,  0.891636,  0.891636 )     2.01695823
          WF centre and spread    3  (  0.891636,  1.807167,  0.891636 )     2.01695823
          WF centre and spread    4  (  0.891636,  0.891636,  1.807167 )     2.01695824
          WF centre and spread    5  (  0.226733,  0.226733,  0.226733 )     2.37014516
          WF centre and spread    6  (  0.226733, -0.226733, -0.226733 )     2.37014508
          WF centre and spread    7  ( -0.226733,  0.226733, -0.226733 )     2.37014515
          WF centre and spread    8  ( -0.226733, -0.226733,  0.226733 )     2.37014514
          Sum of centres and spreads (  5.397608,  5.397608,  5.397608 )    17.54841346

                 Spreads (Ang^2)       Omega I      =    12.335808933
                ================       Omega D      =     0.177593840
                                       Omega OD     =     5.035010692
            Final Spread (Ang^2)       Omega Total  =    17.548413465
         ------------------------------------------------------------------------------
         	
    :::

    ::: {#tab11.2}
      MP mesh             $\Omega$   $\Omega\ifmmode _{\mbox{\scriptsize{I}}} \else$ \_I $~\fi$   $\Omega\ifmmode _{\mbox{\scriptsize{OD}}} \else$ \_OD $~\fi$   $\Omega\ifmmode _{\mbox{\scriptsize{D}}} \else$ \_D $~\fi$
      ------------------- ---------- ------------------------------------------------------------ -------------------------------------------------------------- ------------------------------------------------------------
      $4\times4\times4$   17.54841   12.3358                                                      5.03501                                                        0.17759

      : Converged values of the components of spread functional and
      their sum, given in Å$^2$.
    :::

    []{#tab11.2 label="tab11.2"}

    <figure id="fig11.2">
    <p><br />
    </p>
    <figcaption>Eight MLWFs with <span
    class="math inline"><em>s</em><em>p</em>3</span> character, four on each
    Si atom in the unit cell.</figcaption>
    </figure>

-   *Plot the bandstructure.*

    The interpolated bandstructure is given in
    Fig. [38](#fig11.3){reference-type="ref" reference="fig11.3"}.

    ![Bandstructure of silicon from DFT calculation (solid black) and
    from Wannier interpolation (solid
    red).](figure/example11/silicon_bandstructure.pdf){#fig11.3
    width="0.7\\columnwidth"}

## Further ideas {#further-ideas-1 .unnumbered}

-   *Compare the Wannier-interpolated bandstructure with the full pwscf
    bandstructure with a finer $k$-point grid.*

    Result for a $8\times8\times8$ mesh is shown in
    Fig. [39](#fig11.4){reference-type="ref" reference="fig11.4"}.

    ![Bandstructure of silicon from DFT calculation (solid black) and
    from Wannier interpolation with a $4\times4\times4$ mesh (solid red)
    and $8\times8\times8$ mesh (solid
    blue).](figure/example11/silicon_bs_DFT_vs_W90_finer_grid.pdf){#fig11.4
    width="0.7\\columnwidth"}

-   *Compute four MLWFs spanning the low-lying conduction states.*

    The MLWFs spanning the 4 low-lying conduction states are shown in
    Fig. [40](#fig11.5){reference-type="ref" reference="fig11.5"}. The
    initial projections were 4 $sp3$ on the Si atom at (0,0,0).

    <figure id="fig11.5">

    <figcaption></figcaption>
    </figure>

# Benzene --- valence and low lying conduction states {#sec12:benzene}

![Benzene molecule in periodic cell plotted with the
[XCrySDen]{.smallcaps} program.](figure/example12/benzene.png){#fig12.0
width="0.25\\columnwidth"}

## Valence States {#valence-states-1 .unnumbered}

-   Outline: *Obtain MLWFs for the valence bands of benzene.*

```{=html}
<!-- -->
```
-   *Inspect the output file benzene.wout. The total spread converges to
    its minimum value after just a few iterations.*

    Convergence of total spread $\Omega$ is shown in
    Fig. [43](#fig12.1){reference-type="ref" reference="fig12.1"}. The
    spread converges very quickly, and after only 15 iterations the
    $|\Delta\Omega|$ is already below $10^{-8}$. Below is shown the
    final state of the minimization, after 22 iterations:

    ::: tcolorbox
        	 Final State
          WF centre and spread    1  ( -6.875141,  7.935472,  7.937658 )     0.65233309
          WF centre and spread    2  (  4.748245,  7.935472,  7.937658 )     0.65234300
          WF centre and spread    3  (  5.809324,  6.097678,  7.937658 )     0.65186298
          WF centre and spread    4  (  7.816785, -7.396927,  7.648002 )     1.20135958
          WF centre and spread    5  (  7.816785, -7.396927, -7.648002 )     1.20135948
          WF centre and spread    6  (  7.936083,  7.324218, -7.937658 )     0.60992987
          WF centre and spread    7  (  7.935073, -6.095184,  7.937658 )     0.65161009
          WF centre and spread    8  (  6.874208, -6.712573,  7.937658 )     0.61134798
          WF centre and spread    9  (  6.874224,  6.850489, -7.648754 )     1.20500764
          WF centre and spread   10  ( -7.936214,  6.097679,  7.937658 )     0.65186256
          WF centre and spread   11  (  5.813353, -6.095183,  7.937658 )     0.65160955
          WF centre and spread   12  (  6.874224,  6.850489,  7.648754 )     1.20500766
          WF centre and spread   13  (  5.812342,  7.324217,  7.937658 )     0.60993084
          WF centre and spread   14  (  5.931632, -7.396946,  7.648001 )     1.20136423
          WF centre and spread   15  (  5.931632, -7.396946, -7.648001 )     1.20136424
          Sum of centres and spreads ( 71.362557,  7.925029, 55.563607 )    12.95829277

                 Spreads (Ang^2)       Omega I      =    10.455434168
                ================       Omega D      =     0.000000000
                                       Omega OD     =     2.502858604
            Final Spread (Ang^2)       Omega Total  =    12.958292772
         ------------------------------------------------------------------------------
        	
    :::

-   *Plot the MLWFs 2-4*

    MLWFs are shown in Fig. [42](#fig12.2){reference-type="ref"
    reference="fig12.2"}.

    <figure id="fig12.2">

    <figcaption>MLWFs 2, 3 and 4, with Vesta from
    <span><code>cube</code></span> format.</figcaption>
    </figure>

![Convergence of total spread $\Omega$. The red curve refers to the left
y-axis, i.e. the actual value of the total spread at each iteration. The
blue curve refers to the right y-axis, i.e. the absolute difference
between between the spread functional at iteration $i$ and $i-1$, i.e.
$\Delta\Omega$.](figure/example12/spread_convergence.pdf){#fig12.1
width="0.7\\columnwidth"}

## Valence + Conduction States {#valence-conduction-states-1 .unnumbered}

-   Outline: *Obtain MLWFs for the valence and low-lying conduction
    states of benzene.*

```{=html}
<!-- -->
```
-   *First, we minimise $\Omega_I$. Then, we minimise
    $\Omega_D + \Omega_{OD}$.*

    Extract from the `.wout` output file for the disentanglement
    procedure with initial and final value of $\Omega_I$

    ::: tcolorbox
                           Extraction of optimally-connected subspace
                           ------------------------------------------
         +---------------------------------------------------------------------+<-- DIS
         |  Iter     Omega_I(i-1)      Omega_I(i)      Delta (frac.)    Time   |<-- DIS
         +---------------------------------------------------------------------+<-- DIS
               1      14.77292507      14.36793746       2.819E-02      0.06    <-- DIS
               .		.					.				.			 .
               .		.					.				.			 .
              76      14.26979011      14.26979011       7.234E-11      0.34    <-- DIS

                     <<<      Delta < 1.000E-10  over  3 iterations     >>>
                     <<< Disentanglement convergence criteria satisfied >>>

                Final Omega_I    14.26979011 (Ang^2)

         +----------------------------------------------------------------------------+
    :::

    <figure id="fig12.3">

    <figcaption>MLWFs 1, 3 and 13 with Vesta from
    <span><code>cube</code></span> format.</figcaption>
    </figure>

    Below a snippet from the `.wout` output file, showing the finale
    state of the minimisation of $\Omega_D$ and $\Omega_OD$.

    ::: tcolorbox
         Final State
          WF centre and spread    1  ( -6.872991, -7.937658,  7.937657 )     0.64685210
          WF centre and spread    2  ( -7.937084, -6.094542, -7.937656 )     0.64620663
          WF centre and spread    3  (  5.810194, -6.094542,  7.937654 )     0.64620672
          WF centre and spread    4  (  4.746097, -7.937658,  7.937657 )     0.64685949
          WF centre and spread    5  (  5.810193,  6.094541, -7.937657 )     0.64620720
          WF centre and spread    6  ( -7.937083,  6.094541,  7.937655 )     0.64620743
          WF centre and spread    7  (  6.874209,  6.725491,  7.937658 )     0.58709128
          WF centre and spread    8  (  6.874209, -6.725491,  7.937657 )     0.58709059
          WF centre and spread    9  (  5.824168, -7.331141,  7.937657 )     0.58581741
          WF centre and spread   10  (  5.824168,  7.331141, -7.937658 )     0.58581778
          WF centre and spread   11  (  7.924257,  7.331142,  7.937657 )     0.58581692
          WF centre and spread   12  (  7.924257, -7.331141, -7.937658 )     0.58581674
          WF centre and spread   13  ( -7.514755,  7.937658, -7.937656 )     1.57077288
          WF centre and spread   14  (  7.622957, -6.642249,  7.937655 )     1.58413409
          WF centre and spread   15  (  6.125499, -6.642276, -7.937651 )     1.58406298
          WF centre and spread   16  (  5.387848, -7.937658, -7.937656 )     1.57079972
          WF centre and spread   17  (  6.125498,  6.642275,  7.937656 )     1.58407213
          WF centre and spread   18  (  7.622957,  6.642249, -7.937653 )     1.58414266
          Sum of centres and spreads ( 60.234596,-15.875317, 15.875317 )    16.87397474

                 Spreads (Ang^2)       Omega I      =    14.269790106
                ================       Omega D      =     0.000000000
                                       Omega OD     =     2.604184635
            Final Spread (Ang^2)       Omega Total  =    16.873974742
         ------------------------------------------------------------------------------
    :::

```{=html}
<!-- -->
```
-   *Plot the MLWFs 1, 7 and 13.*

    MLWFs are shown in Fig. [44](#fig12.3){reference-type="ref"
    reference="fig12.3"}.

# (5,5) Carbon Nanotube --- Transport properites {#sec13:cnt}

-   Outline: *Obtain the bandstructure, quantum conductance and density
    of states of a metallic (5,5) carbon nanotube.*

<figure id="fig13.0">

<figcaption>5 unit cells for the carbon nanotube system from a) side
view and b) prospective top view plotted with the <span
class="smallcaps">XCrySDen</span> program.</figcaption>
</figure>

-   *Run pwscf and `wannier90`. Inspect the output file `cnt55.wout`.
    The minimisation of the spread occurs in a two-step proce- dure.
    First, we minimise $\Omega_I$. Then, we minimise
    $\Omega_D + \Omega_{OD}$.*

    Below, an extract from the `.wout` file showing a summary of the
    disentanglement procedure (minimisation of $\Omega_I$)

::: tcolorbox
                       Extraction of optimally-connected subspace
                       ------------------------------------------
     +---------------------------------------------------------------------+<-- DIS
     |  Iter     Omega_I(i-1)      Omega_I(i)      Delta (frac.)    Time   |<-- DIS
     +---------------------------------------------------------------------+<-- DIS
           1      33.96797815      33.91073784       1.688E-03      0.00    <-- DIS
           2      33.92937273      33.90274574       7.854E-04      0.02    <-- DIS
           .		.					.				.			 .
           .		.					.				.			 .

          45      33.89889125      33.89889125       4.172E-11      0.50    <-- DIS
          46      33.89889125      33.89889125       1.626E-11      0.51    <-- DIS

                 <<<      Delta < 1.000E-10  over  3 iterations     >>>
                 <<< Disentanglement convergence criteria satisfied >>>

            Final Omega_I    33.89889125 (Ang^2)

     +----------------------------------------------------------------------------+
    	
:::

Below, an extract from the `.wout` file showing the final state for the
minimisation of $\Omega_D + \Omega_{OD}s$

::: tcolorbox
        	 Final State
      WF centre and spread    1  ( -6.875141,  7.935472,  7.937658 )     0.65233309
      WF centre and spread    2  (  4.748245,  7.935472,  7.937658 )     0.65234300
      WF centre and spread    3  (  5.809324,  6.097678,  7.937658 )     0.65186298
      WF centre and spread    4  (  7.816785, -7.396927,  7.648002 )     1.20135958
      WF centre and spread    5  (  7.816785, -7.396927, -7.648002 )     1.20135948
      WF centre and spread    6  (  7.936083,  7.324218, -7.937658 )     0.60992987
      WF centre and spread    7  (  7.935073, -6.095184,  7.937658 )     0.65161009
      WF centre and spread    8  (  6.874208, -6.712573,  7.937658 )     0.61134798
      WF centre and spread    9  (  6.874224,  6.850489, -7.648754 )     1.20500764
      WF centre and spread   10  ( -7.936214,  6.097679,  7.937658 )     0.65186256
      WF centre and spread   11  (  5.813353, -6.095183,  7.937658 )     0.65160955
      WF centre and spread   12  (  6.874224,  6.850489,  7.648754 )     1.20500766
      WF centre and spread   13  (  5.812342,  7.324217,  7.937658 )     0.60993084
      WF centre and spread   14  (  5.931632, -7.396946,  7.648001 )     1.20136423
      WF centre and spread   15  (  5.931632, -7.396946, -7.648001 )     1.20136424
      Sum of centres and spreads ( 71.362557,  7.925029, 55.563607 )    12.95829277

             Spreads (Ang^2)       Omega I      =    10.455434168
            ================       Omega D      =     0.000000000
                                   Omega OD     =     2.502858604
        Final Spread (Ang^2)       Omega Total  =    12.958292772
     ------------------------------------------------------------------------------
:::

1.  *Note that the initial $p_z$ projections on the carbon atoms are
    oriented in the radial direction with respect to the nanotube axis.*

    ` `

    > Begin Projections
    >
    > Ang
    >
    > c= 3.3780, -0.7128, -0.6157 :pz :z= 3.3780, -0.7128, 0.0000
    > :x=0,0,1

2.  *The interpolated bandstructure is written to `cnt55_band.agr`*

    To plot the interpolated bands, the quantum conductance and the
    Density of States as shown in Fig. 6 in the Wannier90 tutorial, one
    can use the `xmgrace` program.

::: tcolorbox
Run the `xmgrace` plotting program from command line as

` `

> \$ \> xmgrace

Before importing the data to be plotted, we have to reorganize the
layout by selecting

> Edit $\mapsto$ Arrange graphs...

here we can generate a grid of graphs by selecting the number of columns
and rows from the drop menus. For this particular example, we want to
increase the number of columns to 3, i.e. `Cols: 3` and leave the number
of rows to 1 in the `Matrix` section. Moreover, we don't want any gap
between the graphs so we also need to modify the value of `Hgap/width`
in the bottom `Spacing` section, i.e `Hgap/width 0`. Once we have
generated the three graphs we need to import the data. This can be
achieved by

> Data $\mapsto$ Import $\mapsto$ ASCII...

The three files to import are `cnt_band.agr`, `cnt_qc.dat` and
`cnt_dos.dat`, respectively. For each file we need to select the graph
in the `Read to graph:` section, i.e. `G(0), G(1)` and `G(2)`,
respectively.

In order to flip the x-axis with the y-axis, one need to perform the
following

` `

> Data $\mapsto$ Transformations $\mapsto$ Evaluate expressions...

In the `Formula:` section write

` `

> s1.x=s0.y; s1.y=s0.x

and then click `apply`.
:::

![Reproduction of Fig. 6 in the wannier90
tutorial.](figure/example13/cnt55_band.pdf){#fig13.1
width="0.7\\columnwidth"}

# Linear Sodium Chain --- Transport properties {#sec14:nachain}

-   Outline: *Compare the quantum conductance of a periodic linear chain
    of Sodium atoms with that of a defected chain*

<figure id="fig14.1">
<p><br />
<span id="fig14.1" label="fig14.1"></span></p>
<figcaption>Unit cell of a periodic linear Sodium chain (left panel) and
of a defected linear Sodium chain (right panel) plotted with the <span
class="smallcaps">XCrySDen</span> program. The former consists of 2 Na
atom per unit cell (6 unit cells have been drawn for comparison with the
defected system). The latter consists of 13 Na atoms per unit
cell.</figcaption>
</figure>

1.  *Run pwscf and wannier90 for the periodic and defected systems.*

2.  *Compare the quantum conductance of the periodic (bulk) calculation
    with the defected (LCR) calculation.*

    The quantum conductance and the DOS are shown in
    Fig. [48](#fig14.2){reference-type="ref" reference="fig14.2"}.

    ![DOS (left) and quantum conductance (right) of periodic (solid
    black) and defected (solid red) Sodium linear
    chain.](figure/example14/Na_chain_dos_qc.pdf){#fig14.2
    width="1.0\\columnwidth" height="0.5\\columnwidth"}

# (5,0) Carbon Nanotube --- Transport properties {#sec15:CNT}

# Silicon --- Boltzmann transport {#sec16:SiBZT}

-   Outline: *Obtain MLWFs for the valence and low-lying conduction
    states of Si. Calculate the electrical conductivity, the Seebeck
    coefficient and the thermal conductivity in the constant relaxation
    time approximation using the `BoltzWann` module.*

![Unit cell of Silicon crystal plotted with the [XCrySDen]{.smallcaps}
program.](figure/example11/silicon.png){#fig16.0
width="0.25\\columnwidth"}

-   For this example we are only going to show the solutions from point
    6 onwards, as the first 5 steps are the usual steps to obtain MLWFs.

-   *Run `postw90` to calculate the transport coefficients.*

-   *Inspect the output file `Si.wpout`. Check if no warnings are
    issued. Note that if no special flags are passed to BoltzWann, it
    assumes that the *ab initio* calculation did not include
    magnetization effects, and thus it sets to 2 the number of electrons
    per state.*

    Below the section in the `Si.wpout` relative to the Boltzmann
    transport, where it reports the number of electrons per state and
    the relaxation time in fs.

::: tcolorbox
     *---------------------------------------------------------------------------*
     |                   Boltzmann Transport (BoltzWann module)                  |
     *---------------------------------------------------------------------------*
     | Please cite the following paper when publishing results obtained using    |
     | the BoltzWann module:                                                     |
     | G. Pizzi, D. Volja, B. Kozinsky, M. Fornari, and N. Marzari,              |
     | Comp. Phys. Comm. 185, 422 (2014); DOI:10.1016/j.cpc.2013.09.015          |
     *---------------------------------------------------------------------------*

       Calculating Transport Distribution function (TDF) and DOS...
         k-grid used for band interpolation in BoltzWann: 40x40x40
         Number of electrons per state: 2
         Relaxation time (fs):    10.00000000
       TDF and DOS calculated.

       Transport properties calculated.

     *---------------------------------------------------------------------------*
     |                        End of the BoltzWann module                        |
     *---------------------------------------------------------------------------*
:::

*Using your favourite plotting program, plot the `Si_boltzdos.dat` file
to inspect the DOS.*

Plot shown in Fig. [50](#fig16.1){reference-type="ref"
reference="fig16.1"}-(a).

*Using your favourite plotting program, plot columns 1 and 3 of the
`Si_seebeck.dat` file to inspect the $S_{xx}$ component of the Seebeck
coefficient as a function of the chemical potential $\mu$, at $T = 300$
K.*

Plot shown in Fig. [50](#fig16.1){reference-type="ref"
reference="fig16.1"}-(b).

<figure id="fig16.1">

<figcaption>Panel (a) DOS of Silicon computed with
<span>BoltzWann</span>. Panel (b) <span
class="math inline"><em>S</em><sub><em>x</em></sub><em>x</em></span>
component of the Seebeck tensor as function of the chemical potential
<span class="math inline"><em>μ</em></span> computed with
<span>BoltzWann</span> at <span
class="math inline"><em>T</em> = 300</span> K.</figcaption>
</figure>

## Further ideas {#further-ideas-2 .unnumbered}

-   *Change the interpolation to a $60\times60\times60$ mesh and run
    again postw90 to check if the results for the transport properties
    are converged.*

    Plot of the two DOS with $40\times40\times40$ (red line) and
    $60\times60\times60$ (blue line) are shown in
    Fig. [51](#fig16.3){reference-type="ref" reference="fig16.3"}. We
    can see that all the peaks for the valence and conduction states in
    the $60\times60\times60$ DOS are also reproduced in the
    $40\times40\times40$ DOS (even though there is some noise, which
    however does not affect the qualitative description.)

    ![Convergence of
    DOS](figure/example16/Si_boltzdos_convergence.pdf){#fig16.3
    width="0.7\\columnwidth"}

-   *Change the `Si.win` input file so that it calculates the transport
    coefficients for temperatures from $300$ to $700$ K, with steps of
    $200$ K. Rerun `postw90` and verify that the increase in execution
    time is negligible (in fact, most of the time is spent to
    interpolate the band structure on the k mesh). Plot the Seebeck
    coefficient for the three temperatures $T = 300$ K, $T = 500$ K and
    $T = 700$ K. To do this, you have to filter the `Si_seebeck.dat` to
    select only those lines where the second column is equal to the
    required temperature. A possible script to select the $S_{xx}$
    component of the Seebeck coefficient for $T = 500$ K using the
    `awk/gawk` command line program is the following:*

    ` `

    > awk 'if (\$2 == 500) print \$1, \$3;' \< Si_seebeck.dat \>
    > Si_seebeck_xx_500K.dat

    Below is shown the total wall-time for the two calculations done
    with the original set up, i.e. $T_{min} = T_{max} = 300$ K and
    $T_{min} = 300$ K, $T_{max} = 700$ K, $\Delta T = 200$ K.

    ::: tcolorbox
         Total Execution Time          16.356 (sec)
    :::

    ::: tcolorbox
         Total Execution Time          16.108 (sec)
    :::

    The plot of $S_{xx}(\mu)$ for different values of $T$ is shown in
    Fig. [52](#fig16.4){reference-type="ref" reference="fig16.4"}

    ![$S_{xx}$ component of the Seebeck tensor as function of the
    chemical potential $\mu$ for different values of the temperature:
    $T = 300$ K (solid purple), $T = 500$ K (solid green) and $T = 700$
    K (solid blue).](figure/example16/Si_seebeck_T.pdf){#fig16.4
    width="0.7\\columnwidth"}

-   *Try to calculate the Seebeck coefficient as a function of the
    temperature, for a n--doped sample with, e.g., $n = 10^{18}$
    cm$^{-3}$. Note that to this aim, you need to calculate consistently
    the value $\mu(T)$ of the chemical potential as a function of the
    temperature, so as to reproduce the given value of $n$. Then, you
    have to write a small program/script to interpolate the output of
    BoltzWann, that you should have run on a suitable grid of $(\mu,T)$
    points.*

    ASSUMPTIONS: 1) The addition of a $n-$type dopant does not modify
    the electronic structure, it only moves the Fermi level up; 2) The
    density of states is temperature-independent. $\mu(T)$ is a
    decreasing monotonic function of $T$.

    To obtain a $\mu(T)$ in a consistent way we use the above
    assumptions and the following equation:
    $$N_c + N_v = \int_{-\infty}^{+\infty}\! \mathrm{d}\varepsilon \, g(\varepsilon,T=0)\,f(\varepsilon,\mu(T)),\label{eq16.1}$$
    where $N_v=8$, number of valence electrons per unit cell when no
    dopants are considered, $N_c= nV_{cell}$ is the number of carriers
    per unit cell ($V_{cell}$ is the volume of the unit cell in
    cm$^{-3}$). $g(\varepsilon,T=0)$ is the density of states at $T=0$K
    and by assumption it does not change with $T$. Finally,
    $f(\varepsilon,\mu(T))$ is the Fermi-Dirac distribution as a
    function of $\varepsilon$ and $T$
    $$f(\varepsilon,\mu(T)) = \frac{1}{1 + \exp[\frac{\varepsilon - \mu(T)}{\kappa_B T}]}$$
    For each $T$ we find the value of the $\mu(T)$ such as the integral
    is (approximately) $N_c + N_v$. [^2] The values of $\mu$ for $T$ in
    the range \[300 K--700 K\] are shown in
    Tab. [11](#tab16.1){reference-type="ref" reference="tab16.1"}

    ::: {#tab16.1}
       $T$ \[K\]   $\mu$ \[eV\]
      ----------- --------------
          300         6.839
          400         6.798
          500         6.752
          600         6.707
          700         6.677

      : Values of the chemical potential $\mu$ in eV as a function of
      $T$ in K, computed by numerically solving
      Eq. [\[eq16.1\]](#eq16.1){reference-type="ref"
      reference="eq16.1"}.
    :::

    []{#tab16.1 label="tab16.1"}

    In practice we do not perform an interpolation but we run a single
    calculation with $\Delta \mu = 0.001$ eV since these are not
    expensive and then we filter out the result from `Si_seebeck.dat`
    with the following simple script

    ` `

    > mulist='cat mu.dat \| awk 'printf \"i4\" \$1''; i=0; for mu in
    > \$mulist; do i='echo \$i+1\|bc' ; cat Si_seebeck.dat \| awk -v
    > \"mu=\$mu\" 'if(\$1==mu) print \$1,\$2,\$3,\$7,\$11' \| awk -v
    > \"Tcol=\$i\" 'if(NR==Tcol) print \$1, \$2, \$3, \$4, \$5' \>\>
    > Si_seebeck_vs_T.dat;done

    where `mu.dat` is a data file containing the second column of
    Tab. [11](#tab16.1){reference-type="ref" reference="tab16.1"}.
    Fig. [53](#fig16.5){reference-type="ref" reference="fig16.5"} shows
    the plots of the diagonal coefficients of the Seebeck tensor with
    respect to $T$ generated by the above script and stored in
    `Si_seebeck_vs_T.dat`.

    ![](figure/example16/Si_seebeck_vs_T.pdf){#fig16.5
    width="0.7\\columnwidth"}

# Iron --- Spin-orbit-coupled bands and Fermi-surface contours {#sec17:IronSO}

-   Outline: *Plot the spin-orbit-coupled bands of ferromagnetic bcc Fe.
    Plot the Fermi-surface contours on a plane in the Brillouin zone.*

![Unit cell of Iron crystal plotted with the [XCrySDen]{.smallcaps}
program.](figure/example08/iron.png){#fig17.0 width="0.25\\columnwidth"}

-   Compute the MLWFs and compute the energy eigenvalues and spin
    expectation values.

    The final state for all the 18 MLWFs is

    ::: tcolorbox
        	 Final State
          WF centre and spread    1  ( -0.709848,  0.000000,  0.000000 )     1.08973288
          WF centre and spread    2  ( -0.685480, -0.000000,  0.000000 )     1.10285536
          WF centre and spread    3  (  0.709848, -0.000000,  0.000000 )     1.08973288
          WF centre and spread    4  (  0.685480, -0.000000,  0.000000 )     1.10285536
          WF centre and spread    5  ( -0.000000, -0.709848,  0.000000 )     1.08973287
          WF centre and spread    6  (  0.000000, -0.685480,  0.000000 )     1.10285536
          WF centre and spread    7  (  0.000000,  0.709848,  0.000000 )     1.08973288
          WF centre and spread    8  (  0.000000,  0.685480,  0.000000 )     1.10285536
          WF centre and spread    9  (  0.000000, -0.000000, -0.709835 )     1.08977307
          WF centre and spread   10  (  0.000000, -0.000000, -0.685503 )     1.10302800
          WF centre and spread   11  ( -0.000000, -0.000000,  0.709835 )     1.08977304
          WF centre and spread   12  (  0.000000, -0.000000,  0.685503 )     1.10302800
          WF centre and spread   13  ( -0.000000, -0.000000, -0.000000 )     0.43232470
          WF centre and spread   14  ( -0.000000, -0.000000, -0.000000 )     0.41118748
          WF centre and spread   15  (  0.000000,  0.000000, -0.000000 )     0.43232470
          WF centre and spread   16  (  0.000000,  0.000000, -0.000000 )     0.41118748
          WF centre and spread   17  ( -0.000000,  0.000000,  0.000000 )     0.43232866
          WF centre and spread   18  ( -0.000000,  0.000000,  0.000000 )     0.41119649
          Sum of centres and spreads (  0.000000,  0.000000,  0.000000 )    15.68650457

                 Spreads (Ang^2)       Omega I      =    11.898334117
                ================       Omega D      =     0.031570932
                                       Omega OD     =     3.756599523
            Final Spread (Ang^2)       Omega Total  =    15.686504572
         ------------------------------------------------------------------------------
        	
    :::

-   *To plot the bands using `python`*

    ` `

    > \$\> python Fe-bands.py

    The interpolated band structure of Fe with spin-orbit interaction
    using the module `kpath` is shown in
    Fig. [55](#fig17.1){reference-type="ref" reference="fig17.1"}. The
    color scheme is used to show the expectation value of the spin
    operator $\hat{S}_z$ in units of $\hbar/2$.

    ![[wannier90]{.smallcaps} interpolated bands of Fe computed from a
    DFT calculation with spin-orbit interaction. Colour-scheme shows the
    expectation value $\braket{\hat{S}_z}$ in units of
    $\hbar/2$.](figure/example17/Fe_bandstructure.pdf){#fig17.1
    width="0.6\\columnwidth"}

-   *Next we plot the Fermi-surface contours on the (010) plane
    $k_y = 0$, using the `kslice` module.*

    <figure id="fig17.2">

    <figcaption></figcaption>
    </figure>

    ## Further ideas {#further-ideas-3 .unnumbered}

    -   *Redraw the Fermi surface contours on the (010) plane starting
        from a calculation without spin-orbit coupling (SOC), by adding
        to the input files iron\_{up,down}.win in Example 8.*

        The Fermi surface contours on the (010) plane without SOC are
        shown in Fig. [57](#fig17.4){reference-type="ref"
        reference="fig17.4"}-(a).

    -   *For a spinor calculation we can still spin-decompose the DOS.*

        <figure id="fig17.4">

        <figcaption>Spin-decomposed DOS (panel a) with spin-up (red) and
        spin-down (blue) components. Projected DOS on odd-indexed MLWFs (red)
        and even-indexed (blue).</figcaption>
        </figure>

# Iron---Berry curvature, anomalous Hall conductivity and optical conductivity {#sec18:IronBerry}

-   Outline: *Calculate the Berry curvature, anomalous Hall
    conductivity, and (magneto)optical conductivity of ferromagnetic bcc
    Fe with spin-orbit coupling. In preparation for this example it may
    be useful to read Ref.  and Ch. 11 of the User Guide.*

```{=html}
<!-- -->
```
-   *Compute the MLWFs and compute the energy eigenvalues and spin
    expectation values.*

    These are the same six steps of
    Ex. [17](#sec17:IronSO){reference-type="ref"
    reference="sec17:IronSO"} and therefore the results are not going to
    be showed here again.

## Berry curvature plots {#berry-curvature-plots .unnumbered}

-   *The Berry curvature $\Omega_{\alpha\beta}(\mathbf{k})$ of the
    occupied states is defined in Eq. (11.18) of the User Guide.* *Plot
    the Berry curvature component $\Omega_z(\ifmmode  \mathbf{k}  \else$
    $~\fi) = \Omega_{xy}(\ifmmode  \mathbf{k}  \else$ $~\fi)$ along the
    magnetization direction.*

    The Fermi energy should be $12.6283$ eV. With this value we obtain
    the energy bands and the Berry curvature component
    $\Omega_z(\ifmmode  \mathbf{k}  \else$
    $~\fi) = \Omega_{xy}(\ifmmode  \mathbf{k}  \else$ $~\fi)$ along
    high-symmetry points shown in
    Fig. [58](#fig18.1){reference-type="ref" reference="fig18.1"} and
    Fig. [59](#fig18.2){reference-type="ref" reference="fig18.2"}. Eq.
    (11.18) of the User Guide is reported below for completeness.

    $$\Omega_{\alpha\beta}(\mathbf{k}) = \sum_{n}^{occ} f_{n\mathbf{k}}\Omega_{n,\alpha\beta},$$
    with
    $$\Omega_{n,\alpha\beta} = \varepsilon_{\alpha\beta\gamma}\Omega_{n,\gamma} = -2\;\mathrm{Im}\braket{\nabla_{k_\alpha}\ifmmode u_{n\ifmmode  \mathbf{k}  \else $ \mathbf{k} $~\fi} \else $ u_{n\ifmmode  \mathbf{k}  \else $ \mathbf{k} $~\fi} $~\fi\vert \nabla_{k_\beta}\ifmmode u_{n\ifmmode  \mathbf{k}  \else $ \mathbf{k} $~\fi} \else $ u_{n\ifmmode  \mathbf{k}  \else $ \mathbf{k} $~\fi} $~\fi},$$
    where the Greek letters indicate Cartesian coordinates,
    $\varepsilon_{\alpha\beta\gamma}$ is the Levi-Civita antisymmetric
    tensor, and
    $\ket{u_{n\ifmmode  \mathbf{k}  \else $ \mathbf{k} $~\fi}}$s are the
    cell-periodic Bloch functions.

![Band structure of Fe along symmetry lines
$\Gamma$-H-P-N-$\Gamma$-H-N-$\Gamma$-P-N.](figure/example18/Fe_bandstructure.pdf){#fig18.1
width="0.8\\columnwidth"}

![Berry curvature $\Omega_z(\ifmmode  \mathbf{k}  \else$ $~\fi)$ in Fe
along symmetry lines.](figure/example18/Fe_Berry_phase.png){#fig18.2
width="0.7\\columnwidth"}

-   *Combine the plot of the Fermi lines on the $k_y$ plane with a
    heat-map plot of (minus) the Berry curvature*

    The plot of the Fermi lines with a colour-map of
    $-\Omega_z(k_x,0,k_z)$ is shown in
    Fig. [59](#fig18.2){reference-type="ref" reference="fig18.2"}.

![(Colour online) Calculated total Berry curvature
$-\Omega_z(\ifmmode  \mathbf{k}  \else$ $~\fi)$ in the plane $k_y=0$
(note log scale). Intersections of the Fermi surface with this plane are
shown.](figure/example18/Fe_Fermi_surface+Berry_phase.png){#fig18.3
width="0.5\\columnwidth"}

## Anomalous Hall conductivity {#anomalous-hall-conductivity .unnumbered}

-   *AHC converges rather slowly with k-point sampling, and a
    $25 \times 25 \times 25$ does not yield a well-converged value.
    Compare the converged AHC value with those obtained in Refs.  and .*

    The *x,y,z*-components of the AHC for a $25\times25\times25$ BZ mesh
    are shown in the snippet below. The converged result reported in
    Refs.  and for the *z*-component is 756.76
    ($(\Omega \mathrm{cm})^{-1}$). Hence, a $25\times25\times25$ BZ mesh
    clearly gives a very inaccurate value ($\sim 36.4\%$ error). Even
    with adaptive refinement the error is still very large
    ($\sim 31.7\%$). It is worth to note that the adaptive refinement
    slightly breaks the symmetry and gives non-zero values for the
    *x*-component and *y*-component, although these are opposite in
    sign.

    ::: tcolorbox
         Properties calculated in module  b e r r y
         ------------------------------------------

           * Anomalous Hall conductivity

         Interpolation grid: 25 25 25

         Fermi energy (ev):   12.6283

         AHC (S/cm)       x          y          z
         ==========    -0.0000     0.0000   554.6437


         Total Execution Time          59.112 (sec)
    :::

    ::: tcolorbox
         Properties calculated in module  b e r r y
         ------------------------------------------

           * Anomalous Hall conductivity

         Regular interpolation grid: 25 25 25
           Adaptive refinement grid: 5 5 5
               Refinement threshold: Berry curvature >100.00 bohr^2
          Points triggering refinement:   42( 0.27%)

         Fermi energy (ev):   12.6283

         AHC (S/cm)       x          y          z
         ==========     2.4602    -2.4602   574.2950
    :::

    Since these are quite demanding calculations, we only report the
    value of the AHC for a $125\times125\times125$ BZ mesh with a
    $5\times5\times5$ adaptive refinement grid (see snippet below). The
    value for the *z*-component is 729.8276 $(\Omega \mathrm{cm})^{-1}$,
    which is in much closer agreement with the converged result from
    Refs.  and . Also, the magnitude of *x,y*-component is greatly
    reduced as expected.

    ::: tcolorbox
         Properties calculated in module  b e r r y
         ------------------------------------------

           * Anomalous Hall conductivity

         Regular interpolation grid: 125 125 125
           Adaptive refinement grid: 5 5 5
               Refinement threshold: Berry curvature >100.00 Ang^2
          Points triggering refinement: 1818( 0.09%)

         Fermi energy (ev):   12.6283

         AHC (S/cm)       x          y          z
         ==========    -0.2775     0.2775   729.8276
    :::

```{=html}
<!-- -->
```
-   *The Wannier-interpolation formula for the Berry curvature comprises
    three terms, denoted $J0$, $J1$, and $J2$ in Ref.  .*

    From Ref.  $$-2\;\mathrm{Im} G_{\alpha\beta} = J0 + J1 + J2,$$ where
    $$G_{\alpha\beta} = Tr[(\partial_\alpha \hat{P})\hat{Q}\hat{H}\hat{Q}(\partial_\beta\hat{P})]$$

    The three components $J0, J1$ and $J2$ for the $k$-point sampling of
    $125\times125\times125$ and a $5\times5\times5$ adaptive refinement
    grid are shown in the snippet below

::: tcolorbox
     J0 term :      0.0002    -0.0002     2.8479
     J1 term :      0.0004    -0.0004    18.4855
     J2 term :     -0.2782     0.2782   708.4942
     -------------------------------------------
:::

## Optical conductivity {#optical-conductivity .unnumbered}

-   *The optical conductivity tensor of bcc Fe with magnetization along
    $\hat{\mathbf{z}}$ has the form*
    $$\boldsymbol{\sigma} = \boldsymbol{\sigma}_\mathrm{S} + \boldsymbol{\sigma}_{\mathrm{A}} =
    \begin{pmatrix}
    \sigma_{xx} & 0                       & 0 \\
     0          & \sigma_{yy}=\sigma_{xx}  & 0 \\
     0          &  0                     & \sigma_{zz}
    \end{pmatrix} + \begin{pmatrix} 0 & \sigma_{xy} & 0 \\ -\sigma_{yx} & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}$$

-    *The DC AHC calculated earlier corresponds to $\sigma_{xy}$ in the
    limit $\omega \rightarrow 0$. At finite frequency
    $\sigma_{xy} = -\sigma_{yx}$ acquires an imaginary part which
    describes magnetic circular dichroism (MCD). Compute the complex
    optical conductivity for $\hbar\omega$ up to $7$ eV*

    The plot for the ac AHC is shown in
    Fig. [61](#fig18.4){reference-type="ref" reference="fig18.4"}.

![Plot of the real part of the complex optical conductivity with a
$50\times50\times50$ $k$-point mesh (black) and $125\times125\times125$
$k$-point mesh (red). The inset is a magnification of the region
\[0-0.1\] eV.](figure/example18/Fe-kubo_A_xy_125.pdf){#fig18.4
width="0.7\\columnwidth"}

-   *Compare the $\omega \rightarrow 0$ limit of $\sigma_{xy}$ with the
    result obtained earlier by integrating the Berry curvature.*

    The result obtained by integrating the Berry curvature is 729.83
    $(\Omega \mathrm{cm})^{-1}$ and the $\omega \rightarrow 0$ limit of
    the complex optical conductivity is $669.37$
    $(\Omega \mathrm{cm})^{-1}$.

    *Plot the MCD spectrum.*

    The plot of the magnetic circular dichroism is shown in
    Fig. [62](#fig18.5){reference-type="ref" reference="fig18.5"}.

![The magnetic circular dichroism from interpolation of the
Kubo-Greenwood
formula.](figure/example18/Fe_MCD_xy_125_sp3d2_projections.pdf){#fig18.5
width="0.7\\columnwidth"}

## Further ideas {#further-ideas-4 .unnumbered}

-   *Recompute the AHC and optical spectra of bcc Fe using projected s,
    p, and d-type Wannier functions instead of the hybridrised MLWFs
    (see Example 8), and compare the results.*

    First we have to modify the projection block in the input file
    `Fe.win` as did in Ex. [8](#sec8:Iron){reference-type="ref"
    reference="sec8:Iron"}

    ` `

    > begin projections Fe:s;p;d end projections

    Then we need to re-do points 3,4 and 6.

    Below there is the extract from the output file `Fe.wpout`. The
    result obtained from $s,p$ and $d$ projections for the $z$
    component, i.e. $\sigma_{xy}$, of the AHC is exactly the same as the
    one obtained from $sp_3d_2,d_{xy},d_{xz}$, and $d_{yz}$ projections.
    Plot of AHC and MCD are shown in
    Fig. [63](#fig18.6){reference-type="ref" reference="fig18.6"}.

    ::: tcolorbox
         Properties calculated in module  b e r r y
         ------------------------------------------

           * Anomalous Hall conductivity

         Regular interpolation grid: 25 25 25
           Adaptive refinement grid: 5 5 5
               Refinement threshold: Berry curvature >100.00 bohr^2
          Points triggering refinement:   42( 0.27%)

         Fermi energy (ev):   12.6283

         AHC (S/cm)       x          y          z
         ==========
         J0 term :      0.0006    -0.0006    -2.9033
         J1 term :      0.0032    -0.0032    10.0566
         J2 term :      2.4564    -2.4564   567.1417
         -------------------------------------------
         Total   :      2.4602    -2.4602   574.2950

        	
    :::

    <figure id="fig18.6">

    <figcaption>Left panel: Anomalous Hall conductivity. Right panel:
    Magnetic circular dichroism for <span
    class="math inline">ℏ<em>ω</em></span> up to <span
    class="math inline">7</span> eV, starting from <span
    class="math inline"><em>s</em>; <em>p</em>; <em>d</em></span> initial
    projections</figcaption>
    </figure>

# Iron---Orbital magnetization {#sec19:IronOM}

-   Outline: *Calculate the orbital magnetization of ferromagnetic bcc
    Fe by Wannier interpolation.*

```{=html}
<!-- -->
```
-   These are the same steps performed for
    Ex. [17](#sec17:IronSO){reference-type="ref"
    reference="sec17:IronSO"} and
    Ex. [18](#sec18:IronBerry){reference-type="ref"
    reference="sec18:IronBerry"}. Hence, they are not repeated here.

-   *The orbital magnetization is computed as the BZ integral of the
    quantity $\mathbf{M}^{\mathrm{orb}}(\ifmmode  \mathbf{k}  \else$
    $~\fi)$ defined in Eq. (12.20) of the User Guide.*

    Below we report Eq. (11.20) from the User Guide, and the total
    orbital magnetization as the integral of
    $\mathbf{M}^{\mathrm{orb}}(\ifmmode  \mathbf{k}  \else$ $~\fi)$ over
    the BZ $$\begin{aligned}
    \mathbf{M}^{\mathrm{orb}}(\ifmmode  \mathbf{k}  \else $ \mathbf{k} $~\fi) & = \sum_{n}\frac{1}{2}f_{n\ifmmode  \mathbf{k}  \else $ \mathbf{k} $~\fi}\;\mathrm{Im} \braket{\nabla_{\ifmmode  \mathbf{k}  \else $ \mathbf{k} $~\fi}\ifmmode u_{n\ifmmode  \mathbf{k}  \else $ \mathbf{k} $~\fi} \else $ u_{n\ifmmode  \mathbf{k}  \else $ \mathbf{k} $~\fi} $~\fi\vert \times (H_\ifmmode  \mathbf{k}  \else $ \mathbf{k} $~\fi+ \epsilon_\ifmmode  \mathbf{k}  \else $ \mathbf{k} $~\fi- 2\epsilon_{\mathrm{F}}) \vert \nabla_\ifmmode  \mathbf{k}  \else $ \mathbf{k} $~\fi\ifmmode u_{n\ifmmode  \mathbf{k}  \else $ \mathbf{k} $~\fi} \else $ u_{n\ifmmode  \mathbf{k}  \else $ \mathbf{k} $~\fi} $~\fi}
    \label{eq19.1} \\
    \mathbf{M}^{\mathrm{orb}}_{\mathrm{tot}} & = V\int\,\frac{\ifmmode \ifmmode \mathrm{d} \else $ \mathbf{d} $~\fi\ifmmode  \mathbf{k}  \else $ \mathbf{k} $~\fi\else $ \ifmmode \mathrm{d} \else $ \mathbf{d} $~\fi\ifmmode  \mathbf{k}  \else $ \mathbf{k} $~\fi$~\fi}{(2\pi)^3} \mathbf{M}^{\mathrm{orb}}(\ifmmode  \mathbf{k}  \else $ \mathbf{k} $~\fi)
    \label{eq19.2}
    \end{aligned}$$

    The two snippets below show the components of the total orbital
    magnetization computed according to
    Eq. ([\[eq19.2\]](#eq19.2){reference-type="ref"
    reference="eq19.2"}), and the spin magnetisation from the DFT
    calculation respectively

    ::: tcolorbox
         Properties calculated in module  b e r r y
         ------------------------------------------

           * Orbital magnetization

         Interpolation grid: 25 25 25


         Fermi energy (ev) =   12.628300


         M_orb (bohr magn/cell)        x          y          z
         ======================
            Local circulation :      0.0000    -0.0000     0.0935
         Itinerant circulation:      0.0000     0.0000    -0.0180
         --------------------------------------------------------
                      Total   :      0.0000    -0.0000     0.0755
    :::

    ::: tcolorbox
             total magnetization       =     0.00    -0.00    -2.22 Bohr mag/cell
             absolute magnetization    =     2.34 Bohr mag/cell
    :::

-   *Plot $\mathbf{M}^{\mathrm{orb}}(\ifmmode  \mathbf{k}  \else$
    $~\fi)$ along high-symmetry lines and compare the result with Fig. 2
    of Ref.  .*

    Before comparing the result of our calculation with the result in
    Fig. 2 of Ref.  , we need to fix a unit-conversion problem in the
    python script `Fe-bands+morb_z.py`. In fact, the units of
    $\mathbf{M}^{\mathrm{orb}}(\ifmmode  \mathbf{k}  \else$ $~\fi)$ are
    not Ry$\cdot\si{~\AA}^2$ as stated in the python script but
    eV$\cdot\si{~\AA}^2$ instead (as also stated in the User Guide).
    Moreover, in Ref. 
    $\mathbf{M}^{\mathrm{orb}}(\ifmmode  \mathbf{k}  \else$ $~\fi)$ is
    given in atomic units, i.e. Hartree$\cdot$bohr radii$^2$. In order
    to have a meaningful comparison we need to modify the python script
    accordingly. Open `Fe-bands+morb_z.py` and modify the following
    lines

    ` `

    > data = np.loadtxt('Fe-morb.dat')
    >
    > x=data\[:,0\]
    >
    > y=data\[:,3\]

    as

    ` `

    > data = np.loadtxt('Fe-morb.dat')
    >
    > x=data\[:,0\]
    >
    > y=data\[:,3\] \* 0.131234

    where $0.131234$ is the conversion factor from eV$\cdot\si{~\AA}^2$
    to a.u. We also need to modify the label for the y-axis from

    >     	pl.ylabel(r'$M^{\rm{orb}}_z(\mathbf{k})$  [ Ry$\cdot\AA^2$ ]')

    to

    >     	pl.ylabel(r'$M^{\rm{orb}}_z(\mathbf{k})$  [ a.u. ]')

    Now we can run the python script

    ` `

    > \$\> python Fe-bands+morb_z.py

    and look at the plot, here shown in
    Fig. [64](#fig19.1){reference-type="ref" reference="fig19.1"}. The
    difference between the quantities in the two plot is roughly the
    $-\frac{1}{2}$ factor due to the two different definitions of
    $\mathbf{M}^{\mathrm{orb}}$.

    ![Plot of $\mathbf{M}^{\mathrm{orb}}(\ifmmode  \mathbf{k}  \else$
    $~\fi)$ calculated by Wannier interpolation along the path
    $\Gamma$--H--P in the Brillouin
    zone.](figure/example19/Fe-morb_z.pdf){#fig19.1
    width="0.6\\columnwidth"}

    *Plot $\mathbf{M}^{\mathrm{orb}}(\ifmmode  \mathbf{k}  \else$
    $~\fi)$ together with the Fermi contours on the (010) BZ plane*

    ![Plot of $\mathbf{M}^{\mathrm{orb}}(\ifmmode  \mathbf{k}  \else$
    $~\fi)$ together with the Fermi contours on the (010) BZ
    plane](figure/example19/Fe-kslice-morb_z+fermi_lines.pdf){#fig19.3
    width="0.7\\columnwidth"}

# Disentanglement restricted inside spherical regions of *k*-space LaVO$_3$. {#sec20:LaVO3}

-   Outline: *Obtain disentangled MLWFs for strained $\mathrm{LaVO}_3$.*

<figure id="fig20">

<figcaption>Left: atomic structure of epitaxially-strained (tetragonal)
LaVO<span class="math inline"><sub>3</sub></span>. Right: atomic
structure of epitaxially-strained (tetragonal) SrMnO<span
class="math inline"><sub>3</sub></span>. Both structures have been
plotted with the <span class="smallcaps">XCrySDen</span>
program.</figcaption>
</figure>

-   These are the usual steps to generate MLWFs and are not reported
    here.

-   *Inspect the output file `LaVO3.wout`. In the initial summary, you
    will see that the disentanglement was performed only within one
    sphere of radius 0.2 around the point `A = (0.5, 0.5, 0.5)` in
    reciprocal space:*

    ::: tcolorbox
         *------------------------------- DISENTANGLE --------------------------------*
         |  Using band disentanglement                :                 T             |


         	...

         |  Number of spheres in k-space              :                 1             |
         |   center n.   1 :     0.500   0.500   0.500,    radius   =   0.200         |
    :::

-   *Compare the band structure that [wannier90]{.smallcaps} produced
    with the one obtained using Quantum ESPRESSO.*

    To obtain the band structure from the Quantum ESPRESSO calculation
    we can use the `bands.x` program available at
    <http://www.tcm.phy.cam.ac.uk/~jry20/bands.html>, see mini-tutorial
    at the end of Ex. [6](#sec6:copper){reference-type="ref"
    reference="sec6:copper"}. Here, we only report the `.inp` file used
    to generate the $k$-point mesh for the non-scf calculation

    ::: tcolorbox
        7.03 0.00 0.00
        0.00 7.03 0.00
        0.00 0.00 7.6627

        30

        G       0.00000  0.00000  0.00000  M       0.50000  0.50000  0.00000
        M       0.50000  0.50000  0.00000  X       0.50000  0.00000  0.00000
        X       0.50000  0.00000  0.00000  G       0.00000  0.00000  0.00000
        G       0.00000  0.00000  0.00000  Z       0.00000  0.00000  0.50000
        Z       0.00000  0.00000  0.50000  A       0.50000  0.50000  0.50000
        A       0.50000  0.50000  0.50000  R       0.50000  0.00000  0.50000
        R       0.50000  0.00000  0.50000  X       0.50000  0.00000  0.00000
    :::

    Remember to add the following line to the `.bands` file in order to
    show the eigenvalues at each k-point.

    ` `

    > verbosity = 'high'

    Plot of the interpolated band structure is shown in
    Fig. ([67](#fig20.1){reference-type="ref" reference="fig20.1"}). In
    the top panel, the full band structure is shown. In the bottom panel
    a magnification around the Fermi energy is shown (similar to Fig. 9
    in the Tutorial).

## Further ideas {#further-ideas-5 .unnumbered}

-   *Try to obtain the Wannier functions using the standard
    disentanglement procedure ...*

    Plots of the band structure of LaVO$_3$ with full disentanglement
    and no disentanglement are shown in
    Fig. ([68](#fig20.2){reference-type="ref" reference="fig20.2"}).
    These are plotted against the quantum ESPRESSO band structure (solid
    black lines) and the [wannier90]{.smallcaps}-interpolated one with
    disentanglement performed only within a sphere centred in A (red
    dots). We see that the other two methods diverge from the DFT
    calculation in region of $k$-space where the bands of interest are
    not entangled with other unwanted bands. For example, in the zone
    between $\Gamma$ and M and Z and A the interpolated bands with full
    disentanglement and no disentanglement diverge substantially from
    the DFT calculation.

-   *In order to illustrate all possible cases, it is instructive to
    apply this method to SrMnO$_3$ ...*

    Plots of the interpolated bands for the different cases are shown in
    Fig. ([69](#fig20.4){reference-type="ref" reference="fig20.4"}). In
    this case, the disentanglement for all the Mn-3d-derived states
    (empty red circles in Fig. ([69](#fig20.4){reference-type="ref"
    reference="fig20.4"})) is only necessary around the $\Gamma$ point,
    as for all the other points and lines the bands of interest are well
    separated from other bands lower in energy. However, if we only
    consider the $e_g$ states (solid blue circles in
    Fig. ([69](#fig20.4){reference-type="ref" reference="fig20.4"}))
    then the situation is different as these states are entangled with
    the $t_{2g}$ states around X. Of course the $t_{2g}$ states (solid
    green cones in Fig. ([69](#fig20.4){reference-type="ref"
    reference="fig20.4"})) are entangled with $e_{g}$ states around X
    and with lower-lying states at $\Gamma$.

<figure id="fig20.1">
<p><br />
</p>
<figcaption>Top panel: full band structure of epitaxially-strained
(tetragonal) LaVO<span class="math inline"><sub>3</sub></span> along the
<span class="math inline"><em>Γ</em></span>-M-X-<span
class="math inline"><em>Γ</em></span>-Z-A-R-X from DFT calculation
(solid black) and interpolation from <span
class="smallcaps">wannier90</span> (red dots). Bottom panel:
magnification around Fermi energy <span
class="math inline">16.6049</span> (dashed line). The disentanglement
was performed only for <span
class="math inline"><em>k</em></span>-points within a sphere of radius
0.2 <span class="math inline">$\si{~\AA}^{-1}$</span> centred in
A.</figcaption>
</figure>

![Comparison of interpolated band structure of epitaxially-strained
(tetragonal) LaVO$_3$ with disentanglement on a sphere of radius 0.2
$\si{~\AA}^{-1}$ centred in A (red dots), full disentanglement (blue
dots) and no disentanglement (green dots). Fermi energy is shown with a
dashed line.](figure/example20/LaVO3_bandstructure_all.pdf){#fig20.2
width="0.7\\columnwidth"}

![[wannier90]{.smallcaps}-interpolated bands of SrMnO$_3$. From only
$t_{2g}$ states (solid green cones), from only $e_g$ states (solid blue
circles), or all Mn-3d-derived states ($t_{2g} + e_g$) (empty red
circles).](figure/example20/SrMnO3_allbands.pdf){#fig20.4
width="0.7\\columnwidth"}

# Gallium Arsenide---Symmetry-adapted Wannier functions {#sec21:GaAsSA}

-   Outline: *Obtain symmetry-adapted Wannier functions out of four
    valence bands of GaAs. For the theoretical background of the
    symmetry-adapted Wannier functions, see R. Sakuma, Phys. Rev. B
    **87**, 235109 (2013).*

![Unit cell of GaAs crystal plotted with the [XCrySDen]{.smallcaps}
program.](figure/example01/GaAs.png){#fig21 width="0.25\\columnwidth"}

-   These are common to all calculations, and they have already been
    performed in previous examples. Hence, no results are shown here.

    The space group of GaAs is $F{-}43m$ (sequential number 276 in the
    International Tables for Crystallography, Vol. A). In our example
    the Ga atom is placed at the origin, whose Wyckoff letter is $a$ and
    its multiplicity is $4$. The site symmetry group of $a$ is ${-}43m$,
    which is isomorphous to the full point group of the crystal, also
    known as $T_d^2$. This is due to the fact that $F{-}43m$ is
    symmorphic. Hence, $a$ contains 24 symmetry operations (see
    Tab. [\[tab21.1\]](#tab21.1){reference-type="ref"
    reference="tab21.1"}). The As atom is placed at (0.25,0.25,0.25) in
    fractional coordinates, whose Wyckoff letter is $c$ and its
    multiplicity is 4. It also contains 24 symmetry operations (see
    Tab. [13](#tab21.2){reference-type="ref" reference="tab21.2"}).

    ::: {#tab21.2}
      x,y,z     -x, -y, z   -x,y,-z   x,-y,-z   z,x,y     z,-x,-y
      --------- ----------- --------- --------- --------- ---------
      -z,-x,y   -z,x,-y     y,z,x     -y,z,-x   y,-z,-x   -y,-z,x
      y,x,z     -y,-x,z     y,-x,-z   -y,x,-z   x,z,y     -x,z,-y
      -x,-z,y   x,-z,-y     z,y,x     z,-y,-x   -z,y,-x   -z,-y,x

      : 24 symmetry operations for the Wyckoff position $4c$ in $-43m$
      [@bilbaocrystserver].
    :::

    []{#tab21.1 label="tab21.1"}

    ::: {#tab21.2}
      x,y,z   -x+1/2,-y+1/2, z   -x+1/2,y,-z+1/2   x,-y+1/2,-z+1/2
      ------- ------------------ ----------------- -----------------
      z,x,y   z,-x+1/2,-y+1/2    -z+1/2,-x+1/2,y   -z+1/2,x,-y+1/2
      y,z,x   -y+1/2,z,-x+1/2    y,-z+1/2,-x+1/2   -y+1/2,-z+1/2,x
      y,x,z   -y+1/2,-x+1/2,z    y,-x+1/2,-z+1/2   -y+1/2,x,-z+1/2
      x,z,y   -x+1/2,z,-y+1/2    -x+1/2,-z+1/2,y   x,-z+1/2,-y+1/2
      z,y,x   z,-y+1/2,-x+1/2    -z+1/2,y,-x+1/2   -z+1/2,-y+1/2,x

      : 24 symmetry operations for the Wyckoff position $4c$ in $-43m$
      [@bilbaocrystserver].
    :::

    []{#tab21.2 label="tab21.2"}

    The list of site-symmetry operations may be found in the `.sym` file
    and in the output file `pw2wan.out`. In the latter, the list is in
    the section relative to the computation of the $D_{mn}$ matrix (see
    Ref.  ).

## One $s$-like Wannier function centred at Ga {#one-s-like-wannier-function-centred-at-ga .unnumbered}

-   *Compute the symmetry-adapted MLWF.*

    The ${-}43m$ site-symmetry group is isomorphous to $T_d^2$. From the
    table of characters of $T_d^2$ we find 5 irreducible representations
    (*irrep*). The irrep with character $A_1$ is a one-dimensional
    representation, whose eigenfunction is spherically symmetric. Hence,
    a single $s$-like orbital in (0,0,0) may be used. However, this is
    not enough as the choice of the initial guess must also be
    compatible with the symmetry of the bands. In fact, if we tried to
    wannierise only the lowest band, excluding all the other bands (this
    can be done by changing the input file as
    `num_wann = 1, num_bands = 1` and `exclude_bands = 1-5, 7-19`), the
    resulting $1\times1$ $U(\mathbf{k})$ could not fulfill Eq. 19 in
    Ref.  . Similarly, if we tried to wannierise only the three top
    bands.

    ` `

    > begin projections
    >
    > f= 0.0, 0.0, 0.0 : s
    >
    > end projections

::: tcolorbox
      ----------------
      *** Compute DMN
      ----------------

      Number of symmetry operators =    24
    	
:::

<figure id="fig21.1">

<figcaption>One <span class="math inline"><em>s</em></span>-like
symmetry-adapted Wannier function centred on the Gallium atom in
GaAs.</figcaption>
</figure>

## Three $p$-like Wannier functions centred at Ga {#three-p-like-wannier-functions-centred-at-ga .unnumbered}

-   *Compute the symmetry-adapted MLWFs.*

    Another representation of ${-}43m$, namely $T_2$, has dimension 3.
    Its eigenfunctions are linear functions proportional to $x,y,z$.
    Hence, we can use three $p$-like orbitals ($p_x,p_y,p_z$) centred at
    (0,0,0).

    ` `

    > begin projections
    >
    > f= 0.0, 0.0, 0.0 : p
    >
    > end projections

<figure id="fig21.2">

<figcaption>Three <span class="math inline"><em>p</em></span>-like
symmetry-adapted Wannier functions centred on the Gallium atom in
GaAs.</figcaption>
</figure>

## One $s$-like and three $p$-like Wannier functions centred at Ga {#one-s-like-and-three-p-like-wannier-functions-centred-at-ga .unnumbered}

-   *Compute the symmetry-adapted MLWFs.*

    We can construct also construct a representation of dimension
    $4=3+1$ for the 4 valence bands by specifying 1 $s$-like orbital and
    3 $p$-like orbitals on Ga, which corresponds to the irreducible
    representations $A_1$ and $T_2$ respectively. However, it would not
    be possible to

    ` `

    > begin projections
    >
    > f= 0.0, 0.0, 0.0 : s
    >
    > f= 0.0, 0.0, 0.0 : p
    >
    > end projections

<figure id="fig21.3">

<figcaption>One <span class="math inline"><em>s</em></span>-like and
three <span class="math inline"><em>p</em></span>-like Wannier functions
centred on the Gallium atom in GaAs.</figcaption>
</figure>

## One $s$-like and three $p$-like Wannier functions centred at As {#one-s-like-and-three-p-like-wannier-functions-centred-at-as .unnumbered}

The site-symmetry group for the As anion centred at $(0.25,0.25,0.25)$
is ${-}43m$ as well and we can perform the same analysis done for the Ga
cation. Contrary to the Ga case, for the As anion it is possible to
wannierise the bottom band from one $s$-like orbital centred at
$(0.25,0.25,0.25)$ and the top three bands from three $p$-like orbitals
centred at $(0.25,0.25,0.25)$ (see
Fig. [75](#fig21.5){reference-type="ref" reference="fig21.5"}).

-   *Compute the symmetry-adapted MLWFs.*

    ` `

    > begin projections
    >
    > f=0.25,0.25,0.25 : s
    >
    > f=0.25,0.25,0.25 : p
    >
    > end projections

<figure id="fig21.4">

<figcaption>One <span class="math inline"><em>s</em></span>-like and
three <span class="math inline"><em>p</em></span>-like Wannier functions
centred on the Arsenic atom in GaAs.</figcaption>
</figure>

<figure id="fig21.5">
<p><br />
</p>
<figcaption>Interpolated <span class="smallcaps">wannier90</span> bands
of GaAs starting from a) 1 <span
class="math inline"><em>s</em></span>-like centred on the Arsenic anion
and b) three <span class="math inline"><em>p</em></span>-like orbitals
centred on the Arsenic anion, respectively.</figcaption>
</figure>

## Four $s$-like Wannier functions centred on the four Ga-As bonds {#four-s-like-wannier-functions-centred-on-the-four-ga-as-bonds .unnumbered}

From a group-theoretical point of view, the case of four $s$-like
functions centred on four covalent bonds, correspond to the *irrep*
$A_{1g}$ of the site-symmetry group $.3m$ of the Wyckoff position $e$.
There are 6 symmetry operations for each equivalent position (0.125,
0.125, 0.125), (0.125, 0.125, -.375), (-.375, 0.125, 0.125) and (0.125,
-.375, 0.125). The combined 24 symmetry operations turn out to be
exactly that of the full ${-}43m$ group.

-   *Compute the symmetry-adapted MLWFs.*

    ` `

    > begin projections
    >
    > f= 0.125, 0.125, 0.125: s
    >
    > f= 0.125, 0.125, -.375: s
    >
    > f= -.375, 0.125, 0.125: s
    >
    > f= 0.125, -.375, 0.125: s
    >
    > end projections

<figure id="fig21.6">

<figcaption>Four <span
class="math inline"><em>s</em><em>p</em><sub>3</sub></span>-like
symmetry-adapted Wannier functions centred on the Ga-As bonds in
GaAs.</figcaption>
</figure>

# Copper---Symmetry-adapted Wannier functions {#sec22:CopperSA}

-   Outline: *Obtain symmetry-adapted Wannier functions for Cu. By
    symmetry-adapted mode, for example, we can make atomic centered
    s-like Wannier function, which is not possible in the usual
    procedure to create maximally localized Wannier functions. For the
    theoretical background of the symmetry-adapted Wannier functions,
    see R. Sakuma, Phys. Rev. **B** 87, 235109 (2013).*

![Unit cell of Copper
crystal.](figure/example04/copper_crystal.png){#fig22.0
width="0.25\\columnwidth"}

*Each directory creates $s$-like symmetry-adapted Wannier function
centered at different position on top of atomic centered $d$-like
Wannier functions.*

Below it is reported the README file from the example directory

::: tcolorbox
    # Symmetry-adapted mode

    Additional input in Cu.win file
      site_symmetry = .true.   (default value is .false.)
      symmetrize_eps = 1d-9    (default value is 1d-3   )


    Additional input in Cu.pw2wan file
      write_dmn = .true.


    Working directories
      s_at_0.00 : s-like Wannier function centered at (0,0,0)       + atomic-centered d-like WFs
      s_at_0.25 : s-like Wannier function centered at (1/4,1/4,1/4) + atomic-centered d-like WFs
      s_at_0.50 : s-like Wannier function centered at (1/2,1/2,1/2) + atomic-centered d-like WFs



    In s_at_0.25, we use an additional flag "read_sym = .true." to customize the symmetry operations
    to be used.
    We exclude the inversion symmetry to create s-like Wannier function centered at (1/4,1/4,1/4).
    Information on symmetry operations without inversion symmetry is taken from GaAs calculation.
    See more detailed discussion in R. Sakuma, Phys. Rev. B 87, 235109 (2013).
:::

The space group of Cu is $Fm{-}3m$ (sequential number 225 in the
International Tables for Crystallography, Vol. A).

## $s$-like Wannier function centred at the origin {#s-like-wannier-function-centred-at-the-origin .unnumbered}

-   *Compute the symmetry-adapted MLWFs.*

    In this example both the $s$-orbital and $d$-orbitals is placed at
    the origin, whose Wyckoff letter is $a$ and its multiplicity is $4$.
    The site symmetry group of $a$ is $m{-}3m$, which is isomorphous to
    the full point group of the crystal, known as $O_h$ and it contains
    48 symmetry operations. The six MLWFs obtained by placing both the
    initial $s$-orbital and the $d$-orbitals on the Cu atom are shown in
    Fig. [78](#fig22.1){reference-type="ref" reference="fig22.1"}.

<figure id="fig22.1">
<p><br />
</p>
<figcaption>Six symmetry-adapted MLWFs in Cu. The initial <span
class="math inline"><em>s</em></span>-orbital is placed at the
origin.</figcaption>
</figure>

## $s$-like Wannier function centred at (0.25,0.25,0.25) {#s-like-wannier-function-centred-at-0.250.250.25 .unnumbered}

-   *Compute the symmetry-adapted MLWFs.* In this example the
    $s$-orbital is placed at (0.25,0.25,0.25), whose Wyckoff letter is
    $c$ and its multiplicity is $8$. The site symmetry group of $c$ is
    ${-}43m$, which is not isomorphous to the full point group of the
    crystal ($O_h$). This site symmetry group contains only 24 symmetry
    operations, i.e. it comes from $O_h$ when inversion is removed. This
    is the reason why the flag `read_sym =.true. ` in the `.pw2wan` file
    and an additional input is required, namely `Cu.sym`. In fact, this
    file is not automatically generated by `pw2wannier90.x` but is given
    as input to be read. The six MLWFs obtained by placing the initial
    $s$-orbital at (0.25,0.25,0.25) and the $d$-orbitals on the Cu atom
    are shown in Fig. [79](#fig22.2){reference-type="ref"
    reference="fig22.2"}.

<figure id="fig22.2">
<p><br />
</p>
<figcaption>Six symmetry-adapted MLWFs in Cu. The initial <span
class="math inline"><em>s</em></span>-orbital is placed at
(0.25,0.25,0/25).</figcaption>
</figure>

## $s$-like Wannier function centred at (0.5,0.5,0.5) {#s-like-wannier-function-centred-at-0.50.50.5 .unnumbered}

-   *Compute the symmetry-adapted MLWFs.* In this example the $s$
    orbital is placed at (0.5,0.5,0.5), whose Wyckoff letter is $b$ and
    its multiplicity is $4$. The site symmetry group of $c$ is $m{-}3m$,
    which is isomorphous to the full point group of the crystal ($O_h$).
    Hence, no additional input file is needed in this case. The six
    MLWFs obtained by placing the initial $s$-orbital at (0.5,0.5,0.5)
    and the $d$-orbitals on the Cu atom are shown in
    Fig. [80](#fig22.3){reference-type="ref" reference="fig22.3"}.

<figure id="fig22.3">
<p><br />
</p>
<figcaption>Six symmetry-adapted MLWFs in Cu. The initial <span
class="math inline"><em>s</em></span>-orbital is placed at
(0.5,0.5,0.5).</figcaption>
</figure>

# Platinum---Spin Hall conductivity {#sec29:PtSHC}

-   Outline: *Calculate spin Hall conductivity (SHC) and plot Berry
    curvature-like term of fcc Pt considering spin-orbit coupling. To
    gain a better understanding of this example, it is suggested to read
    Ref.  for a detailed description of the theory and Ch. 12.5 of the
    User Guide.*

```{=html}
<!-- -->
```
-   *Compute the MLWFs, spin Hall conductivity and `kpath`, `kslice`
    plots.*

## Spin Hall conductivity {#spin-hall-conductivity .unnumbered}

-   *SHC converges rather slowly with k-point sampling, and a
    $25 \times 25 \times 25$ kmesh does not yield a well-converged
    value. To get a converged SHC value, increase the density of kmesh
    and then compare the converged result with those obtained in Refs. 
    and .*

    The file `Pt-shc-fermiscan.dat` contains the calculated SHC. The SHC
    for a $25\times25\times25$ kmesh are shown in the snippet below.

    ::: tcolorbox
        #No.   Fermi energy(eV)   SHC((hbar/e)*S/cm)
           1     6.000000    0.00000000E+00
        ...
         120    17.900000    0.17230482E+04
         121    18.000000    0.17054542E+04
        ...
         201    26.000000    0.22665760E+03
    :::

    The calculated Fermi energy obtained from `Quantum ESPRESSO` is
    $17.9919$ eV. It may vary among different calculations due to the
    differences between versions of `Quantum ESPRESSO` or compilers, and
    these may lead to deviations from the following results. However,
    the difference should be acceptable and the calculated SHC should be
    essentially the same.

    The SHC at the Fermi energy is 1705 $(\hbar/e)\mathrm{S/cm}$. The
    converged results reported in Refs.  and are around 2200
    $(\hbar/e)\mathrm{S/cm}$. Hence, a $25\times25\times25$ kmesh
    clearly gives an inaccurate value ($\sim 22.5\%$ error).

    Since these are quite demanding calculations, we only report the
    value of the SHC for a $100\times100\times100$ kmesh (see snippet
    below). The value for the SHC at Fermi energy is 2207
    $(\hbar/e)\mathrm{S/cm}$, which is in much closer agreement with the
    converged result from Refs.  and .

    ::: tcolorbox
        #No.   Fermi energy(eV)   SHC((hbar/e)*S/cm)
           1     6.000000    0.00000000E+00
        ...
         120    17.900000    0.21899191E+04
         121    18.000000    0.22066678E+04
        ...
         201    26.000000    0.24919920E+03
    :::

-   To complete the previous discussions, we also compare the Fermi
    energy scan plots of the two calculations as shown in the
    Fig. [81](#fig29.3){reference-type="ref" reference="fig29.3"}.

    ![Fermi energy scan plots for calculations with $25\times25\times25$
    kmesh and $100\times100\times100$
    kmesh.](figure/example29/pt_shc_kmesh.pdf){#fig29.3
    width=".8\\columnwidth"}

-   The `seedname.wpout` will print the percentage of $k$-points which
    have been calculated at the moment, as well as the corresponding
    calculation time, as shown in the following snippet.

    ::: tcolorbox
         Properties calculated in module  b e r r y
         ------------------------------------------

           * Spin Hall Conductivity

             Fermi energy scan

         Calculation started
         -------------------------------
           k-points       wall      diff
          calculated      time      time
          ----------      ----      ----
               0%          0.0       0.0
              10%         22.7      22.7
              20%         36.5      13.8
              30%         50.4      14.0
              40%         64.4      14.0
              50%         78.4      14.0
              60%         92.5      14.1
              70%        106.5      14.0
              80%        120.4      13.9
              90%        134.2      13.8
             100%        147.9      13.7


         Interpolation grid: 25 25 25

         Using adaptive smearing
               adptive smearing prefactor    1.414
               adptive smearing max width    1.000 eV
    :::

    This might be helpful as you can roughly estimate the total
    computational time of your calculation, or it might give credence to
    the code that it is actually functioning :). Note this report is
    merely based on the "root" computation node. It is accurate if the
    `postw90` is run in serial, or the load on each node is balanced if
    running in parallel. However, the estimation is rough if loads are
    not balanced among nodes. This may happen if the performance of
    nodes in your cluster are not identical, or adaptive kmesh
    refinements are triggered so some nodes may compute much more
    $k$-points than others. Besides, if you are careful enough, you may
    find the diff time of 10% is much larger than later ones. This is
    caused by some done-once-and-for-all computations carried out at the
    beginning, thus later computations are much faster.

## Berry curvature-like term plots {#berry-curvature-like-term-plots .unnumbered}

-   *The band-projected Berry curvature-like term
    $\Omega_{n,\alpha\beta}^{\text{spin} \gamma}({\bm k})$ is defined as
    Eq. (12.22) in the User Guide.* *Plot the band structure of Pt and
    color it by the magnitude of its band-projected Berry curvature-like
    term $\Omega_{n,xy}^{\text{spin}z}(\bm k)$, and plot the k-resolved
    Berry curvature-like term $\Omega_{xy}^{\text{spin}z}(\bm k)$ along
    the same path in the BZ.*

    With Fermi energy set as 17.9919 eV we obtain the energy bands
    colored by the
    $\Omega_{n,\alpha\beta}^{\text{spin} \gamma}({\bm k})$ and the
    $k$-resolved Berry curvature-like term
    $\Omega_{xy}^{\text{spin}z}(\bm k)$ along high-symmetry lines as
    shown in Fig. [82](#fig29.1){reference-type="ref"
    reference="fig29.1"}, which contains two plots calculated with
    different fixed smearing width.

<figure id="fig29.1">

<figcaption>Top panels: Band structure of Pt along symmetry lines
W-L-<span class="math inline"><em>Γ</em></span>-X-W-<span
class="math inline"><em>Γ</em></span>, colored by the <span
class="math inline"><em>Ω</em><sub><em>n</em>, <em>x</em><em>y</em></sub><sup>spin<em>z</em></sup>(<strong>k</strong>)</span>.
Bottom panels: <span class="math inline"><em>k</em></span>-resolved
Berry curvature-like term <span
class="math inline"><em>Ω</em><sub><em>x</em><em>y</em></sub><sup>spin<em>z</em></sup>(<strong>k</strong>)</span>
along the symmetry lines.</figcaption>
</figure>

-   *Combine the plot of the Fermi lines on the $(k_x,k_y)$ plane with a
    heatmap plot of the Berry curvature-like term of spin Hall
    conductivity.*

    The plots of the Fermi lines with a heatmap of
    $\Omega_{xy}^{\text{spin}z}(k_x,k_y,0)$ are shown in
    Fig. [83](#fig29.2){reference-type="ref" reference="fig29.2"}.

<figure id="fig29.2">

<figcaption>Calculated <span
class="math inline"><em>k</em></span>-resolved Berry curvature-like term
<span
class="math inline"><em>Ω</em><sub><em>x</em><em>y</em></sub><sup>spin<em>z</em></sup>(<strong>k</strong>)</span>
in the plane <span
class="math inline"><em>k</em><sub><em>z</em></sub> = 0</span> (note the
magnitude of <span
class="math inline"><em>Ω</em><sub><em>x</em><em>y</em></sub><sup>spin<em>z</em></sup>(<strong>k</strong>)</span>
is in log scale). Intersections of the Fermi surface with this plane are
shown as black lines.</figcaption>
</figure>

# Gallium Arsenide---Frequency-dependent spin Hall conductivity {#sec30:GaAsSHC}

-   Outline: *Calculate the alternating current (ac) spin Hall
    conductivity of gallium arsenide considering spin-orbit coupling. To
    gain a better understanding of this example, it is suggested to read
    Ref.  for a detailed description of the theory and Ch. 12.5 of the
    User Guide.*

```{=html}
<!-- -->
```
-   *Compute the MLWFs and compute the ac spin Hall conductivity.*

## ac spin Hall conductivity {#ac-spin-hall-conductivity .unnumbered}

-   *The ac SHC of GaAs converges rather slowly with $k$-point sampling,
    and a $100 \times 100 \times 100$ kmesh does not yield a
    well-converged value. To get a converged SHC value, increase the
    density of kmesh and then compare the converged result with those
    obtained in Refs.  .*

    The file `GaAs-shc-freqscan.dat` contains the calculated ac SHC. The
    snippet below shows a calculated result with $100\times100\times100$
    kmesh, a fixed smearing width of 0.05 eV and no scissors shift
    applied.

    ::: tcolorbox
        #No.   Frequency(eV)   Re(sigma)((hbar/e)*S/cm)   Im(sigma)((hbar/e)*S/cm)
           1     0.000000    -0.68114601E+00    0.00000000E+00
        ...
         801     8.000000    -0.39471936E+01   -0.29928198E+02
    :::

    The ac SHC is plotted as Fig. [84](#fig30.1){reference-type="ref"
    reference="fig30.1"}.

    ![Frequency scan plot for GaAs ac SHC, using a low kmesh of
    $100\times100\times100$.](figure/example30/gaas_freqscan_100kpt.pdf){#fig30.1
    width=".8\\columnwidth"}

-   If further increasing the density of kmesh to
    $250\times250\times250$, and using the adaptive smearing, a nice
    converged plot could be produced as
    Fig. [85](#fig30.2){reference-type="ref" reference="fig30.2"}. Note
    that by using keywords a scissors shift of 1.117 eV is applied.
    Fig. [85](#fig30.2){reference-type="ref" reference="fig30.2"} can be
    viewed as Fig. [84](#fig30.1){reference-type="ref"
    reference="fig30.1"} translated by $\sim1$ eV along the horizontal
    axis.

    ![Frequency scan plots for GaAs ac SHC, using a dense kmesh of
    $250\times250\times250$. Two kinds of smearing are
    compared.](figure/example30/gaas_freqscan.pdf){#fig30.2
    width="0.8\\columnwidth"}

[^1]: To install [wannier90]{.smallcaps} you can follow the instructions
    in the [readme]{.smallcaps} file of the [wannier90]{.smallcaps}
    distribution. For an introduction to the theory, you can look at the
    [wannier90]{.smallcaps} User guide
    <http://www.wannier.org/user_guide.html>, the
    [wannier90]{.smallcaps} Tutorial and references therein.

[^2]: This can be easily achieved with any code, e.g. Python, MATLAB or
    even bash.
