# Tutorial solutions

## Preliminaries

Welcome to `wannier90`! This is the tutorial solutions booklet for the examples in the `wannier90` v3.1.0 tutorial (<http://www.wannier.org/doc/tutorial.pdf>). Info on the installation process and the theory of Maximally Localized Wannier Functions (MLWFs) is not reported here as they can be found elsewhere.[^1] The solutions in this booklet are for v3.1.0 only! The following (open-source) programs are required to reproduce the plots and figures in this booklet:

- `gnuplot` is used to plot bandstructures. It is available for many operating systems and is often installed by default on Unix/Linux distributions. In particular, we used gnuplot 4.6 patchlevel 6.
    <br><http://www.gnuplot.info>

- `Grace` is another plotting tool to visualise bandstructures.
    <br><http://plasma-gate.weizmann.ac.il/Grace/>

- `Vesta` is the default 3D visualisation program [@vesta] adopted in this booklet. It is used to visualise crystal structures, volumetric data (such as WFs and densities). Download is available for several OS here: <http://jp-minerals.org/vesta/en/>.

- `XCrySDen` is also used to visualise crystal structures and Fermi surfaces in particular. It is available for Unix/Linux, Windows (using cygwin), and OSX. To correctly display files from wannier90, version 1.4 or later must be used.
    <br><http://www.xcrysden.org>

- `VMD` may also be used to visualise crystal structures and 3D-fields. It can also read a great variety of input formats and it comes with handy postprocessing tools. <http://www.ks.uiuc.edu/Research/vmd>

**Disclaimer:** All the band structure interpolations have been carried out with `ws_distance = .false.`, which is the default value for the version 2.1. However, in the new `wannier90` release, corresponding to version 3.0, the default value of `ws_distance` has been changed to `.true.`, as, to the best of our knowledge, the Wigner-Seitz interpolation scheme never lowers the quality of the interpolation and it is often superior to the default scheme.

## About these tutorial solutions

This solution manual consists of 24 sections, each containing the solutions, in the form of plots, tables, and text, to the corresponding example in the `wannier90` v3.1.0 tutorial! For each example, only the outline and key questions from the tutorial are reported here. All of the `wannier90` input files have been provided. From example 5 onwards, input files for the pwscf interface (<http://www.quantum-espresso.org>) to `wannier90` have also been provided. You will need a recent working version of the Quantum ESPRESSO package (`v6.2` and above), to run these examples. In particular, you will need `pw.x` and `pw2wannier90.x`, as explained in the `wannier90` v3.1.0 tutorial. Please visit <http://www.quantum-espresso.org> to download the package and follow the instruction on the website for installation. Further details on how to run the calculations for each example may be found in the corresponding section of the `wannier90` v3.1.0 tutorial. There are interfaces to a number of other electronic structure codes including: ABINIT (<http://www.abinit.org>), fleur (<http://www.flapw.de>), OpenMX (<http://www.openmx-square.org/>), GPAW (<https://wiki.fysik.dtu.dk/gpaw/>), VASP (<http://www.vasp.at>), and Wien2K (<http://www.wien2k.at>).

All the tests were performed on an x86_64 octa-core Intel(R) Xeon(R) CPU E5620 @2.40GHz. Two packages `intel-suite/2015.3.187` and `mkl/2015.3.187` were used for the compilation of `wannier90`. We expect some of the numerical results to depend on the architecture, the compiler distribution, e.g. `gcc` vs `gfortran`, the version of the compiler and the libraries. However, general trends should not be affected by these parameters and, in principle, it should be safe to ignore the differences between different set-ups. If you find that your results are significantly different from the one given in this booklet, please report it to the `wannier90` developers team by opening an issue on the [GitHub repository](https://github.com/wannier-developers/wannier90) or by writing an email to the forum at `wannier@quantum-espresso.org` (we strongly recommend the first option). Moreover, if you know how to solve the issue and you have a fix for it, you can open a pull-request on the GitHub repo.

## Contact us

If you have any suggestions on how this solution manual may be improved and for any other issue, open an issue on the official repository of the `wannier90` code on [GitHub](https://github.com/wannier-developers/wannier90) or send an email on the forum at `wannier@quantum-espresso.org` (we strongly recommend the former). For the forum note that you will need to be registered. Emails from non-registered users will be deleted automatically. You can register by following the links at <http://www.wannier.org/forum.html>.

[^1]: To install `wannier90` you can follow the instructions in the `readme` file of the `wannier90` distribution. For an introduction to the theory, you can look at the `wannier90` [User guide](../user_guide/introduction.md), the `wannier90` Tutorial and references therein.
