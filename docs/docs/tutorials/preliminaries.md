---
author:
- Version 3.1
bibliography:
- ../wannier90.bib
title: ": Tutorials"
---

# Preliminaries {#preliminaries .unnumbered}

Welcome to `wannier90`! The examples contained in these tutorials are
designed to help you become familiar with the procedure of generating,
analysing and using maximally-localised Wannier functions (MLWFs). As a
first step, install `wannier90` following the instructions in the
` README` file of the `wannier90` distribution. For an introduction to
the theory underlying MLWFs, you are encouraged to refer to the brief
overview given in the `wannier90` User Guide [@UserGuide], to the two
seminal papers of Refs. [@marzari-prb97; @souza-prb01], a recent review
article [@marzari-rmp12] and to a paper [@mostofi-cpc08] describing
`wannier90`.

The following additional programs may be installed in order to visualise
the output of `wannier90` (they are optional, not all of them are
necessary)

-   `gnuplot` is used to plot bandstructures. It is available for many
    operating systems and is often installed by default on Unix/Linux
    distributions <br>
    <http://www.gnuplot.info>

-   `xmgrace` may also be used to plot bandstructures.<br>
    <http://plasma-gate.weizmann.ac.il/Grace>

-   `XCrySDen` is used to visualise crystal structures, MLWFs, and Fermi
    surfaces. It is available for Unix/Linux, Windows (using cygwin),
    and OSX. To correctly display files from `wannier90`, version 1.4 or
    later must be used.<br>
    <http://www.xcrysden.org>

-   `vmd` can also be used to visualise crystal structures and MLWFs.<br>
    <http://www.ks.uiuc.edu/Research/vmd>

-   `python` with the `numpy` and `matplotlib` modules is used in
    tutorials 17 &#151; 19<br>
    <http://www.python.org><br>
    <http://www.numpy.org><br>
    <http://matplotlib.org>

# Parallel execution {#parallel-execution .unnumbered}

`postw90.x` and `wannier90.x` can
be run in parallel to speed up the calculations, using the MPI
libraries.

To enable the parallel version to be built, you must specify some flags
in the `make.inc` file of `wannier90` and `postw90`; for further
information, please refer to the `README.install` file in the top
directory of the `wannier90` distribution.

Then, to run e.g. with 8 processors, you typically need to run a command
similar to `postw90` as follows:

```bash title="Terminal"
mpirun -np 8 postw90.x seedname
```

(the `mpirun` command and its flags may differ depending on the MPI
libraries installed on your system: refer to your MPI manual and/or to
your system administrator for further information).

# About these tutorials {#about-this-tutorials .unnumbered}

The first part of this collection of tutorials comprises four examples taken from
Refs. [@marzari-prb97; @souza-prb01]: gallium arsenide, lead, silicon
and copper. All of the `wannier90` input files have been provided.

The second part of this collection of tutorials covers the generation of
`wannier90` input files starting from a full electronic structure
calculation. We have provided input files for the
`pwscf` interface (<http://www.quantum-espresso.org>) to
`wannier90`. Therefore, you will need to install and compile elements of
the ` quantum-espresso` package, namely `pw.x` and ` pw2wannier90.x`, in
order to run these tutorials. Please visit
<http://www.quantum-espresso.org> to download the package, and for
installation instructions. The tutorials work with
`pwscf v5.1.x` and `v6.0.x`. The exception are the tutorials on
symmetry adapted Wannier functions which require v6.0.x together with
the very latest version of `pw2wannier90.f90`. This can be found in the
directory `pwscf/v6.0` in the wannier distribution. It should be moved
to `PP/src` in the `pwscf` distribution and compiled using
`make pp`. Later versions v6.x.x should have the most up-to-date version
of pw2wannier90.f90 already included in the Quantum ESPRESSO
distribution.

There are interfaces to a number of other electronic structure codes
including abinit (<http://www.abinit.org>), fleur
(<http://www.flapw.de>), OpenMX (<http://www.openmx-square.org/>), GPAW
(<https://wiki.fysik.dtu.dk/gpaw/>), VASP (<http://www.vasp.at>), and
Wien2k (<http://www.wien2k.at>)

# Contact us {#contact-us .unnumbered}

If you have any suggestions regarding ways in which these tutorials may be
improved, then send us an email.

For other questions, email the `wannier90` forum at
` wannier@quantum-espresso.org`. Note that first you will need to
register in order to post emails. Emails from non-registered users are
deleted automatically. You can register by following the links at
<http://www.wannier.org/forum.html>.
