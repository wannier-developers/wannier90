# Frequently Asked Questions

## General Questions

### What is `wannier90`?

`wannier90` is a computer package, written in Fortran90, for obtaining
maximally-localised Wannier functions, using them to calculate
bandstructures, Fermi surfaces, dielectric properties, sparse
Hamiltonians and many things besides.

### Where can I get `wannier90`?

The most recent release of `wannier90` is always available on our
website <http://www.wannier.org>.

### Where can I get the most recent information about `wannier90`?

The latest news about `wannier90` can be followed on our website
<http://www.wannier.org>.

### Is `wannier90` free?

Yes! `wannier90` is available for use free-of-charge under the GNU
General Public Licence. See the file `LICENSE` that comes with the
`wannier90` distribution or the GNU hopepage at <http://www.gnu.org>.

## Getting Help

### Is there a Tutorial available for `wannier90`?

Yes! The `examples` directory of the `wannier90` distribution contains
input files for a number of tutorial calculations. The ` doc` directory
contains the accompanying tutorial handout.

### Where do I get support for `wannier90`?

There are a number of options:

1.  The `wannier90` User Guide, available in the `doc` directory of the
    distribution, and from the webpage
    (<http://www.wannier.org/user_guide.html>)

2.  The `wannier90` webpage for the most recent announcements
    (<http://www.wannier.org>)

3.  The `wannier90` mailing list (see
    <http://www.wannier.org/forum.html>)

### Is there a mailing list for `wannier90`?

Yes! You need to register: go to <http://www.wannier.org/forum.html> and
follow the instructions.

## Providing Help: Finding and Reporting Bugs

### I think I found a bug. How do I report it?

-   Check and double-check. Make sure it's a bug.

-   Check that it is a bug in `wannier90` and not a bug in the software
    interfaced to `wannier90`.

-   Check that you're using the latest version of `wannier90`.

-   Send us an email. Make sure to describe the problem and to attach
    all input and output files relating to the problem that you have
    found.

### I have got an idea! How do I report a wish?

We're always happy to listen to suggestions. Email your idea to the
`wannier90` developers.

### I want to help! How can I contribute to `wannier90`?

Great! There's always plenty of functionality to add. Email us to let us
know about the functionality you'd like to contribute.

### I like `wannier90`! Should I donate anything to its authors?

Our Swiss bank account number is\... just kidding! There is no need to
donate anything, please just cite our paper in any publications that
arise from your use of `wannier90`:

> \[ref\] G. Pizzi, V. Vitale, R. Arita, S. Blügel, F. Freimuth, G.
> Géranton, M. Gibertini, D. Gresch, C. Johnson, T. Koretsune, J.
> Ibañez-Azpiroz, H. Lee, J.M. Lihm, D. Marchand, A. Marrazzo, Y.
> Mokrousov, J.I. Mustafa, Y. Nohara, Y. Nomura, L. Paulatto, S. Poncé,
> T. Ponweiser, J. Qiao, F. Thöle, S.S. Tsirkin, M. Wierzbowska, N.
> Marzari, D. Vanderbilt, I. Souza, A.A. Mostofi, J.R. Yates,\
> Wannier90 as a community code: new features and applications, *J.
> Phys. Cond. Matt.* **32**, 165902 (2020)\
> <https://doi.org/10.1088/1361-648X/ab51ff>

If you are using versions 2.x of the code, cite instead:

> \[ref\] A. A. Mostofi, J. R. Yates, G. Pizzi, Y.-S. Lee, I. Souza,
> D. Vanderbilt and N. Marzari,\
> An updated version of `wannier90`: A Tool for Obtaining
> Maximally-Localised Wannier Functions, *Comput. Phys. Commun.*
> **185**, 2309 (2014)\
> <http://doi.org/10.1016/j.cpc.2014.05.003>

## Installation

### How do I install `wannier90`?

Follow the instructions in the file `README.install` in the main
directory of the `wannier90` distribution.

### Are there `wannier90` binaries available?

Not at present.

### Is there anything else I need?

Yes. `wannier90` works on top of an electronic structure calculation.

At the time of writing there are public, fully functioning, interfaces
between `wannier90` and [pwscf]{.smallcaps}, abinit
(<http://www.abinit.org>), siesta (<http://www.icmab.es/siesta/>), VASP
(<https://www.vasp.at>), Wien2k (<http://www.wien2k.at>), fleur
(<http://www.fleur.de>), OpenMX (<http://www.openmx-square.org/>), GPAW
(<https://wiki.fysik.dtu.dk/gpaw/>).

To use `wannier90` in combination with [pwscf]{.smallcaps} code (a
plane-wave, pseudopotential, density-functional theory code, which is
part of the `quantum-espresso` package) you will need to download
[pwscf]{.smallcaps} from the webpage <http://www.quantum-espresso.org>.
Then compile [pwscf]{.smallcaps} and the `wannier90` interface program
`pw2wannier90`. For instructions, please refer to the documentation that
comes with the `quantum-espresso` distribution.

For examples of how to use [pwscf]{.smallcaps} and `wannier90` in
conjunction with each other, see the `wannier90` Tutorial.
