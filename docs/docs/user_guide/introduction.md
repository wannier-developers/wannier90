# Introduction

## Getting Help

The latest release of `wannier90` and documentation can always be found
at <https://www.wannier.org>.

The development version may be cloned/downloaded from the official
repository of the `wannier90` code on GitHub (see
<https://github.com/wannier-developers/wannier90>).

There is a `wannier90` mailing list for discussing issues in the
development, theory, coding and algorithms pertinent to MLWF. You can
register for this mailing list by following the links at
<https://www.wannier.org/forum.html>. Alternatively, for technical issues
about the `wannier90` code, check the official repository of
`wannier90` on GitHub where you may raise issues or ask questions about
about its functionalities.

Finally, many frequently asked questions are answered in
Appendix [2](#chap:faq){reference-type="ref" reference="chap:faq"}. An
expanded FAQ session may be found on the Wiki page of the GitHub
repository at
<https://github.com/wannier-developers/wannier90/wiki/FAQ>.

## Citations

We ask that you acknowledge the use of `wannier90` in any publications
arising from the use of this code through the following reference

> \[ref\] G. Pizzi, V. Vitale, R. Arita, S. Blügel, F. Freimuth, G.
> Géranton, M. Gibertini, D. Gresch, C. Johnson, T. Koretsune, J.
> Ibañez-Azpiroz, H. Lee, J.M. Lihm, D. Marchand, A. Marrazzo, Y.
> Mokrousov, J.I. Mustafa, Y. Nohara, Y. Nomura, L. Paulatto, S. Poncé,
> T. Ponweiser, J. Qiao, F. Thöle, S.S. Tsirkin, M. Wierzbowska, N.
> Marzari, D. Vanderbilt, I. Souza, A.A. Mostofi, J.R. Yates,
> Wannier90 as a community code: new features and applications, *J.
> Phys. Cond. Matt.* **32**, 165902 (2020)
> <https://doi.org/10.1088/1361-648X/ab51ff>

If you are using versions 2.x of the code, cite instead:

> \[ref\] A. A. Mostofi, J. R. Yates, G. Pizzi, Y.-S. Lee, I. Souza,
> D. Vanderbilt and N. Marzari,
> An updated version of `wannier90`: A Tool for Obtaining
> Maximally-Localised Wannier Functions, *Comput. Phys. Commun.*
> **185**, 2309 (2014)
> <http://dx.doi.org/10.1016/j.cpc.2014.05.003>

It would also be appropriate to cite the original articles:

> Maximally localized generalized Wannier functions for composite energy bands,
> N. Marzari and D. Vanderbilt, *Phys. Rev. B* **56**, 12847 (1997)

and

> Maximally localized Wannier functions for entangled energy bands,
> I. Souza, N. Marzari and D. Vanderbilt, *Phys. Rev. B* **65**, 035109 (2001)

## Credits

The Wannier90 Developer Group includes Giovanni Pizzi (EPFL, CH),
Valerio Vitale (Cambridge, GB), David Vanderbilt (Rutgers University,
US), Nicola Marzari (EPFL, CH), Ivo Souza (Universidad del Pais Vasco,
ES), Arash A. Mostofi (Imperial College London, GB), and Jonathan R.
Yates (University of Oxford, GB).

The present release of `wannier90` was written by the Wannier90
Developer Group together with Ryotaro Arita (Riken and U. Tokyo, JP),
Stefan Blügel (FZ Jülich, DE), Frank Freimuth (FZ Jülich, DE), Guillame
Géranton (FZ Jülich, DE), Marco Gibertini (EPFL and University of
Geneva, CH), Dominik Gresch (ETHZ, CH), Charles Johnson (Imperial
College London, GB), Takashi Koretsune (Tohoku University and JST
PRESTO, JP), Julen Ibañez-Azpiroz (Universidad del Pais Vasco, ES),
Hyungjun Lee (EPFL, CH), Jae-Mo Lihm (Seoul National University, KR),
Daniel Marchand (EPFL, CH), Antimo Marrazzo (EPFL, CH), Yuriy Mokrousov
(FZ Jülich, DE), Jamal I. Mustafa (UC Berkeley, USA), Yoshiro Nohara
(Tokyo, JP), Yusuke Nomura (U. Tokyo, JP), Lorenzo Paulatto (Sorbonne
Paris, FR), Samuel Poncé (Oxford University, GB), Thomas Ponweiser (RISC
Software GmbH, AT), Florian Thöle (ETHZ, CH), Stepan Tsirkin
(Universidad del Pais Vasco, ES), Małgorzata Wierzbowska (Polish Academy
of Science, PL).

Contributors to the code include: Daniel Aberg (w90pov code), Lampros
Andrinopoulos (w90vdw code), Pablo Aguado Puente (gyrotropic routines),
Raffaello Bianco (k-slice plotting), Marco Buongiorno Nardelli (dosqc
v1.0 subroutines upon which transport.f90 is based), Stefano De
Gironcoli (pw2wannier90.x interface to Quantum ESPRESSO), Pablo Garcia
Fernandez (matrix elements of the position operator), Nicholas D. M.
Hine (w90vdw code), Young-Su Lee (specialised Gamma point routines and
transport), Antoine Levitt (preconditioning), Graham Lopez (extension of
pw2wannier90 to add terms needed for orbital magnetisation), Radu Miron
(constrained centres), Nicolas Poilvert (transport routines), Michel
Posternak (original plotting routines), Rei Sakuma (Symmetry-adapted
Wannier functions), Gabriele Sclauzero (disentanglement in spheres in
k-space), Matthew Shelley (transport routines), Christian Stieger
(routine to print the U matrices), David Strubbe (various
bugfixes/improvements), Timo Thonhauser (extension of pw2wannier90 to
add terms needed for orbital magnetisation), Junfeng Qiao (spin Hall
conductivity, projectability-disentangled Wannier functions).

We also acknowledge individuals not already mentioned above who
participated in the first Wannier90 community meeting (San Sebastian,
2016) for useful discussions: Daniel Fritsch, Victor Garcia Suarez,
Jan-Philipp Hanke, Ji Hoon Ryoo, Jürg Hutter, Javier Junquera, Liang
Liang, Michael Obermeyer, Gianluca Prandini, Paolo Umari.

`wannier90` Version 2.x was written by: Arash A. Mostofi, Giovanni
Pizzi, Ivo Souza, Jonathan R. Yates. `wannier90` Version 1.0 was written
by: Arash A. Mostofi, Jonathan R. Yates, Young-Su Lee. `wannier90` is
based on the Wannier Fortran 77 code written for isolated bands by
Nicola Marzari and David Vanderbilt and for entangled bands by Ivo
Souza, Nicola Marzari, and David Vanderbilt.

`wannier90` © 2007-2020 The Wannier Developer Group and individual
contributors

## Licence

All the material in this distribution is free software; you can
redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
