=========
Wannier90
=========

The Maximally-Localised Generalised Wannier Functions Code 
----------------------------------------------------------

The homepage of the Wannier90 code is http://www.wannier.org

The code is hosted on 
GitHub_.

.. _GitHub: https://github.com/wannier-developers/wannier90

How to contribute
+++++++++++++++++

Discussions about development, roadmap, and further information of interest
to developers and contributors can be found on the 
`GitHub wiki`_.
In particular, if you want to contribute, you should first read the 
`guide to contributors`_ section.

.. _GitHub wiki: https://github.com/wannier-developers/wannier90/wiki/ContributorsGuide
.. _guide to contributors: https://github.com/wannier-developers/wannier90/wiki/ContributorsGuide

How to cite
+++++++++++
Please cite the following paper in any publications arising from the use of 
this code:
                                                         
  A.A. Mostofi, J.R. Yates, G. Pizzi, Y.S. Lee, I. Souza, 
  D Vanderbilt, N Marzari, *An updated version of wannier90: A tool for 
  obtaining maximally-localised Wannier functions*, 
  `Comput. Phys. Commun. 185, 2309 (2014)`_ 

.. _Comput. Phys. Commun. 185, 2309 (2014): http://dx.doi.org/10.1016/j.cpc.2014.05.003

For the method please cite:

  N. Marzari and D. Vanderbilt,
  *Maximally Localized Generalised Wannier Functions for Composite Energy Bands*,    
  `Phys. Rev. B 56 12847 (1997)`_
                                                  
  I. Souza, N. Marzari and D. Vanderbilt,
  *Maximally Localized Wannier Functions for Entangled Energy Bands*, 
  `Phys. Rev. B 65 035109 (2001)`_

  Nicola Marzari, Arash A. Mostofi, Jonathan R. Yates, Ivo Souza, 
  David Vanderbilt,
  *Maximally localized Wannier functions: Theory and applications*, 
  `Rev. Mod. Phys. 84, 1419 (2012)`_ 

.. _Phys. Rev. B 56 12847 (1997): http://dx.doi.org/10.1103/PhysRevB.56.12847
.. _Phys. Rev. B 65 035109 (2001): http://dx.doi.org/10.1103/PhysRevB.65.035109
.. _Rev. Mod. Phys. 84, 1419 (2012): http://dx.doi.org/10.1103/RevModPhys.84.1419


Licence
+++++++

The Wannier90 code is licensed under GPLv2. 
You can read the licence text in the LICENSE file in the root directory 
of the Wannier90 distribution.

Authors and contributors
++++++++++++++++++++++++

Wannier90 Versions 2.x have been written by: 

* Arash A. Mostofi   (Imperial College London, UK)
* Giovanni Pizzi     (EPFL, Switzerland)
* Ivo Souza          (Universidad del Pais Vasco, Spain)
* Jonathan R. Yates  (University of Oxford, UK)

Contributors to the code include:

* Young-Su Lee (KIST, S. Korea): specialised Gamma point routines & transport
* Matthew Shelley (Imperial College London, UK): transport
* Nicolas Poilvert (Penn State University, USA): transport
* Raffaello Bianco (Univ. Pierre et Marie Curie Paris 6 and CNRS):  k-slice plotting
* Gabriele Sclauzero (ETH, Zurich, Switzerland): disentanglement in spheres in k-space
* David Strubbe (MIT, USA): various bugfixes/improvements
* Rei Sakuma (Lund University, Sweden): Symmetry-adapted Wannier functions
* Yusuke Nomura (U. Tokyo, JP): Symmetry-adapted Wannier functions, non-collinear spin with ultrasoft in pw2wannier90
* Takashi Koretsune (Riken, JP): Symmetry-adapted Wannier functions, non-collinear spin with ultrasoft in pw2wannier90
* Yoshiro Nohara (Atomic-scale material simulations, Co., Ltd.): Symmetry-adapted Wannier functions
* Ryotaro Arita (Riken, JP): Symmetry-adapted Wannier functions
* Lorenzo Paulatto (UPMC Paris, FR): Improvements to the interpolation routines, non-collinear spin with ultrasoft in pw2wannier90
* Florian Thole (ETHZ, CH): non-collinear spin with ultrasoft in pw2wannier90
* Pablo Garcia Fernandez (Unican, ES): Matrix elements of the position operator
* Dominik Gresch (ETHZ, CH): FORD infrastructure for code documentation
* Samuel Ponce (Oxford University, UK): Test suite for Wannier90
* Marco Gibertini (EPFL, CH): Improvements to the interpolation routines
* Christian Stieger (ETHZ, CH): Routine to print the U matrices
* Stepan Tsirkin (Universidad del Pais Vasco, Spain): bug fixes in the berry module

Moreover:

* ``w90vdw`` was written by:
  Lampros Andrinopoulos, Nicholas D. M. Hine and Arash A. Mostofi
* ``w90pov`` (povray plotting) was written by:
  Daniel Aberg (LLNL, USA)

We gratefully acknowledge Marco Buongiorno Nardelli for the `dosqc
v1.0` subroutines for the calculation of the density of states and the
quantum transport, upon which `transport.f90` is based. 

Wannier90 Version 1.0 was written by:

* Arash A. Mostofi   (Imperial College London, UK)
* Jonathan R. Yates  (University of Oxford, UK)
* Young-Su Lee       (KIST, S. Korea)

Wannier90 is based on Fortran 77 codes written by:

* Nicola Marzari (EPFL, Switzerland)
* Ivo Souza (Universidad del Pais Vasco, Spain)
* David Vanderbilt (Rutgers University, USA)

