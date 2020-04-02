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

  G. Pizzi, V. Vitale, R. Arita, S. Blügel, F. Freimuth, G. Géranton, M. Gibertini, 
  D. Gresch, C. Johnson, T. Koretsune, J Ibañez-Azpiroz, H. Lee, J.M. Lihm, 
  D. Marchand, A. Marrazzo, Y. Mokrousov, J.I. Mustafa, Y. Nohara, Y. Nomura, 
  L. Paulatto, S. Poncé, T. Ponweiser, J. Qiao, F. Thöle, S.S. Tsirkin, 
  M. Wierzbowska, N. Marzari, D. Vanderbilt, I. Souza, A.A. Mostofi, J.R. Yates, 
  Wannier90 as a community code: new features and applications, 
  `J. Phys. Cond. Matt. 32, 165902`_ (2020)

.. _J. Phys. Cond. Matt. 32, 165902: https://doi.org/10.1088/1361-648X/ab51ff

If you are using versions 2.x of the code, cite instead:
                                                         
  A.A. Mostofi, J.R. Yates, G. Pizzi, Y.S. Lee, I. Souza, 
  D Vanderbilt, N Marzari, *An updated version of wannier90: A tool for 
  obtaining maximally-localised Wannier functions*, 
  `Comput. Phys. Commun. 185, 2309 (2014)`_ 

.. _Comput. Phys. Commun. 185, 2309 (2014): http://doi.org/10.1016/j.cpc.2014.05.003

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

FAQ
++++

Detailed information on how to install, contribute and report issues/bugs may be found on the `FAQ page`_.  Common issues encountered by the Wannier90-user community are also addressed there.

.. _FAQ page: https://github.com/wannier-developers/wannier90/wiki/FAQ

Authors and contributors
++++++++++++++++++++++++

The Wannier90 Developer Group includes:

* Giovanni Pizzi    (EPFL, CH)
* Valerio Vitale    (Cambridge, GB)
* David Vanderbilt  (Rutgers University, US)
* Nicola Marzari    (EPFL, CH)
* Ivo Souza         (Universidad del Pais Vasco, ES)
* Arash A. Mostofi  (Imperial College London, GB)
* Jonathan R. Yates (University of Oxford, GB)

In addition to the Wannier90 Developer Group, the other authors of Wannier90 v.3.x are:

* Ryotaro Arita (Riken and U. Tokyo, JP): Symmetry-adapted Wannier functions
* Stefan Blügel (FZ  Jülich, DE): Parallelization of the core routines
* Frank Freimuth (FZ  Jülich, DE): Parallelization of the core routines
* Guillame Géranton (FZ  Jülich, DE): Parallelization of the core routines
* Marco Gibertini (EPFL and University of Geneva, CH): Improvements to the interpolation routines
* Dominik Gresch (ETHZ, CH): FORD infrastructure for code documentation, AiiDA-Wannier90 interface
* Charles Johnson (Imperial College London, GB): Selectively-localised WFs and constrained centres
* Takashi Koretsune (Tohoku University and JST PRESTO, JP): Symmetry-adapted Wannier functions, non-collinear spin with ultrasoft in pw2wannier90
* Julen Ibañez-Azpiroz (Universidad del Pais Vasco, ES): shift-current calculation
* Hyungjun Lee (EPFL, CH): Spinor-valued WFs, parallelisation of the core Wannier90 routines
* Jae-Mo Lihm (Seoul National University, KR): SCDM-k implementation for non-collinear spin in pw2wannier90
* Daniel Marchand (EPFL, CH): AiiDA-Wannier90 interface
* Antimo Marrazzo (EPFL, CH): GW bands interpolation, AiiDA-Wannier90 interface
* Yuriy Mokrousov (FZ  Jülich, DE): Parallelization of the core routines
* Jamal I. Mustafa (UC Berkeley, USA): Subselection of projections
* Yoshiro Nohara (Tokyo, JP): Symmetry-adapted Wannier functions
* Yusuke Nomura (U. Tokyo, JP): Symmetry-adapted Wannier functions, non-collinear spin with ultrasoft in pw2wannier90
* Lorenzo Paulatto (Sorbonne Paris, FR): Improvements to the interpolation routines, non-collinear spin with ultrasoft in pw2wannier90
* Samuel Poncé (Oxford University, GB): Test suite for Wannier90
* Thomas Ponweiser (RISC Software GmbH, AT): performance optimizations for postw90
* Florian Thöle (ETHZ, CH): non-collinear spin with ultrasoft in pw2wannier90
* Stepan Tsirkin (Universidad del Pais Vasco, ES): GW bands interpolation, gyrotropic module, shift-current calculation, bug fixes in the berry module
* Małgorzata Wierzbowska (Polish Academy of Science, PL): performance optimizations for postw90

Contributors to the code include:

* Daniel Aberg: w90pov code
* Lampros Andrinopoulos: w90vdw code
* Pablo Aguado Puente: gyrotropic routines
* Raffaello Bianco: k-slice plotting
* Marco Buongiorno Nardelli: `dosqc` v1.0 subroutines upon which `transport.f90` is based.
* Stefano De Gironcoli: `pw2wannier90.x` interface to Quantum ESPRESSO
* Pablo Garcia Fernandez: Matrix elements of the position operator
* Nicholas D. M. Hine: w90vdw code
* Young-Su Lee: specialised Gamma point routines & transport
* Antoine Levitt: preconditioning
* Graham Lopez: extension of pw2wannier90 to add terms needed for orbital magnetisation
* Radu Miron: constrained centres
* Nicolas Poilvert: transport routines
* Michel Posternak: original plotting routines
* Rei Sakuma: Symmetry-adapted Wannier functions
* Gabriele Sclauzero: disentanglement in spheres in k-space
* Matthew Shelley: transport routines
* Christian Stieger: routine to print the U matrices
* David Strubbe: various bugfixes/improvements
* Timo Thonhauser: extension of pw2wannier90 to add terms needed for orbital magnetisation

We also acknowledge individuals not already mentioned above who participated in the first Wannier90 community meeting (San Sebastian, 2016) for useful discussions:

* Daniel Fritsch
* Victor Garcia Suarez
* Jan-Philipp Hanke
* Ji Hoon Ryoo
* Jürg Hutter
* Javier Junquera
* Liang Liang
* Michael Obermeyer
* Gianluca Prandini
* Paolo Umari

Wannier90 Version 2.x was written by:

* Arash A. Mostofi   (Imperial College London, GB)
* Giovanni Pizzi     (EPFL, CH)
* Ivo Souza          (Universidad del Pais Vasco, ES)
* Jonathan R. Yates  (University of Oxford, GB)

Wannier90 Version 1.0 was written by:

* Arash A. Mostofi   (Imperial College London, GB)
* Jonathan R. Yates  (University of Oxford, GB)
* Young-Su Lee       (KIST, KR)

Wannier90 is based on the [Wannier Fortran 77 code](http://www.wannier.org/history/) by:

* Nicola Marzari (EPFL, CH)
* Ivo Souza (Universidad del Pais Vasco, ES)
* David Vanderbilt (Rutgers University, US)
