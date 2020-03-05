# CHANGELOG of Wannier90

## v3.1.0 (5th March 2020)

### New features

- Calculation of spin Hall conductivity according to the formalism given in Junfeng Qiao, Jiaqi Zhou, Zhe Yuan and Weisheng Zhao, PRB 98, 214402 (2018) + examples 29 and 30 and tests [[#264]](https://github.com/wannier-developers/wannier90/pull/264)

- Implementation of the SCDM method in Wannier90 for spinor wavefunctions and added example31 for the tutorial [[#277]](https://github.com/wannier-developers/wannier90/pull/277)

- Utility `w90spn2spn.x` for conversion between formatted and unformatted `.spn` files [[#274]](https://github.com/wannier-developers/wannier90/pull/274)

### Various improvements and bugfixes

- Add `install` target for `make`, now `make install` works [[#307]](https://github.com/wannier-developers/wannier90/pull/307)

- Added new documentation and fixed typos [[#308]](https://github.com/wannier-developers/wannier90/pull/308) [[#256]](https://github.com/wannier-developers/wannier90/pull/256) [[#306]](https://github.com/wannier-developers/wannier90/pull/306) [[#300]](https://github.com/wannier-developers/wannier90/pull/300)

- small restructuring of the routines for the AHC calculation with multiple Fermi energies and adaptive refinement [[#289]](https://github.com/wannier-developers/wannier90/pull/289)

- Fix in restart mode (now gathering matrices before dumping `.chk` file) [[#283]](https://github.com/wannier-developers/wannier90/pull/283)

- Fix `pw90common_fourier_R_to_k` subroutines bugs when `use_ws_distance = .true.` [[#273]](https://github.com/wannier-developers/wannier90/pull/273)[[#271]](https://github.com/wannier-developers/wannier90/issues/271)

- Fixes and speed improvements in `pw2wannier90.f90` [[#295]](https://github.com/wannier-developers/wannier90/pull/295) [[#255]](https://github.com/wannier-developers/wannier90/pull/255)

- Bug fixes of `gw2wannier90.py` and `k_mapper.py` [[#294]](https://github.com/wannier-developers/wannier90/pull/294)

- Various bugfixes in input reading and output generation [[#309]](https://github.com/wannier-developers/wannier90/pull/309) [[#290]](https://github.com/wannier-developers/wannier90/pull/290) [[#270]](https://github.com/wannier-developers/wannier90/pull/270) [[#266]](https://github.com/wannier-developers/wannier90/pull/266) [[#258]](https://github.com/wannier-developers/wannier90/pull/258)

- Adding tests for BoltzWann (at least for the electrical conductivity) [[#282]](https://github.com/wannier-developers/wannier90/pull/282)

- Bug fixed in the disentanglement step with SAWFs, especially when running in parallel [[#315]](https://github.com/wannier-developers/wannier90/pull/315)

## v3.0.0 (27th February 2019)

### New Wannier90 features

- Selective localization and constrained centres from Marianetti et al. [PRB 90, 165125 (2014)] [[#187]](https://github.com/wannier-developers/wannier90/pull/187) [[#207]](https://github.com/wannier-developers/wannier90/pull/207) [[#206]](https://github.com/wannier-developers/wannier90/pull/206) [[#205]](https://github.com/wannier-developers/wannier90/pull/205)

- Implementation of the SCDM method in Wannier90 [[#167]](https://github.com/wannier-developers/wannier90/pull/167) [[#202]](https://github.com/wannier-developers/wannier90/pull/202) and added an example for the tutorial [[#194]](https://github.com/wannier-developers/wannier90/pull/194). Note that during development we have adapted the interface of the SCDM input flags; Quantum ESPRESSO v6.4 has the current implementation. [[#238]](https://github.com/wannier-developers/wannier90/pull/238)


- Parallelisation of the core Wannier90 routines [[#193]](https://github.com/wannier-developers/wannier90/pull/193) [[#185]](https://github.com/wannier-developers/wannier90/pull/185) [[#173]](https://github.com/wannier-developers/wannier90/pull/173) [[#151]](https://github.com/wannier-developers/wannier90/pull/151) [[#149]](https://github.com/wannier-developers/wannier90/pull/149) [[#125]](https://github.com/wannier-developers/wannier90/pull/125) [[#120]](https://github.com/wannier-developers/wannier90/pull/120) [[#76]](https://github.com/wannier-developers/wannier90/pull/76)

- Added feature to select projections from a longer list [[#200]](https://github.com/wannier-developers/wannier90/pull/200)
  [[#215]](https://github.com/wannier-developers/wannier90/pull/215)

- Preconditioning of the SD/CG algorithm to optimize the spread [[#98]](https://github.com/wannier-developers/wannier90/pull/98)

- Allow to generate cube files also for non orthogonal unit cells [[#162]](https://github.com/wannier-developers/wannier90/pull/162), documentation [[#179]](https://github.com/wannier-developers/wannier90/pull/179) and examples [[#198]](https://github.com/wannier-developers/wannier90/pull/198)

- Enable plotting of spinor Wannier functions [[#157]](https://github.com/wannier-developers/wannier90/pull/157) [[#155]](https://github.com/wannier-developers/wannier90/pull/155)

- Improve how band paths are generated when discontinuities are present [[#199]](https://github.com/wannier-developers/wannier90/pull/199)
Moreover, it generates a new file `seedname_band.labelinfo.dat` with information on where the high-symmetry points were located.

- Add possibility to print b-vectors on file. [[#176]](https://github.com/wannier-developers/wannier90/pull/176)

- Add various command line arguments (e.g. for dry-run) and command-line help [[#158]](https://github.com/wannier-developers/wannier90/pull/158)

- Change of some variable defaults [[#218]](https://github.com/wannier-developers/wannier90/pull/218)

### New postw90 features, optimizations and new post-processing codes

- Calculation of nonlinear shift currents according to the formalism given in J. Ibañez-Azpiroz, S. S. Tsirkin and I. Souza, arXiv:1804.04030 (2018) + example 25 [[#181]](https://github.com/wannier-developers/wannier90/pull/181) [[#180]](https://github.com/wannier-developers/wannier90/pull/180)
  
- New gyrotropic module implementing the method described in S. S. Tsirkin, P. Aguado Puente and I. Souza, arXiv:1710.03204 [[#160]](https://github.com/wannier-developers/wannier90/pull/160)

- G0W0 interface for Wannier90 (tested with QE and Yambo) and documentation [[#102]](https://github.com/wannier-developers/wannier90/pull/102) [[#96]](https://github.com/wannier-developers/wannier90/pull/96)

- Performance improvements of parallel postw90 routines [[#153]](https://github.com/wannier-developers/wannier90/pull/153) [[#108]](https://github.com/wannier-developers/wannier90/pull/108)

### Various improvements and bugfixes

- New extensive solution booklet, contributed by V. Vitale, with solutions for examples 1-22 [[#233]](https://github.com/wannier-developers/wannier90/pull/233)

- Fix guiding center gamma point bug [[#223]](https://github.com/wannier-developers/wannier90/pull/223)

- Fixed issue in the search for points inside the Wigner-Seitz cell, making the code more robust for very anisotropic/stretched cells [[#216]](https://github.com/wannier-developers/wannier90/pull/216)

- Fixes to bugs in library mode [[#196]](https://github.com/wannier-developers/wannier90/pull/196) [[#224]](https://github.com/wannier-developers/wannier90/pull/224) [[#170]](https://github.com/wannier-developers/wannier90/pull/170)

- Fix U matrices not being written [[#177]](https://github.com/wannier-developers/wannier90/pull/177)

- Fixed the `hamiltonian_write_rmn` routine to correspond to `hamiltonian_write_tb` and use the proper formula for diagonal elements [[#168]](https://github.com/wannier-developers/wannier90/pull/168)

- Bugfix to `read_sym` in pw2wannier90 [[#166]](https://github.com/wannier-developers/wannier90/pull/166)

- Fix python script generated by kslice [[#154]](https://github.com/wannier-developers/wannier90/pull/154)

- Bugfixes for reading and writing uHu ad uIu matrices [[#146]](https://github.com/wannier-developers/wannier90/pull/146) [[#140]](https://github.com/wannier-developers/wannier90/pull/140) [[#97]](https://github.com/wannier-developers/wannier90/pull/97)

- Fixed index range when plotting U matrices to file [[#145]](https://github.com/wannier-developers/wannier90/pull/145)

- Fix a bug in Gamma-only routines [[#132]](https://github.com/wannier-developers/wannier90/pull/132)

- Fix to geninterp when run with more processors than kpoints [[#129]](https://github.com/wannier-developers/wannier90/pull/129)

- Update of `get_oper.F90` to correctly pass matrix elements [[#126]](https://github.com/wannier-developers/wannier90/pull/126)

- Fix Fermi energy initialisation in parameters [[#121]](https://github.com/wannier-developers/wannier90/pull/121)

- Improvements to the Wigner-Seitz detection routines [[#117]](https://github.com/wannier-developers/wannier90/pull/117) [[#109]](https://github.com/wannier-developers/wannier90/pull/109)

- Fix berry_task check for morb, and add check for kpoint_path block in parameters [[#258]](https://github.com/wannier-developers/wannier90/pull/258)

- Use 64 bit integer in io_wallclocktime [[#266]](https://github.com/wannier-developers/wannier90/pull/266)



### Improvement for compilation and for developers
- Improvements to the testing infrastructure and addding code coverage check [[#220]](https://github.com/wannier-developers/wannier90/pull/220) [[#217]](https://github.com/wannier-developers/wannier90/pull/217) [[#184]](https://github.com/wannier-developers/wannier90/pull/184) [[#182]](https://github.com/wannier-developers/wannier90/pull/182) [[#163]](https://github.com/wannier-developers/wannier90/pull/163) [[#204]](https://github.com/wannier-developers/wannier90/pull/204) [[#156]](https://github.com/wannier-developers/wannier90/pull/156) [[#141]](https://github.com/wannier-developers/wannier90/pull/141) [[#127]](https://github.com/wannier-developers/wannier90/pull/127)

- Added pre-commit hooks to fix indentation and trailing whitespace [[#203]](https://github.com/wannier-developers/wannier90/pull/203)

- Add support for dynamic library compilation [[#188]](https://github.com/wannier-developers/wannier90/pull/188)
  
- Added config file for Mac OS X using HomeBrew [[#118]](https://github.com/wannier-developers/wannier90/pull/118)

- Small fix to have some compilers (NAG) not to complain [[#93]](https://github.com/wannier-developers/wannier90/pull/93)




## v2.1.0 (13th January 2017)

### New features

- Implementation of the symmetry-adapted Wannier functions
  (see R. Sakuma, Phys. Rev. B 87, 235109 (2013), courtesy
   of R. Sakuma (Lund University, Sweden), T. Koretsune (Riken, JP),
   Y. Nomura (U. Tokyo, JP), Y. Nohara (Atomic-Scale Material
   Simulations, Co., Ltd.), R. Arita (Riken, JP))

- Streamlined the interface between wannier90 and tight-binding
  codes such as pythtb (new input variable: `write_tb`). Also,
  matrix elements of the position operator can now be printed
  (courtesy P. Garcia Fernandez, Unican, ES)

- Non-collinear sping with ultrasoft pseudos now implemented
  in the pw2wannier90 interface with Quantum ESPRESSO, working
  also in parallel (courtesy F. Thoele (ETHZ, CH), T. Koretsune
  (Riken, JP), L. Paulatto (UPMC Paris))

- postw90 prints an error message when trying to evaluate a quantity
  that requires a Fermi level, but none has been provided by the user
  (input variable fermi_level now has no default value)

- Fixed bug in kslice.F90 (trying to plot isoenergy contours colored
  by the spin gave an error message at runtime).

- Adaptively-refined mesh now implemented correctly for even sizes
  (e.g., 4x4).

- Fixed bug in the calculation of orbital magnetization in berry.F90

- Fixed the `skip_B1_tests` flag, that was ignored

- Fix to documentation of geninterp: k units were 1/ang and not
  2pi/ang as previously stated; also made names for absolute/fractional
  coordinates consistent with the rest of the code:now 'crystal' or
  'frac' mean fractional coordinates for k-points, 'cart' and
  'abs' mean absolute.

- Improvements to various Makefiles

- Renamed make.sys to make.inc to cope with possible issues when this
  file is sent over email (some providers think it's a virus and
  remove it from attachments). If you have already compiled wannier90,
  after updating the version, just rename make.sys to make.inc.

- Added the `use_ws_distance` flag to improve the interpolation of
  band structures (courtesy of L. Paulatto, UPMC Paris).

- Improved the interface with the Z2pack code (courtesy D. Gresch, ETHZ)

### New features for developers

- Added a test-suite, and integrated with GitHub and Travis-CI for
  continuous integration. A number of tests have been added
  (contributed mainly by S. Ponce, Oxford). Also compilation on the
  buildbot test farm at the Oxford Materials Modelling Laboratory
  has been activated.
- First implementation of the FORD infrastructure for code documentation
  (courtesy D. Gresch, ETHZ).

## v2.0.1 (2nd April 2015)

- Added the possibility to disentangle only in small spherical regions
  in k space around a selected number of k points - see also Example 20.
  (contributed by Gabriele Sclauzero, ETH Zurich)

- Added the `skip_B1_tests` flag to be compatible with the Z2PACK code
  (http://www.physics.rutgers.edu/z2pack/)

- Reintroduced the `boltz_bandshift` flag for a rigid shift of the bands
  in BoltzWann; updated reference of BoltzWann paper; added a version
  of BoltzWann tutorial example (example 16) that works also without
  Quantum ESPRESSO.

- Added documentation on pw2wannier90.f90 for QE 5.1.x

- bug fix in w90chk2chk.F90 (missing allocation if there were excluded bands)

- bug fix in kslice.F90 (missing transpose when plotting Fermi lines in
  the python scripts)

## v2.0 (14th October 2013)

- Updated interface to PWscf v5.0 (including a minor bug fix to the
  definition of the projection functions).

- Enabled very general specification of spinor projections

- Generalized the definition of the supercell for real-space plotting of
  the Wannier functions: now the supercell size (`wannier_plot_supercell`)
  can be different along the three directions. Moreover, the home cell
  is left approximately at the center of the supercell, rather than near
  one of its edges.

- Now it is possible to provide both 'seedname' and 'seedname.win' on the
  command line, both will work

- Tabs can now be present in the input file without problems

- Added the new (parallel) postw90.x executable for calculations which use
  the MLWFs calculated by wannier90.x as an input.

- Added the following postw90.x modules:

  - dos
  - berry
  - kpath
  - kslice
  - BoltzWann
  - geninterp

- Added the w90chk2chk.x utility to convert the checkpoint file between
  the formatted/unformatted formats, to move it between different computers
  with different architectures

- Added the w90vdw utility to calculate van der Waals energies with MLWFs

- Added the w90pov utility to render ray-traced isosurfaces using POV-Ray

- A few bugfixes

## v1.2 (14th Jan 2010)

- The information written to the file `seedname_hr.dat` (the Hamiltonian
  in the WF basis) has been extended to include the number of WF (`num_wann`),
  the number of Wigner-Seitz points (nrpts) and the degeneracy of each point
  (ndegen).

- The information contained in seedname.chk has been extended to
  include the number of bands (`num_bands`), the number of excluded
  bands (`num_exclude_bands`), the excluded band indices (`exclude_bands`)
  and the Monkhorst-Pack grid dimensions (`mp_grid`). As a result
  v1.2 is not compatible with checkfiles written with older versions of
  the code.

- Automated lcr transport calculations from a single supercell.
  Includes robust sorting and Wannier function parity determination
  algorithms (new input variables: `tran_num_cell_ll`, `tran_num_cell_rr`,
  `tran_group_threshold`, `easy_fix`)

- Added examples 14, 15 to tutorial to displaying new lcr transport
  functionality

- Updated interface to PWscf v4.1.2 (new input variable `write_unkg`)

- New utility: PL_assessment. Investigates principal layer size
  from bulk transport, many k-point lead calculations.

- Altered `<seedname>_centres.xyz` output file to include atomic
  positions in transport calculations

- Fixes for `tran_num_cell_ll=1` in `tran_hr_one_dim_cut`

- Minor bug fixes

- Addition of Matthew Shelley and Nicolas Poilvert as contributors.

- Report of estimated memory usage

- Improved memory efficiency

## v1.1 (21st Dec 2007)

- Addition of specific algorithms for when only Gamma-point
  sampling is used (new input variable: `gamma_only`)

- Addition of routines for quantum transport and DoS calculations
  (new input variables: `transport`, `transport_mode`, etc.)

- Option to write out hamiltonian matrix elements in the Wannier
  function basis (new input variable: `hr_plot`)

- Option to set a convergence threshold for localisation procedure
  (new input variables: `conv_tol`, `conv_window`)

- Improved minimisation algorithms for localisation routines
  (new input variables: `conv_noise_amp`, `conv_noise_num`)

- Option to specify the number of shells that are searched to find
  nearest neighbour b-vectors (new input variable: `search_shells`)

- Option to plot bandstructures in xmgrace format (`bands_plot_format=xmgrace`)

- Option to plot Wannier functions in cube format (`wannier_plot_format=cube`)
  -- works for isolated molecules, further testing for periodic systems is
  required -- significantly reduces WF file-size
  (new input variable: `wannier_plot_radius`)

- Optional capability to specify some projections in input file and
  have the remaining centres chosen randomly by the code

- Checkpointing and restarts all done via the .chk file (`_um.dat` file
  now obsolete)

- Further enhancements to the way projections are specified

- Option to map Wannier functions onto bandstructure
  (new input variable: `bands_plot_project`)

- Option to have spinor Wannier functions
  (new input variable: `spinors`)

- A few new tutorial examples

- Improvements to "library mode" functionality

## v1.0.2 (1st Dec 2006)

- Addition of "library mode" functionality

- Introduction of "range vectors" for specifying `exclude_bands` and
  `wannier_plot_list` in input file

- Option to specify random projections

- Option to use Bloch phases for initial projections

- Addition of `timing_level` input flag to control how much timing
  information is outputted

- Option to translate final centres to the home unit cell

- Option to write final centres in xyz format

- Acceleration of disentanglement procedure

- Speed-up of localisation routines

- Improved robustness of plotting routines

## v1.0.1 (17th May 2006)

- Bug fix in wannierise minimiser - caused poor convergence in large systems

- Increase precision of k-points in `*.nnkp` file

- More robust selection of eigenvectors in disentanglement (`dis_proj_froz`)

- Addition of `write_proj` keyword -- outputs projection of original bands
  on final Wannier functions

- Longer strings for atom labels in `*.win`

- Minor format change to `*_bands.dat`

- Check restart keyword in `*.win`
