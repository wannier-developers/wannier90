!-*- mode: F90 -*-!
!------------------------------------------------------------!
! This file is distributed as part of the Wannier90 code and !
! under the terms of the GNU General Public License. See the !
! file `LICENSE' in the root directory of the Wannier90      !
! distribution, or http://www.gnu.org/copyleft/gpl.txt       !
!                                                            !
! The webpage of the Wannier90 code is www.wannier.org       !
!                                                            !
! The Wannier90 code is hosted on GitHub:                    !
!                                                            !
! https://github.com/wannier-developers/wannier90            !
!------------------------------------------------------------!
!                                                            !
!  w90_berry: computes Berry phase and related properties    !
!                                                            !
!------------------------------------------------------------!

module w90_berry

  !! This module computes various "Berry phase" related properties
  !!
  !! Key REFERENCES
  !!
  !! *  WYSV06 = PRB 74, 195118 (2006)  (anomalous Hall conductivity - AHC)
  !! *  YWVS07 = PRB 75, 195121 (2007)  (Kubo frequency-dependent conductivity)
  !! *  LVTS12 = PRB 85, 014435 (2012)  (orbital magnetization and AHC)
  !! *  CTVR06 = PRB 74, 024408 (2006)  (  "          "       )
  !! *  IATS18 = PRB 97, 245143 (2018)  (nonlinear shift current)
  !! *  QZYZ18 = PRB 98, 214402 (2018)  (spin Hall conductivity - SHC)
  !! *  RPS19  = PRB 99, 235113 (2019)  (spin Hall conductivity - SHC)
  !! *  IAdJS19 = arXiv:1910.06172 (2019) (quasi-degenerate k.p)
  ! ---------------------------------------------------------------
  !
  ! * Undocumented, works for limited purposes only:
  !                                 reading k-points and weights from file

  use w90_constants, only: dp
  use w90_error, only: w90_error_type, set_error_alloc, set_error_dealloc, set_error_fatal, &
    set_error_input, set_error_fatal, set_error_file

  implicit none

  private

  public :: berry_get_imfgh_klist
  public :: berry_get_imf_klist
  public :: berry_get_kdotp
  public :: berry_get_shc_klist
  public :: berry_main

  ! Pseudovector <--> Antisymmetric tensor
  !
  ! x <--> (y,z)
  ! y <--> (z,x)
  ! z <--> (x,y)
  !
  integer, dimension(3), parameter :: alpha_A = (/2, 3, 1/)
  integer, dimension(3), parameter :: beta_A = (/3, 1, 2/)

  ! Independent components of a symmetric tensor
  !
  ! 1 <--> xx
  ! 2 <--> yy
  ! 3 <--> zz
  ! 4 <--> xy
  ! 5 <--> xz
  ! 6 <--> yz
  !
  integer, dimension(6), parameter :: alpha_S = (/1, 2, 3, 1, 1, 2/)
  integer, dimension(6), parameter :: beta_S = (/1, 2, 3, 2, 3, 3/)
  integer, dimension(6), parameter, public :: berry_alpha_S = alpha_S
  integer, dimension(6), parameter, public::  berry_beta_S = beta_S
  integer, parameter, public:: berry_alpha_beta_S(3, 3) = &
                               reshape((/1, 4, 5, 4, 2, 6, 5, 6, 3/), (/3, 3/))

contains

  !================================================!
  !                   PUBLIC PROCEDURES
  !================================================!
  subroutine berry_main(pw90_berry, dis_manifold, fermi_energy_list, kmesh_info, kpoint_dist, &
                        kpt_latt, pw90_band_deriv_degen, pw90_oper_read, pw90_spin, physics, &
                        ws_region, pw90_spin_hall, wannier_data, ws_distance, wigner_seitz, &
                        print_output, AA_R, BB_R, CC_R, HH_R, SH_R, SHR_R, SR_R, SS_R, SAA_R, &
                        SBB_R, u_matrix, v_matrix, eigval, real_lattice, scissors_shift, mp_grid, &
                        fermi_n, num_wann, num_kpts, num_bands, num_valence_bands, &
                        effective_model, have_disentangled, spin_decomp, seedname, stdout, timer, &
                        error, comm)
    !================================================!
    !
    !! Computes the following quantities:
    !!   (i) Anomalous Hall conductivity (from Berry curvature)
    !!  (ii) Complex optical conductivity (Kubo-Greenwood) & JDOS
    !! (iii) Orbital magnetization
    !!  (iv) Nonlinear shift current
    !!   (v) Spin Hall conductivity
    !
    !================================================!

    use w90_comms, only: comms_reduce, w90comm_type, mpirank, mpisize
    use w90_constants, only: dp, cmplx_0, pi, pw90_physical_constants_type
    use w90_utility, only: utility_recip_lattice_base
    use w90_get_oper, only: get_HH_R, get_AA_R, get_BB_R, get_CC_R, get_SS_R, get_SHC_R, &
      get_SAA_R, get_SBB_R
    use w90_io, only: io_file_unit, io_stopwatch_start, io_stopwatch_stop
    use w90_types, only: print_output_type, wannier_data_type, &
      dis_manifold_type, kmesh_info_type, ws_region_type, ws_distance_type, timer_list_type
    use w90_postw90_types, only: pw90_berry_mod_type, pw90_spin_mod_type, &
      pw90_spin_hall_type, pw90_band_deriv_degen_type, pw90_oper_read_type, wigner_seitz_type, &
      kpoint_dist_type

    implicit none

    ! arguments
    type(pw90_berry_mod_type), intent(inout) :: pw90_berry
    type(dis_manifold_type), intent(in) :: dis_manifold
    type(kmesh_info_type), intent(in) :: kmesh_info
    type(kpoint_dist_type), intent(in) :: kpoint_dist
    type(pw90_band_deriv_degen_type), intent(in) :: pw90_band_deriv_degen
    type(pw90_oper_read_type), intent(in) :: pw90_oper_read
    type(pw90_spin_mod_type), intent(in) :: pw90_spin
    type(print_output_type), intent(in) :: print_output
    type(pw90_physical_constants_type), intent(in) :: physics
    type(ws_region_type), intent(in) :: ws_region
    type(pw90_spin_hall_type), intent(in) :: pw90_spin_hall
    type(w90comm_type), intent(in) :: comm
    type(wannier_data_type), intent(in) :: wannier_data
    type(wigner_seitz_type), intent(inout) :: wigner_seitz
    type(ws_distance_type), intent(inout) :: ws_distance
    type(timer_list_type), intent(inout) :: timer
    type(w90_error_type), allocatable, intent(out) :: error

    complex(kind=dp), allocatable, intent(inout) :: AA_R(:, :, :, :) ! <0n|r|Rm>
    complex(kind=dp), allocatable, intent(inout) :: BB_R(:, :, :, :) ! <0|H(r-R)|R>
    complex(kind=dp), allocatable, intent(inout) :: CC_R(:, :, :, :, :) ! <0|r_alpha.H(r-R)_beta|R>
    complex(kind=dp), allocatable, intent(inout) :: HH_R(:, :, :) !  <0n|r|Rm>
    complex(kind=dp), allocatable, intent(inout) :: SH_R(:, :, :, :) ! <0n|sigma_x,y,z.H|Rm>
    complex(kind=dp), allocatable, intent(inout) :: SHR_R(:, :, :, :, :) ! <0n|sigma_x,y,z.H.(r-R)_alpha|Rm>
    complex(kind=dp), allocatable, intent(inout) :: SR_R(:, :, :, :, :) ! <0n|sigma_x,y,z.(r-R)_alpha|Rm>
    complex(kind=dp), allocatable, intent(inout) :: SS_R(:, :, :, :) ! <0n|sigma_x,y,z|Rm>
    !spin Hall using Ryoo's method
    complex(kind=dp), allocatable, intent(inout) :: SAA_R(:, :, :, :, :) ! <0n|sigma_x,y,z.(r-R)_alpha|Rm>
    !! $$\langle 0n | \sigma_{x,y,z}.(\hat{r}-R)_{\alpha}  | Rm \rangle$$
    complex(kind=dp), allocatable, intent(inout) :: SBB_R(:, :, :, :, :) ! <0n|sigma_x,y,z.H.(r-R)_alpha|Rm>
    !! $$\langle 0n | \sigma_{x,y,z}.H.(\hat{r}-R)_{\alpha}  | Rm \rangle$$
    complex(kind=dp), intent(in) :: u_matrix(:, :, :), v_matrix(:, :, :)

    real(kind=dp), intent(in) :: eigval(:, :)
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    real(kind=dp), intent(in) :: scissors_shift
    real(kind=dp), allocatable, intent(in) :: fermi_energy_list(:)
    real(kind=dp), intent(in) :: kpt_latt(:, :)

    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: num_wann, num_kpts, num_bands, num_valence_bands, fermi_n
    integer, intent(in) :: stdout

    character(len=50), intent(in) :: seedname
    logical, intent(in) :: have_disentangled
    logical, intent(in) :: spin_decomp
    logical, intent(in) :: effective_model

    ! local variables
    real(kind=dp), allocatable :: adkpt(:, :)

    ! AHC and orbital magnetization, calculated for a list of Fermi levels
    !
    ! First index labels J0,J1,J2 terms, second labels the Cartesian component
    !
    real(kind=dp) :: imf_k_list(3, 3, fermi_n), imf_list(3, 3, fermi_n), imf_list2(3, 3, fermi_n)
    real(kind=dp) :: img_k_list(3, 3, fermi_n), img_list(3, 3, fermi_n)
    real(kind=dp) :: imh_k_list(3, 3, fermi_n), imh_list(3, 3, fermi_n)
    real(kind=dp) :: ahc_list(3, 3, fermi_n)
    real(kind=dp) :: LCtil_list(3, 3, fermi_n), ICtil_list(3, 3, fermi_n), Morb_list(3, 3, fermi_n)
    real(kind=dp) :: imf_k_list_dummy(3, 3, fermi_n) ! adaptive refinement of AHC
    ! shift current
    real(kind=dp), allocatable :: sc_k_list(:, :, :)
    real(kind=dp), allocatable :: sc_list(:, :, :)
    ! kdotp
    complex(kind=dp), allocatable :: kdotp(:, :, :, :, :)
    ! Complex optical conductivity, dividided into Hermitean and
    ! anti-Hermitean parts
    !
    complex(kind=dp), allocatable :: kubo_H_k(:, :, :)
    complex(kind=dp), allocatable :: kubo_H(:, :, :)
    complex(kind=dp), allocatable :: kubo_AH_k(:, :, :)
    complex(kind=dp), allocatable :: kubo_AH(:, :, :)
    ! decomposition into up-up, down-down and spin-flip transitions
    complex(kind=dp), allocatable :: kubo_H_k_spn(:, :, :, :)
    complex(kind=dp), allocatable :: kubo_H_spn(:, :, :, :)
    complex(kind=dp), allocatable :: kubo_AH_k_spn(:, :, :, :)
    complex(kind=dp), allocatable :: kubo_AH_spn(:, :, :, :)

    ! Joint density of states
    !
    real(kind=dp), allocatable :: jdos_k(:)
    real(kind=dp), allocatable :: jdos(:)
    ! decomposition into up-up, down-down and spin-flip transitions
    real(kind=dp), allocatable :: jdos_k_spn(:, :)
    real(kind=dp), allocatable :: jdos_spn(:, :)

    ! Spin Hall conductivity
    real(kind=dp), allocatable :: shc_fermi(:), shc_k_fermi(:)
    complex(kind=dp), allocatable :: shc_freq(:), shc_k_freq(:)
    ! for fermi energy scan, adaptive kmesh
    real(kind=dp), allocatable :: shc_k_fermi_dummy(:)

    real(kind=dp) :: cell_volume
    real(kind=dp) :: kweight, kweight_adpt, kpt(3), db1, db2, db3, fac, rdum, vdum(3)

    integer :: n, i, j, k, jk, ikpt, if, ierr, loop_x, loop_y, loop_z, kdotp_nbands
    integer :: loop_xyz, loop_adpt, adpt_counter_list(fermi_n), ifreq, file_unit
    integer :: my_node_id, num_nodes

    character(len=120) :: file_name

    logical :: eval_ahc, eval_morb, eval_kubo, not_scannable, eval_sc, eval_shc, eval_kdotp
    logical :: ladpt_kmesh
    logical :: ladpt(fermi_n)

    my_node_id = mpirank(comm)
    num_nodes = mpisize(comm)

    if (fermi_n == 0) then
      call set_error_input(error, 'Must specify one or more Fermi levels when berry=true', comm)
      return
    endif

    if (print_output%timing_level > 1 .and. print_output%iprint > 0) &
      call io_stopwatch_start('berry: prelims', timer)

    cell_volume = real_lattice(1, 1)*(real_lattice(2, 2)*real_lattice(3, 3) - &
                                      real_lattice(3, 2)*real_lattice(2, 3)) + &
                  real_lattice(1, 2)*(real_lattice(2, 3)*real_lattice(3, 1) - &
                                      real_lattice(3, 3)*real_lattice(2, 1)) + &
                  real_lattice(1, 3)*(real_lattice(2, 1)*real_lattice(3, 2) - &
                                      real_lattice(3, 1)*real_lattice(2, 2))

    ! Mesh spacing in reduced coordinates
    !
    db1 = 1.0_dp/real(pw90_berry%kmesh%mesh(1), dp)
    db2 = 1.0_dp/real(pw90_berry%kmesh%mesh(2), dp)
    db3 = 1.0_dp/real(pw90_berry%kmesh%mesh(3), dp)

    eval_ahc = .false.
    eval_morb = .false.
    eval_kubo = .false.
    eval_sc = .false.
    eval_shc = .false.
    eval_kdotp = .false.

    if (index(pw90_berry%task, 'ahc') > 0) eval_ahc = .true.
    if (index(pw90_berry%task, 'morb') > 0) eval_morb = .true.
    if (index(pw90_berry%task, 'kubo') > 0) eval_kubo = .true.
    if (index(pw90_berry%task, 'sc') > 0) eval_sc = .true.
    if (index(pw90_berry%task, 'shc') > 0) eval_shc = .true.
    if (index(pw90_berry%task, 'kdotp') > 0) eval_kdotp = .true.

    ! Wannier matrix elements, allocations and initializations
    !
    if (eval_ahc) then
      call get_HH_R(dis_manifold, kpt_latt, print_output, wigner_seitz, HH_R, u_matrix, v_matrix, &
                    eigval, real_lattice, scissors_shift, num_bands, num_kpts, num_wann, &
                    num_valence_bands, effective_model, have_disentangled, seedname, stdout, &
                    timer, error, comm)
      if (allocated(error)) return
      call get_AA_R(pw90_berry, dis_manifold, kmesh_info, kpt_latt, print_output, AA_R, HH_R, &
                    v_matrix, eigval, wigner_seitz%irvec, wigner_seitz%nrpts, num_bands, num_kpts, &
                    num_wann, effective_model, have_disentangled, seedname, stdout, timer, error, &
                    comm)
      if (allocated(error)) return
      imf_list = 0.0_dp
      adpt_counter_list = 0
    endif

    if (eval_morb) then
      call get_HH_R(dis_manifold, kpt_latt, print_output, wigner_seitz, HH_R, u_matrix, v_matrix, &
                    eigval, real_lattice, scissors_shift, num_bands, num_kpts, num_wann, &
                    num_valence_bands, effective_model, have_disentangled, seedname, stdout, &
                    timer, error, comm)
      if (allocated(error)) return
      call get_AA_R(pw90_berry, dis_manifold, kmesh_info, kpt_latt, print_output, AA_R, HH_R, &
                    v_matrix, eigval, wigner_seitz%irvec, wigner_seitz%nrpts, num_bands, num_kpts, &
                    num_wann, effective_model, have_disentangled, seedname, stdout, timer, error, &
                    comm)
      if (allocated(error)) return
      call get_BB_R(dis_manifold, kmesh_info, kpt_latt, print_output, BB_R, v_matrix, eigval, &
                    scissors_shift, wigner_seitz%irvec, wigner_seitz%nrpts, num_bands, num_kpts, &
                    num_wann, have_disentangled, seedname, stdout, timer, error, comm)
      if (allocated(error)) return
      call get_CC_R(dis_manifold, kmesh_info, kpt_latt, print_output, pw90_oper_read, CC_R, &
                    v_matrix, eigval, scissors_shift, wigner_seitz%irvec, wigner_seitz%nrpts, &
                    num_bands, num_kpts, num_wann, have_disentangled, seedname, stdout, timer, &
                    error, comm)
      if (allocated(error)) return

      imf_list2 = 0.0_dp
      img_list = 0.0_dp
      imh_list = 0.0_dp
    endif

    ! List here berry_tasks that assume nfermi=1
    !
    not_scannable = eval_kubo .or. (eval_shc .and. pw90_spin_hall%freq_scan)
    if (not_scannable .and. fermi_n .ne. 1) then
      call set_error_input(error, 'The berry_task(s, comm, comm) you chose require that you specify a single ' &
                           //'Fermi energy: scanning the Fermi energy is not implemented', comm)
      return
    endif

    if (eval_kubo) then
      call get_HH_R(dis_manifold, kpt_latt, print_output, wigner_seitz, HH_R, u_matrix, v_matrix, &
                    eigval, real_lattice, scissors_shift, num_bands, num_kpts, num_wann, &
                    num_valence_bands, effective_model, have_disentangled, seedname, stdout, &
                    timer, error, comm)
      if (allocated(error)) return
      call get_AA_R(pw90_berry, dis_manifold, kmesh_info, kpt_latt, print_output, AA_R, HH_R, &
                    v_matrix, eigval, wigner_seitz%irvec, wigner_seitz%nrpts, num_bands, num_kpts, &
                    num_wann, effective_model, have_disentangled, seedname, stdout, timer, &
                    error, comm)
      if (allocated(error)) return
      allocate (kubo_H_k(3, 3, pw90_berry%kubo_nfreq))
      allocate (kubo_H(3, 3, pw90_berry%kubo_nfreq))
      allocate (kubo_AH_k(3, 3, pw90_berry%kubo_nfreq))
      allocate (kubo_AH(3, 3, pw90_berry%kubo_nfreq))
      allocate (jdos_k(pw90_berry%kubo_nfreq))
      allocate (jdos(pw90_berry%kubo_nfreq))
      kubo_H = cmplx_0
      kubo_AH = cmplx_0
      jdos = 0.0_dp
      if (spin_decomp) then

        call get_SS_R(dis_manifold, kpt_latt, print_output, pw90_oper_read, SS_R, v_matrix, &
                      eigval, wigner_seitz%irvec, wigner_seitz%nrpts, num_bands, num_kpts, &
                      num_wann, have_disentangled, seedname, stdout, timer, error, comm)
        if (allocated(error)) return
        allocate (kubo_H_k_spn(3, 3, 3, pw90_berry%kubo_nfreq))
        allocate (kubo_H_spn(3, 3, 3, pw90_berry%kubo_nfreq))
        allocate (kubo_AH_k_spn(3, 3, 3, pw90_berry%kubo_nfreq))
        allocate (kubo_AH_spn(3, 3, 3, pw90_berry%kubo_nfreq))
        allocate (jdos_k_spn(3, pw90_berry%kubo_nfreq))
        allocate (jdos_spn(3, pw90_berry%kubo_nfreq))
        ! fixme, shouldn't we pedantically check these allocs also JJ
        kubo_H_spn = cmplx_0
        kubo_AH_spn = cmplx_0
        jdos_spn = 0.0_dp
      endif
    endif

    if (eval_sc) then
      call get_HH_R(dis_manifold, kpt_latt, print_output, wigner_seitz, HH_R, u_matrix, v_matrix, &
                    eigval, real_lattice, scissors_shift, num_bands, num_kpts, num_wann, &
                    num_valence_bands, effective_model, have_disentangled, seedname, stdout, &
                    timer, error, comm)
      if (allocated(error)) return
      call get_AA_R(pw90_berry, dis_manifold, kmesh_info, kpt_latt, print_output, AA_R, HH_R, &
                    v_matrix, eigval, wigner_seitz%irvec, wigner_seitz%nrpts, num_bands, num_kpts, &
                    num_wann, effective_model, have_disentangled, seedname, stdout, timer, error, &
                    comm)
      if (allocated(error)) return
      allocate (sc_k_list(3, 6, pw90_berry%kubo_nfreq))
      allocate (sc_list(3, 6, pw90_berry%kubo_nfreq))
      sc_k_list = 0.0_dp
      sc_list = 0.0_dp
    endif

    if (eval_shc) then

      call get_HH_R(dis_manifold, kpt_latt, print_output, wigner_seitz, HH_R, u_matrix, v_matrix, &
                    eigval, real_lattice, scissors_shift, num_bands, num_kpts, num_wann, &
                    num_valence_bands, effective_model, have_disentangled, seedname, stdout, &
                    timer, error, comm)
      if (allocated(error)) return
      call get_AA_R(pw90_berry, dis_manifold, kmesh_info, kpt_latt, print_output, AA_R, HH_R, &
                    v_matrix, eigval, wigner_seitz%irvec, wigner_seitz%nrpts, num_bands, num_kpts, &
                    num_wann, effective_model, have_disentangled, seedname, stdout, timer, error, &
                    comm)
      if (allocated(error)) return
      call get_SS_R(dis_manifold, kpt_latt, print_output, pw90_oper_read, SS_R, v_matrix, eigval, &
                    wigner_seitz%irvec, wigner_seitz%nrpts, num_bands, num_kpts, num_wann, &
                    have_disentangled, seedname, stdout, timer, error, comm)
      if (allocated(error)) return

      if (index(pw90_spin_hall%method, 'qiao') > 0) then
        call get_SHC_R(dis_manifold, kmesh_info, kpt_latt, print_output, pw90_oper_read, &
                       pw90_spin_hall, SH_R, SHR_R, SR_R, v_matrix, eigval, scissors_shift, &
                       wigner_seitz%irvec, wigner_seitz%nrpts, num_bands, num_kpts, num_wann, &
                       num_valence_bands, have_disentangled, seedname, stdout, timer, error, comm)
        if (allocated(error)) return
      else
        call get_SAA_R(dis_manifold, kmesh_info, kpt_latt, print_output, SAA_R, v_matrix, &
                       scissors_shift, wigner_seitz%irvec, wigner_seitz%nrpts, num_bands, &
                       num_kpts, num_wann, have_disentangled, seedname, stdout, timer, error, comm)
        if (allocated(error)) return
        call get_SBB_R(dis_manifold, kmesh_info, kpt_latt, print_output, SBB_R, v_matrix, &
                       scissors_shift, wigner_seitz%irvec, wigner_seitz%nrpts, num_bands, &
                       num_kpts, num_wann, have_disentangled, seedname, stdout, timer, error, comm)
        if (allocated(error)) return
      endif

      if (pw90_spin_hall%freq_scan) then
        allocate (shc_freq(pw90_berry%kubo_nfreq))
        allocate (shc_k_freq(pw90_berry%kubo_nfreq))
        shc_freq = 0.0_dp
        shc_k_freq = 0.0_dp
      else
        allocate (shc_fermi(fermi_n))
        allocate (shc_k_fermi(fermi_n))
        allocate (shc_k_fermi_dummy(fermi_n))
        shc_fermi = 0.0_dp
        shc_k_fermi = 0.0_dp
        !only used for fermiscan & adpt kmesh
        shc_k_fermi_dummy = 0.0_dp
        adpt_counter_list = 0
      endif

    endif

    if (eval_kdotp) then
      call get_HH_R(dis_manifold, kpt_latt, print_output, wigner_seitz, HH_R, u_matrix, v_matrix, &
                    eigval, real_lattice, scissors_shift, num_bands, num_kpts, num_wann, &
                    num_valence_bands, effective_model, have_disentangled, seedname, stdout, &
                    timer, error, comm)
      if (allocated(error)) return
      kdotp_nbands = size(pw90_berry%kdotp_bands); 
      allocate (kdotp(kdotp_nbands, kdotp_nbands, 3, 3, 3))
      kdotp = cmplx_0
    endif

    if (print_output%iprint > 0) then

      write (stdout, '(/,/,1x,a)') &
        'Properties calculated in module  b e r r y'
      write (stdout, '(1x,a)') &
        '------------------------------------------'

      if (eval_ahc) write (stdout, '(/,3x,a)') &
        '* Anomalous Hall conductivity'

      if (eval_morb) write (stdout, '(/,3x,a)') '* Orbital magnetization'

      if (eval_kubo) then
        if (spin_decomp) then
          write (stdout, '(/,3x,a)') &
            '* Complex optical conductivity and its spin-decomposition'
          write (stdout, '(/,3x,a)') &
            '* Joint density of states and its spin-decomposition'
        else
          write (stdout, '(/,3x,a)') '* Complex optical conductivity'
          write (stdout, '(/,3x,a)') '* Joint density of states'
        endif
      endif

      if (eval_sc) write (stdout, '(/,3x,a)') &
        '* Shift current'

      if (eval_shc) then
        write (stdout, '(/,3x,a)') '* Spin Hall Conductivity'
        if (index(pw90_spin_hall%method, 'qiao') > 0) then
          write (stdout, '(/,3x,a)') '  Qiao''s SHC (Phys.Rev.B 98.214402)'
        else
          write (stdout, '(/,3x,a)') '  Ryoo''s SHC (Phys.Rev.B 99.235113)'
        endif
        if (pw90_spin_hall%freq_scan) then
          write (stdout, '(/,3x,a)') '  Frequency scan'
        else
          write (stdout, '(/,3x,a)') '  Fermi energy scan'
        endif
      endif

      if (eval_kdotp) write (stdout, '(/,3x,a)') '* k.p expansion coefficients'

      if (pw90_berry%transl_inv) then
        if (eval_morb) then
          call set_error_input(error, 'transl_inv=T disabled for morb', comm)
          return
        endif
        write (stdout, '(/,1x,a)') 'Using a translationally-invariant discretization for the'
        write (stdout, '(1x,a)') 'band-diagonal Wannier matrix elements of r, etc.'
      endif

      if (print_output%timing_level > 1) then
        call io_stopwatch_stop('berry: prelims', timer)
        call io_stopwatch_start('berry: k-interpolation', timer)
      endif

      if (eval_kdotp) then
        ! JJ pw90_berry%kdotp_bands is only allocated on process 0
        ! this causes a segfault at 2895 (accessing nonexistent element zero)
        ! moving to process 0 (on_root) only
        call berry_get_kdotp(kdotp, dis_manifold, kpt_latt, print_output, pw90_berry, &
                             pw90_band_deriv_degen, wannier_data, ws_distance, wigner_seitz, &
                             ws_region, HH_R, u_matrix, v_matrix, eigval, real_lattice, &
                             scissors_shift, mp_grid, num_bands, num_kpts, num_wann, &
                             num_valence_bands, effective_model, have_disentangled, seedname, &
                             stdout, timer, error, comm)
        if (allocated(error)) return
      end if
    end if ! print_output%iprint > 0, aka "on_root"

    ! Set up adaptive refinement mesh
    !
    allocate (adkpt(3, pw90_berry%curv_adpt_kmesh**3), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating adkpt in berry', comm)
      return
    endif
    ikpt = 0
    !
    ! OLD VERSION (only works correctly for odd grids including original point)
    !
    ! do i=-(pw90_berry_curv_adpt_kmesh-1)/2,(pw90_berry_curv_adpt_kmesh-1)/2
    !    do j=-(pw90_berry_curv_adpt_kmesh-1)/2,(pw90_berry_curv_adpt_kmesh-1)/2
    !       do k=-(pw90_berry_curv_adpt_kmesh-1)/2,(pw90_berry_curv_adpt_kmesh-1)/2
    !          ikpt=ikpt+1
    !          adkpt(1,ikpt)=i*db1/pw90_berry_curv_adpt_kmesh
    !          adkpt(2,ikpt)=j*db2/pw90_berry_curv_adpt_kmesh
    !          adkpt(3,ikpt)=k*db3/pw90_berry_curv_adpt_kmesh
    !       end do
    !    end do
    ! end do
    !
    ! NEW VERSION (both even and odd grids)
    !
    do i = 0, pw90_berry%curv_adpt_kmesh - 1
      do j = 0, pw90_berry%curv_adpt_kmesh - 1
        do k = 0, pw90_berry%curv_adpt_kmesh - 1
          ikpt = ikpt + 1
          adkpt(1, ikpt) = db1*((i + 0.5_dp)/pw90_berry%curv_adpt_kmesh - 0.5_dp)
          adkpt(2, ikpt) = db2*((j + 0.5_dp)/pw90_berry%curv_adpt_kmesh - 0.5_dp)
          adkpt(3, ikpt) = db3*((k + 0.5_dp)/pw90_berry%curv_adpt_kmesh - 0.5_dp)
        end do
      end do
    end do

    ! Loop over interpolation k-points
    !
    if (pw90_berry%wanint_kpoint_file) then

      ! NOTE: still need to specify pw90_berry_kmesh in the input file
      !
      !        - Must use the correct nominal value in order to
      !          correctly set up adaptive smearing in kubo

      if (print_output%iprint > 0) write (stdout, '(/,1x,a,i10,a)') &
        'Reading interpolation grid from file kpoint.dat: ', &
        sum(kpoint_dist%num_int_kpts_on_node), ' points'

      ! Loop over k-points on the irreducible wedge of the Brillouin
      ! zone, read from file 'kpoint.dat'
      !
      do loop_xyz = 1, kpoint_dist%num_int_kpts_on_node(my_node_id)
        kpt(:) = kpoint_dist%int_kpts(:, loop_xyz)
        kweight = kpoint_dist%weight(loop_xyz)
        kweight_adpt = kweight/pw90_berry%curv_adpt_kmesh**3
        !               .
        ! ***BEGIN COPY OF CODE BLOCK 1***
        !
        if (eval_ahc) then
          call berry_get_imf_klist(dis_manifold, fermi_energy_list, kpt_latt, ws_region, &
                                   print_output, wannier_data, ws_distance, wigner_seitz, AA_R, &
                                   BB_R, CC_R, HH_R, u_matrix, v_matrix, eigval, kpt, &
                                   real_lattice, imf_k_list, scissors_shift, mp_grid, num_bands, &
                                   num_kpts, num_wann, num_valence_bands, effective_model, &
                                   have_disentangled, seedname, stdout, timer, error, comm)
          if (allocated(error)) return

          ladpt = .false.
          do if = 1, fermi_n
            vdum(1) = sum(imf_k_list(:, 1, if))
            vdum(2) = sum(imf_k_list(:, 2, if))
            vdum(3) = sum(imf_k_list(:, 3, if))
            if (pw90_berry%curv_unit == 'bohr2') vdum = vdum/physics%bohr**2
            rdum = sqrt(dot_product(vdum, vdum))
            if (rdum > pw90_berry%curv_adpt_kmesh_thresh) then
              adpt_counter_list(if) = adpt_counter_list(if) + 1
              ladpt(if) = .true.
            else
              imf_list(:, :, if) = imf_list(:, :, if) + imf_k_list(:, :, if)*kweight
            endif
          enddo
          if (any(ladpt)) then
            do loop_adpt = 1, pw90_berry%curv_adpt_kmesh**3
              ! Using imf_k_list here would corrupt values for other
              ! frequencies, hence dummy. Only if-th element is used
              call berry_get_imf_klist(dis_manifold, fermi_energy_list, kpt_latt, ws_region, &
                                       print_output, wannier_data, ws_distance, wigner_seitz, &
                                       AA_R, BB_R, CC_R, HH_R, u_matrix, v_matrix, eigval, &
                                       kpt(:) + adkpt(:, loop_adpt), real_lattice, &
                                       imf_k_list_dummy, scissors_shift, mp_grid, num_bands, &
                                       num_kpts, num_wann, num_valence_bands, effective_model, &
                                       have_disentangled, seedname, stdout, timer, error, comm, &
                                       ladpt=ladpt)
              if (allocated(error)) return

              do if = 1, fermi_n
                if (ladpt(if)) then
                  imf_list(:, :, if) = imf_list(:, :, if) &
                                       + imf_k_list_dummy(:, :, if)*kweight_adpt
                endif
              enddo
            end do
          endif
        end if

        if (eval_morb) then
          call berry_get_imfgh_klist(dis_manifold, fermi_energy_list, kpt_latt, ws_region, &
                                     print_output, wannier_data, ws_distance, wigner_seitz, AA_R, &
                                     BB_R, CC_R, HH_R, u_matrix, v_matrix, eigval, kpt, &
                                     real_lattice, scissors_shift, mp_grid, fermi_n, num_bands, &
                                     num_kpts, num_wann, num_valence_bands, effective_model, &
                                     have_disentangled, seedname, stdout, timer, error, comm, &
                                     imf_k_list, img_k_list, imh_k_list)
          if (allocated(error)) return

          imf_list2 = imf_list2 + imf_k_list*kweight
          img_list = img_list + img_k_list*kweight
          imh_list = imh_list + imh_k_List*kweight
        endif

        if (eval_kubo) then
          if (spin_decomp) then
            call berry_get_kubo_k(pw90_berry, dis_manifold, fermi_energy_list, kpt_latt, &
                                  pw90_band_deriv_degen, pw90_spin, ws_region, print_output, &
                                  wannier_data, ws_distance, wigner_seitz, AA_R, HH_R, kubo_AH_k, &
                                  kubo_H_k, SS_R, u_matrix, v_matrix, eigval, kpt, real_lattice, &
                                  jdos_k, scissors_shift, mp_grid, num_bands, num_kpts, num_wann, &
                                  num_valence_bands, effective_model, have_disentangled, &
                                  spin_decomp, seedname, stdout, timer, error, comm, &
                                  kubo_AH_k_spn, kubo_H_k_spn, jdos_k_spn)
            if (allocated(error)) return
          else
            call berry_get_kubo_k(pw90_berry, dis_manifold, fermi_energy_list, kpt_latt, &
                                  pw90_band_deriv_degen, pw90_spin, ws_region, print_output, &
                                  wannier_data, ws_distance, wigner_seitz, AA_R, HH_R, kubo_AH_k, &
                                  kubo_H_k, SS_R, u_matrix, v_matrix, eigval, kpt, real_lattice, &
                                  jdos_k, scissors_shift, mp_grid, num_bands, num_kpts, num_wann, &
                                  num_valence_bands, effective_model, have_disentangled, &
                                  spin_decomp, seedname, stdout, timer, error, comm)
            if (allocated(error)) return
          endif
          kubo_H = kubo_H + kubo_H_k*kweight
          kubo_AH = kubo_AH + kubo_AH_k*kweight
          jdos = jdos + jdos_k*kweight
          if (spin_decomp) then
            kubo_H_spn = kubo_H_spn + kubo_H_k_spn*kweight
            kubo_AH_spn = kubo_AH_spn + kubo_AH_k_spn*kweight
            jdos_spn = jdos_spn + jdos_k_spn*kweight
          endif
        endif

        if (eval_sc) then
          call berry_get_sc_klist(pw90_berry, dis_manifold, fermi_energy_list, kmesh_info, &
                                  kpt_latt, ws_region, print_output, pw90_band_deriv_degen, &
                                  wannier_data, ws_distance, wigner_seitz, AA_R, HH_R, u_matrix, &
                                  v_matrix, eigval, kpt, real_lattice, sc_k_list, scissors_shift, &
                                  mp_grid, num_bands, num_kpts, num_wann, num_valence_bands, &
                                  effective_model, have_disentangled, seedname, stdout, timer, &
                                  error, comm)
          if (allocated(error)) return
          sc_list = sc_list + sc_k_list*kweight
        end if

        ! ***END COPY OF CODE BLOCK 1***

        if (eval_shc) then
          ! print calculation progress, from 0%, 10%, ... to 100%
          ! Note the 1st call to berry_get_shc_klist will be much longer
          ! than later calls due to the time spent on
          !   berry_get_shc_klist -> wham_get_eig_deleig ->
          !   pw90common_fourier_R_to_k -> ws_translate_dist
          if (print_output%iprint > 0) then
            call berry_print_progress(kpoint_dist%num_int_kpts_on_node(my_node_id), loop_xyz, &
                                      1, 1, stdout)
          endif
          if (.not. pw90_spin_hall%freq_scan) then
            call berry_get_shc_klist(pw90_berry, dis_manifold, fermi_energy_list, kpt_latt, &
                                     pw90_band_deriv_degen, ws_region, pw90_spin_hall, &
                                     print_output, wannier_data, ws_distance, wigner_seitz, AA_R, &
                                     HH_R, SH_R, SHR_R, SR_R, SS_R, SAA_R, SBB_R, u_matrix, &
                                     v_matrix, eigval, kpt, real_lattice, scissors_shift, mp_grid, &
                                     fermi_n, num_bands, num_kpts, num_wann, num_valence_bands, &
                                     effective_model, have_disentangled, seedname, stdout, timer, &
                                     error, comm, shc_k_fermi=shc_k_fermi)
            if (allocated(error)) return

            !check whether needs to tigger adpt kmesh or not.
            !Since the calculated shc_k at one Fermi energy can be reused
            !by all the Fermi energies, if we find out that at a specific
            !Fermi energy shc_k(if) > thresh, then we will update shc_k at
            !all the Fermi energies as well.
            !This also avoids repeated calculation if shc_k(if) > thresh
            !is satisfied at more than one Fermi energy.
            ladpt_kmesh = .false.
            !if adpt_kmesh==1, no need to calculate on the same kpt again.
            !This happens if adpt_kmesh==1 while adpt_kmesh_thresh is low.
            if (pw90_berry%curv_adpt_kmesh > 1) then
              do if = 1, fermi_n
                rdum = abs(shc_k_fermi(if))
                if (pw90_berry%curv_unit == 'bohr2') rdum = rdum/physics%bohr**2
                if (rdum > pw90_berry%curv_adpt_kmesh_thresh) then
                  adpt_counter_list(1) = adpt_counter_list(1) + 1
                  ladpt_kmesh = .true.
                  exit
                endif
              enddo
            else
              ladpt_kmesh = .false.
            end if
            if (ladpt_kmesh) then
              do loop_adpt = 1, pw90_berry%curv_adpt_kmesh**3
                !Using shc_k here would corrupt values for other
                !kpt, hence dummy. Only if-th element is used.
                call berry_get_shc_klist(pw90_berry, dis_manifold, fermi_energy_list, kpt_latt, &
                                         pw90_band_deriv_degen, ws_region, pw90_spin_hall, &
                                         print_output, wannier_data, ws_distance, wigner_seitz, &
                                         AA_R, HH_R, SH_R, SHR_R, SR_R, SS_R, SAA_R, SBB_R, &
                                         u_matrix, v_matrix, eigval, kpt(:) + adkpt(:, loop_adpt), &
                                         real_lattice, scissors_shift, mp_grid, fermi_n, &
                                         num_bands, num_kpts, num_wann, num_valence_bands, &
                                         effective_model, have_disentangled, seedname, stdout, &
                                         timer, error, comm, shc_k_fermi=shc_k_fermi_dummy)
                if (allocated(error)) return

                shc_fermi = shc_fermi + kweight_adpt*shc_k_fermi_dummy
              end do
            else
              shc_fermi = shc_fermi + kweight*shc_k_fermi
            end if
          else ! freq_scan, no adaptive kmesh
            call berry_get_shc_klist(pw90_berry, dis_manifold, fermi_energy_list, kpt_latt, &
                                     pw90_band_deriv_degen, ws_region, pw90_spin_hall, &
                                     print_output, wannier_data, ws_distance, wigner_seitz, AA_R, &
                                     HH_R, SH_R, SHR_R, SR_R, SS_R, SAA_R, SBB_R, u_matrix, &
                                     v_matrix, eigval, kpt, real_lattice, scissors_shift, mp_grid, &
                                     fermi_n, num_bands, num_kpts, num_wann, num_valence_bands, &
                                     effective_model, have_disentangled, seedname, stdout, timer, &
                                     error, comm, shc_k_freq=shc_k_freq)
            if (allocated(error)) return

            shc_freq = shc_freq + kweight*shc_k_freq
          end if
        end if

      end do !loop_xyz

    else! Do not read 'kpoint.dat'. Loop over a regular grid in the full BZ

      kweight = db1*db2*db3
      kweight_adpt = kweight/pw90_berry%curv_adpt_kmesh**3

      do loop_xyz = my_node_id, PRODUCT(pw90_berry%kmesh%mesh) - 1, num_nodes
        loop_x = loop_xyz/(pw90_berry%kmesh%mesh(2)*pw90_berry%kmesh%mesh(3))
        loop_y = (loop_xyz - loop_x*(pw90_berry%kmesh%mesh(2) &
                                     *pw90_berry%kmesh%mesh(3)))/pw90_berry%kmesh%mesh(3)
        loop_z = loop_xyz - loop_x*(pw90_berry%kmesh%mesh(2)*pw90_berry%kmesh%mesh(3)) &
                 - loop_y*pw90_berry%kmesh%mesh(3)
        kpt(1) = loop_x*db1
        kpt(2) = loop_y*db2
        kpt(3) = loop_z*db3

        ! ***BEGIN CODE BLOCK 1***
        if (eval_ahc) then

          call berry_get_imf_klist(dis_manifold, fermi_energy_list, kpt_latt, ws_region, &
                                   print_output, wannier_data, ws_distance, wigner_seitz, AA_R, &
                                   BB_R, CC_R, HH_R, u_matrix, v_matrix, eigval, kpt, &
                                   real_lattice, imf_k_list, scissors_shift, mp_grid, num_bands, &
                                   num_kpts, num_wann, num_valence_bands, effective_model, &
                                   have_disentangled, seedname, stdout, timer, error, comm)
          if (allocated(error)) return

          ladpt = .false.
          do if = 1, fermi_n
            vdum(1) = sum(imf_k_list(:, 1, if))
            vdum(2) = sum(imf_k_list(:, 2, if))
            vdum(3) = sum(imf_k_list(:, 3, if))
            if (pw90_berry%curv_unit == 'bohr2') vdum = vdum/physics%bohr**2
            rdum = sqrt(dot_product(vdum, vdum))
            if (rdum > pw90_berry%curv_adpt_kmesh_thresh) then
              adpt_counter_list(if) = adpt_counter_list(if) + 1
              ladpt(if) = .true.
            else
              imf_list(:, :, if) = imf_list(:, :, if) + imf_k_list(:, :, if)*kweight
            endif
          enddo
          if (any(ladpt)) then
            do loop_adpt = 1, pw90_berry%curv_adpt_kmesh**3
              ! Using imf_k_list here would corrupt values for other
              ! frequencies, hence dummy. Only if-th element is used
              call berry_get_imf_klist(dis_manifold, fermi_energy_list, kpt_latt, ws_region, &
                                       print_output, wannier_data, ws_distance, wigner_seitz, &
                                       AA_R, BB_R, CC_R, HH_R, u_matrix, v_matrix, eigval, &
                                       kpt(:) + adkpt(:, loop_adpt), real_lattice, &
                                       imf_k_list_dummy, scissors_shift, mp_grid, num_bands, &
                                       num_kpts, num_wann, num_valence_bands, effective_model, &
                                       have_disentangled, seedname, stdout, timer, error, comm, &
                                       ladpt=ladpt)
              if (allocated(error)) return

              do if = 1, fermi_n
                if (ladpt(if)) then
                  imf_list(:, :, if) = imf_list(:, :, if) &
                                       + imf_k_list_dummy(:, :, if)*kweight_adpt
                endif
              enddo
            end do
          endif
        end if

        if (eval_morb) then
          call berry_get_imfgh_klist(dis_manifold, fermi_energy_list, kpt_latt, ws_region, &
                                     print_output, wannier_data, ws_distance, wigner_seitz, AA_R, &
                                     BB_R, CC_R, HH_R, u_matrix, v_matrix, eigval, kpt, &
                                     real_lattice, scissors_shift, mp_grid, fermi_n, num_bands, &
                                     num_kpts, num_wann, num_valence_bands, effective_model, &
                                     have_disentangled, seedname, stdout, timer, error, comm, &
                                     imf_k_list, img_k_list, imh_k_list)
          if (allocated(error)) return

          imf_list2 = imf_list2 + imf_k_list*kweight
          img_list = img_list + img_k_list*kweight
          imh_list = imh_list + imh_k_List*kweight
        endif

        if (eval_kubo) then
          if (spin_decomp) then
            call berry_get_kubo_k(pw90_berry, dis_manifold, fermi_energy_list, kpt_latt, &
                                  pw90_band_deriv_degen, pw90_spin, ws_region, print_output, &
                                  wannier_data, ws_distance, wigner_seitz, AA_R, HH_R, kubo_AH_k, &
                                  kubo_H_k, SS_R, u_matrix, v_matrix, eigval, kpt, real_lattice, &
                                  jdos_k, scissors_shift, mp_grid, num_bands, num_kpts, num_wann, &
                                  num_valence_bands, effective_model, have_disentangled, &
                                  spin_decomp, seedname, stdout, timer, error, comm, &
                                  kubo_AH_k_spn, kubo_H_k_spn, jdos_k_spn)
            if (allocated(error)) return

          else
            call berry_get_kubo_k(pw90_berry, dis_manifold, fermi_energy_list, kpt_latt, &
                                  pw90_band_deriv_degen, pw90_spin, ws_region, print_output, &
                                  wannier_data, ws_distance, wigner_seitz, AA_R, HH_R, kubo_AH_k, &
                                  kubo_H_k, SS_R, u_matrix, v_matrix, eigval, kpt, real_lattice, &
                                  jdos_k, scissors_shift, mp_grid, num_bands, num_kpts, num_wann, &
                                  num_valence_bands, effective_model, have_disentangled, &
                                  spin_decomp, seedname, stdout, timer, error, comm)
            if (allocated(error)) return

          endif
          kubo_H = kubo_H + kubo_H_k*kweight
          kubo_AH = kubo_AH + kubo_AH_k*kweight
          jdos = jdos + jdos_k*kweight
          if (spin_decomp) then
            kubo_H_spn = kubo_H_spn + kubo_H_k_spn*kweight
            kubo_AH_spn = kubo_AH_spn + kubo_AH_k_spn*kweight
            jdos_spn = jdos_spn + jdos_k_spn*kweight
          endif
        endif

        if (eval_sc) then
          call berry_get_sc_klist(pw90_berry, dis_manifold, fermi_energy_list, kmesh_info, &
                                  kpt_latt, ws_region, print_output, pw90_band_deriv_degen, &
                                  wannier_data, ws_distance, wigner_seitz, AA_R, HH_R, u_matrix, &
                                  v_matrix, eigval, kpt, real_lattice, sc_k_list, scissors_shift, &
                                  mp_grid, num_bands, num_kpts, num_wann, num_valence_bands, &
                                  effective_model, have_disentangled, seedname, stdout, timer, &
                                  error, comm)
          if (allocated(error)) return

          sc_list = sc_list + sc_k_list*kweight
        end if

        ! ***END CODE BLOCK 1***

        if (eval_shc) then
          ! print calculation progress, from 0%, 10%, ... to 100%
          ! Note the 1st call to berry_get_shc_klist will be much longer
          ! than later calls due to the time spent on
          !   berry_get_shc_klist -> wham_get_eig_deleig ->
          !   pw90common_fourier_R_to_k -> ws_translate_dist
          if (print_output%iprint > 0) then
            call berry_print_progress(PRODUCT(pw90_berry%kmesh%mesh) - 1, loop_xyz, my_node_id, &
                                      num_nodes, stdout)
          endif
          if (.not. pw90_spin_hall%freq_scan) then
            call berry_get_shc_klist(pw90_berry, dis_manifold, fermi_energy_list, kpt_latt, &
                                     pw90_band_deriv_degen, ws_region, pw90_spin_hall, &
                                     print_output, wannier_data, ws_distance, wigner_seitz, AA_R, &
                                     HH_R, SH_R, SHR_R, SR_R, SS_R, SAA_R, SBB_R, u_matrix, &
                                     v_matrix, eigval, kpt, real_lattice, scissors_shift, mp_grid, &
                                     fermi_n, num_bands, num_kpts, num_wann, num_valence_bands, &
                                     effective_model, have_disentangled, seedname, stdout, timer, &
                                     error, comm, shc_k_fermi=shc_k_fermi)
            if (allocated(error)) return

            !check whether needs to tigger adpt kmesh or not.
            !Since the calculated shc_k at one Fermi energy can be reused
            !by all the Fermi energies, if we find out that at a specific
            !Fermi energy shc_k(if) > thresh, then we will update shc_k at
            !all the Fermi energies as well.
            !This also avoids repeated calculation if shc_k(if) > thresh
            !is satisfied at more than one Fermi energy.
            ladpt_kmesh = .false.
            !if adpt_kmesh==1, no need to calculate on the same kpt again.
            !This happens if adpt_kmesh==1 while adpt_kmesh_thresh is low.
            if (pw90_berry%curv_adpt_kmesh > 1) then
              do if = 1, fermi_n
                rdum = abs(shc_k_fermi(if))
                if (pw90_berry%curv_unit == 'bohr2') rdum = rdum/physics%bohr**2
                if (rdum > pw90_berry%curv_adpt_kmesh_thresh) then
                  adpt_counter_list(1) = adpt_counter_list(1) + 1
                  ladpt_kmesh = .true.
                  exit
                endif
              enddo
            else
              ladpt_kmesh = .false.
            end if
            if (ladpt_kmesh) then
              do loop_adpt = 1, pw90_berry%curv_adpt_kmesh**3
                !Using shc_k here would corrupt values for other
                !kpt, hence dummy. Only if-th element is used.
                call berry_get_shc_klist(pw90_berry, dis_manifold, fermi_energy_list, kpt_latt, &
                                         pw90_band_deriv_degen, ws_region, pw90_spin_hall, &
                                         print_output, wannier_data, ws_distance, wigner_seitz, &
                                         AA_R, HH_R, SH_R, SHR_R, SR_R, SS_R, SAA_R, SBB_R, &
                                         u_matrix, v_matrix, eigval, kpt(:) + adkpt(:, loop_adpt), &
                                         real_lattice, scissors_shift, mp_grid, fermi_n, &
                                         num_bands, num_kpts, num_wann, num_valence_bands, &
                                         effective_model, have_disentangled, seedname, stdout, &
                                         timer, error, comm, shc_k_fermi=shc_k_fermi_dummy)
                if (allocated(error)) return

                shc_fermi = shc_fermi + kweight_adpt*shc_k_fermi_dummy
              end do
            else
              shc_fermi = shc_fermi + kweight*shc_k_fermi
            end if
          else ! freq_scan, no adaptive kmesh
            call berry_get_shc_klist(pw90_berry, dis_manifold, fermi_energy_list, kpt_latt, &
                                     pw90_band_deriv_degen, ws_region, pw90_spin_hall, &
                                     print_output, wannier_data, ws_distance, wigner_seitz, AA_R, &
                                     HH_R, SH_R, SHR_R, SR_R, SS_R, SAA_R, SBB_R, u_matrix, &
                                     v_matrix, eigval, kpt, real_lattice, scissors_shift, mp_grid, &
                                     fermi_n, num_bands, num_kpts, num_wann, num_valence_bands, &
                                     effective_model, have_disentangled, seedname, stdout, timer, &
                                     error, comm, shc_k_freq=shc_k_freq)
            if (allocated(error)) return

            shc_freq = shc_freq + kweight*shc_k_freq
          end if
        end if

      end do !loop_xyz

    end if !wanint_kpoint_file

    ! Collect contributions from all nodes
    if (eval_ahc) then
      call comms_reduce(imf_list(1, 1, 1), 3*3*fermi_n, 'SUM', error, comm)
      if (allocated(error)) return
      call comms_reduce(adpt_counter_list(1), fermi_n, 'SUM', error, comm)
      if (allocated(error)) return
    endif

    if (eval_morb) then
      call comms_reduce(imf_list2(1, 1, 1), 3*3*fermi_n, 'SUM', error, comm)
      if (allocated(error)) return
      call comms_reduce(img_list(1, 1, 1), 3*3*fermi_n, 'SUM', error, comm)
      if (allocated(error)) return
      call comms_reduce(imh_list(1, 1, 1), 3*3*fermi_n, 'SUM', error, comm)
      if (allocated(error)) return
    end if

    if (eval_kubo) then
      call comms_reduce(kubo_H(1, 1, 1), 3*3*pw90_berry%kubo_nfreq, 'SUM', error, comm)
      if (allocated(error)) return
      call comms_reduce(kubo_AH(1, 1, 1), 3*3*pw90_berry%kubo_nfreq, 'SUM', error, comm)
      if (allocated(error)) return
      call comms_reduce(jdos(1), pw90_berry%kubo_nfreq, 'SUM', error, comm)
      if (allocated(error)) return
      if (spin_decomp) then
        call comms_reduce(kubo_H_spn(1, 1, 1, 1), 3*3*3*pw90_berry%kubo_nfreq, 'SUM', error, comm)
        if (allocated(error)) return
        call comms_reduce(kubo_AH_spn(1, 1, 1, 1), 3*3*3*pw90_berry%kubo_nfreq, 'SUM', error, comm)
        if (allocated(error)) return
        call comms_reduce(jdos_spn(1, 1), 3*pw90_berry%kubo_nfreq, 'SUM', error, comm)
        if (allocated(error)) return
      endif
    endif

    if (eval_sc) then
      call comms_reduce(sc_list(1, 1, 1), 3*6*pw90_berry%kubo_nfreq, 'SUM', error, comm)
      if (allocated(error)) return
    end if

    if (eval_shc) then
      if (pw90_spin_hall%freq_scan) then
        call comms_reduce(shc_freq(1), pw90_berry%kubo_nfreq, 'SUM', error, comm)
        if (allocated(error)) return
      else
        call comms_reduce(shc_fermi(1), fermi_n, 'SUM', error, comm)
        if (allocated(error)) return
        call comms_reduce(adpt_counter_list(1), fermi_n, 'SUM', error, comm)
        if (allocated(error)) return
      end if
    end if

    if (print_output%iprint > 0) then

      if (print_output%timing_level > 1) call io_stopwatch_stop('berry: k-interpolation', timer)
      write (stdout, '(1x,a)') ' '
      if (eval_ahc .and. pw90_berry%curv_adpt_kmesh .ne. 1) then
        if (.not. pw90_berry%wanint_kpoint_file) write (stdout, '(1x,a28,3(i0,1x))') &
          'Regular interpolation grid: ', pw90_berry%kmesh%mesh
        write (stdout, '(1x,a28,3(i0,1x))') 'Adaptive refinement grid: ', &
          pw90_berry%curv_adpt_kmesh, pw90_berry%curv_adpt_kmesh, pw90_berry%curv_adpt_kmesh
        if (pw90_berry%curv_unit == 'ang2') then
          write (stdout, '(1x,a28,a17,f6.2,a)') &
            'Refinement threshold: ', 'Berry curvature >', &
            pw90_berry%curv_adpt_kmesh_thresh, ' Ang^2'
        elseif (pw90_berry%curv_unit == 'bohr2') then
          write (stdout, '(1x,a28,a17,f6.2,a)') &
            'Refinement threshold: ', 'Berry curvature >', &
            pw90_berry%curv_adpt_kmesh_thresh, ' bohr^2'
        endif
        if (fermi_n == 1) then
          if (pw90_berry%wanint_kpoint_file) then
            write (stdout, '(1x,a30,i5,a,f5.2,a)') &
              ' Points triggering refinement: ', &
              adpt_counter_list(1), '(', &
              100*real(adpt_counter_list(1), dp) &
              /sum(kpoint_dist%num_int_kpts_on_node), '%)'
          else
            write (stdout, '(1x,a30,i5,a,f5.2,a)') &
              ' Points triggering refinement: ', &
              adpt_counter_list(1), '(', &
              100*real(adpt_counter_list(1), dp)/product(pw90_berry%kmesh%mesh), '%)'
          endif
        endif
      elseif (eval_shc) then
        if (pw90_berry%curv_adpt_kmesh .ne. 1) then
          if (.not. pw90_berry%wanint_kpoint_file) write (stdout, '(1x,a28,3(i0,1x))') &
            'Regular interpolation grid: ', pw90_berry%kmesh%mesh
          if (.not. pw90_spin_hall%freq_scan) then
            write (stdout, '(1x,a28,3(i0,1x))') &
              'Adaptive refinement grid: ', &
              pw90_berry%curv_adpt_kmesh, pw90_berry%curv_adpt_kmesh, pw90_berry%curv_adpt_kmesh
            if (pw90_berry%curv_unit == 'ang2') then
              write (stdout, '(1x,a28,f12.2,a)') &
                'Refinement threshold: ', &
                pw90_berry%curv_adpt_kmesh_thresh, ' Ang^2'
            elseif (pw90_berry%curv_unit == 'bohr2') then
              write (stdout, '(1x,a28,f12.2,a)') &
                'Refinement threshold: ', &
                pw90_berry%curv_adpt_kmesh_thresh, ' bohr^2'
            endif
            if (pw90_berry%wanint_kpoint_file) then
              write (stdout, '(1x,a30,i8,a,f6.2,a)') &
                ' Points triggering refinement: ', adpt_counter_list(1), '(', &
                100*real(adpt_counter_list(1), dp)/sum(kpoint_dist%num_int_kpts_on_node), '%)'
            else
              write (stdout, '(1x,a30,i8,a,f6.2,a)') &
                ' Points triggering refinement: ', adpt_counter_list(1), '(', &
                100*real(adpt_counter_list(1), dp)/product(pw90_berry%kmesh%mesh), '%)'
            endif
          endif
        else
          if (.not. pw90_berry%wanint_kpoint_file) write (stdout, &
                                                          '(1x,a20,3(i0,1x))') 'Interpolation grid: ', pw90_berry%kmesh%mesh(1:3)
        endif
        write (stdout, '(a)') ''
        if (pw90_berry%kubo_smearing%use_adaptive) then
          write (stdout, '(1x,a)') 'Using adaptive smearing'
          write (stdout, '(7x,a,f8.3)') 'adaptive smearing prefactor ', &
            pw90_berry%kubo_smearing%adaptive_prefactor
          write (stdout, '(7x,a,f8.3,a)') 'adaptive smearing max width ', &
            pw90_berry%kubo_smearing%adaptive_max_width, ' eV'
        else
          write (stdout, '(1x,a)') 'Using fixed smearing'
          write (stdout, '(7x,a,f8.3,a)') 'fixed smearing width ', &
            pw90_berry%kubo_smearing%fixed_width, ' eV'
        endif
        write (stdout, '(a)') ''
        if (abs(scissors_shift) > 1.0e-7_dp) then
          write (stdout, '(1X,A,I0,A,G18.10,A)') "Using scissors_shift to shift energy bands with index > ", &
            num_valence_bands, " by ", scissors_shift, " eV."
        endif
        if (pw90_spin_hall%bandshift) then
          write (stdout, '(1X,A,I0,A,G18.10,A)') "Using shc_bandshift to shift energy bands with index >= ", &
            pw90_spin_hall%bandshift_firstband, " by ", pw90_spin_hall%bandshift_energyshift, " eV."
        endif
      else
        if (.not. pw90_berry%wanint_kpoint_file) write (stdout, &
                                                        '(1x,a20,3(i0,1x))') 'Interpolation grid: ', pw90_berry%kmesh%mesh(1:3)
      endif

      if (eval_ahc) then
        !
        ! --------------------------------------------------------------------
        ! At this point imf contains
        !
        ! (1/N) sum_k Omega_{alpha beta}(k),
        !
        ! an approximation to
        !
        ! V_c.int dk/(2.pi)^3 Omega_{alpha beta}(k) dk
        !
        ! (V_c is the cell volume). We want
        !
        ! sigma_{alpha beta}=-(e^2/hbar) int dk/(2.pi)^3 Omega(k) dk
        !
        ! Hence need to multiply by -(e^2/hbar.V_c).
        ! To get a conductivity in units of S/cm,
        !
        ! (i)   Divide by V_c to obtain (1/N) sum_k omega(k)/V_c, with units
        !       of [L]^{-1} (Berry curvature Omega(k) has units of [L]^2)
        ! (ii)  [L] = Angstrom. Multiply by 10^8 to convert to (cm)^{-1}
        ! (iii) Multiply by -e^2/hbar in SI, with has units ofconductance,
        !       (Ohm)^{-1}, or Siemens (S), to get the final result in S/cm
        !
        !==================================================
        ! fac = -e^2/(hbar.V_c*10^-8)
        !==================================================
        !
        ! with 'V_c' in Angstroms^3, and 'e', 'hbar' in SI units
        ! --------------------------------------------------------------------
        !
        fac = -1.0e8_dp*physics%elem_charge_SI**2/(physics%hbar_SI*cell_volume)
        ahc_list(:, :, :) = imf_list(:, :, :)*fac
        if (fermi_n > 1) then
          write (stdout, '(/,1x,a)') &
            '---------------------------------'
          write (stdout, '(1x,a)') &
            'Output data files related to AHC:'
          write (stdout, '(1x,a)') &
            '---------------------------------'
          file_name = trim(seedname)//'-ahc-fermiscan.dat'
          write (stdout, '(/,3x,a)') '* '//file_name
          file_unit = io_file_unit()
          open (file_unit, FILE=file_name, STATUS='UNKNOWN', FORM='FORMATTED')
        endif
        do if = 1, fermi_n
          if (fermi_n > 1) write (file_unit, '(4(F12.6,1x))') &
            fermi_energy_list(if), sum(ahc_list(:, 1, if)), &
            sum(ahc_list(:, 2, if)), sum(ahc_list(:, 3, if))
          write (stdout, '(/,1x,a18,F10.4)') 'Fermi energy (ev):', &
            fermi_energy_list(if)
          if (fermi_n > 1) then
            if (pw90_berry%wanint_kpoint_file) then
              write (stdout, '(1x,a30,i5,a,f5.2,a)') &
                ' Points triggering refinement: ', &
                adpt_counter_list(if), '(', &
                100*real(adpt_counter_list(if), dp) &
                /sum(kpoint_dist%num_int_kpts_on_node), '%)'
            else
              write (stdout, '(1x,a30,i5,a,f5.2,a)') &
                ' Points triggering refinement: ', &
                adpt_counter_list(if), '(', &
                100*real(adpt_counter_list(if), dp) &
                /product(pw90_berry%kmesh%mesh), '%)'
            endif
          endif
          write (stdout, '(/,1x,a)') &
            'AHC (S/cm)       x          y          z'
          if (print_output%iprint > 1) then
            write (stdout, '(1x,a)') &
              '=========='
            write (stdout, '(1x,a9,2x,3(f10.4,1x))') 'J0 term :', &
              ahc_list(1, 1, if), ahc_list(1, 2, if), ahc_list(1, 3, if)
            write (stdout, '(1x,a9,2x,3(f10.4,1x))') 'J1 term :', &
              ahc_list(2, 1, if), ahc_list(2, 2, if), ahc_list(2, 3, if)
            write (stdout, '(1x,a9,2x,3(f10.4,1x))') 'J2 term :', &
              ahc_list(3, 1, if), ahc_list(3, 2, if), ahc_list(3, 3, if)
            write (stdout, '(1x,a)') &
              '-------------------------------------------'
            write (stdout, '(1x,a9,2x,3(f10.4,1x),/)') 'Total   :', &
              sum(ahc_list(:, 1, if)), sum(ahc_list(:, 2, if)), &
              sum(ahc_list(:, 3, if))
          else
            write (stdout, '(1x,a10,1x,3(f10.4,1x),/)') '==========', &
              sum(ahc_list(:, 1, if)), sum(ahc_list(:, 2, if)), &
              sum(ahc_list(:, 3, if))
          endif
        enddo
        if (fermi_n > 1) close (file_unit)
      endif

      if (eval_morb) then
        !
        ! --------------------------------------------------------------------
        ! At this point X=img_ab(:)-fermi_energy*imf_ab(:) and
        !               Y=imh_ab(:)-fermi_energy*imf_ab(:)
        ! contain, eg,
        !
        ! (1/N) sum_k X(k), where X(k)=-2*Im[g(k)-E_F.f(k)]
        !
        ! This is an approximation to
        !
        ! V_c.int dk/(2.pi)^3 X(k) dk
        !
        ! (V_c is the cell volume). We want a magnetic moment per cell,
        ! in units of the Bohr magneton. The magnetization-like quantity is
        !
        ! \tilde{M}^LC=-(e/2.hbar) int dk/(2.pi)^3 X(k) dk
        !
        ! So we take X and
        !
        !  (i)  The summand is an energy in eV times a Berry curvature in
        !       Ang^2. To convert to a.u., divide by 27.2 and by 0.529^2
        !  (ii) Multiply by -(e/2.hbar)=-1/2 in atomic units
        ! (iii) At this point we have a magnetic moment (per cell) in atomic
        !       units. 1 Bohr magneton = 1/2 atomic unit, so need to multiply
        !       by 2 to convert it to Bohr magnetons
        ! --------------------------------------------------------------------
        !
        fac = -physics%eV_au/physics%bohr**2
        if (fermi_n > 1) then
          write (stdout, '(/,1x,a)') &
            '---------------------------------'
          write (stdout, '(1x,a)') &
            'Output data files related to the orbital magnetization:'
          write (stdout, '(1x,a)') &
            '---------------------------------'
          file_name = trim(seedname)//'-morb-fermiscan.dat'
          write (stdout, '(/,3x,a)') '* '//file_name
          file_unit = io_file_unit()
          open (file_unit, FILE=file_name, STATUS='UNKNOWN', FORM='FORMATTED')
        endif
        do if = 1, fermi_n
          LCtil_list(:, :, if) = (img_list(:, :, if) &
                                  - fermi_energy_list(if)*imf_list2(:, :, if))*fac
          ICtil_list(:, :, if) = (imh_list(:, :, if) &
                                  - fermi_energy_list(if)*imf_list2(:, :, if))*fac
          Morb_list(:, :, if) = LCtil_list(:, :, if) + ICtil_list(:, :, if)
          if (fermi_n > 1) write (file_unit, '(4(F12.6,1x))') &
            fermi_energy_list(if), sum(Morb_list(1:3, 1, if)), &
            sum(Morb_list(1:3, 2, if)), sum(Morb_list(1:3, 3, if))
          write (stdout, '(/,/,1x,a,F12.6)') 'Fermi energy (ev) =', &
            fermi_energy_list(if)
          write (stdout, '(/,/,1x,a)') &
            'M_orb (bohr magn/cell)        x          y          z'
          if (print_output%iprint > 1) then
            write (stdout, '(1x,a)') &
              '======================'
            write (stdout, '(1x,a22,2x,3(f10.4,1x))') 'Local circulation :', &
              sum(LCtil_list(1:3, 1, if)), sum(LCtil_list(1:3, 2, if)), &
              sum(LCtil_list(1:3, 3, if))
            write (stdout, '(1x,a22,2x,3(f10.4,1x))') &
              'Itinerant circulation:', &
              sum(ICtil_list(1:3, 1, if)), sum(ICtil_list(1:3, 2, if)), &
              sum(ICtil_list(1:3, 3, if))
            write (stdout, '(1x,a)') &
              '--------------------------------------------------------'
            write (stdout, '(1x,a22,2x,3(f10.4,1x),/)') 'Total   :', &
              sum(Morb_list(1:3, 1, if)), sum(Morb_list(1:3, 2, if)), &
              sum(Morb_list(1:3, 3, if))
          else
            write (stdout, '(1x,a22,2x,3(f10.4,1x),/)') &
              '======================', &
              sum(Morb_list(1:3, 1, if)), sum(Morb_list(1:3, 2, if)), &
              sum(Morb_list(1:3, 3, if))
          endif
        enddo
        if (fermi_n > 1) close (file_unit)
      endif

      ! -----------------------------!
      ! Complex optical conductivity !
      ! -----------------------------!
      !
      if (eval_kubo) then
        !
        ! Convert to S/cm
        fac = 1.0e8_dp*physics%elem_charge_SI**2/(physics%hbar_SI*cell_volume)
        kubo_H = kubo_H*fac
        kubo_AH = kubo_AH*fac
        if (spin_decomp) then
          kubo_H_spn = kubo_H_spn*fac
          kubo_AH_spn = kubo_AH_spn*fac
        endif
        !
        write (stdout, '(/,1x,a)') &
          '----------------------------------------------------------'
        write (stdout, '(1x,a)') &
          'Output data files related to complex optical conductivity:'
        write (stdout, '(1x,a)') &
          '----------------------------------------------------------'
        !
        ! Symmetric: real (imaginary) part is Hermitean (anti-Hermitean)
        !
        do n = 1, 6
          i = alpha_S(n)
          j = beta_S(n)
          file_name = trim(seedname)//'-kubo_S_'// &
                      achar(119 + i)//achar(119 + j)//'.dat'
          file_name = trim(file_name)
          file_unit = io_file_unit()
          write (stdout, '(/,3x,a)') '* '//file_name
          open (file_unit, FILE=file_name, STATUS='UNKNOWN', FORM='FORMATTED')
          do ifreq = 1, pw90_berry%kubo_nfreq
            if (spin_decomp) then
              write (file_unit, '(9E16.8)') real(pw90_berry%kubo_freq_list(ifreq), dp), &
                real(0.5_dp*(kubo_H(i, j, ifreq) + kubo_H(j, i, ifreq)), dp), &
                aimag(0.5_dp*(kubo_AH(i, j, ifreq) + kubo_AH(j, i, ifreq))), &
                real(0.5_dp*(kubo_H_spn(i, j, 1, ifreq) &
                             + kubo_H_spn(j, i, 1, ifreq)), dp), &
                aimag(0.5_dp*(kubo_AH_spn(i, j, 1, ifreq) &
                              + kubo_AH_spn(j, i, 1, ifreq))), &
                real(0.5_dp*(kubo_H_spn(i, j, 2, ifreq) &
                             + kubo_H_spn(j, i, 2, ifreq)), dp), &
                aimag(0.5_dp*(kubo_AH_spn(i, j, 2, ifreq) &
                              + kubo_AH_spn(j, i, 2, ifreq))), &
                real(0.5_dp*(kubo_H_spn(i, j, 3, ifreq) &
                             + kubo_H_spn(j, i, 3, ifreq)), dp), &
                aimag(0.5_dp*(kubo_AH_spn(i, j, 3, ifreq) &
                              + kubo_AH_spn(j, i, 3, ifreq)))
            else
              write (file_unit, '(3E16.8)') real(pw90_berry%kubo_freq_list(ifreq), dp), &
                real(0.5_dp*(kubo_H(i, j, ifreq) + kubo_H(j, i, ifreq)), dp), &
                aimag(0.5_dp*(kubo_AH(i, j, ifreq) + kubo_AH(j, i, ifreq)))
            endif
          enddo
          close (file_unit)
        enddo
        !
        ! Antisymmetric: real (imaginary) part is anti-Hermitean (Hermitean)
        !
        do n = 1, 3
          i = alpha_A(n)
          j = beta_A(n)
          file_name = trim(seedname)//'-kubo_A_'// &
                      achar(119 + i)//achar(119 + j)//'.dat'
          file_name = trim(file_name)
          file_unit = io_file_unit()
          write (stdout, '(/,3x,a)') '* '//file_name
          open (file_unit, FILE=file_name, STATUS='UNKNOWN', FORM='FORMATTED')
          do ifreq = 1, pw90_berry%kubo_nfreq
            if (spin_decomp) then
              write (file_unit, '(9E16.8)') real(pw90_berry%kubo_freq_list(ifreq), dp), &
                real(0.5_dp*(kubo_AH(i, j, ifreq) - kubo_AH(j, i, ifreq)), dp), &
                aimag(0.5_dp*(kubo_H(i, j, ifreq) - kubo_H(j, i, ifreq))), &
                real(0.5_dp*(kubo_AH_spn(i, j, 1, ifreq) &
                             - kubo_AH_spn(j, i, 1, ifreq)), dp), &
                aimag(0.5_dp*(kubo_H_spn(i, j, 1, ifreq) &
                              - kubo_H_spn(j, i, 1, ifreq))), &
                real(0.5_dp*(kubo_AH_spn(i, j, 2, ifreq) &
                             - kubo_AH_spn(j, i, 2, ifreq)), dp), &
                aimag(0.5_dp*(kubo_H_spn(i, j, 2, ifreq) &
                              - kubo_H_spn(j, i, 2, ifreq))), &
                real(0.5_dp*(kubo_AH_spn(i, j, 3, ifreq) &
                             - kubo_AH_spn(j, i, 3, ifreq)), dp), &
                aimag(0.5_dp*(kubo_H_spn(i, j, 3, ifreq) &
                              - kubo_H_spn(j, i, 3, ifreq)))
            else
              write (file_unit, '(3E16.8)') real(pw90_berry%kubo_freq_list(ifreq), dp), &
                real(0.5_dp*(kubo_AH(i, j, ifreq) - kubo_AH(j, i, ifreq)), dp), &
                aimag(0.5_dp*(kubo_H(i, j, ifreq) - kubo_H(j, i, ifreq)))
            endif
          enddo
          close (file_unit)
        enddo
        !
        ! Joint density of states
        !
        file_name = trim(seedname)//'-jdos.dat'
        write (stdout, '(/,3x,a)') '* '//file_name
        file_unit = io_file_unit()
        open (file_unit, FILE=file_name, STATUS='UNKNOWN', FORM='FORMATTED')
        do ifreq = 1, pw90_berry%kubo_nfreq
          if (spin_decomp) then
            write (file_unit, '(5E16.8)') real(pw90_berry%kubo_freq_list(ifreq), dp), &
              jdos(ifreq), jdos_spn(:, ifreq)
          else
            write (file_unit, '(2E16.8)') real(pw90_berry%kubo_freq_list(ifreq), dp), &
              jdos(ifreq)
          endif
        enddo
        close (file_unit)
      endif

      if (eval_sc) then
        ! -----------------------------!
        ! Nonlinear shift current
        ! -----------------------------!

        ! --------------------------------------------------------------------
        ! At this point sc_list contains
        !
        ! (1/N) sum_k (r_^{b}r^{c}_{a}+r_^{c}r^{b}_{a})(k) delta(w),
        !
        ! an approximation to
        !
        ! V_c.int dk/(2.pi)^3 (r_^{b}r^{c}_{a}+r_^{c}r^{b}_{a})(k) delta(w) dk
        !
        ! (V_c is the cell volume). We want
        !
        ! sigma_{abc}=( pi.e^3/(4.hbar^2) ) int dk/(2.pi)^3 Im[ (r_^{b}r^{c}_{a}+r_^{c}r^{b}_{a})(k) delta(w) ] dk
        !
        ! Note factor 1/4 instead of 1/2 as compared to SS PRB 61 5337 (2000) (Eq. 57),
        ! because we introduce 2 delta functions instead of 1.
        ! Hence we need to multiply by  pi.e^3/(4.hbar^2.V_c).
        ! To get the nonlinear response in units of A/V^2,
        !
        ! (i)   Divide by V_c to obtain (1/N) sum_k (r_^{b}r^{c}_{a}+r_^{c}r^{b}_{a})delta(w)/V_c, with units
        !       of [T] (integrand terms r_^{b}r^{c}_{a} delta(w) have units of [T].[L]^3)
        ! (ii)  Multiply by eV_seconds to convert the units of [T] from eV to seconds (coming from delta function)
        ! (iii) Multiply by ( pi.e^3/(4.hbar^2) ) in SI, which multiplied by [T] in seconds from (ii), gives final
        !       units of A/V^2
        !
        !==================================================
        ! fac = eV_seconds.( pi.e^3/(4.hbar^2.V_c) )
        !==================================================
        !
        ! with 'V_c' in Angstroms^3, and 'e', 'hbar' in SI units
        ! --------------------------------------------------------------------

        fac = physics%eV_seconds*pi*physics%elem_charge_SI**3/(4*physics%hbar_SI**(2)*cell_volume)
        write (stdout, '(/,1x,a)') &
          '----------------------------------------------------------'
        write (stdout, '(1x,a)') &
          'Output data files related to shift current:               '
        write (stdout, '(1x,a)') &
          '----------------------------------------------------------'

        do i = 1, 3
          do jk = 1, 6
            j = alpha_S(jk)
            k = beta_S(jk)
            file_name = trim(seedname)//'-sc_'// &
                        achar(119 + i)//achar(119 + j)//achar(119 + k)//'.dat'
            file_name = trim(file_name)
            file_unit = io_file_unit()
            write (stdout, '(/,3x,a)') '* '//file_name
            open (file_unit, FILE=file_name, STATUS='UNKNOWN', FORM='FORMATTED')
            do ifreq = 1, pw90_berry%kubo_nfreq
              write (file_unit, '(2E18.8E3)') real(pw90_berry%kubo_freq_list(ifreq), dp), &
                fac*sc_list(i, jk, ifreq)
            enddo
            close (file_unit)
          enddo
        enddo

      endif

      ! -----------------------!
      ! Spin Hall conductivity !
      ! -----------------------!
      !
      if (eval_shc) then
        !
        ! Convert to the unit: (hbar/e) S/cm
        ! at this point, we need to
        ! (i)   multiply -e^2/hbar/(V*N_k) as in the QZYZ18 Eq.(5),
        !       note 1/N_k has already been applied by the kweight
        ! (ii)  convert charge current to spin current:
        !       divide the result by -e and multiply hbar/2 to
        !       recover the spin current, so the overall
        !       effect is -hbar/2/e
        ! (iii) multiply 1e8 to convert it to the unit S/cm
        ! So, the overall factor is
        !   fac = 1.0e8 * e^2 / hbar / V / 2.0
        ! and the final unit of spin Hall conductivity is (hbar/e)S/cm
        !
        fac = 1.0e8_dp*physics%elem_charge_SI**2/(physics%hbar_SI*cell_volume)/2.0_dp
        if (pw90_spin_hall%freq_scan) then
          shc_freq = shc_freq*fac
        else
          shc_fermi = shc_fermi*fac
        endif
        !
        write (stdout, '(/,1x,a)') &
          '----------------------------------------------------------'
        write (stdout, '(1x,a)') &
          'Output data files related to Spin Hall conductivity:'
        write (stdout, '(1x,a)') &
          '----------------------------------------------------------'
        !
        if (.not. pw90_spin_hall%freq_scan) then
          file_name = trim(seedname)//'-shc-fermiscan'//'.dat'
        else
          file_name = trim(seedname)//'-shc-freqscan'//'.dat'
        endif
        file_name = trim(file_name)
        file_unit = io_file_unit()
        write (stdout, '(/,3x,a)') '* '//file_name
        open (file_unit, FILE=file_name, STATUS='UNKNOWN', FORM='FORMATTED')
        if (.not. pw90_spin_hall%freq_scan) then
          write (file_unit, '(a,3x,a,3x,a)') &
            '#No.', 'Fermi energy(eV)', 'SHC((hbar/e)*S/cm)'
          do n = 1, fermi_n
            write (file_unit, '(I4,1x,F12.6,1x,E17.8)') &
              n, fermi_energy_list(n), shc_fermi(n)
          enddo
        else
          write (file_unit, '(a,3x,a,3x,a,3x,a)') '#No.', 'Frequency(eV)', &
            'Re(sigma)((hbar/e)*S/cm)', 'Im(sigma)((hbar/e)*S/cm)'
          do n = 1, pw90_berry%kubo_nfreq
            write (file_unit, '(I4,1x,F12.6,1x,1x,2(E17.8,1x))') n, &
              real(pw90_berry%kubo_freq_list(n), dp), real(shc_freq(n), dp), aimag(shc_freq(n))
          enddo
        endif
        close (file_unit)

      endif

      if (eval_kdotp) then
        ! -----------------------------!
        ! k.p expansion coefficients
        ! -----------------------------!

        write (stdout, '(/,1x,a)') &
          '----------------------------------------------------------'
        write (stdout, '(1x,a)') &
          'Output data files related to k.p:                         '
        write (stdout, '(1x,a)') &
          '----------------------------------------------------------'
        ! zeroth order in k
        file_name = trim(seedname)//'-kdotp_0.dat'
        file_name = trim(file_name)
        file_unit = io_file_unit()
        write (stdout, '(/,3x,a)') '* '//file_name
        open (file_unit, FILE=file_name, STATUS='UNKNOWN', FORM='FORMATTED')
        write (file_unit, '(2E18.8E3)') kdotp(:, :, 1, 1, 1)
        close (file_unit)

        ! first order in k
        file_name = trim(seedname)//'-kdotp_1.dat'
        write (stdout, '(/,3x,a)') '* '//file_name
        open (file_unit, FILE=file_name, STATUS='UNKNOWN', FORM='FORMATTED')
        do i = 1, 3
          write (file_unit, '(2E18.8E3)') kdotp(:, :, 2, i, 1)
        end do
        close (file_unit)

        ! second order in k
        file_name = trim(seedname)//'-kdotp_2.dat'
        write (stdout, '(/,3x,a)') '* '//file_name
        open (file_unit, FILE=file_name, STATUS='UNKNOWN', FORM='FORMATTED')
        do i = 1, 3
          do j = 1, 3
            write (file_unit, '(2E18.8E3)') kdotp(:, :, 3, i, j)
          end do
        end do
        close (file_unit)

      end if

    end if !print_output%iprint >0, aka "on_root"

  end subroutine berry_main

  !================================================!
  subroutine berry_get_imf_klist(dis_manifold, fermi_energy_list, kpt_latt, ws_region, &
                                 print_output, wannier_data, ws_distance, wigner_seitz, AA_R, &
                                 BB_R, CC_R, HH_R, u_matrix, v_matrix, eigval, kpt, real_lattice, &
                                 imf_k_list, scissors_shift, mp_grid, num_bands, num_kpts, &
                                 num_wann, num_valence_bands, effective_model, have_disentangled, &
                                 seedname, stdout, timer, error, comm, occ, ladpt)
    !================================================!
    !
    !! Calculates the Berry curvature traced over the occupied
    !! states, -2Im[f(k)] [Eq.33 CTVR06, Eq.6 LVTS12] for a list
    !! of Fermi energies, and stores it in axial-vector form
    !
    !================================================!

    use w90_types, only: print_output_type, wannier_data_type, &
      dis_manifold_type, ws_region_type, ws_distance_type, timer_list_type
    use w90_comms, only: w90comm_type
    use w90_postw90_types, only: wigner_seitz_type

    implicit none

    ! arguments
    type(dis_manifold_type), intent(in) :: dis_manifold
    real(kind=dp), allocatable, intent(in) :: fermi_energy_list(:)
    real(kind=dp), intent(in) :: kpt_latt(:, :)
    type(print_output_type), intent(in) :: print_output
    type(ws_region_type), intent(in) :: ws_region
    type(w90comm_type), intent(in) :: comm
    type(wannier_data_type), intent(in) :: wannier_data
    type(wigner_seitz_type), intent(inout) :: wigner_seitz
    type(ws_distance_type), intent(inout) :: ws_distance
    type(timer_list_type), intent(inout) :: timer
    type(w90_error_type), allocatable, intent(out) :: error

    integer, intent(in) :: num_wann, num_bands, num_kpts, num_valence_bands
    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: stdout

    real(kind=dp), intent(in) :: eigval(:, :)
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    real(kind=dp), intent(in) :: kpt(3)
    real(kind=dp), intent(out) :: imf_k_list(:, :, :)
    real(kind=dp), intent(in) :: scissors_shift

    complex(kind=dp), intent(in) :: u_matrix(:, :, :), v_matrix(:, :, :)
    complex(kind=dp), allocatable, intent(inout) :: AA_R(:, :, :, :) ! <0n|r|Rm>
    complex(kind=dp), allocatable, intent(inout) :: BB_R(:, :, :, :) ! <0|H(r-R)|R>
    complex(kind=dp), allocatable, intent(inout) :: CC_R(:, :, :, :, :) ! <0|r_alpha.H(r-R)_beta|R>
    complex(kind=dp), allocatable, intent(inout) :: HH_R(:, :, :) !  <0n|r|Rm>

    character(len=50), intent(in) :: seedname
    logical, intent(in) :: have_disentangled
    logical, intent(in) :: effective_model

    real(kind=dp), intent(in), optional :: occ(:)
    logical, intent(in), optional :: ladpt(:)

    ! local variables
    integer :: fermi_n

    fermi_n = 0
    if (allocated(fermi_energy_list)) fermi_n = size(fermi_energy_list)
    if (present(occ)) then
      call berry_get_imfgh_klist(dis_manifold, fermi_energy_list, kpt_latt, ws_region, print_output, &
                                 wannier_data, ws_distance, wigner_seitz, AA_R, BB_R, CC_R, HH_R, u_matrix, &
                                 v_matrix, eigval, kpt, real_lattice, scissors_shift, mp_grid, &
                                 fermi_n, num_bands, num_kpts, num_wann, num_valence_bands, &
                                 effective_model, have_disentangled, seedname, &
                                 stdout, timer, error, comm, imf_k_list, occ=occ)
      if (allocated(error)) return

    else
      if (present(ladpt)) then
        call berry_get_imfgh_klist(dis_manifold, fermi_energy_list, kpt_latt, ws_region, print_output, &
                                   wannier_data, ws_distance, wigner_seitz, AA_R, BB_R, CC_R, HH_R, &
                                   u_matrix, v_matrix, eigval, kpt, real_lattice, scissors_shift, &
                                   mp_grid, fermi_n, num_bands, num_kpts, num_wann, &
                                   num_valence_bands, effective_model, have_disentangled, &
                                   seedname, stdout, timer, error, comm, imf_k_list, ladpt=ladpt)
        if (allocated(error)) return

      else
        call berry_get_imfgh_klist(dis_manifold, fermi_energy_list, kpt_latt, ws_region, print_output, &
                                   wannier_data, ws_distance, wigner_seitz, AA_R, BB_R, CC_R, HH_R, &
                                   u_matrix, v_matrix, eigval, kpt, real_lattice, scissors_shift, &
                                   mp_grid, fermi_n, num_bands, num_kpts, num_wann, &
                                   num_valence_bands, effective_model, have_disentangled, &
                                   seedname, stdout, timer, error, comm, imf_k_list)
        if (allocated(error)) return

      endif
    endif

  end subroutine berry_get_imf_klist

  !================================================!
  subroutine berry_get_imfgh_klist(dis_manifold, fermi_energy_list, kpt_latt, ws_region, print_output, &
                                   wannier_data, ws_distance, wigner_seitz, AA_R, BB_R, CC_R, HH_R, &
                                   u_matrix, v_matrix, eigval, kpt, real_lattice, &
                                   scissors_shift, mp_grid, fermi_n, num_bands, num_kpts, &
                                   num_wann, num_valence_bands, effective_model, &
                                   have_disentangled, seedname, stdout, timer, error, comm, &
                                   imf_k_list, img_k_list, imh_k_list, occ, ladpt)
    !================================================!
    !
    !! Calculates the three quantities needed for the orbital
    !! magnetization:
    !!
    !! * -2Im[f(k)] [Eq.33 CTVR06, Eq.6 LVTS12]
    !! * -2Im[g(k)] [Eq.34 CTVR06, Eq.7 LVTS12]
    !! * -2Im[h(k)] [Eq.35 CTVR06, Eq.8 LVTS12]
    !! They are calculated together (to reduce the number of
    !! Fourier calls) for a list of Fermi energies, and stored
    !! in axial-vector form.
    !
    ! The two optional output parameters 'imh_k_list' and
    ! 'img_k_list' are only calculated if both of them are
    ! present.
    !
    !================================================!

    use w90_comms, only: w90comm_type, mpirank
    use w90_constants, only: dp, cmplx_i
    use w90_types, only: print_output_type, wannier_data_type, &
      dis_manifold_type, kmesh_info_type, ws_region_type, ws_distance_type, timer_list_type
    use w90_postw90_common, only: pw90common_fourier_R_to_k_vec, pw90common_fourier_R_to_k
    use w90_postw90_types, only: wigner_seitz_type
    use w90_utility, only: utility_re_tr_prod, utility_im_tr_prod, utility_zgemm_new
    use w90_wan_ham, only: wham_get_eig_UU_HH_JJlist, wham_get_occ_mat_list

    implicit none

    ! arguments
    type(dis_manifold_type), intent(in) :: dis_manifold
    real(kind=dp), allocatable, intent(in) :: fermi_energy_list(:)
    real(kind=dp), intent(in) :: kpt_latt(:, :)
    type(print_output_type), intent(in) :: print_output
    type(ws_region_type), intent(in) :: ws_region
    type(w90comm_type), intent(in) :: comm
    type(wannier_data_type), intent(in) :: wannier_data
    type(wigner_seitz_type), intent(inout) :: wigner_seitz
    type(ws_distance_type), intent(inout) :: ws_distance
    type(timer_list_type), intent(inout) :: timer
    type(w90_error_type), allocatable, intent(out) :: error

    integer, intent(in) :: num_wann, num_bands, num_kpts, num_valence_bands, fermi_n
    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: stdout

    real(kind=dp), intent(in) :: kpt(3)
    real(kind=dp), intent(in) :: eigval(:, :)
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    real(kind=dp), intent(in) :: scissors_shift

    complex(kind=dp), intent(in) :: u_matrix(:, :, :), v_matrix(:, :, :)
    complex(kind=dp), allocatable, intent(inout) :: HH_R(:, :, :) !  <0n|r|Rm>
    complex(kind=dp), allocatable, intent(inout) :: AA_R(:, :, :, :) ! <0n|r|Rm>
    complex(kind=dp), allocatable, intent(inout) :: BB_R(:, :, :, :) ! <0|H(r-R)|R>
    complex(kind=dp), allocatable, intent(inout) :: CC_R(:, :, :, :, :) ! <0|r_alpha.H(r-R)_beta|R>

    character(len=50), intent(in) :: seedname
    logical, intent(in) :: have_disentangled

    real(kind=dp), intent(out), optional :: imf_k_list(:, :, :)
    real(kind=dp), intent(out), optional :: img_k_list(:, :, :)
    real(kind=dp), intent(out), optional :: imh_k_list(:, :, :)
    real(kind=dp), intent(in), optional :: occ(:)
    logical, intent(in) :: effective_model
    logical, intent(in), optional :: ladpt(:)

    ! local variables
    complex(kind=dp), allocatable :: HH(:, :)
    complex(kind=dp), allocatable :: UU(:, :)
    complex(kind=dp), allocatable :: f_list(:, :, :)
    complex(kind=dp), allocatable :: g_list(:, :, :)
    complex(kind=dp), allocatable :: AA(:, :, :)
    complex(kind=dp), allocatable :: BB(:, :, :)
    complex(kind=dp), allocatable :: CC(:, :, :, :)
    complex(kind=dp), allocatable :: OOmega(:, :, :)
    complex(kind=dp), allocatable :: JJp_list(:, :, :, :)
    complex(kind=dp), allocatable :: JJm_list(:, :, :, :)
    ! Temporary space for matrix products
    complex(kind=dp), allocatable :: tmp(:, :, :)

    real(kind=dp) :: eig(num_wann)
    real(kind=dp) :: s

    integer :: i, j, ife, nfermi_loc

    logical :: todo(fermi_n)

    if (present(occ)) then
      nfermi_loc = 1
    else
      nfermi_loc = fermi_n
    endif

    if (present(ladpt)) then
      todo = ladpt
    else
      todo = .true.
    endif

    allocate (HH(num_wann, num_wann))
    allocate (UU(num_wann, num_wann))
    allocate (f_list(num_wann, num_wann, nfermi_loc))
    allocate (g_list(num_wann, num_wann, nfermi_loc))
    allocate (JJp_list(num_wann, num_wann, nfermi_loc, 3))
    allocate (JJm_list(num_wann, num_wann, nfermi_loc, 3))
    allocate (AA(num_wann, num_wann, 3))
    allocate (OOmega(num_wann, num_wann, 3))

    ! Gather W-gauge matrix objects
    !

    if (present(occ)) then
      call wham_get_eig_UU_HH_JJlist(dis_manifold, fermi_energy_list, kpt_latt, ws_region, &
                                     print_output, wannier_data, ws_distance, wigner_seitz, HH, &
                                     HH_R, JJm_list, JJp_list, u_matrix, UU, v_matrix, eig, &
                                     eigval, kpt, real_lattice, scissors_shift, mp_grid, &
                                     num_bands, num_kpts, num_wann, num_valence_bands, &
                                     effective_model, have_disentangled, seedname, stdout, timer, &
                                     error, comm, occ=occ)
      if (allocated(error)) return
      call wham_get_occ_mat_list(fermi_energy_list, f_list, g_list, UU, num_wann, error, comm, &
                                 occ=occ)
      if (allocated(error)) return

    else
      call wham_get_eig_UU_HH_JJlist(dis_manifold, fermi_energy_list, kpt_latt, ws_region, &
                                     print_output, wannier_data, ws_distance, wigner_seitz, HH, &
                                     HH_R, JJm_list, JJp_list, u_matrix, UU, v_matrix, eig, &
                                     eigval, kpt, real_lattice, scissors_shift, mp_grid, &
                                     num_bands, num_kpts, num_wann, num_valence_bands, &
                                     effective_model, have_disentangled, seedname, stdout, timer, &
                                     error, comm)
      if (allocated(error)) return
      call wham_get_occ_mat_list(fermi_energy_list, f_list, g_list, UU, num_wann, error, comm, &
                                 eig=eig)
      if (allocated(error)) return

    endif

    call pw90common_fourier_R_to_k_vec(ws_region, wannier_data, ws_distance, wigner_seitz, AA_R, &
                                       kpt, real_lattice, mp_grid, num_wann, error, comm, &
                                       OO_true=AA, OO_pseudo=OOmega)
    if (allocated(error)) return

    if (present(imf_k_list)) then
      ! Trace formula for -2Im[f], Eq.(51) LVTS12
      !
      do ife = 1, nfermi_loc
        if (todo(ife)) then
          do i = 1, 3
            !
            ! J0 term (Omega_bar term of WYSV06)
            imf_k_list(1, i, ife) = &
              utility_re_tr_prod(f_list(:, :, ife), OOmega(:, :, i))
            !
            ! J1 term (DA term of WYSV06)
            imf_k_list(2, i, ife) = -2.0_dp* &
                                    ( &
                                    utility_im_tr_prod(AA(:, :, alpha_A(i)), JJp_list(:, :, ife, beta_A(i))) &
                                    + utility_im_tr_prod(JJm_list(:, :, ife, alpha_A(i)), AA(:, :, beta_A(i))) &
                                    )
            !
            ! J2 term (DD of WYSV06)
            imf_k_list(3, i, ife) = -2.0_dp* &
                                    utility_im_tr_prod(JJm_list(:, :, ife, alpha_A(i)), JJp_list(:, :, ife, beta_A(i)))
          end do
        endif
      end do
    end if

    if (present(img_k_list)) img_k_list = 0.0_dp
    if (present(imh_k_list)) imh_k_list = 0.0_dp

    if (present(img_k_list) .and. present(imh_k_list)) then
      allocate (BB(num_wann, num_wann, 3))
      allocate (CC(num_wann, num_wann, 3, 3))

      allocate (tmp(num_wann, num_wann, 5))
      ! tmp(:,:,1:3) ... not dependent on inner loop variables
      ! tmp(:,:,1) ..... HH . AA(:,:,alpha_A(i))
      ! tmp(:,:,2) ..... LLambda_ij [Eq. (37) LVTS12] expressed as a pseudovector
      ! tmp(:,:,3) ..... HH . OOmega(:,:,i)
      ! tmp(:,:,4:5) ... working matrices for matrix products of inner loop

      call pw90common_fourier_R_to_k_vec(ws_region, wannier_data, ws_distance, wigner_seitz, BB_R, &
                                         kpt, real_lattice, mp_grid, num_wann, error, comm, &
                                         OO_true=BB)
      if (allocated(error)) return

      do j = 1, 3
        do i = 1, j
          call pw90common_fourier_R_to_k(ws_region, wannier_data, ws_distance, wigner_seitz, &
                                         CC(:, :, i, j), CC_R(:, :, :, i, j), kpt, real_lattice, &
                                         mp_grid, 0, num_wann, error, comm)
          if (allocated(error)) return

          CC(:, :, j, i) = conjg(transpose(CC(:, :, i, j)))
        end do
      end do

      ! Trace formula for -2Im[g], Eq.(66) LVTS12
      ! Trace formula for -2Im[h], Eq.(56) LVTS12
      !
      do i = 1, 3
        call utility_zgemm_new(HH, AA(:, :, alpha_A(i)), tmp(:, :, 1))
        call utility_zgemm_new(HH, OOmega(:, :, i), tmp(:, :, 3))
        !
        ! LLambda_ij [Eq. (37) LVTS12] expressed as a pseudovector
        tmp(:, :, 2) = cmplx_i*(CC(:, :, alpha_A(i), beta_A(i)) &
                                - conjg(transpose(CC(:, :, alpha_A(i), beta_A(i)))))

        do ife = 1, nfermi_loc
          !
          ! J0 terms for -2Im[g] and -2Im[h]
          !
          ! tmp(:,:,5) = HH . AA(:,:,alpha_A(i)) . f_list(:,:,ife) . AA(:,:,beta_A(i))
          call utility_zgemm_new(tmp(:, :, 1), f_list(:, :, ife), tmp(:, :, 4))
          call utility_zgemm_new(tmp(:, :, 4), AA(:, :, beta_A(i)), tmp(:, :, 5))

          s = 2.0_dp*utility_im_tr_prod(f_list(:, :, ife), tmp(:, :, 5)); 
          img_k_list(1, i, ife) = utility_re_tr_prod(f_list(:, :, ife), tmp(:, :, 2)) - s
          imh_k_list(1, i, ife) = utility_re_tr_prod(f_list(:, :, ife), tmp(:, :, 3)) + s

          !
          ! J1 terms for -2Im[g] and -2Im[h]
          !
          ! tmp(:,:,1) = HH . AA(:,:,alpha_A(i))
          ! tmp(:,:,4) = HH . JJm_list(:,:,ife,alpha_A(i))
          call utility_zgemm_new(HH, JJm_list(:, :, ife, alpha_A(i)), tmp(:, :, 4))

          img_k_list(2, i, ife) = -2.0_dp* &
                                  ( &
                                  utility_im_tr_prod(JJm_list(:, :, ife, alpha_A(i)), BB(:, :, beta_A(i))) &
                                  - utility_im_tr_prod(JJm_list(:, :, ife, beta_A(i)), BB(:, :, alpha_A(i))) &
                                  )
          imh_k_list(2, i, ife) = -2.0_dp* &
                                  ( &
                                  utility_im_tr_prod(tmp(:, :, 1), JJp_list(:, :, ife, beta_A(i))) &
                                  + utility_im_tr_prod(tmp(:, :, 4), AA(:, :, beta_A(i))) &
                                  )

          !
          ! J2 terms for -2Im[g] and -2Im[h]
          !
          ! tmp(:,:,4) = JJm_list(:,:,ife,alpha_A(i)) . HH
          ! tmp(:,:,5) = HH . JJm_list(:,:,ife,alpha_A(i))
          call utility_zgemm_new(JJm_list(:, :, ife, alpha_A(i)), HH, tmp(:, :, 4))
          call utility_zgemm_new(HH, JJm_list(:, :, ife, alpha_A(i)), tmp(:, :, 5))

          img_k_list(3, i, ife) = -2.0_dp* &
                                  utility_im_tr_prod(tmp(:, :, 4), JJp_list(:, :, ife, beta_A(i)))
          imh_k_list(3, i, ife) = -2.0_dp* &
                                  utility_im_tr_prod(tmp(:, :, 5), JJp_list(:, :, ife, beta_A(i)))
        end do
      end do
      deallocate (tmp)
    end if

  end subroutine berry_get_imfgh_klist

  !================================================!
  !                   PRIVATE PROCEDURES
  !================================================!

  subroutine berry_get_kubo_k(pw90_berry, dis_manifold, fermi_energy_list, kpt_latt, &
                              pw90_band_deriv_degen, pw90_spin, ws_region, print_output, &
                              wannier_data, ws_distance, wigner_seitz, AA_R, HH_R, kubo_AH_k, &
                              kubo_H_k, SS_R, u_matrix, v_matrix, eigval, kpt, real_lattice, &
                              jdos_k, scissors_shift, mp_grid, num_bands, num_kpts, num_wann, &
                              num_valence_bands, effective_model, have_disentangled, spin_decomp, &
                              seedname, stdout, timer, error, comm, kubo_AH_k_spn, kubo_H_k_spn, &
                              jdos_k_spn)
    !================================================!
    !
    !! Contribution from point k to the complex interband optical
    !! conductivity, separated into Hermitian (H) and anti-Hermitian (AH)
    !! parts. Also returns the joint density of states
    !
    !================================================!

    use w90_constants, only: dp, cmplx_0, cmplx_i, pi
    use w90_comms, only: w90comm_type
    use w90_utility, only: utility_diagonalize, utility_rotate, utility_w0gauss, &
      utility_recip_lattice_base
    use w90_types, only: print_output_type, wannier_data_type, &
      dis_manifold_type, ws_region_type, ws_distance_type, timer_list_type
    use w90_postw90_types, only: pw90_berry_mod_type, pw90_spin_mod_type, &
      pw90_band_deriv_degen_type, wigner_seitz_type
    use w90_postw90_common, only: pw90common_get_occ, pw90common_fourier_R_to_k_new, &
      pw90common_fourier_R_to_k_vec, pw90common_kmesh_spacing
    use w90_spin, only: spin_get_nk
    use w90_wan_ham, only: wham_get_D_h, wham_get_eig_deleig

    implicit none

    ! arguments
    type(pw90_berry_mod_type), intent(inout) :: pw90_berry
    type(dis_manifold_type), intent(in) :: dis_manifold
    real(kind=dp), allocatable, intent(in) :: fermi_energy_list(:)
    real(kind=dp), intent(in) :: kpt_latt(:, :)
    type(pw90_band_deriv_degen_type), intent(in) :: pw90_band_deriv_degen
    type(pw90_spin_mod_type), intent(in) :: pw90_spin
    type(print_output_type), intent(in) :: print_output
    type(ws_region_type), intent(in) :: ws_region
    type(w90comm_type), intent(in) :: comm
    type(wannier_data_type), intent(in) :: wannier_data
    type(wigner_seitz_type), intent(inout) :: wigner_seitz
    type(ws_distance_type), intent(inout) :: ws_distance
    type(timer_list_type), intent(inout) :: timer
    type(w90_error_type), allocatable, intent(out) :: error

    integer, intent(in) :: num_wann, num_bands, num_kpts, num_valence_bands
    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: stdout

    real(kind=dp), intent(in) :: kpt(3)
    real(kind=dp), intent(out) :: jdos_k(:)
    real(kind=dp), intent(in) :: eigval(:, :)
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    real(kind=dp), intent(in) :: scissors_shift

    complex(kind=dp), intent(out) :: kubo_H_k(:, :, :)
    complex(kind=dp), intent(out) :: kubo_AH_k(:, :, :)
    complex(kind=dp), intent(in) :: u_matrix(:, :, :), v_matrix(:, :, :)
    complex(kind=dp), allocatable, intent(inout) :: HH_R(:, :, :) !  <0n|r|Rm>
    complex(kind=dp), allocatable, intent(inout) :: AA_R(:, :, :, :) ! <0n|r|Rm>
    complex(kind=dp), allocatable, intent(inout) :: SS_R(:, :, :, :) ! <0n|sigma_x,y,z|Rm>

    character(len=50), intent(in) :: seedname
    logical, intent(in) :: have_disentangled
    logical, intent(in) :: spin_decomp
    logical, intent(in) :: effective_model

    ! Last three arguments should be present iff spin_decomp=T (but
    ! this is not checked: do it?)
    real(kind=dp), optional, intent(out) :: jdos_k_spn(:, :)
    complex(kind=dp), optional, intent(out) :: kubo_AH_k_spn(:, :, :, :)
    complex(kind=dp), optional, intent(out) :: kubo_H_k_spn(:, :, :, :)

    ! local variables
    complex(kind=dp), allocatable :: HH(:, :)
    complex(kind=dp), allocatable :: delHH(:, :, :)
    complex(kind=dp), allocatable :: UU(:, :)
    complex(kind=dp), allocatable :: D_h(:, :, :)
    complex(kind=dp), allocatable :: AA(:, :, :)

    real(kind=dp) :: recip_lattice(3, 3), volume
    ! Adaptive smearing
    !
    real(kind=dp) :: del_eig(num_wann, 3), joint_level_spacing, eta_smr, Delta_k, arg, vdum(3)
    real(kind=dp) :: eig(num_wann), occ(num_wann), delta, rfac1, rfac2, occ_prod, spn_nk(num_wann)
    complex(kind=dp) :: cfac, omega
    integer :: i, j, n, m, ifreq, ispn

    allocate (HH(num_wann, num_wann))
    allocate (delHH(num_wann, num_wann, 3))
    allocate (UU(num_wann, num_wann))
    allocate (D_h(num_wann, num_wann, 3))
    allocate (AA(num_wann, num_wann, 3))

    if (pw90_berry%kubo_smearing%use_adaptive) then
      call wham_get_eig_deleig(dis_manifold, kpt_latt, pw90_band_deriv_degen, ws_region, print_output, wannier_data, &
                               ws_distance, wigner_seitz, delHH, HH, HH_R, u_matrix, UU, v_matrix, &
                               del_eig, eig, eigval, kpt, real_lattice, scissors_shift, mp_grid, &
                               num_bands, num_kpts, num_wann, num_valence_bands, effective_model, &
                               have_disentangled, seedname, stdout, timer, error, comm)
      if (allocated(error)) return

      call utility_recip_lattice_base(real_lattice, recip_lattice, volume)
      Delta_k = pw90common_kmesh_spacing(pw90_berry%kmesh%mesh, recip_lattice)
    else
      call pw90common_fourier_R_to_k_new(ws_region, wannier_data, ws_distance, wigner_seitz, HH_R, &
                                         kpt, real_lattice, mp_grid, num_wann, error, comm, OO=HH, &
                                         OO_dx=delHH(:, :, 1), OO_dy=delHH(:, :, 2), &
                                         OO_dz=delHH(:, :, 3))
      if (allocated(error)) return

      call utility_diagonalize(HH, num_wann, eig, UU, error, comm)
      if (allocated(error)) return
    endif
    call pw90common_get_occ(fermi_energy_list(1), eig, occ, num_wann)

    call wham_get_D_h(delHH, D_h, UU, eig, num_wann)

    call pw90common_fourier_R_to_k_vec(ws_region, wannier_data, ws_distance, wigner_seitz, AA_R, &
                                       kpt, real_lattice, mp_grid, num_wann, error, comm, &
                                       OO_true=AA)
    if (allocated(error)) return

    do i = 1, 3
      AA(:, :, i) = utility_rotate(AA(:, :, i), UU, num_wann)
    enddo
    AA = AA + cmplx_i*D_h ! Eq.(25) WYSV06

    ! Replace imaginary part of frequency with a fixed value
    if (.not. pw90_berry%kubo_smearing%use_adaptive .and. pw90_berry%kubo_smearing%fixed_width /= 0.0_dp) &
      pw90_berry%kubo_freq_list = real(pw90_berry%kubo_freq_list, dp) &
                                  + cmplx_i*pw90_berry%kubo_smearing%fixed_width

    kubo_H_k = cmplx_0
    kubo_AH_k = cmplx_0
    jdos_k = 0.0_dp
    if (spin_decomp) then
      call spin_get_nk(ws_region, pw90_spin, wannier_data, ws_distance, wigner_seitz, HH_R, SS_R, &
                       kpt, real_lattice, spn_nk, mp_grid, num_wann, error, comm)
      if (allocated(error)) return

      kubo_H_k_spn = cmplx_0
      kubo_AH_k_spn = cmplx_0
      jdos_k_spn = 0.0_dp
    end if
    do m = 1, num_wann
      do n = 1, num_wann
        if (n == m) cycle
        if (eig(m) > pw90_berry%kubo_eigval_max .or. eig(n) > pw90_berry%kubo_eigval_max) cycle
        if (spin_decomp) then
          if (spn_nk(n) >= 0 .and. spn_nk(m) >= 0) then
            ispn = 1 ! up --> up transition
          elseif (spn_nk(n) < 0 .and. spn_nk(m) < 0) then
            ispn = 2 ! down --> down
          else
            ispn = 3 ! spin-flip
          end if
        end if
        if (pw90_berry%kubo_smearing%use_adaptive) then
          ! Eq.(35) YWVS07
          vdum(:) = del_eig(m, :) - del_eig(n, :)
          joint_level_spacing = sqrt(dot_product(vdum(:), vdum(:)))*Delta_k
          eta_smr = min(joint_level_spacing*pw90_berry%kubo_smearing%adaptive_prefactor, &
                        pw90_berry%kubo_smearing%adaptive_max_width)
        else
          eta_smr = pw90_berry%kubo_smearing%fixed_width
        endif
        rfac1 = (occ(m) - occ(n))*(eig(m) - eig(n))
        occ_prod = occ(n)*(1.0_dp - occ(m))
        do ifreq = 1, pw90_berry%kubo_nfreq
          !
          ! Complex frequency for the anti-Hermitian conductivity
          !
          if (pw90_berry%kubo_smearing%use_adaptive) then
            omega = real(pw90_berry%kubo_freq_list(ifreq), dp) + cmplx_i*eta_smr
          else
            omega = pw90_berry%kubo_freq_list(ifreq)
          endif
          !
          ! Broadened delta function for the Hermitian conductivity and JDOS
          !
          arg = (eig(m) - eig(n) - real(omega, dp))/eta_smr
          ! If only Hermitean part were computed, could speed up
          ! by inserting here 'if(abs(arg)>10.0_dp) cycle'
          delta = utility_w0gauss(arg, pw90_berry%kubo_smearing%type_index, error, comm)/eta_smr
          if (allocated(error)) return
          !
          ! Lorentzian shape (for testing purposes)
          ! delta=1.0_dp/(1.0_dp+arg*arg)/pi
          ! delta=delta/eta_smr
          !
          jdos_k(ifreq) = jdos_k(ifreq) + occ_prod*delta
          if (spin_decomp) &
            jdos_k_spn(ispn, ifreq) = jdos_k_spn(ispn, ifreq) + occ_prod*delta
          cfac = cmplx_i*rfac1/(eig(m) - eig(n) - omega)
          rfac2 = -pi*rfac1*delta
          do j = 1, 3
            do i = 1, 3
              kubo_H_k(i, j, ifreq) = kubo_H_k(i, j, ifreq) &
                                      + rfac2*AA(n, m, i)*AA(m, n, j)
              kubo_AH_k(i, j, ifreq) = kubo_AH_k(i, j, ifreq) &
                                       + cfac*AA(n, m, i)*AA(m, n, j)
              if (spin_decomp) then
                kubo_H_k_spn(i, j, ispn, ifreq) = &
                  kubo_H_k_spn(i, j, ispn, ifreq) &
                  + rfac2*AA(n, m, i)*AA(m, n, j)
                kubo_AH_k_spn(i, j, ispn, ifreq) = &
                  kubo_AH_k_spn(i, j, ispn, ifreq) &
                  + cfac*AA(n, m, i)*AA(m, n, j)
              endif
            enddo
          enddo
        enddo
      enddo
    enddo

  end subroutine berry_get_kubo_k

  subroutine berry_get_sc_klist(pw90_berry, dis_manifold, fermi_energy_list, kmesh_info, kpt_latt, &
                                ws_region, print_output, pw90_band_deriv_degen, wannier_data, &
                                ws_distance, wigner_seitz, AA_R, HH_R, u_matrix, v_matrix, eigval, &
                                kpt, real_lattice, sc_k_list, scissors_shift, mp_grid, num_bands, &
                                num_kpts, num_wann, num_valence_bands, effective_model, &
                                have_disentangled, seedname, stdout, timer, error, comm)
    !================================================!
    !
    !  Contribution from point k to the nonlinear shift current
    !  [integrand of Eq.8 IATS18]
    !  Notation correspondence with IATS18:
    !  AA_da_bar              <-->   \mathbbm{b}
    !  AA_bar                 <-->   \mathbbm{a}
    !  HH_da_bar              <-->   \mathbbm{v}
    !  HH_dadb_bar            <-->   \mathbbm{w}
    !  D_h(n,m)               <-->   \mathbbm{v}_{nm} * Re[1/(E_{m}-E_{n}+i*sc_eta)]
    !  D_h_no_eta(n,m)        <-->   \mathbbm{v}_{nm} / (E_{m}-E_{n})
    !  sum_AD                 <-->   summatory of Eq. 32 IATS18
    !  sum_HD                 <-->   summatory of Eq. 30 IATS18
    !  eig_da(n)-eig_da(m)    <-->   \mathbbm{Delta}_{nm}
    !
    !================================================!

    use w90_constants, only: dp, cmplx_0, cmplx_i
    use w90_utility, only: utility_re_tr, utility_im_tr, utility_w0gauss, utility_w0gauss_vec
    use w90_types, only: print_output_type, wannier_data_type, &
      dis_manifold_type, kmesh_info_type, ws_region_type, ws_distance_type, timer_list_type
    use w90_postw90_types, only: pw90_berry_mod_type, pw90_band_deriv_degen_type, wigner_seitz_type
    use w90_postw90_common, only: pw90common_fourier_R_to_k_vec_dadb, &
      pw90common_fourier_R_to_k_new_second_d, pw90common_get_occ, &
      pw90common_kmesh_spacing, pw90common_fourier_R_to_k_vec_dadb_TB_conv
    use w90_wan_ham, only: wham_get_D_h, &
      wham_get_eig_UU_HH_AA_sc, wham_get_eig_deleig, wham_get_D_h_P_value, &
      wham_get_eig_deleig_TB_conv, wham_get_eig_UU_HH_AA_sc_TB_conv
    use w90_comms, only: w90comm_type
    use w90_utility, only: utility_rotate, utility_zdotu, utility_recip_lattice_base

    implicit none

    ! arguments
    type(pw90_berry_mod_type), intent(in) :: pw90_berry
    type(dis_manifold_type), intent(in) :: dis_manifold
    type(kmesh_info_type), intent(in) :: kmesh_info
    type(pw90_band_deriv_degen_type), intent(in) :: pw90_band_deriv_degen
    type(print_output_type), intent(in) :: print_output
    type(ws_region_type), intent(in) :: ws_region
    type(w90comm_type), intent(in) :: comm
    type(wannier_data_type), intent(in) :: wannier_data
    type(wigner_seitz_type), intent(inout) :: wigner_seitz
    type(ws_distance_type), intent(inout) :: ws_distance
    type(timer_list_type), intent(inout) :: timer
    type(w90_error_type), allocatable, intent(out) :: error

    integer, intent(in) :: num_wann, num_bands, num_kpts, num_valence_bands
    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: stdout

    real(kind=dp), intent(in) :: kpt(3)
    real(kind=dp), intent(out) :: sc_k_list(:, :, :)
    real(kind=dp), intent(in) :: eigval(:, :)
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    real(kind=dp), intent(in) :: scissors_shift
    real(kind=dp), allocatable, intent(in) :: fermi_energy_list(:)
    real(kind=dp), intent(in) :: kpt_latt(:, :)

    complex(kind=dp), intent(in) :: u_matrix(:, :, :), v_matrix(:, :, :)
    complex(kind=dp), allocatable, intent(inout) :: AA_R(:, :, :, :) ! <0n|r|Rm>
    complex(kind=dp), allocatable, intent(inout) :: HH_R(:, :, :) !  <0n|r|Rm>

    character(len=50), intent(in) :: seedname
    logical, intent(in) :: have_disentangled
    logical, intent(in) :: effective_model

    ! local variables
    complex(kind=dp), allocatable :: UU(:, :)
    complex(kind=dp), allocatable :: AA(:, :, :), AA_bar(:, :, :)
    complex(kind=dp), allocatable :: AA_da(:, :, :, :), AA_da_bar(:, :, :, :)
    complex(kind=dp), allocatable :: HH_da(:, :, :), HH_da_bar(:, :, :)
    complex(kind=dp), allocatable :: HH_dadb(:, :, :, :), HH_dadb_bar(:, :, :, :)
    complex(kind=dp), allocatable :: HH(:, :)
    complex(kind=dp), allocatable :: D_h(:, :, :), D_h_no_eta(:, :, :)

    real(kind=dp), allocatable :: eig(:)
    real(kind=dp), allocatable :: eig_da(:, :)
    real(kind=dp), allocatable :: occ(:)

    real(kind=dp) :: recip_lattice(3, 3), volume
    complex(kind=dp) :: sum_AD(3, 3), sum_HD(3, 3), r_mn(3), gen_r_nm(3)
    integer :: a, b, c, bc, n, m, istart, iend
    integer :: p ! i, if, r, ifreq
    real(kind=dp) :: I_nm(3, 6)
    real(kind=dp) :: omega(pw90_berry%kubo_nfreq), delta(pw90_berry%kubo_nfreq), joint_level_spacing
    real(kind=dp) :: eta_smr, Delta_k, vdum(3), occ_fac, wstep, wmin, wmax

    allocate (UU(num_wann, num_wann))
    allocate (AA(num_wann, num_wann, 3))
    allocate (AA_bar(num_wann, num_wann, 3))
    allocate (AA_da(num_wann, num_wann, 3, 3))
    allocate (AA_da_bar(num_wann, num_wann, 3, 3))
    allocate (HH_da(num_wann, num_wann, 3))
    allocate (HH_da_bar(num_wann, num_wann, 3))
    allocate (HH_dadb(num_wann, num_wann, 3, 3))
    allocate (HH_dadb_bar(num_wann, num_wann, 3, 3))
    allocate (HH(num_wann, num_wann))
    allocate (D_h(num_wann, num_wann, 3))
    allocate (D_h_no_eta(num_wann, num_wann, 3))
    allocate (eig(num_wann))
    allocate (occ(num_wann))
    allocate (eig_da(num_wann, 3))

    ! Initialize shift current array at point k
    sc_k_list = 0.d0

    ! Gather W-gauge matrix objects !

    ! choose the convention for the FT sums
    if (pw90_berry%sc_phase_conv .eq. 1) then ! use Wannier centres in the FT exponentials (so called TB convention)
      ! get Hamiltonian and its first and second derivatives
      ! Note that below we calculate the UU matrix--> we have to use the same UU from here on for
      ! maintaining the gauge-covariance of the whole matrix element
      call wham_get_eig_UU_HH_AA_sc_TB_conv(pw90_berry, dis_manifold, kmesh_info, kpt_latt, &
                                            ws_region, print_output, wannier_data, ws_distance, &
                                            wigner_seitz, AA_R, HH, HH_da, HH_dadb, HH_R, &
                                            u_matrix, UU, v_matrix, eig, eigval, kpt, &
                                            real_lattice, scissors_shift, mp_grid, num_bands, &
                                            num_kpts, num_wann, num_valence_bands, &
                                            effective_model, have_disentangled, seedname, stdout, &
                                            timer, error, comm)
      if (allocated(error)) return
      ! get position operator and its derivative
      ! note that AA_da(:,:,a,b) \propto \sum_R exp(iRk)*iR_{b}*<0|r_{a}|R>
      call pw90common_fourier_R_to_k_vec_dadb_TB_conv(ws_region, wannier_data, ws_distance, &
                                                      wigner_seitz, AA_R, kpt, real_lattice, &
                                                      mp_grid, num_wann, error, comm, OO_da=AA, &
                                                      OO_dadb=AA_da)
      if (allocated(error)) return

      ! get eigenvalues and their k-derivatives
      call wham_get_eig_deleig_TB_conv(pw90_band_deriv_degen, HH_da, UU, eig, eig_da, num_wann, &
                                       error, comm)
      if (allocated(error)) return
    elseif (pw90_berry%sc_phase_conv .eq. 2) then ! do not use Wannier centres in the FT exponentials (usual W90 convention)
      ! same as above
      call wham_get_eig_UU_HH_AA_sc(dis_manifold, kpt_latt, ws_region, print_output, wannier_data, &
                                    ws_distance, wigner_seitz, HH, HH_da, HH_dadb, HH_R, u_matrix, UU, &
                                    v_matrix, eig, eigval, kpt, real_lattice, scissors_shift, &
                                    mp_grid, num_bands, num_kpts, num_wann, num_valence_bands, &
                                    effective_model, have_disentangled, seedname, stdout, timer, &
                                    error, comm)
      if (allocated(error)) return

      call pw90common_fourier_R_to_k_vec_dadb(ws_region, wannier_data, ws_distance, wigner_seitz, &
                                              AA_R, kpt, real_lattice, mp_grid, num_wann, error, &
                                              comm, OO_da=AA, OO_dadb=AA_da)
      if (allocated(error)) return

      call wham_get_eig_deleig(dis_manifold, kpt_latt, pw90_band_deriv_degen, ws_region, print_output, wannier_data, &
                               ws_distance, wigner_seitz, HH_da, HH, HH_R, u_matrix, UU, v_matrix, &
                               eig_da, eig, eigval, kpt, real_lattice, scissors_shift, mp_grid, &
                               num_bands, num_kpts, num_wann, num_valence_bands, effective_model, &
                               have_disentangled, seedname, stdout, timer, error, comm)
      if (allocated(error)) return

    end if

    ! get electronic occupations
    call pw90common_get_occ(fermi_energy_list(1), eig, occ, num_wann)

    ! get D_h (Eq. (24) WYSV06)
    call wham_get_D_h_P_value(pw90_berry, HH_da, D_h, UU, eig, num_wann)
    call wham_get_D_h(HH_da, D_h_no_eta, UU, eig, num_wann)

    ! calculate k-spacing in case of adaptive smearing
    if (pw90_berry%kubo_smearing%use_adaptive) then
      call utility_recip_lattice_base(real_lattice, recip_lattice, volume)
      Delta_k = pw90common_kmesh_spacing(pw90_berry%kmesh%mesh, recip_lattice)
    endif

    ! rotate quantities from W to H gauge (we follow wham_get_D_h for delHH_bar_i)
    do a = 1, 3
      ! Berry connection A
      AA_bar(:, :, a) = utility_rotate(AA(:, :, a), UU, num_wann)
      ! first derivative of Hamiltonian dH_da
      HH_da_bar(:, :, a) = utility_rotate(HH_da(:, :, a), UU, num_wann)
      do b = 1, 3
        ! derivative of Berry connection dA_da
        AA_da_bar(:, :, a, b) = utility_rotate(AA_da(:, :, a, b), UU, num_wann)
        ! second derivative of Hamiltonian d^{2}H_dadb
        HH_dadb_bar(:, :, a, b) = utility_rotate(HH_dadb(:, :, a, b), UU, num_wann)
      enddo
    enddo

    ! setup for frequency-related quantities
    omega = real(pw90_berry%kubo_freq_list(:), dp)
    wmin = omega(1)
    wmax = omega(pw90_berry%kubo_nfreq)
    wstep = omega(2) - omega(1)

    ! loop on initial and final bands
    do n = 1, num_wann
      do m = 1, num_wann
        ! cycle diagonal matrix elements and bands above the maximum
        if (n == m) cycle
        if (eig(m) > pw90_berry%kubo_eigval_max .or. eig(n) > pw90_berry%kubo_eigval_max) cycle
        ! setup T=0 occupation factors
        occ_fac = (occ(n) - occ(m))
        if (abs(occ_fac) < 1e-10) cycle

        ! set delta function smearing
        if (pw90_berry%kubo_smearing%use_adaptive) then
          vdum(:) = eig_da(m, :) - eig_da(n, :)
          joint_level_spacing = sqrt(dot_product(vdum(:), vdum(:)))*Delta_k
          eta_smr = min(joint_level_spacing*pw90_berry%kubo_smearing%adaptive_prefactor, &
                        pw90_berry%kubo_smearing%adaptive_max_width)
        else
          eta_smr = pw90_berry%kubo_smearing%fixed_width
        endif

        ! restrict to energy window spanning [-sc_w_thr*eta_smr,+sc_w_thr*eta_smr]
        ! outside this range, the two delta functions are virtually zero
        if (((eig(n) - eig(m) + pw90_berry%sc_w_thr*eta_smr < wmin) .or. &
             (eig(n) - eig(m) - pw90_berry%sc_w_thr*eta_smr > wmax)) .and. &
            ((eig(m) - eig(n) + pw90_berry%sc_w_thr*eta_smr < wmin) .or. &
             (eig(m) - eig(n) - pw90_berry%sc_w_thr*eta_smr > wmax))) cycle

        ! first compute the two sums over intermediate states between AA_bar and HH_da_bar with D_h
        ! appearing in Eqs. (30) and (32) of IATS18
        sum_AD = cmplx_0
        sum_HD = cmplx_0
        do a = 1, 3
          do c = 1, 3
            ! Note that we substract diagonal elements in AA_bar and
            ! HH_da_bar to match the convention in IATS18
            ! (diagonals in D_h are automatically zero, so we do not substract them)
            sum_AD(c, a) = (utility_zdotu(AA_bar(n, :, c), D_h(:, m, a)) - AA_bar(n, n, c)*D_h(n, m, a)) &
                           - (utility_zdotu(D_h(n, :, a), AA_bar(:, m, c)) - D_h(n, m, a)*AA_bar(m, m, c))
            sum_HD(c, a) = (utility_zdotu(HH_da_bar(n, :, c), D_h(:, m, a)) - HH_da_bar(n, n, c)*D_h(n, m, a)) &
                           - (utility_zdotu(D_h(n, :, a), HH_da_bar(:, m, c)) - D_h(n, m, a)*HH_da_bar(m, m, c))
          enddo
        enddo

        ! dipole matrix element
        r_mn(:) = AA_bar(m, n, :) + cmplx_i*D_h_no_eta(m, n, :)

        ! loop over direction of generalized derivative
        do a = 1, 3
          ! store generalized derivative as an array on the additional spatial index,
          ! its composed of 8 terms in total, see Eq (34) combined with (30) and
          ! (32) of IATS18
          gen_r_nm(:) = (AA_da_bar(n, m, :, a) &
                         + ((AA_bar(n, n, :) - AA_bar(m, m, :))*D_h_no_eta(n, m, a) + &
                            (AA_bar(n, n, a) - AA_bar(m, m, a))*D_h_no_eta(n, m, :)) &
                         - cmplx_i*AA_bar(n, m, :)*(AA_bar(n, n, a) - AA_bar(m, m, a)) &
                         + sum_AD(:, a) &
                         + cmplx_i*(HH_dadb_bar(n, m, :, a) &
                                    + sum_HD(:, a) &
                                    + (D_h_no_eta(n, m, :)*(eig_da(n, a) - eig_da(m, a)) + &
                                       D_h_no_eta(n, m, a)*(eig_da(n, :) - eig_da(m, :)))) &
                         /(eig(m) - eig(n)))

          ! Correction term due to finite sc_eta
          ! See Eq. (19) of Phys. Rev. B 103, 247101 (2021)
          if (pw90_berry%sc_use_eta_corr) then
            do p = 1, num_wann
              if (p == n .or. p == m) cycle
              gen_r_nm(:) = gen_r_nm(:) &
                            - pw90_berry%sc_eta**2/((eig(p) - eig(m))**2 &
                                                    + pw90_berry%sc_eta**2)/(eig(n) - eig(m)) &
                            *(AA_bar(n, p, :)*HH_da_bar(p, m, a) &
                              - (HH_da_bar(n, p, :) + cmplx_i*(eig(n) &
                                                               - eig(p))*AA_bar(n, p, :))*AA_bar(p, m, a)) &
                            + pw90_berry%sc_eta**2/((eig(n) - eig(p))**2 &
                                                    + pw90_berry%sc_eta**2)/(eig(n) - eig(m)) &
                            *(HH_da_bar(n, p, a)*AA_bar(p, m, :) &
                              - AA_bar(n, p, a)*(HH_da_bar(p, m, :) + cmplx_i*(eig(p) - eig(m))*AA_bar(p, m, :)))
            enddo
          endif

          ! loop over the remaining two indexes of the matrix product.
          ! Note that shift current is symmetric under b <--> c exchange,
          ! so we avoid computing all combinations using alpha_S and beta_S
          do bc = 1, 6
            b = alpha_S(bc)
            c = beta_S(bc)
            I_nm(a, bc) = aimag(r_mn(b)*gen_r_nm(c) + r_mn(c)*gen_r_nm(b))
          enddo ! bc
        enddo ! a

        ! compute delta(E_nm-w)
        ! choose energy window spanning [-sc_w_thr*eta_smr,+sc_w_thr*eta_smr]
        istart = max(int((eig(n) - eig(m) - pw90_berry%sc_w_thr*eta_smr - wmin)/wstep + 1), 1)
        iend = min(int((eig(n) - eig(m) + pw90_berry%sc_w_thr*eta_smr - wmin)/wstep + 1), pw90_berry%kubo_nfreq)
        ! multiply matrix elements with delta function for the relevant frequencies
        if (istart <= iend) then
          delta = 0.0
          delta(istart:iend) = &
            utility_w0gauss_vec((eig(m) - eig(n) + omega(istart:iend))/eta_smr, &
                                pw90_berry%kubo_smearing%type_index, error, comm)/eta_smr
          if (allocated(error)) return
          call DGER(18, iend - istart + 1, occ_fac, I_nm, 1, delta(istart:iend), 1, sc_k_list(:, :, istart:iend), 18)
        endif
        ! same for delta(E_mn-w)
        istart = max(int((eig(m) - eig(n) - pw90_berry%sc_w_thr*eta_smr - wmin)/wstep + 1), 1)
        iend = min(int((eig(m) - eig(n) + pw90_berry%sc_w_thr*eta_smr - wmin)/wstep + 1), pw90_berry%kubo_nfreq)
        if (istart <= iend) then
          delta = 0.0
          delta(istart:iend) = &
            utility_w0gauss_vec((eig(n) - eig(m) + omega(istart:iend))/eta_smr, &
                                pw90_berry%kubo_smearing%type_index, error, comm)/eta_smr
          if (allocated(error)) return
          call DGER(18, iend - istart + 1, occ_fac, I_nm, 1, delta(istart:iend), 1, sc_k_list(:, :, istart:iend), 18)
        endif

      enddo ! bands
    enddo ! bands

  end subroutine berry_get_sc_klist

  !================================================!
  subroutine berry_get_shc_klist(pw90_berry, dis_manifold, fermi_energy_list, kpt_latt, &
                                 pw90_band_deriv_degen, ws_region, pw90_spin_hall, print_output, &
                                 wannier_data, ws_distance, wigner_seitz, AA_R, HH_R, SH_R, SHR_R, &
                                 SR_R, SS_R, SAA_R, SBB_R, u_matrix, v_matrix, eigval, kpt, &
                                 real_lattice, scissors_shift, mp_grid, fermi_n, num_bands, &
                                 num_kpts, num_wann, num_valence_bands, effective_model, &
                                 have_disentangled, seedname, stdout, timer, error, comm, &
                                 shc_k_fermi, shc_k_freq, shc_k_band)
    !================================================!
    !
    ! Contribution from a k-point to the spin Hall conductivity on a list
    ! of Fermi energies or a list of frequencies or a list of energy bands
    !   sigma_{alpha,beta}^{gamma}(k), alpha, beta, gamma = 1, 2, 3
    !                                                      (x, y, z, respectively)
    ! i.e. the Berry curvature-like term of QZYZ18 Eq.(3) & (4).
    ! The unit is angstrom^2, similar to that of Berry curvature of AHC.
    !
    !  Note the berry_get_js_k() has not been multiplied by hbar/2 (as
    !  required by spin operator) and not been divided by hbar (as required
    !  by the velocity operator). The second velocity operator has not been
    !  divided by hbar as well. But these two hbar required by velocity
    !  operators are canceled by the preceding hbar^2 of QZYZ18 Eq.(3).
    !
    !    shc_k_fermi: return a list for different Fermi energies
    !    shc_k_freq:  return a list for different frequencies
    !    shc_k_band:  return a list for each energy band
    !
    !   Junfeng Qiao (18/8/2018)
    !================================================!

    use w90_constants, only: dp, cmplx_0, cmplx_i
    use w90_utility, only: utility_rotate, utility_recip_lattice_base
    use w90_comms, only: w90comm_type
    use w90_types, only: print_output_type, wannier_data_type, &
      dis_manifold_type, kmesh_info_type, ws_region_type, ws_distance_type, timer_list_type
    use w90_postw90_types, only: pw90_berry_mod_type, pw90_spin_hall_type, &
      pw90_band_deriv_degen_type, wigner_seitz_type
    use w90_postw90_common, only: pw90common_get_occ, pw90common_fourier_R_to_k_vec, &
      pw90common_kmesh_spacing
    use w90_wan_ham, only: wham_get_D_h, wham_get_eig_deleig

    implicit none

    ! arguments
    type(pw90_berry_mod_type), intent(in) :: pw90_berry
    type(dis_manifold_type), intent(in) :: dis_manifold
    real(kind=dp), allocatable, intent(in) :: fermi_energy_list(:)
    real(kind=dp), intent(in) :: kpt_latt(:, :)
    type(pw90_band_deriv_degen_type), intent(in) :: pw90_band_deriv_degen
    type(print_output_type), intent(in) :: print_output
    type(ws_region_type), intent(in) :: ws_region
    type(pw90_spin_hall_type), intent(in) :: pw90_spin_hall
    type(w90comm_type), intent(in) :: comm
    type(wannier_data_type), intent(in) :: wannier_data
    type(wigner_seitz_type), intent(inout) :: wigner_seitz
    type(ws_distance_type), intent(inout) :: ws_distance
    type(timer_list_type), intent(inout) :: timer
    type(w90_error_type), allocatable, intent(out) :: error

    integer, intent(in) :: num_wann, num_bands, num_kpts, num_valence_bands, fermi_n
    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: stdout

    real(kind=dp), intent(in) :: kpt(3)
    real(kind=dp), intent(in) :: eigval(:, :)
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    real(kind=dp), intent(in) :: scissors_shift

    complex(kind=dp), intent(in) :: u_matrix(:, :, :), v_matrix(:, :, :)
    complex(kind=dp), allocatable, intent(inout) :: AA_R(:, :, :, :) ! <0n|r|Rm>
    complex(kind=dp), allocatable, intent(inout) :: HH_R(:, :, :) !  <0n|r|Rm>
    complex(kind=dp), allocatable, intent(inout) :: SR_R(:, :, :, :, :) ! <0n|sigma_x,y,z.(r-R)_alpha|Rm>
    complex(kind=dp), allocatable, intent(inout) :: SHR_R(:, :, :, :, :) ! <0n|sigma_x,y,z.H.(r-R)_alpha|Rm>
    complex(kind=dp), allocatable, intent(inout) :: SH_R(:, :, :, :) ! <0n|sigma_x,y,z.H|Rm>
    complex(kind=dp), allocatable, intent(inout) :: SS_R(:, :, :, :) ! <0n|sigma_x,y,z|Rm>
    complex(kind=dp), allocatable, intent(inout) :: SAA_R(:, :, :, :, :)
    complex(kind=dp), allocatable, intent(inout) :: SBB_R(:, :, :, :, :)

    character(len=50), intent(in) :: seedname
    logical, intent(in) :: have_disentangled
    logical, intent(in) :: effective_model

    complex(kind=dp), optional, intent(out) :: shc_k_freq(pw90_berry%kubo_nfreq)
    real(kind=dp), optional, intent(out) :: shc_k_fermi(fermi_n)
    real(kind=dp), optional, intent(out) :: shc_k_band(num_wann)

    ! internal vars
    complex(kind=dp), allocatable :: HH(:, :)
    complex(kind=dp), allocatable :: delHH(:, :, :)
    complex(kind=dp), allocatable :: UU(:, :)
    complex(kind=dp), allocatable :: D_h(:, :, :)
    complex(kind=dp), allocatable :: AA(:, :, :)
    complex(kind=dp) :: js_k(num_wann, num_wann)

    logical :: lfreq, lfermi, lband

    real(kind=dp) :: recip_lattice(3, 3), volume
    integer :: n, m, i, ifreq

    ! Adaptive smearing
    real(kind=dp) :: del_eig(num_wann, 3), joint_level_spacing, eta_smr, Delta_k, vdum(3)
    real(kind=dp) :: eig(num_wann)
    real(kind=dp) :: occ_fermi(num_wann, fermi_n), occ_freq(num_wann)
    real(kind=dp) :: omega, rfac

    complex(kind=dp) :: omega_list(pw90_berry%kubo_nfreq)
    complex(kind=dp) :: prod, cdum, cfac

    allocate (HH(num_wann, num_wann))
    allocate (delHH(num_wann, num_wann, 3))
    allocate (UU(num_wann, num_wann))
    allocate (D_h(num_wann, num_wann, 3))
    allocate (AA(num_wann, num_wann, 3))

    lfreq = .false.
    lfermi = .false.
    lband = .false.
    if (present(shc_k_freq)) then
      shc_k_freq = 0.0_dp
      lfreq = .true.
    endif
    if (present(shc_k_fermi)) then
      shc_k_fermi = 0.0_dp
      lfermi = .true.
    endif
    if (present(shc_k_band)) then
      shc_k_band = 0.0_dp
      lband = .true.
    endif

    call wham_get_eig_deleig(dis_manifold, kpt_latt, pw90_band_deriv_degen, ws_region, &
                             print_output, wannier_data, ws_distance, wigner_seitz, delHH, HH, &
                             HH_R, u_matrix, UU, v_matrix, del_eig, eig, eigval, kpt, &
                             real_lattice, scissors_shift, mp_grid, num_bands, num_kpts, num_wann, &
                             num_valence_bands, effective_model, have_disentangled, seedname, &
                             stdout, timer, error, comm)
    if (allocated(error)) return

    call wham_get_D_h(delHH, D_h, UU, eig, num_wann)

    ! Here I apply a scissor operator to the conduction bands, if required in the input
    if (pw90_spin_hall%bandshift) then
      eig(pw90_spin_hall%bandshift_firstband:) = eig(pw90_spin_hall%bandshift_firstband:) + pw90_spin_hall%bandshift_energyshift
    end if

    call pw90common_fourier_R_to_k_vec(ws_region, wannier_data, ws_distance, wigner_seitz, AA_R, &
                                       kpt, real_lattice, mp_grid, num_wann, error, comm, &
                                       OO_true=AA)
    if (allocated(error)) return

    do i = 1, 3
      AA(:, :, i) = utility_rotate(AA(:, :, i), UU, num_wann)
    enddo
    AA = AA + cmplx_i*D_h ! Eq.(25) WYSV06

    call berry_get_js_k(ws_region, pw90_spin_hall, wannier_data, ws_distance, wigner_seitz, &
                        D_h(:, :, pw90_spin_hall%alpha), js_k, SH_R, SHR_R, SR_R, SS_R, SAA_R, &
                        SBB_R, UU, eig, del_eig(:, pw90_spin_hall%alpha), &
                        delHH(:, :, pw90_spin_hall%alpha), kpt, real_lattice, mp_grid, num_wann)

    ! adpt_smr only works with pw90_berry_kmesh, so do not use
    ! adpt_smr in kpath or kslice plots.
    if (pw90_berry%kubo_smearing%use_adaptive) then
      call utility_recip_lattice_base(real_lattice, recip_lattice, volume)
      Delta_k = pw90common_kmesh_spacing(pw90_berry%kmesh%mesh, recip_lattice)
    endif
    if (lfreq) then
      call pw90common_get_occ(fermi_energy_list(1), eig, occ_freq, num_wann)
    elseif (lfermi) then
      ! get occ for different fermi_energy
      do i = 1, fermi_n
        call pw90common_get_occ(fermi_energy_list(i), eig, occ_fermi(:, i), num_wann)
      end do
    end if
    do n = 1, num_wann
      ! get Omega_{n,alpha beta}^{gamma}
      if (lfreq) then
        omega_list = cmplx_0
      else if (lfermi .or. lband) then
        omega = 0.0_dp
      end if
      do m = 1, num_wann
        if (m == n) cycle
        if (eig(m) > pw90_berry%kubo_eigval_max .or. eig(n) > pw90_berry%kubo_eigval_max) cycle

        rfac = eig(m) - eig(n)
        !this will calculate AHC
        !prod = -rfac*cmplx_i*AA(n, m, shc_alpha) * rfac*cmplx_i*AA(m, n, shc_beta)
        prod = js_k(n, m)*cmplx_i*rfac*AA(m, n, pw90_spin_hall%beta)
        if (pw90_berry%kubo_smearing%use_adaptive) then
          ! Eq.(35) YWVS07
          vdum(:) = del_eig(m, :) - del_eig(n, :)
          joint_level_spacing = sqrt(dot_product(vdum(:), vdum(:)))*Delta_k
          eta_smr = min(joint_level_spacing*pw90_berry%kubo_smearing%adaptive_prefactor, &
                        pw90_berry%kubo_smearing%adaptive_max_width)
        else
          eta_smr = pw90_berry%kubo_smearing%fixed_width
        endif
        if (lfreq) then
          do ifreq = 1, pw90_berry%kubo_nfreq
            cdum = real(pw90_berry%kubo_freq_list(ifreq), dp) + cmplx_i*eta_smr
            cfac = -2.0_dp/(rfac**2 - cdum**2)
            omega_list(ifreq) = omega_list(ifreq) + cfac*aimag(prod)
          end do
        else if (lfermi .or. lband) then
          rfac = -2.0_dp/(rfac**2 + eta_smr**2)
          omega = omega + rfac*aimag(prod)
        end if
      enddo

      if (lfermi) then
        do i = 1, fermi_n
          shc_k_fermi(i) = shc_k_fermi(i) + occ_fermi(n, i)*omega
        end do
      else if (lfreq) then
        shc_k_freq = shc_k_freq + occ_freq(n)*omega_list
      else if (lband) then
        shc_k_band(n) = omega
      end if
    enddo

    !if (lfermi) then
    !  write (*, '(3(f9.6,1x),f16.8,1x,1E16.8)') &
    !    kpt(1), kpt(2), kpt(3), fermi_energy_list(1), shc_k_fermi(1)
    !end if

    return

  contains

    !================================================!
    !                   PRIVATE PROCEDURES
    !================================================!
    subroutine berry_get_js_k(ws_region, pw90_spin_hall, wannier_data, ws_distance, wigner_seitz, &
                              D_alpha_h, js_k, SH_R, SHR_R, SR_R, SS_R, SAA_R, SBB_R, UU, eig, &
                              del_alpha_eig, delHH_alpha, kpt, real_lattice, mp_grid, num_wann)
      !================================================!
      !
      ! Contribution from point k to the
      !   <psi_k | 1/2*(sigma_gamma*v_alpha + v_alpha*sigma_gamma) | psi_k>
      !
      !  QZYZ18 Eq.(23) without hbar/2 (required by spin operator) and
      !  not divided by hbar (required by velocity operator)
      !
      !  Junfeng Qiao (8/7/2018)
      !
      !================================================!
      use w90_constants, only: dp, cmplx_0, cmplx_i
      use w90_utility, only: utility_rotate
      use w90_types, only: print_output_type, wannier_data_type, ws_region_type, &
        ws_distance_type
      use w90_postw90_types, only: pw90_spin_hall_type, wigner_seitz_type
      use w90_postw90_common, only: pw90common_fourier_R_to_k_new, pw90common_fourier_R_to_k_vec

      implicit none

      ! arguments
      type(ws_region_type), intent(in) :: ws_region
      type(pw90_spin_hall_type), intent(in) :: pw90_spin_hall
      type(wannier_data_type), intent(in) :: wannier_data
      type(wigner_seitz_type), intent(in) :: wigner_seitz
      type(ws_distance_type), intent(inout) :: ws_distance

      integer, intent(in) :: mp_grid(3)
      integer, intent(in) :: num_wann

      real(kind=dp), intent(in) :: kpt(3)
      real(kind=dp), intent(in) :: eig(:)
      real(kind=dp), intent(in) :: del_alpha_eig(:)
      real(kind=dp), intent(in) :: real_lattice(3, 3)

      complex(kind=dp), dimension(:, :), intent(in)  :: delHH_alpha
      complex(kind=dp), intent(in) :: D_alpha_h(:, :)
      complex(kind=dp), intent(in) :: UU(:, :)
      complex(kind=dp), intent(out) :: js_k(:, :)
      complex(kind=dp), allocatable, intent(inout) :: SR_R(:, :, :, :, :) ! <0n|sigma_x,y,z.(r-R)_alpha|Rm>
      complex(kind=dp), allocatable, intent(inout) :: SHR_R(:, :, :, :, :) ! <0n|sigma_x,y,z.H.(r-R)_alpha|Rm>
      complex(kind=dp), allocatable, intent(inout) :: SH_R(:, :, :, :) ! <0n|sigma_x,y,z.H|Rm>
      complex(kind=dp), allocatable, intent(inout) :: SS_R(:, :, :, :) ! <0n|sigma_x,y,z|Rm>
      complex(kind=dp), allocatable, intent(inout) :: SAA_R(:, :, :, :, :)
      complex(kind=dp), allocatable, intent(inout) :: SBB_R(:, :, :, :, :)

      ! internal vars
      complex(kind=dp) :: B_k(num_wann, num_wann)
      complex(kind=dp) :: K_k(num_wann, num_wann)
      complex(kind=dp) :: L_k(num_wann, num_wann)
      complex(kind=dp) :: S_w(num_wann, num_wann)
      complex(kind=dp) :: S_k(num_wann, num_wann)
      complex(kind=dp) :: SR_w(num_wann, num_wann, 3)
      complex(kind=dp) :: SR_alpha_k(num_wann, num_wann)
      complex(kind=dp) :: SHR_w(num_wann, num_wann, 3)
      complex(kind=dp) :: SHR_alpha_k(num_wann, num_wann)
      complex(kind=dp) :: SH_w(num_wann, num_wann, 3)
      complex(kind=dp) :: SH_k(num_wann, num_wann)
      complex(kind=dp) :: eig_mat(num_wann, num_wann)
      complex(kind=dp) :: del_eig_mat(num_wann, num_wann)

      !ryoo
      complex(kind=dp)    :: SAA(num_wann, num_wann, 3, 3)
      complex(kind=dp)    :: SBB(num_wann, num_wann, 3, 3)
      complex(kind=dp)    :: VV0(num_wann, num_wann)
      complex(kind=dp)    :: spinvel0(num_wann, num_wann)
      integer :: i

      !================================================
      js_k = cmplx_0

      !================================================ S_k ===========
      ! < u_k | sigma_gamma | u_k >, QZYZ18 Eq.(25)
      ! QZYZ18 Eq.(36)
      call pw90common_fourier_R_to_k_new(ws_region, wannier_data, ws_distance, wigner_seitz, &
                                         SS_R(:, :, :, pw90_spin_hall%gamma), kpt, real_lattice, &
                                         mp_grid, num_wann, error, comm, OO=S_w)
      if (allocated(error)) return

      ! QZYZ18 Eq.(30)
      S_k = utility_rotate(S_w, UU, num_wann)

      if (index(pw90_spin_hall%method, 'qiao') > 0) then !if Qiao
        !================================================ K_k ===========
        ! < u_k | sigma_gamma | \partial_alpha u_k >, QZYZ18 Eq.(26)
        ! QZYZ18 Eq.(37)
        call pw90common_fourier_R_to_k_vec(ws_region, wannier_data, ws_distance, wigner_seitz, &
                                           SR_R(:, :, :, pw90_spin_hall%gamma, :), kpt, &
                                           real_lattice, mp_grid, num_wann, error, comm, &
                                           OO_true=SR_w)
        if (allocated(error)) return

        ! QZYZ18 Eq.(31)
        SR_alpha_k = -cmplx_i*utility_rotate(SR_w(:, :, pw90_spin_hall%alpha), UU, num_wann)
        K_k = SR_alpha_k + matmul(S_k, D_alpha_h)

        !================================================ L_k ===========
        ! < u_k | sigma_gamma.H | \partial_alpha u_k >, QZYZ18 Eq.(27)
        ! QZYZ18 Eq.(38)
        call pw90common_fourier_R_to_k_vec(ws_region, wannier_data, ws_distance, wigner_seitz, &
                                           SHR_R(:, :, :, pw90_spin_hall%gamma, :), kpt, &
                                           real_lattice, mp_grid, num_wann, error, comm, &
                                           OO_true=SHR_w)
        if (allocated(error)) return
        ! QZYZ18 Eq.(32)
        SHR_alpha_k = -cmplx_i*utility_rotate(SHR_w(:, :, pw90_spin_hall%alpha), UU, num_wann)
        ! QZYZ18 Eq.(39)
        call pw90common_fourier_R_to_k_vec(ws_region, wannier_data, ws_distance, wigner_seitz, &
                                           SH_R, kpt, real_lattice, mp_grid, num_wann, error, &
                                           comm, OO_true=SH_w)
        if (allocated(error)) return

        ! QZYZ18 Eq.(32)
        SH_k = utility_rotate(SH_w(:, :, pw90_spin_hall%gamma), UU, num_wann)
        L_k = SHR_alpha_k + matmul(SH_k, D_alpha_h)

        !================================================ B_k ===========
        ! < \psi_nk | sigma_gamma v_alpha | \psi_mk >, QZYZ18 Eq.(24)
        B_k = cmplx_0
        do i = 1, num_wann
          eig_mat(i, :) = eig(:)
          del_eig_mat(i, :) = del_alpha_eig(:)
        end do
        ! note * is not matmul
        B_k = del_eig_mat*S_k + eig_mat*K_k - L_k

        !================================================ js_k ===========
        ! QZYZ18 Eq.(23)
        ! note the S in SR_R,SHR_R,SH_R of get_SHC_R is sigma,
        ! to get spin current, we need to multiply it by hbar/2,
        ! also we need to divide it by hbar to recover the velocity
        ! operator, these are done outside of this subroutine
        js_k = 1.0_dp/2.0_dp*(B_k + conjg(transpose(B_k)))

      else !if Ryoo  (PRB RPS19 Eq.(21))
        !RPS19 Eqs.(37)-(40)
        call pw90common_fourier_R_to_k_new(ws_region, wannier_data, ws_distance, wigner_seitz, &
                                           SAA_R(:, :, :, pw90_spin_hall%gamma, &
                                                 pw90_spin_hall%alpha), kpt, real_lattice, mp_grid, &
                                           num_wann, error, comm, &
                                           OO=SAA(:, :, pw90_spin_hall%gamma, pw90_spin_hall%alpha))
        if (allocated(error)) return
        call pw90common_fourier_R_to_k_new(ws_region, wannier_data, ws_distance, wigner_seitz, &
                                           SBB_R(:, :, :, pw90_spin_hall%gamma, &
                                                 pw90_spin_hall%alpha), kpt, real_lattice, mp_grid, &
                                           num_wann, error, comm, &
                                           OO=SBB(:, :, pw90_spin_hall%gamma, pw90_spin_hall%alpha))

        if (allocated(error)) return
        call pw90common_fourier_R_to_k_new(ws_region, wannier_data, ws_distance, wigner_seitz, &
                                           HH_R, kpt, real_lattice, mp_grid, num_wann, error, &
                                           comm, OO=HH, OO_dx=delHH(:, :, 1), &
                                           OO_dy=delHH(:, :, 2), OO_dz=delHH(:, :, 3))
        if (allocated(error)) return

        VV0(:, :) = utility_rotate(delHH_alpha(:, :), UU, num_wann)
        SAA(:, :, pw90_spin_hall%gamma, pw90_spin_hall%alpha) = &
          utility_rotate(SAA(:, :, pw90_spin_hall%gamma, pw90_spin_hall%alpha), UU, num_wann)
        SBB(:, :, pw90_spin_hall%gamma, pw90_spin_hall%alpha) = &
          utility_rotate(SBB(:, :, pw90_spin_hall%gamma, pw90_spin_hall%alpha), UU, num_wann)

        spinVel0(:, :) = matmul(VV0(:, :), S_k(:, :)) + &
                         matmul(S_k(:, :), VV0(:, :))

        do n = 1, num_wann
          do m = 1, num_wann !RPS19 Eq.(21) and Eq.(26)
            js_k(n, m) = spinVel0(n, m) &
                         - cmplx_i*(eig(m)*SAA(n, m, pw90_spin_hall%gamma, pw90_spin_hall%alpha) &
                                    - SBB(n, m, pw90_spin_hall%gamma, pw90_spin_hall%alpha))
            js_k(n, m) = js_k(n, m) &
                         + cmplx_i*(eig(n)*conjg(SAA(m, n, pw90_spin_hall%gamma, pw90_spin_hall%alpha)) &
                                    - conjg(SBB(m, n, pw90_spin_hall%gamma, pw90_spin_hall%alpha)))
          enddo
        enddo
        js_k = js_k/2.0_dp
      endif

    end subroutine berry_get_js_k

  end subroutine berry_get_shc_klist

  !================================================!
  subroutine berry_print_progress(end_k, loop_k, start_k, step_k, stdout)
    !================================================!
    ! print k-points calculation progress, seperated into 11 points,
    ! from 0%, 10%, ... to 100%
    ! start_k, end_k are inclusive
    ! loop_k should in the array start_k to end_k with step step_k
    !
    ! only call from root MPI process!
    !================================================!

    use w90_io, only: io_wallclocktime

    implicit none

    ! arguments
    integer, intent(in) :: loop_k, start_k, end_k, step_k, stdout

    ! local variables
    real(kind=dp) :: cur_time, finished
    real(kind=dp), save :: prev_time
    integer :: i, j, n, last_k
    logical, dimension(9) :: kmesh_processed = (/(.false., i=1, 9)/)

    ! The last loop_k in the array start:step:end
    ! e.g. 4 of 0:4:7 = [0, 4], 11 of 3:4:11 = [3, 7, 11]
    last_k = (CEILING((end_k - start_k + 1)/real(step_k)) - 1)*step_k + start_k

    if (loop_k == start_k) then
      write (stdout, '(1x,a)') ''
      write (stdout, '(1x,a)') 'Calculation started'
      write (stdout, '(1x,a)') '-------------------------------'
      write (stdout, '(1x,a)') '  k-points       wall      diff'
      write (stdout, '(1x,a)') ' calculated      time      time'
      write (stdout, '(1x,a)') ' ----------      ----      ----'
      cur_time = io_wallclocktime()
      prev_time = cur_time
      write (stdout, '(5x,a,3x,f10.1,f10.1)') '  0%', cur_time, cur_time - prev_time
    else if (loop_k == last_k) then
      cur_time = io_wallclocktime()
      write (stdout, '(5x,a,3x,f10.1,f10.1)') '100%', cur_time, cur_time - prev_time
      write (stdout, '(1x,a)') ''
    else
      finished = 10.0_dp*real(loop_k - start_k + 1)/real(end_k - start_k + 1)
      do n = 1, size(kmesh_processed)
        if ((.not. kmesh_processed(n)) .and. (finished >= n)) then
          do i = n, size(kmesh_processed)
            if (i <= finished) then
              j = i
              kmesh_processed(i) = .true.
            end if
          end do
          cur_time = io_wallclocktime()
          write (stdout, '(5x,i2,a,3x,f10.1,f10.1)') j, '0%', cur_time, cur_time - prev_time
          prev_time = cur_time
          exit
        end if
      end do
    end if

  end subroutine berry_print_progress

  subroutine berry_get_kdotp(kdotp, dis_manifold, kpt_latt, print_output, pw90_berry, &
                             pw90_band_deriv_degen, wannier_data, ws_distance, wigner_seitz, &
                             ws_region, HH_R, u_matrix, v_matrix, eigval, real_lattice, &
                             scissors_shift, mp_grid, num_bands, num_kpts, num_wann, &
                             num_valence_bands, effective_model, have_disentangled, seedname, &
                             stdout, timer, error, comm)
    !================================================!
    !  Extracts k.p expansion coefficients using quasi-degenerate
    !  (Lowdin) perturbation theory, adapted to the Wannier formalism,
    !  see Appendix in IAdJS19 for details
    !================================================!

    use w90_constants, only: dp, cmplx_0, cmplx_i
    use w90_wan_ham, only: wham_get_D_h, wham_get_eig_UU_HH_AA_sc, wham_get_eig_deleig, &
      wham_get_D_h_P_value
    use w90_utility, only: utility_rotate
    use w90_types, only: print_output_type, wannier_data_type, &
      dis_manifold_type, kmesh_info_type, ws_region_type, ws_distance_type, timer_list_type
    use w90_postw90_types, only: pw90_berry_mod_type, pw90_spin_mod_type, &
      pw90_spin_hall_type, pw90_band_deriv_degen_type, pw90_oper_read_type, wigner_seitz_type, &
      kpoint_dist_type
    use w90_comms, only: w90comm_type

    implicit none

    ! Arguments
    type(pw90_berry_mod_type), intent(in) :: pw90_berry
    type(dis_manifold_type), intent(in) :: dis_manifold
    type(pw90_band_deriv_degen_type), intent(in) :: pw90_band_deriv_degen
    type(print_output_type), intent(in) :: print_output
    type(ws_region_type), intent(in) :: ws_region
    type(wannier_data_type), intent(in) :: wannier_data
    type(wigner_seitz_type), intent(inout) :: wigner_seitz
    type(ws_distance_type), intent(inout) :: ws_distance
    type(timer_list_type), intent(inout) :: timer
    type(w90comm_type), intent(in) :: comm
    type(w90_error_type), allocatable, intent(out) :: error

    complex(kind=dp), intent(out), dimension(:, :, :, :, :)     :: kdotp
    !complex(kind=dp), allocatable, intent(inout) :: AA_R(:, :, :, :) ! <0n|r|Rm>
    complex(kind=dp), allocatable, intent(inout) :: HH_R(:, :, :) !  <0n|r|Rm>
    complex(kind=dp), intent(in) :: u_matrix(:, :, :), v_matrix(:, :, :)

    real(kind=dp), intent(in) :: eigval(:, :)
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    real(kind=dp), intent(in) :: scissors_shift
    real(kind=dp), intent(in) :: kpt_latt(:, :)

    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: num_wann, num_kpts, num_bands, num_valence_bands
    integer, intent(in) :: stdout

    character(len=50), intent(in) :: seedname
    logical, intent(in) :: have_disentangled
    logical, intent(in) :: effective_model

    complex(kind=dp), allocatable :: UU(:, :)
    complex(kind=dp), allocatable :: HH_da(:, :, :), HH_da_bar(:, :, :)
    complex(kind=dp), allocatable :: HH_dadb(:, :, :, :), HH_dadb_bar(:, :, :, :)
    complex(kind=dp), allocatable :: HH(:, :), HH_bar(:, :)
    real(kind=dp), allocatable    :: eig(:)
    real(kind=dp), allocatable    :: eig_da(:, :)
    complex(kind=dp), allocatable :: D_h(:, :, :)

    ! local variables
    !real(kind=dp) :: DeltaE_n, DeltaE_m
    integer :: kdotp_num_bands
    integer :: i, a, b, n, m, r !, c, bc,if, ifreq, istart, iend
    logical :: break_loop

    allocate (UU(num_wann, num_wann))
    allocate (HH_da(num_wann, num_wann, 3))
    allocate (HH_da_bar(num_wann, num_wann, 3))
    allocate (HH_dadb(num_wann, num_wann, 3, 3))
    allocate (HH_dadb_bar(num_wann, num_wann, 3, 3))
    allocate (HH(num_wann, num_wann))
    allocate (HH_bar(num_wann, num_wann))
    allocate (eig(num_wann))
    allocate (eig_da(num_wann, 3))
    allocate (D_h(num_wann, num_wann, 3))

    ! Gather W-gauge matrix objects !

    ! get Hamiltonian and its first and second derivatives
    call wham_get_eig_UU_HH_AA_sc(dis_manifold, kpt_latt, ws_region, &
                                  print_output, wannier_data, ws_distance, wigner_seitz, HH, &
                                  HH_da, HH_dadb, HH_R, u_matrix, UU, v_matrix, eig, eigval, &
                                  pw90_berry%kdotp_kpoint, real_lattice, scissors_shift, mp_grid, &
                                  num_bands, num_kpts, num_wann, num_valence_bands, &
                                  effective_model, have_disentangled, seedname, stdout, timer, &
                                  error, comm)
    if (allocated(error)) return

    ! get eigenvalues and their k-derivatives
    call wham_get_eig_deleig(dis_manifold, kpt_latt, pw90_band_deriv_degen, ws_region, &
                             print_output, wannier_data, ws_distance, wigner_seitz, HH_da, HH, &
                             HH_R, u_matrix, UU, v_matrix, eig_da, eig, eigval, &
                             pw90_berry%kdotp_kpoint, real_lattice, scissors_shift, mp_grid, &
                             num_bands, num_kpts, num_wann, num_valence_bands, &
                             effective_model, have_disentangled, seedname, stdout, timer, &
                             error, comm)
    if (allocated(error)) return

    ! get D_h (Eq. (24) WYSV06)
    call wham_get_D_h_P_value(pw90_berry, HH_da, D_h, UU, eig, num_wann)

    ! rotate quantities from W to H gauge
    HH_bar(:, :) = utility_rotate(HH(:, :), UU, num_wann)
    do a = 1, 3
      ! first derivative of Hamiltonian dH_da
      HH_da_bar(:, :, a) = utility_rotate(HH_da(:, :, a), UU, num_wann)
      do b = 1, 3
        ! second derivative of Hamiltonian d^{2}H_dadb
        HH_dadb_bar(:, :, a, b) = utility_rotate(HH_dadb(:, :, a, b), UU, num_wann)
      enddo
    enddo

    kdotp_num_bands = size(pw90_berry%kdotp_bands)
    ! loop on initial and final bands in k.p set (subset A in IAdJS19)
    do n = 1, kdotp_num_bands
      do m = 1, kdotp_num_bands

        ! zeroth order term
        if (n == m) kdotp(n, m, 1, 1, 1) = eig(pw90_berry%kdotp_bands(n))
        ! first order term
        do a = 1, 3
          kdotp(n, m, 2, a, 1) = HH_da_bar(pw90_berry%kdotp_bands(n), pw90_berry%kdotp_bands(m), a)
        end do
        ! second order term
        do a = 1, 3
          do b = 1, 3
            ! add contribution independent of other states
            kdotp(n, m, 3, a, b) = 0.5*(HH_dadb_bar(pw90_berry%kdotp_bands(n), &
                                                    pw90_berry%kdotp_bands(m), a, b))

            ! add contribution dependent on other states (subset B in IAdJS19)
            do r = 1, num_wann

              ! cycle for bands in the k.p set (subset A)
              break_loop = .false.
              do i = 1, kdotp_num_bands
                if (r == pw90_berry%kdotp_bands(i)) break_loop = .true.
              end do
              if (break_loop) cycle

              kdotp(n, m, 3, a, b) = kdotp(n, m, 3, a, b) + &
                                     0.5*HH_da_bar(pw90_berry%kdotp_bands(n), r, a) &
                                     *HH_da_bar(r, pw90_berry%kdotp_bands(m), b) &
                                     *((eig(pw90_berry%kdotp_bands(n)) - eig(r))**(-1) &
                                       + (eig(pw90_berry%kdotp_bands(m)) - eig(r))**(-1))

            end do
          end do
        end do

      enddo ! bands
    enddo ! bands

  end subroutine berry_get_kdotp

end module w90_berry
