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

! ---------------------------------------------------------------

module w90_berry
  !! This module computes various "Berry phase" related properties
  !!
  !! Key REFERENCES
  !!
  !! *  WYSV06 = PRB 74, 195118 (2006)  (anomalous Hall conductivity - AHC)
  !! *  YWVS07 = PRB 75, 195121 (2007)  (Kubo frequency-dependent conductivity)
  !! *  LVTS12 = PRB 85, 014435 (2012)  (orbital magnetization and AHC)
  !! *  CTVR06 = PRB 74, 024408 (2006)  (  "          "       )
  !! *  IATS18 = arXiv:1804.04030 (2018) (nonlinear shift current)
  !! *  QZYZ18 = PRB 98, 214402 (2018)  (spin Hall conductivity - SHC)
  ! ---------------------------------------------------------------
  !
  ! * Undocumented, works for limited purposes only:
  !                                 reading k-points and weights from file

  use w90_constants, only: dp

  implicit none

  private

  public :: berry_main, berry_get_imf_klist, berry_get_imfgh_klist, berry_get_sc_klist, &
            berry_get_shc_klist!, berry_alpha_S, berry_alpha_beta_S, berry_beta_S

  ! Pseudovector <--> Antisymmetric tensor
  !
  ! x <--> (y,z)
  ! y <--> (z,x)
  ! z <--> (x,y)
  !
  integer, dimension(3), parameter :: alpha_A = (/2, 3, 1/)
  integer, dimension(3), parameter ::  beta_A = (/3, 1, 2/)

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
  integer, dimension(6), parameter ::  beta_S = (/1, 2, 3, 2, 3, 3/)
  integer, dimension(6), parameter, public :: berry_alpha_S = alpha_S
  integer, dimension(6), parameter, public::  berry_beta_S = beta_S
!  integer,   dimension(3,3) , parameter, public::  berry_alpha_beta_S=  (/  (/1,4,5/), (/ 4,2,6 /)  , (/ 5,6,3 /)   /)
  integer, parameter, public:: berry_alpha_beta_S(3, 3) = reshape((/1, 4, 5, 4, 2, 6, 5, 6, 3/), (/3, 3/))
!(/  (/1,4,5/), (/ 4,2,6 /)  , (/ 5,6,3 /)   /)

contains

  !===========================================================!
  !                   PUBLIC PROCEDURES                       !
  !===========================================================!

  subroutine berry_main
    !============================================================!
    !                                                            !
    !! Computes the following quantities:
    !!   (i) Anomalous Hall conductivity (from Berry curvature)
    !!  (ii) Complex optical conductivity (Kubo-Greenwood) & JDOS
    !! (iii) Orbital magnetization
    !!  (iv) Nonlinear shift current
    !!   (v) Spin Hall conductivity
    !                                                            !
    !============================================================!

    use w90_constants, only: dp, cmplx_0, cmplx_i, elem_charge_SI, hbar_SI, &
      eV_au, bohr, pi, eV_seconds
    use w90_comms, only: on_root, num_nodes, my_node_id, comms_reduce
    use w90_io, only: io_error, stdout, io_file_unit, seedname, &
      io_stopwatch
    use w90_postw90_common, only: nrpts, irvec, num_int_kpts_on_node, int_kpts, &
      weight
    use w90_parameters, only: timing_level, iprint, num_wann, berry_kmesh, &
      berry_curv_adpt_kmesh, &
      berry_curv_adpt_kmesh_thresh, &
      wanint_kpoint_file, cell_volume, transl_inv, &
      berry_task, berry_curv_unit, spin_decomp, &
      kubo_nfreq, kubo_freq_list, nfermi, &
      fermi_energy_list, shc_freq_scan, &
      kubo_adpt_smr, kubo_adpt_smr_fac, &
      kubo_adpt_smr_max, kubo_smr_fixed_en_width, &
      scissors_shift, num_valence_bands, &
      shc_bandshift, shc_bandshift_firstband, shc_bandshift_energyshift
    use w90_get_oper, only: get_HH_R, get_AA_R, get_BB_R, get_CC_R, &
      get_SS_R, get_SHC_R

    real(kind=dp), allocatable    :: adkpt(:, :)

    ! AHC and orbital magnetization, calculated for a list of Fermi levels
    !
    ! First index labels J0,J1,J2 terms, second labels the Cartesian component
    !
    real(kind=dp) :: imf_k_list(3, 3, nfermi), imf_list(3, 3, nfermi), imf_list2(3, 3, nfermi)
    real(kind=dp) :: img_k_list(3, 3, nfermi), img_list(3, 3, nfermi)
    real(kind=dp) :: imh_k_list(3, 3, nfermi), imh_list(3, 3, nfermi)
    real(kind=dp) :: ahc_list(3, 3, nfermi)
    real(kind=dp) :: LCtil_list(3, 3, nfermi), ICtil_list(3, 3, nfermi), &
                     Morb_list(3, 3, nfermi)
    real(kind=dp) :: imf_k_list_dummy(3, 3, nfermi) ! adaptive refinement of AHC
    ! shift current
    real(kind=dp), allocatable :: sc_k_list(:, :, :)
    real(kind=dp), allocatable :: sc_list(:, :, :)
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
    real(kind=dp), allocatable    :: shc_fermi(:), shc_k_fermi(:)
    complex(kind=dp), allocatable :: shc_freq(:), shc_k_freq(:)
    ! for fermi energy scan, adaptive kmesh
    real(kind=dp), allocatable    :: shc_k_fermi_dummy(:)

    real(kind=dp)     :: kweight, kweight_adpt, kpt(3), kpt_ad(3), &
                         db1, db2, db3, fac, freq, rdum, vdum(3)
    integer           :: n, i, j, k, jk, ikpt, if, ispn, ierr, loop_x, loop_y, loop_z, &
                         loop_xyz, loop_adpt, adpt_counter_list(nfermi), ifreq, &
                         file_unit
    character(len=120) :: file_name
    logical           :: eval_ahc, eval_morb, eval_kubo, not_scannable, eval_sc, eval_shc
    logical           :: ladpt_kmesh
    logical           :: ladpt(nfermi)

    if (nfermi == 0) call io_error( &
      'Must specify one or more Fermi levels when berry=true')

    if (timing_level > 1 .and. on_root) call io_stopwatch('berry: prelims', 1)

    ! Mesh spacing in reduced coordinates
    !
    db1 = 1.0_dp/real(berry_kmesh(1), dp)
    db2 = 1.0_dp/real(berry_kmesh(2), dp)
    db3 = 1.0_dp/real(berry_kmesh(3), dp)

    eval_ahc = .false.
    eval_morb = .false.
    eval_kubo = .false.
    eval_sc = .false.
    eval_shc = .false.
    if (index(berry_task, 'ahc') > 0) eval_ahc = .true.
    if (index(berry_task, 'morb') > 0) eval_morb = .true.
    if (index(berry_task, 'kubo') > 0) eval_kubo = .true.
    if (index(berry_task, 'sc') > 0) eval_sc = .true.
    if (index(berry_task, 'shc') > 0) eval_shc = .true.

    ! Wannier matrix elements, allocations and initializations
    !
    if (eval_ahc) then
      call get_HH_R
      call get_AA_R
      imf_list = 0.0_dp
      adpt_counter_list = 0
    endif

    if (eval_morb) then
      call get_HH_R
      call get_AA_R
      call get_BB_R
      call get_CC_R
      imf_list2 = 0.0_dp
      img_list = 0.0_dp
      imh_list = 0.0_dp
    endif

    ! List here berry_tasks that assume nfermi=1
    !
    not_scannable = eval_kubo .or. (eval_shc .and. shc_freq_scan)
    if (not_scannable .and. nfermi .ne. 1) call io_error( &
      'The berry_task(s) you chose require that you specify a single ' &
      //'Fermi energy: scanning the Fermi energy is not implemented')

    if (eval_kubo) then
      call get_HH_R
      call get_AA_R
      allocate (kubo_H_k(3, 3, kubo_nfreq))
      allocate (kubo_H(3, 3, kubo_nfreq))
      allocate (kubo_AH_k(3, 3, kubo_nfreq))
      allocate (kubo_AH(3, 3, kubo_nfreq))
      allocate (jdos_k(kubo_nfreq))
      allocate (jdos(kubo_nfreq))
      kubo_H = cmplx_0
      kubo_AH = cmplx_0
      jdos = 0.0_dp
      if (spin_decomp) then
        call get_SS_R
        allocate (kubo_H_k_spn(3, 3, 3, kubo_nfreq))
        allocate (kubo_H_spn(3, 3, 3, kubo_nfreq))
        allocate (kubo_AH_k_spn(3, 3, 3, kubo_nfreq))
        allocate (kubo_AH_spn(3, 3, 3, kubo_nfreq))
        allocate (jdos_k_spn(3, kubo_nfreq))
        allocate (jdos_spn(3, kubo_nfreq))
        kubo_H_spn = cmplx_0
        kubo_AH_spn = cmplx_0
        jdos_spn = 0.0_dp
      endif
    endif

    if (eval_sc) then
      call get_HH_R
      call get_AA_R
      allocate (sc_k_list(3, 6, kubo_nfreq))
      allocate (sc_list(3, 6, kubo_nfreq))
      sc_k_list = 0.0_dp
      sc_list = 0.0_dp
    endif

    if (eval_shc) then
      call get_HH_R
      call get_AA_R
      call get_SS_R
      call get_SHC_R

      if (shc_freq_scan) then
        allocate (shc_freq(kubo_nfreq))
        allocate (shc_k_freq(kubo_nfreq))
        shc_freq = 0.0_dp
        shc_k_freq = 0.0_dp
      else
        allocate (shc_fermi(nfermi))
        allocate (shc_k_fermi(nfermi))
        allocate (shc_k_fermi_dummy(nfermi))
        shc_fermi = 0.0_dp
        shc_k_fermi = 0.0_dp
        !only used for fermiscan & adpt kmesh
        shc_k_fermi_dummy = 0.0_dp
        adpt_counter_list = 0
      endif
    endif

    if (on_root) then

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
        if (shc_freq_scan) then
          write (stdout, '(/,3x,a)') '  Frequency scan'
        else
          write (stdout, '(/,3x,a)') '  Fermi energy scan'
        endif
      endif

      if (transl_inv) then
        if (eval_morb) &
          call io_error('transl_inv=T disabled for morb')
        write (stdout, '(/,1x,a)') &
          'Using a translationally-invariant discretization for the'
        write (stdout, '(1x,a)') &
          'band-diagonal Wannier matrix elements of r, etc.'
      endif

      if (timing_level > 1) then
        call io_stopwatch('berry: prelims', 2)
        call io_stopwatch('berry: k-interpolation', 1)
      endif

    end if !on_root

    ! Set up adaptive refinement mesh
    !
    allocate (adkpt(3, berry_curv_adpt_kmesh**3), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating adkpt in berry')
    ikpt = 0
    !
    ! OLD VERSION (only works correctly for odd grids including original point)
    !
    ! do i=-(berry_curv_adpt_kmesh-1)/2,(berry_curv_adpt_kmesh-1)/2
    !    do j=-(berry_curv_adpt_kmesh-1)/2,(berry_curv_adpt_kmesh-1)/2
    !       do k=-(berry_curv_adpt_kmesh-1)/2,(berry_curv_adpt_kmesh-1)/2
    !          ikpt=ikpt+1
    !          adkpt(1,ikpt)=i*db1/berry_curv_adpt_kmesh
    !          adkpt(2,ikpt)=j*db2/berry_curv_adpt_kmesh
    !          adkpt(3,ikpt)=k*db3/berry_curv_adpt_kmesh
    !       end do
    !    end do
    ! end do
    !
    ! NEW VERSION (both even and odd grids)
    !
    do i = 0, berry_curv_adpt_kmesh - 1
      do j = 0, berry_curv_adpt_kmesh - 1
        do k = 0, berry_curv_adpt_kmesh - 1
          ikpt = ikpt + 1
          adkpt(1, ikpt) = db1*((i + 0.5_dp)/berry_curv_adpt_kmesh - 0.5_dp)
          adkpt(2, ikpt) = db2*((j + 0.5_dp)/berry_curv_adpt_kmesh - 0.5_dp)
          adkpt(3, ikpt) = db3*((k + 0.5_dp)/berry_curv_adpt_kmesh - 0.5_dp)
        end do
      end do
    end do

    ! Loop over interpolation k-points
    !
    if (wanint_kpoint_file) then

      ! NOTE: still need to specify berry_kmesh in the input file
      !
      !        - Must use the correct nominal value in order to
      !          correctly set up adaptive smearing in kubo

      if (on_root) write (stdout, '(/,1x,a,i10,a)') &
        'Reading interpolation grid from file kpoint.dat: ', &
        sum(num_int_kpts_on_node), ' points'

      ! Loop over k-points on the irreducible wedge of the Brillouin
      ! zone, read from file 'kpoint.dat'
      !
      do loop_xyz = 1, num_int_kpts_on_node(my_node_id)
        kpt(:) = int_kpts(:, loop_xyz)
        kweight = weight(loop_xyz)
        kweight_adpt = kweight/berry_curv_adpt_kmesh**3
        !               .
        ! ***BEGIN COPY OF CODE BLOCK 1***
        !
        if (eval_ahc) then
          call berry_get_imf_klist(kpt, imf_k_list)
          ladpt = .false.
          do if = 1, nfermi
            vdum(1) = sum(imf_k_list(:, 1, if))
            vdum(2) = sum(imf_k_list(:, 2, if))
            vdum(3) = sum(imf_k_list(:, 3, if))
            if (berry_curv_unit == 'bohr2') vdum = vdum/bohr**2
            rdum = sqrt(dot_product(vdum, vdum))
            if (rdum > berry_curv_adpt_kmesh_thresh) then
              adpt_counter_list(if) = adpt_counter_list(if) + 1
              ladpt(if) = .true.
            else
              imf_list(:, :, if) = imf_list(:, :, if) + imf_k_list(:, :, if)*kweight
            endif
          enddo
          if (any(ladpt)) then
            do loop_adpt = 1, berry_curv_adpt_kmesh**3
              ! Using imf_k_list here would corrupt values for other
              ! frequencies, hence dummy. Only if-th element is used
              call berry_get_imf_klist(kpt(:) + adkpt(:, loop_adpt), &
                                       imf_k_list_dummy, ladpt=ladpt)
              do if = 1, nfermi
                if (ladpt(if)) then
                  imf_list(:, :, if) = imf_list(:, :, if) &
                                       + imf_k_list_dummy(:, :, if)*kweight_adpt
                endif
              enddo
            end do
          endif
        end if

        if (eval_morb) then
          call berry_get_imfgh_klist(kpt, imf_k_list, img_k_list, imh_k_list)
          imf_list2 = imf_list2 + imf_k_list*kweight
          img_list = img_list + img_k_list*kweight
          imh_list = imh_list + imh_k_List*kweight
        endif

        if (eval_kubo) then
          if (spin_decomp) then
            call berry_get_kubo_k(kpt, kubo_H_k, kubo_AH_k, jdos_k, &
                                  kubo_H_k_spn, kubo_AH_k_spn, jdos_k_spn)
          else
            call berry_get_kubo_k(kpt, kubo_H_k, kubo_AH_k, jdos_k)
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
          call berry_get_sc_klist(kpt, sc_k_list)
          sc_list = sc_list + sc_k_list*kweight
        end if

        !
        ! ***END COPY OF CODE BLOCK 1***

        if (eval_shc) then
          ! print calculation progress, from 0%, 10%, ... to 100%
          ! Note the 1st call to berry_get_shc_klist will be much longer
          ! than later calls due to the time spent on
          !   berry_get_shc_klist -> wham_get_eig_deleig ->
          !   pw90common_fourier_R_to_k -> ws_translate_dist
          call berry_print_progress(loop_xyz, 1, num_int_kpts_on_node(my_node_id), 1)
          if (.not. shc_freq_scan) then
            call berry_get_shc_klist(kpt, shc_k_fermi=shc_k_fermi)
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
            if (berry_curv_adpt_kmesh > 1) then
              do if = 1, nfermi
                rdum = abs(shc_k_fermi(if))
                if (berry_curv_unit == 'bohr2') rdum = rdum/bohr**2
                if (rdum > berry_curv_adpt_kmesh_thresh) then
                  adpt_counter_list(1) = adpt_counter_list(1) + 1
                  ladpt_kmesh = .true.
                  exit
                endif
              enddo
            else
              ladpt_kmesh = .false.
            end if
            if (ladpt_kmesh) then
              do loop_adpt = 1, berry_curv_adpt_kmesh**3
                !Using shc_k here would corrupt values for other
                !kpt, hence dummy. Only if-th element is used.
                call berry_get_shc_klist(kpt(:) + adkpt(:, loop_adpt), &
                                         shc_k_fermi=shc_k_fermi_dummy)
                shc_fermi = shc_fermi + kweight_adpt*shc_k_fermi_dummy
              end do
            else
              shc_fermi = shc_fermi + kweight*shc_k_fermi
            end if
          else ! freq_scan, no adaptive kmesh
            call berry_get_shc_klist(kpt, shc_k_freq=shc_k_freq)
            shc_freq = shc_freq + kweight*shc_k_freq
          end if
        end if

      end do !loop_xyz

    else ! Do not read 'kpoint.dat'. Loop over a regular grid in the full BZ

      kweight = db1*db2*db3
      kweight_adpt = kweight/berry_curv_adpt_kmesh**3

      do loop_xyz = my_node_id, PRODUCT(berry_kmesh) - 1, num_nodes
        loop_x = loop_xyz/(berry_kmesh(2)*berry_kmesh(3))
        loop_y = (loop_xyz - loop_x*(berry_kmesh(2) &
                                     *berry_kmesh(3)))/berry_kmesh(3)
        loop_z = loop_xyz - loop_x*(berry_kmesh(2)*berry_kmesh(3)) &
                 - loop_y*berry_kmesh(3)
        kpt(1) = loop_x*db1
        kpt(2) = loop_y*db2
        kpt(3) = loop_z*db3

        ! ***BEGIN CODE BLOCK 1***
        !
        if (eval_ahc) then
          call berry_get_imf_klist(kpt, imf_k_list)
          ladpt = .false.
          do if = 1, nfermi
            vdum(1) = sum(imf_k_list(:, 1, if))
            vdum(2) = sum(imf_k_list(:, 2, if))
            vdum(3) = sum(imf_k_list(:, 3, if))
            if (berry_curv_unit == 'bohr2') vdum = vdum/bohr**2
            rdum = sqrt(dot_product(vdum, vdum))
            if (rdum > berry_curv_adpt_kmesh_thresh) then
              adpt_counter_list(if) = adpt_counter_list(if) + 1
              ladpt(if) = .true.
            else
              imf_list(:, :, if) = imf_list(:, :, if) + imf_k_list(:, :, if)*kweight
            endif
          enddo
          if (any(ladpt)) then
            do loop_adpt = 1, berry_curv_adpt_kmesh**3
              ! Using imf_k_list here would corrupt values for other
              ! frequencies, hence dummy. Only if-th element is used
              call berry_get_imf_klist(kpt(:) + adkpt(:, loop_adpt), &
                                       imf_k_list_dummy, ladpt=ladpt)
              do if = 1, nfermi
                if (ladpt(if)) then
                  imf_list(:, :, if) = imf_list(:, :, if) &
                                       + imf_k_list_dummy(:, :, if)*kweight_adpt
                endif
              enddo
            end do
          endif
        end if

        if (eval_morb) then
          call berry_get_imfgh_klist(kpt, imf_k_list, img_k_list, imh_k_list)
          imf_list2 = imf_list2 + imf_k_list*kweight
          img_list = img_list + img_k_list*kweight
          imh_list = imh_list + imh_k_List*kweight
        endif

        if (eval_kubo) then
          if (spin_decomp) then
            call berry_get_kubo_k(kpt, kubo_H_k, kubo_AH_k, jdos_k, &
                                  kubo_H_k_spn, kubo_AH_k_spn, jdos_k_spn)
          else
            call berry_get_kubo_k(kpt, kubo_H_k, kubo_AH_k, jdos_k)
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
          call berry_get_sc_klist(kpt, sc_k_list)
          sc_list = sc_list + sc_k_list*kweight
        end if

        !
        ! ***END CODE BLOCK 1***

        if (eval_shc) then
          ! print calculation progress, from 0%, 10%, ... to 100%
          ! Note the 1st call to berry_get_shc_klist will be much longer
          ! than later calls due to the time spent on
          !   berry_get_shc_klist -> wham_get_eig_deleig ->
          !   pw90common_fourier_R_to_k -> ws_translate_dist
          call berry_print_progress(loop_xyz, my_node_id, PRODUCT(berry_kmesh) - 1, num_nodes)
          if (.not. shc_freq_scan) then
            call berry_get_shc_klist(kpt, shc_k_fermi=shc_k_fermi)
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
            if (berry_curv_adpt_kmesh > 1) then
              do if = 1, nfermi
                rdum = abs(shc_k_fermi(if))
                if (berry_curv_unit == 'bohr2') rdum = rdum/bohr**2
                if (rdum > berry_curv_adpt_kmesh_thresh) then
                  adpt_counter_list(1) = adpt_counter_list(1) + 1
                  ladpt_kmesh = .true.
                  exit
                endif
              enddo
            else
              ladpt_kmesh = .false.
            end if
            if (ladpt_kmesh) then
              do loop_adpt = 1, berry_curv_adpt_kmesh**3
                !Using shc_k here would corrupt values for other
                !kpt, hence dummy. Only if-th element is used.
                call berry_get_shc_klist(kpt(:) + adkpt(:, loop_adpt), &
                                         shc_k_fermi=shc_k_fermi_dummy)
                shc_fermi = shc_fermi + kweight_adpt*shc_k_fermi_dummy
              end do
            else
              shc_fermi = shc_fermi + kweight*shc_k_fermi
            end if
          else ! freq_scan, no adaptive kmesh
            call berry_get_shc_klist(kpt, shc_k_freq=shc_k_freq)
            shc_freq = shc_freq + kweight*shc_k_freq
          end if
        end if

      end do !loop_xyz

    end if !wanint_kpoint_file

    ! Collect contributions from all nodes
    !
    if (eval_ahc) then
      call comms_reduce(imf_list(1, 1, 1), 3*3*nfermi, 'SUM')
      call comms_reduce(adpt_counter_list(1), nfermi, 'SUM')
    endif

    if (eval_morb) then
      call comms_reduce(imf_list2(1, 1, 1), 3*3*nfermi, 'SUM')
      call comms_reduce(img_list(1, 1, 1), 3*3*nfermi, 'SUM')
      call comms_reduce(imh_list(1, 1, 1), 3*3*nfermi, 'SUM')
    end if

    if (eval_kubo) then
      call comms_reduce(kubo_H(1, 1, 1), 3*3*kubo_nfreq, 'SUM')
      call comms_reduce(kubo_AH(1, 1, 1), 3*3*kubo_nfreq, 'SUM')
      call comms_reduce(jdos(1), kubo_nfreq, 'SUM')
      if (spin_decomp) then
        call comms_reduce(kubo_H_spn(1, 1, 1, 1), 3*3*3*kubo_nfreq, 'SUM')
        call comms_reduce(kubo_AH_spn(1, 1, 1, 1), 3*3*3*kubo_nfreq, 'SUM')
        call comms_reduce(jdos_spn(1, 1), 3*kubo_nfreq, 'SUM')
      endif
    endif

    if (eval_sc) then
      call comms_reduce(sc_list(1, 1, 1), 3*6*kubo_nfreq, 'SUM')
    end if

    if (eval_shc) then
      if (shc_freq_scan) then
        call comms_reduce(shc_freq(1), kubo_nfreq, 'SUM')
      else
        call comms_reduce(shc_fermi(1), nfermi, 'SUM')
        call comms_reduce(adpt_counter_list(1), nfermi, 'SUM')
      end if
    end if

    if (on_root) then

      if (timing_level > 1) call io_stopwatch('berry: k-interpolation', 2)
      write (stdout, '(1x,a)') ' '
      if (eval_ahc .and. berry_curv_adpt_kmesh .ne. 1) then
        if (.not. wanint_kpoint_file) write (stdout, '(1x,a28,3(i0,1x))') &
          'Regular interpolation grid: ', berry_kmesh
        write (stdout, '(1x,a28,3(i0,1x))') 'Adaptive refinement grid: ', &
          berry_curv_adpt_kmesh, berry_curv_adpt_kmesh, berry_curv_adpt_kmesh
        if (berry_curv_unit == 'ang2') then
          write (stdout, '(1x,a28,a17,f6.2,a)') &
            'Refinement threshold: ', 'Berry curvature >', &
            berry_curv_adpt_kmesh_thresh, ' Ang^2'
        elseif (berry_curv_unit == 'bohr2') then
          write (stdout, '(1x,a28,a17,f6.2,a)') &
            'Refinement threshold: ', 'Berry curvature >', &
            berry_curv_adpt_kmesh_thresh, ' bohr^2'
        endif
        if (nfermi == 1) then
          if (wanint_kpoint_file) then
            write (stdout, '(1x,a30,i5,a,f5.2,a)') &
              ' Points triggering refinement: ', &
              adpt_counter_list(1), '(', &
              100*real(adpt_counter_list(1), dp) &
              /sum(num_int_kpts_on_node), '%)'
          else
            write (stdout, '(1x,a30,i5,a,f5.2,a)') &
              ' Points triggering refinement: ', &
              adpt_counter_list(1), '(', &
              100*real(adpt_counter_list(1), dp)/product(berry_kmesh), '%)'
          endif
        endif
      elseif (eval_shc) then
        if (berry_curv_adpt_kmesh .ne. 1) then
          if (.not. wanint_kpoint_file) write (stdout, '(1x,a28,3(i0,1x))') &
            'Regular interpolation grid: ', berry_kmesh
          if (.not. shc_freq_scan) then
            write (stdout, '(1x,a28,3(i0,1x))') &
              'Adaptive refinement grid: ', &
              berry_curv_adpt_kmesh, berry_curv_adpt_kmesh, berry_curv_adpt_kmesh
            if (berry_curv_unit == 'ang2') then
              write (stdout, '(1x,a28,f12.2,a)') &
                'Refinement threshold: ', &
                berry_curv_adpt_kmesh_thresh, ' Ang^2'
            elseif (berry_curv_unit == 'bohr2') then
              write (stdout, '(1x,a28,f12.2,a)') &
                'Refinement threshold: ', &
                berry_curv_adpt_kmesh_thresh, ' bohr^2'
            endif
            if (wanint_kpoint_file) then
              write (stdout, '(1x,a30,i8,a,f6.2,a)') &
                ' Points triggering refinement: ', adpt_counter_list(1), '(', &
                100*real(adpt_counter_list(1), dp)/sum(num_int_kpts_on_node), '%)'
            else
              write (stdout, '(1x,a30,i8,a,f6.2,a)') &
                ' Points triggering refinement: ', adpt_counter_list(1), '(', &
                100*real(adpt_counter_list(1), dp)/product(berry_kmesh), '%)'
            endif
          endif
        else
          if (.not. wanint_kpoint_file) write (stdout, '(1x,a20,3(i0,1x))') &
            'Interpolation grid: ', berry_kmesh(1:3)
        endif
        write (stdout, '(a)') ''
        if (kubo_adpt_smr) then
          write (stdout, '(1x,a)') 'Using adaptive smearing'
          write (stdout, '(7x,a,f8.3)') 'adaptive smearing prefactor ', kubo_adpt_smr_fac
          write (stdout, '(7x,a,f8.3,a)') 'adaptive smearing max width ', kubo_adpt_smr_max, ' eV'
        else
          write (stdout, '(1x,a)') 'Using fixed smearing'
          write (stdout, '(7x,a,f8.3,a)') 'fixed smearing width ', &
            kubo_smr_fixed_en_width, ' eV'
        endif
        write (stdout, '(a)') ''
        if (abs(scissors_shift) > 1.0e-7_dp) then
          write (stdout, '(1X,A,I0,A,G18.10,A)') "Using scissors_shift to shift energy bands with index > ", &
            num_valence_bands, " by ", scissors_shift, " eV."
        endif
        if (shc_bandshift) then
          write (stdout, '(1X,A,I0,A,G18.10,A)') "Using shc_bandshift to shift energy bands with index >= ", &
            shc_bandshift_firstband, " by ", shc_bandshift_energyshift, " eV."
        endif
      else
        if (.not. wanint_kpoint_file) write (stdout, '(1x,a20,3(i0,1x))') &
          'Interpolation grid: ', berry_kmesh(1:3)
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
        ! ===========================
        ! fac = -e^2/(hbar.V_c*10^-8)
        ! ===========================
        !
        ! with 'V_c' in Angstroms^3, and 'e', 'hbar' in SI units
        ! --------------------------------------------------------------------
        !
        fac = -1.0e8_dp*elem_charge_SI**2/(hbar_SI*cell_volume)
        ahc_list(:, :, :) = imf_list(:, :, :)*fac
        if (nfermi > 1) then
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
        do if = 1, nfermi
          if (nfermi > 1) write (file_unit, '(4(F12.6,1x))') &
            fermi_energy_list(if), sum(ahc_list(:, 1, if)), &
            sum(ahc_list(:, 2, if)), sum(ahc_list(:, 3, if))
          write (stdout, '(/,1x,a18,F10.4)') 'Fermi energy (ev):', &
            fermi_energy_list(if)
          if (nfermi > 1) then
            if (wanint_kpoint_file) then
              write (stdout, '(1x,a30,i5,a,f5.2,a)') &
                ' Points triggering refinement: ', &
                adpt_counter_list(if), '(', &
                100*real(adpt_counter_list(if), dp) &
                /sum(num_int_kpts_on_node), '%)'
            else
              write (stdout, '(1x,a30,i5,a,f5.2,a)') &
                ' Points triggering refinement: ', &
                adpt_counter_list(if), '(', &
                100*real(adpt_counter_list(if), dp) &
                /product(berry_kmesh), '%)'
            endif
          endif
          write (stdout, '(/,1x,a)') &
            'AHC (S/cm)       x          y          z'
          if (iprint > 1) then
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
        if (nfermi > 1) close (file_unit)
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
        fac = -eV_au/bohr**2
        if (nfermi > 1) then
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
        do if = 1, nfermi
          LCtil_list(:, :, if) = (img_list(:, :, if) &
                                  - fermi_energy_list(if)*imf_list2(:, :, if))*fac
          ICtil_list(:, :, if) = (imh_list(:, :, if) &
                                  - fermi_energy_list(if)*imf_list2(:, :, if))*fac
          Morb_list(:, :, if) = LCtil_list(:, :, if) + ICtil_list(:, :, if)
          if (nfermi > 1) write (file_unit, '(4(F12.6,1x))') &
            fermi_energy_list(if), sum(Morb_list(1:3, 1, if)), &
            sum(Morb_list(1:3, 2, if)), sum(Morb_list(1:3, 3, if))
          write (stdout, '(/,/,1x,a,F12.6)') 'Fermi energy (ev) =', &
            fermi_energy_list(if)
          write (stdout, '(/,/,1x,a)') &
            'M_orb (bohr magn/cell)        x          y          z'
          if (iprint > 1) then
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
        if (nfermi > 1) close (file_unit)
      endif

      ! -----------------------------!
      ! Complex optical conductivity !
      ! -----------------------------!
      !
      if (eval_kubo) then
        !
        ! Convert to S/cm
        fac = 1.0e8_dp*elem_charge_SI**2/(hbar_SI*cell_volume)
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
          do ifreq = 1, kubo_nfreq
            if (spin_decomp) then
              write (file_unit, '(9E16.8)') real(kubo_freq_list(ifreq), dp), &
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
              write (file_unit, '(3E16.8)') real(kubo_freq_list(ifreq), dp), &
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
          do ifreq = 1, kubo_nfreq
            if (spin_decomp) then
              write (file_unit, '(9E16.8)') real(kubo_freq_list(ifreq), dp), &
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
              write (file_unit, '(3E16.8)') real(kubo_freq_list(ifreq), dp), &
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
        do ifreq = 1, kubo_nfreq
          if (spin_decomp) then
            write (file_unit, '(5E16.8)') real(kubo_freq_list(ifreq), dp), &
              jdos(ifreq), jdos_spn(:, ifreq)
          else
            write (file_unit, '(2E16.8)') real(kubo_freq_list(ifreq), dp), &
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
        ! ===========================
        ! fac = eV_seconds.( pi.e^3/(4.hbar^2.V_c) )
        ! ===========================
        !
        ! with 'V_c' in Angstroms^3, and 'e', 'hbar' in SI units
        ! --------------------------------------------------------------------

        fac = eV_seconds*pi*elem_charge_SI**3/(4*hbar_SI**(2)*cell_volume)
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
            do ifreq = 1, kubo_nfreq
              write (file_unit, '(2E18.8E3)') real(kubo_freq_list(ifreq), dp), &
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
        fac = 1.0e8_dp*elem_charge_SI**2/(hbar_SI*cell_volume)/2.0_dp
        if (shc_freq_scan) then
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
        if (.not. shc_freq_scan) then
          file_name = trim(seedname)//'-shc-fermiscan'//'.dat'
        else
          file_name = trim(seedname)//'-shc-freqscan'//'.dat'
        endif
        file_name = trim(file_name)
        file_unit = io_file_unit()
        write (stdout, '(/,3x,a)') '* '//file_name
        open (file_unit, FILE=file_name, STATUS='UNKNOWN', FORM='FORMATTED')
        if (.not. shc_freq_scan) then
          write (file_unit, '(a,3x,a,3x,a)') &
            '#No.', 'Fermi energy(eV)', 'SHC((hbar/e)*S/cm)'
          do n = 1, nfermi
            write (file_unit, '(I4,1x,F12.6,1x,E17.8)') &
              n, fermi_energy_list(n), shc_fermi(n)
          enddo
        else
          write (file_unit, '(a,3x,a,3x,a,3x,a)') '#No.', 'Frequency(eV)', &
            'Re(sigma)((hbar/e)*S/cm)', 'Im(sigma)((hbar/e)*S/cm)'
          do n = 1, kubo_nfreq
            write (file_unit, '(I4,1x,F12.6,1x,1x,2(E17.8,1x))') n, &
              real(kubo_freq_list(n), dp), real(shc_freq(n), dp), aimag(shc_freq(n))
          enddo
        endif
        close (file_unit)

      endif

    end if !on_root

  end subroutine berry_main

  subroutine berry_get_imf_klist(kpt, imf_k_list, occ, ladpt)
    !============================================================!
    !                                                            !
    !! Calculates the Berry curvature traced over the occupied
    !! states, -2Im[f(k)] [Eq.33 CTVR06, Eq.6 LVTS12] for a list
    !! of Fermi energies, and stores it in axial-vector form
    !                                                            !
    !============================================================!
    ! Arguments
    !
    real(kind=dp), intent(in)                    :: kpt(3)
    real(kind=dp), intent(out), dimension(:, :, :) :: imf_k_list
    real(kind=dp), intent(in), optional, dimension(:) :: occ
    logical, intent(in), optional, dimension(:) :: ladpt

    if (present(occ)) then
      call berry_get_imfgh_klist(kpt, imf_k_list, occ=occ)
    else
      if (present(ladpt)) then
        call berry_get_imfgh_klist(kpt, imf_k_list, ladpt=ladpt)
      else
        call berry_get_imfgh_klist(kpt, imf_k_list)
      endif
    endif

  end subroutine berry_get_imf_klist

  subroutine berry_get_imfgh_klist(kpt, imf_k_list, img_k_list, imh_k_list, occ, ladpt)
    !=========================================================!
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
    !=========================================================!

    use w90_constants, only: dp, cmplx_0, cmplx_i
    use w90_utility, only: utility_re_tr_prod, utility_im_tr_prod
    use w90_parameters, only: num_wann, nfermi
    use w90_postw90_common, only: pw90common_fourier_R_to_k_vec, pw90common_fourier_R_to_k
    use w90_wan_ham, only: wham_get_eig_UU_HH_JJlist, wham_get_occ_mat_list
    use w90_get_oper, only: AA_R, BB_R, CC_R
    use w90_utility, only: utility_zgemm_new

    ! Arguments
    !
    real(kind=dp), intent(in)     :: kpt(3)
    real(kind=dp), intent(out), dimension(:, :, :), optional &
      :: imf_k_list, img_k_list, imh_k_list
    real(kind=dp), intent(in), optional, dimension(:) :: occ
    logical, intent(in), optional, dimension(:) :: ladpt

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
    real(kind=dp)                 :: eig(num_wann)
    integer                       :: i, j, ife, nfermi_loc
    real(kind=dp)                 :: s
    logical                       :: todo(nfermi)

    ! Temporary space for matrix products
    complex(kind=dp), allocatable, dimension(:, :, :) :: tmp

    if (present(occ)) then
      nfermi_loc = 1
    else
      nfermi_loc = nfermi
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
      call wham_get_eig_UU_HH_JJlist(kpt, eig, UU, HH, JJp_list, JJm_list, occ=occ)
      call wham_get_occ_mat_list(UU, f_list, g_list, occ=occ)
    else
      call wham_get_eig_UU_HH_JJlist(kpt, eig, UU, HH, JJp_list, JJm_list)
      call wham_get_occ_mat_list(UU, f_list, g_list, eig=eig)
    endif

    call pw90common_fourier_R_to_k_vec(kpt, AA_R, OO_true=AA, OO_pseudo=OOmega)

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

      call pw90common_fourier_R_to_k_vec(kpt, BB_R, OO_true=BB)
      do j = 1, 3
        do i = 1, j
          call pw90common_fourier_R_to_k(kpt, CC_R(:, :, :, i, j), CC(:, :, i, j), 0)
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

  !===========================================================!
  !                   PRIVATE PROCEDURES                      !
  !===========================================================!

  subroutine berry_get_kubo_k(kpt, kubo_H_k, kubo_AH_k, jdos_k, &
                              kubo_H_k_spn, kubo_AH_k_spn, jdos_k_spn)
    !====================================================================!
    !                                                                    !
    !! Contribution from point k to the complex interband optical
    !! conductivity, separated into Hermitian (H) and anti-Hermitian (AH)
    !! parts. Also returns the joint density of states
    !                                                                    !
    !====================================================================!

    use w90_constants, only: dp, cmplx_0, cmplx_i, pi
    use w90_utility, only: utility_diagonalize, utility_rotate, utility_w0gauss
    use w90_parameters, only: num_wann, kubo_nfreq, kubo_freq_list, &
      fermi_energy_list, kubo_eigval_max, &
      kubo_adpt_smr, kubo_smr_fixed_en_width, &
      kubo_adpt_smr_max, kubo_adpt_smr_fac, &
      kubo_smr_index, berry_kmesh, spin_decomp
    use w90_postw90_common, only: pw90common_get_occ, pw90common_fourier_R_to_k_new, &
      pw90common_fourier_R_to_k_vec, pw90common_kmesh_spacing
    use w90_wan_ham, only: wham_get_D_h, wham_get_eig_deleig
    use w90_get_oper, only: HH_R, AA_R
    use w90_spin, only: spin_get_nk

    ! Arguments
    !
    ! Last three arguments should be present iff spin_decomp=T (but
    ! this is not checked: do it?)
    !
    real(kind=dp), intent(in)  :: kpt(3)
    complex(kind=dp), dimension(:, :, :), intent(out) :: kubo_H_k
    complex(kind=dp), dimension(:, :, :), intent(out) :: kubo_AH_k
    real(kind=dp), dimension(:), intent(out) :: jdos_k
    complex(kind=dp), optional, dimension(:, :, :, :), intent(out) :: kubo_H_k_spn
    complex(kind=dp), optional, dimension(:, :, :, :), intent(out) :: kubo_AH_k_spn
    real(kind=dp), optional, dimension(:, :), intent(out) :: jdos_k_spn

    complex(kind=dp), allocatable :: HH(:, :)
    complex(kind=dp), allocatable :: delHH(:, :, :)
    complex(kind=dp), allocatable :: UU(:, :)
    complex(kind=dp), allocatable :: D_h(:, :, :)
    complex(kind=dp), allocatable :: AA(:, :, :)

    ! Adaptive smearing
    !
    real(kind=dp)    :: del_eig(num_wann, 3), joint_level_spacing, &
                        eta_smr, Delta_k, arg, vdum(3)

    integer          :: i, j, n, m, ifreq, ispn
    real(kind=dp)    :: eig(num_wann), occ(num_wann), delta, &
                        rfac1, rfac2, occ_prod, spn_nk(num_wann)
    complex(kind=dp) :: cfac, omega

    allocate (HH(num_wann, num_wann))
    allocate (delHH(num_wann, num_wann, 3))
    allocate (UU(num_wann, num_wann))
    allocate (D_h(num_wann, num_wann, 3))
    allocate (AA(num_wann, num_wann, 3))

    if (kubo_adpt_smr) then
      call wham_get_eig_deleig(kpt, eig, del_eig, HH, delHH, UU)
      Delta_k = pw90common_kmesh_spacing(berry_kmesh)
    else
      call pw90common_fourier_R_to_k_new(kpt, HH_R, OO=HH, &
                                         OO_dx=delHH(:, :, 1), &
                                         OO_dy=delHH(:, :, 2), &
                                         OO_dz=delHH(:, :, 3))
      call utility_diagonalize(HH, num_wann, eig, UU)
    endif
    call pw90common_get_occ(eig, occ, fermi_energy_list(1))
    call wham_get_D_h(delHH, UU, eig, D_h)

    call pw90common_fourier_R_to_k_vec(kpt, AA_R, OO_true=AA)
    do i = 1, 3
      AA(:, :, i) = utility_rotate(AA(:, :, i), UU, num_wann)
    enddo
    AA = AA + cmplx_i*D_h ! Eq.(25) WYSV06

    ! Replace imaginary part of frequency with a fixed value
    if (.not. kubo_adpt_smr .and. kubo_smr_fixed_en_width /= 0.0_dp) &
      kubo_freq_list = real(kubo_freq_list, dp) &
                       + cmplx_i*kubo_smr_fixed_en_width

    kubo_H_k = cmplx_0
    kubo_AH_k = cmplx_0
    jdos_k = 0.0_dp
    if (spin_decomp) then
      call spin_get_nk(kpt, spn_nk)
      kubo_H_k_spn = cmplx_0
      kubo_AH_k_spn = cmplx_0
      jdos_k_spn = 0.0_dp
    end if
    do m = 1, num_wann
      do n = 1, num_wann
        if (n == m) cycle
        if (eig(m) > kubo_eigval_max .or. eig(n) > kubo_eigval_max) cycle
        if (spin_decomp) then
          if (spn_nk(n) >= 0 .and. spn_nk(m) >= 0) then
            ispn = 1 ! up --> up transition
          elseif (spn_nk(n) < 0 .and. spn_nk(m) < 0) then
            ispn = 2 ! down --> down
          else
            ispn = 3 ! spin-flip
          end if
        end if
        if (kubo_adpt_smr) then
          ! Eq.(35) YWVS07
          vdum(:) = del_eig(m, :) - del_eig(n, :)
          joint_level_spacing = sqrt(dot_product(vdum(:), vdum(:)))*Delta_k
          eta_smr = min(joint_level_spacing*kubo_adpt_smr_fac, &
                        kubo_adpt_smr_max)
        else
          eta_smr = kubo_smr_fixed_en_width
        endif
        rfac1 = (occ(m) - occ(n))*(eig(m) - eig(n))
        occ_prod = occ(n)*(1.0_dp - occ(m))
        do ifreq = 1, kubo_nfreq
          !
          ! Complex frequency for the anti-Hermitian conductivity
          !
          if (kubo_adpt_smr) then
            omega = real(kubo_freq_list(ifreq), dp) + cmplx_i*eta_smr
          else
            omega = kubo_freq_list(ifreq)
          endif
          !
          ! Broadened delta function for the Hermitian conductivity and JDOS
          !
          arg = (eig(m) - eig(n) - real(omega, dp))/eta_smr
          ! If only Hermitean part were computed, could speed up
          ! by inserting here 'if(abs(arg)>10.0_dp) cycle'
          delta = utility_w0gauss(arg, kubo_smr_index)/eta_smr
          !
          ! Lorentzian shape (for testing purposes)
!             delta=1.0_dp/(1.0_dp+arg*arg)/pi
!             delta=delta/eta_smr
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

  subroutine berry_get_sc_klist(kpt, sc_k_list)
    !====================================================================!
    !                                                                    !
    !  Contribution from point k to the nonlinear shift current
    !  [integrand of Eq.8 IATS18]
    !  Notation correspondence with IATS18:
    !  AA_da_bar              <-->   \mathbbm{b}
    !  AA_bar                 <-->   \mathbbm{a}
    !  HH_dadb_bar            <-->   \mathbbm{w}
    !  D_h(n,m)               <-->   \mathbbm{v}_{nm}/(E_{m}-E_{n})
    !  sum_AD                 <-->   summatory of Eq. 32 IATS18
    !  sum_HD                 <-->   summatory of Eq. 30 IATS18
    !  eig_da(n)-eig_da(m)    <-->   \mathbbm{Delta}_{nm}
    !                                                                    !
    !====================================================================!

    ! Arguments
    !
    use w90_constants, only: dp, cmplx_0, cmplx_i
    use w90_utility, only: utility_re_tr, utility_im_tr, utility_w0gauss, utility_w0gauss_vec
    use w90_parameters, only: num_wann, nfermi, kubo_nfreq, kubo_freq_list, fermi_energy_list, &
      kubo_smr_index, berry_kmesh, kubo_adpt_smr_fac, &
      kubo_adpt_smr_max, kubo_adpt_smr, kubo_eigval_max, &
      kubo_smr_fixed_en_width, sc_phase_conv, sc_w_thr
    use w90_postw90_common, only: pw90common_fourier_R_to_k_vec_dadb, &
      pw90common_fourier_R_to_k_new_second_d, pw90common_get_occ, &
      pw90common_kmesh_spacing, pw90common_fourier_R_to_k_vec_dadb_TB_conv
    use w90_wan_ham, only: wham_get_eig_UU_HH_JJlist, wham_get_occ_mat_list, wham_get_D_h, &
      wham_get_eig_UU_HH_AA_sc, wham_get_eig_deleig, wham_get_D_h_P_value, &
      wham_get_eig_deleig_TB_conv, wham_get_eig_UU_HH_AA_sc_TB_conv
    use w90_get_oper, only: AA_R
    use w90_utility, only: utility_rotate, utility_zdotu
    ! Arguments
    !
    real(kind=dp), intent(in)                        :: kpt(3)
    real(kind=dp), intent(out), dimension(:, :, :)     :: sc_k_list

    complex(kind=dp), allocatable :: UU(:, :)
    complex(kind=dp), allocatable :: AA(:, :, :), AA_bar(:, :, :)
    complex(kind=dp), allocatable :: AA_da(:, :, :, :), AA_da_bar(:, :, :, :)
    complex(kind=dp), allocatable :: HH_da(:, :, :), HH_da_bar(:, :, :)
    complex(kind=dp), allocatable :: HH_dadb(:, :, :, :), HH_dadb_bar(:, :, :, :)
    complex(kind=dp), allocatable :: HH(:, :)
    complex(kind=dp), allocatable :: D_h(:, :, :)
    real(kind=dp), allocatable    :: eig(:)
    real(kind=dp), allocatable    :: eig_da(:, :)
    real(kind=dp), allocatable    :: occ(:)

    complex(kind=dp)              :: sum_AD(3, 3), sum_HD(3, 3), r_mn(3), gen_r_nm(3)
    integer                       :: i, if, a, b, c, bc, n, m, r, ifreq, istart, iend
    real(kind=dp)                 :: I_nm(3, 6), &
                                     omega(kubo_nfreq), delta(kubo_nfreq), joint_level_spacing, &
                                     eta_smr, Delta_k, arg, vdum(3), occ_fac, wstep, wmin, wmax

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
    allocate (eig(num_wann))
    allocate (occ(num_wann))
    allocate (eig_da(num_wann, 3))

    ! Initialize shift current array at point k
    sc_k_list = 0.d0

    ! Gather W-gauge matrix objects !

    ! choose the convention for the FT sums
    if (sc_phase_conv .eq. 1) then ! use Wannier centres in the FT exponentials (so called TB convention)
      ! get Hamiltonian and its first and second derivatives
      ! Note that below we calculate the UU matrix--> we have to use the same UU from here on for
      ! maintaining the gauge-covariance of the whole matrix element
      call wham_get_eig_UU_HH_AA_sc_TB_conv(kpt, eig, UU, HH, HH_da, HH_dadb)
      ! get position operator and its derivative
      ! note that AA_da(:,:,a,b) \propto \sum_R exp(iRk)*iR_{b}*<0|r_{a}|R>
      call pw90common_fourier_R_to_k_vec_dadb_TB_conv(kpt, AA_R, OO_da=AA, OO_dadb=AA_da)
      ! get eigenvalues and their k-derivatives
      call wham_get_eig_deleig_TB_conv(kpt, eig, eig_da, HH_da, UU)
    elseif (sc_phase_conv .eq. 2) then ! do not use Wannier centres in the FT exponentials (usual W90 convention)
      ! same as above
      call wham_get_eig_UU_HH_AA_sc(kpt, eig, UU, HH, HH_da, HH_dadb)
      call pw90common_fourier_R_to_k_vec_dadb(kpt, AA_R, OO_da=AA, OO_dadb=AA_da)
      call wham_get_eig_deleig(kpt, eig, eig_da, HH, HH_da, UU)
    end if

    ! get electronic occupations
    call pw90common_get_occ(eig, occ, fermi_energy_list(1))

    ! get D_h (Eq. (24) WYSV06)
    call wham_get_D_h_P_value(HH_da, UU, eig, D_h)

    ! calculate k-spacing in case of adaptive smearing
    if (kubo_adpt_smr) Delta_k = pw90common_kmesh_spacing(berry_kmesh)

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
    omega = real(kubo_freq_list(:), dp)
    wmin = omega(1)
    wmax = omega(kubo_nfreq)
    wstep = omega(2) - omega(1)

    ! loop on initial and final bands
    do n = 1, num_wann
      do m = 1, num_wann
        ! cycle diagonal matrix elements and bands above the maximum
        if (n == m) cycle
        if (eig(m) > kubo_eigval_max .or. eig(n) > kubo_eigval_max) cycle
        ! setup T=0 occupation factors
        occ_fac = (occ(n) - occ(m))
        if (abs(occ_fac) < 1e-10) cycle

        ! set delta function smearing
        if (kubo_adpt_smr) then
          vdum(:) = eig_da(m, :) - eig_da(n, :)
          joint_level_spacing = sqrt(dot_product(vdum(:), vdum(:)))*Delta_k
          eta_smr = min(joint_level_spacing*kubo_adpt_smr_fac, &
                        kubo_adpt_smr_max)
        else
          eta_smr = kubo_smr_fixed_en_width
        endif

        ! restrict to energy window spanning [-sc_w_thr*eta_smr,+sc_w_thr*eta_smr]
        ! outside this range, the two delta functions are virtually zero
        if (((eig(n) - eig(m) + sc_w_thr*eta_smr < wmin) .or. (eig(n) - eig(m) - sc_w_thr*eta_smr > wmax)) .and. &
            ((eig(m) - eig(n) + sc_w_thr*eta_smr < wmin) .or. (eig(m) - eig(n) - sc_w_thr*eta_smr > wmax))) cycle

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
        r_mn(:) = AA_bar(m, n, :) + cmplx_i*D_h(m, n, :)

        ! loop over direction of generalized derivative
        do a = 1, 3
          ! store generalized derivative as an array on the additional spatial index,
          ! its composed of 8 terms in total, see Eq (34) combined with (30) and
          ! (32) of IATS18
          gen_r_nm(:) = (AA_da_bar(n, m, :, a) &
                         + ((AA_bar(n, n, :) - AA_bar(m, m, :))*D_h(n, m, a) + &
                            (AA_bar(n, n, a) - AA_bar(m, m, a))*D_h(n, m, :)) &
                         - cmplx_i*AA_bar(n, m, :)*(AA_bar(n, n, a) - AA_bar(m, m, a)) &
                         + sum_AD(:, a) &
                         + cmplx_i*(HH_dadb_bar(n, m, :, a) &
                                    + sum_HD(:, a) &
                                    + (D_h(n, m, :)*(eig_da(n, a) - eig_da(m, a)) + &
                                       D_h(n, m, a)*(eig_da(n, :) - eig_da(m, :)))) &
                         /(eig(m) - eig(n)))

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
        istart = max(int((eig(n) - eig(m) - sc_w_thr*eta_smr - wmin)/wstep + 1), 1)
        iend = min(int((eig(n) - eig(m) + sc_w_thr*eta_smr - wmin)/wstep + 1), kubo_nfreq)
        ! multiply matrix elements with delta function for the relevant frequencies
        if (istart <= iend) then
          delta = 0.0
          delta(istart:iend) = &
            utility_w0gauss_vec((eig(m) - eig(n) + omega(istart:iend))/eta_smr, kubo_smr_index)/eta_smr
          call DGER(18, iend - istart + 1, occ_fac, I_nm, 1, delta(istart:iend), 1, sc_k_list(:, :, istart:iend), 18)
        endif
        ! same for delta(E_mn-w)
        istart = max(int((eig(m) - eig(n) - sc_w_thr*eta_smr - wmin)/wstep + 1), 1)
        iend = min(int((eig(m) - eig(n) + sc_w_thr*eta_smr - wmin)/wstep + 1), kubo_nfreq)
        if (istart <= iend) then
          delta = 0.0
          delta(istart:iend) = &
            utility_w0gauss_vec((eig(n) - eig(m) + omega(istart:iend))/eta_smr, kubo_smr_index)/eta_smr
          call DGER(18, iend - istart + 1, occ_fac, I_nm, 1, delta(istart:iend), 1, sc_k_list(:, :, istart:iend), 18)
        endif

      enddo ! bands
    enddo ! bands

  end subroutine berry_get_sc_klist

  subroutine berry_get_shc_klist(kpt, shc_k_fermi, shc_k_freq, shc_k_band)
    !====================================================================!
    !                                                                    !
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
    !   Junfeng Qiao (18/8/2018)                                         !
    !====================================================================!

    use w90_constants, only: dp, cmplx_0, cmplx_i
    use w90_utility, only: utility_rotate
    use w90_parameters, only: num_wann, kubo_eigval_max, kubo_nfreq, &
      kubo_freq_list, kubo_adpt_smr, kubo_smr_fixed_en_width, &
      kubo_adpt_smr_max, kubo_adpt_smr_fac, berry_kmesh, &
      fermi_energy_list, nfermi, shc_alpha, shc_beta, shc_gamma, &
      shc_bandshift, shc_bandshift_firstband, shc_bandshift_energyshift
    use w90_postw90_common, only: pw90common_get_occ, &
      pw90common_fourier_R_to_k_vec, pw90common_kmesh_spacing
    use w90_wan_ham, only: wham_get_D_h, wham_get_eig_deleig
    use w90_get_oper, only: AA_R
    !use w90_comms, only: my_node_id
    !!!

    ! args
    real(kind=dp), intent(in)  :: kpt(3)
    real(kind=dp), optional, intent(out) :: shc_k_fermi(nfermi)
    complex(kind=dp), optional, intent(out) :: shc_k_freq(kubo_nfreq)
    real(kind=dp), optional, intent(out) :: shc_k_band(num_wann)

    ! internal vars
    logical                       :: lfreq, lfermi, lband
    complex(kind=dp), allocatable :: HH(:, :)
    complex(kind=dp), allocatable :: delHH(:, :, :)
    complex(kind=dp), allocatable :: UU(:, :)
    complex(kind=dp), allocatable :: D_h(:, :, :)
    complex(kind=dp), allocatable :: AA(:, :, :)

    complex(kind=dp)              :: js_k(num_wann, num_wann)

    ! Adaptive smearing
    !
    real(kind=dp)    :: del_eig(num_wann, 3), joint_level_spacing, &
                        eta_smr, Delta_k, vdum(3)

    integer          :: n, m, i, ifreq
    real(kind=dp)    :: eig(num_wann)
    real(kind=dp)    :: occ_fermi(num_wann, nfermi), occ_freq(num_wann)
    complex(kind=dp) :: omega_list(kubo_nfreq)
    real(kind=dp)    :: omega, rfac
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

    call wham_get_eig_deleig(kpt, eig, del_eig, HH, delHH, UU)
    call wham_get_D_h(delHH, UU, eig, D_h)

    ! Here I apply a scissor operator to the conduction bands, if required in the input
    if (shc_bandshift) then
      eig(shc_bandshift_firstband:) = eig(shc_bandshift_firstband:) + shc_bandshift_energyshift
    end if

    call pw90common_fourier_R_to_k_vec(kpt, AA_R, OO_true=AA)
    do i = 1, 3
      AA(:, :, i) = utility_rotate(AA(:, :, i), UU, num_wann)
    enddo
    AA = AA + cmplx_i*D_h ! Eq.(25) WYSV06

    call berry_get_js_k(kpt, eig, del_eig(:, shc_alpha), &
                        D_h(:, :, shc_alpha), UU, js_k)

    ! adpt_smr only works with berry_kmesh, so do not use
    ! adpt_smr in kpath or kslice plots.
    if (kubo_adpt_smr) then
      Delta_k = pw90common_kmesh_spacing(berry_kmesh)
    endif
    if (lfreq) then
      call pw90common_get_occ(eig, occ_freq, fermi_energy_list(1))
    elseif (lfermi) then
      ! get occ for different fermi_energy
      do i = 1, nfermi
        call pw90common_get_occ(eig, occ_fermi(:, i), fermi_energy_list(i))
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
        if (eig(m) > kubo_eigval_max .or. eig(n) > kubo_eigval_max) cycle

        rfac = eig(m) - eig(n)
        !this will calculate AHC
        !prod = -rfac*cmplx_i*AA(n, m, shc_alpha) * rfac*cmplx_i*AA(m, n, shc_beta)
        prod = js_k(n, m)*cmplx_i*rfac*AA(m, n, shc_beta)
        if (kubo_adpt_smr) then
          ! Eq.(35) YWVS07
          vdum(:) = del_eig(m, :) - del_eig(n, :)
          joint_level_spacing = sqrt(dot_product(vdum(:), vdum(:)))*Delta_k
          eta_smr = min(joint_level_spacing*kubo_adpt_smr_fac, &
                        kubo_adpt_smr_max)
        else
          eta_smr = kubo_smr_fixed_en_width
        endif
        if (lfreq) then
          do ifreq = 1, kubo_nfreq
            cdum = real(kubo_freq_list(ifreq), dp) + cmplx_i*eta_smr
            cfac = -2.0_dp/(rfac**2 - cdum**2)
            omega_list(ifreq) = omega_list(ifreq) + cfac*aimag(prod)
          end do
        else if (lfermi .or. lband) then
          rfac = -2.0_dp/(rfac**2 + eta_smr**2)
          omega = omega + rfac*aimag(prod)
        end if
      enddo

      if (lfermi) then
        do i = 1, nfermi
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

    !===========================================================!
    !                   PRIVATE PROCEDURES                      !
    !===========================================================!
    subroutine berry_get_js_k(kpt, eig, del_alpha_eig, D_alpha_h, UU, js_k)
      !====================================================================!
      !                                                                    !
      ! Contribution from point k to the
      !   <psi_k | 1/2*(sigma_gamma*v_alpha + v_alpha*sigma_gamma) | psi_k>
      !
      !  QZYZ18 Eq.(23) without hbar/2 (required by spin operator) and
      !  not divided by hbar (required by velocity operator)
      !
      !  Junfeng Qiao (8/7/2018)
      !                                                                    !
      !====================================================================!

      use w90_constants, only: dp, cmplx_0, cmplx_i
      use w90_utility, only: utility_rotate
      use w90_parameters, only: num_wann, shc_alpha, shc_gamma
      use w90_postw90_common, only: pw90common_fourier_R_to_k_new, &
        pw90common_fourier_R_to_k_vec
      use w90_get_oper, only: SS_R, SR_R, SHR_R, SH_R

      ! args
      real(kind=dp), intent(in)  :: kpt(3)
      real(kind=dp), dimension(:), intent(in)  :: eig
      real(kind=dp), dimension(:), intent(in)  :: del_alpha_eig
      complex(kind=dp), dimension(:, :), intent(in)  :: D_alpha_h
      complex(kind=dp), dimension(:, :), intent(in)  :: UU
      complex(kind=dp), dimension(:, :), intent(out) :: js_k

      ! internal vars
      complex(kind=dp)    :: B_k(num_wann, num_wann)
      complex(kind=dp)    :: K_k(num_wann, num_wann)
      complex(kind=dp)    :: L_k(num_wann, num_wann)
      complex(kind=dp)    :: S_w(num_wann, num_wann)
      complex(kind=dp)    :: S_k(num_wann, num_wann)
      complex(kind=dp)    :: SR_w(num_wann, num_wann, 3)
      complex(kind=dp)    :: SR_alpha_k(num_wann, num_wann)
      complex(kind=dp)    :: SHR_w(num_wann, num_wann, 3)
      complex(kind=dp)    :: SHR_alpha_k(num_wann, num_wann)
      complex(kind=dp)    :: SH_w(num_wann, num_wann, 3)
      complex(kind=dp)    :: SH_k(num_wann, num_wann)
      complex(kind=dp)    :: eig_mat(num_wann, num_wann)
      complex(kind=dp)    :: del_eig_mat(num_wann, num_wann)

      !===========
      js_k = cmplx_0

      !=========== S_k ===========
      ! < u_k | sigma_gamma | u_k >, QZYZ18 Eq.(25)
      ! QZYZ18 Eq.(36)
      call pw90common_fourier_R_to_k_new(kpt, SS_R(:, :, :, shc_gamma), OO=S_w)
      ! QZYZ18 Eq.(30)
      S_k = utility_rotate(S_w, UU, num_wann)

      !=========== K_k ===========
      ! < u_k | sigma_gamma | \partial_alpha u_k >, QZYZ18 Eq.(26)
      ! QZYZ18 Eq.(37)
      call pw90common_fourier_R_to_k_vec(kpt, SR_R(:, :, :, shc_gamma, :), OO_true=SR_w)
      ! QZYZ18 Eq.(31)
      SR_alpha_k = -cmplx_i*utility_rotate(SR_w(:, :, shc_alpha), UU, num_wann)
      K_k = SR_alpha_k + matmul(S_k, D_alpha_h)

      !=========== L_k ===========
      ! < u_k | sigma_gamma.H | \partial_alpha u_k >, QZYZ18 Eq.(27)
      ! QZYZ18 Eq.(38)
      call pw90common_fourier_R_to_k_vec(kpt, SHR_R(:, :, :, shc_gamma, :), OO_true=SHR_w)
      ! QZYZ18 Eq.(32)
      SHR_alpha_k = -cmplx_i*utility_rotate(SHR_w(:, :, shc_alpha), UU, num_wann)
      ! QZYZ18 Eq.(39)
      call pw90common_fourier_R_to_k_vec(kpt, SH_R, OO_true=SH_w)
      ! QZYZ18 Eq.(32)
      SH_k = utility_rotate(SH_w(:, :, shc_gamma), UU, num_wann)
      L_k = SHR_alpha_k + matmul(SH_k, D_alpha_h)

      !=========== B_k ===========
      ! < \psi_nk | sigma_gamma v_alpha | \psi_mk >, QZYZ18 Eq.(24)
      B_k = cmplx_0
      do i = 1, num_wann
        eig_mat(i, :) = eig(:)
        del_eig_mat(i, :) = del_alpha_eig(:)
      end do
      ! note * is not matmul
      B_k = del_eig_mat*S_k + eig_mat*K_k - L_k

      !=========== js_k ===========
      ! QZYZ18 Eq.(23)
      ! note the S in SR_R,SHR_R,SH_R of get_SHC_R is sigma,
      ! to get spin current, we need to multiply it by hbar/2,
      ! also we need to divide it by hbar to recover the velocity
      ! operator, these are done outside of this subroutine
      js_k = 1.0_dp/2.0_dp*(B_k + conjg(transpose(B_k)))

    end subroutine berry_get_js_k

  end subroutine berry_get_shc_klist

  subroutine berry_print_progress(loop_k, start_k, end_k, step_k)
    !============================================================!
    ! print k-points calculation progress, seperated into 11 points,
    ! from 0%, 10%, ... to 100%
    ! start_k, end_k are inclusive
    ! loop_k should in the array start_k to end_k with step step_k
    !============================================================!
    use w90_comms, only: on_root
    use w90_io, only: stdout, io_wallclocktime

    integer, intent(in) :: loop_k, start_k, end_k, step_k

    real(kind=dp) :: cur_time, finished
    real(kind=dp), save :: prev_time
    integer :: i, j, n, last_k
    logical, dimension(9) :: kmesh_processed = (/(.false., i=1, 9)/)

    if (on_root) then
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
    end if ! on_root

  end subroutine berry_print_progress

end module w90_berry
