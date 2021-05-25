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

module w90_gyrotropic
  !! This module computes various "gyrotropic" effects
  !! as described in :
  !!    TAS17 =  arXiv:1710.03204 (2017) Gyrotropic effects in trigonal tellurium studied from first principles
  !!                   S.S.Tsirkin, P. Aguado Puente, I. Souza
  ! ---------------------------------------------------------------

  use w90_constants, only: dp
  use w90_berry, only: berry_get_imf_klist, berry_get_imfgh_klist
  implicit none

  private

  public :: gyrotropic_main

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
  !integer, dimension(6), parameter :: alpha_S = (/1, 2, 3, 1, 1, 2/)
  !integer, dimension(6), parameter ::  beta_S = (/1, 2, 3, 2, 3, 3/)

contains

  !===========================================================!
  !                   PUBLIC PROCEDURES                       !
  !===========================================================!

  subroutine gyrotropic_main(num_bands, num_kpts, num_wann, nrpts, irvec, ndegen, rpt_origin, &
                             eigval, real_lattice, crvec, u_matrix, v_matrix, dis_data, &
                             kmesh_info, k_points, param_input, pw90_common, postw90_oper, &
                             pw90_ham, comm, physics, stdout, seedname, gyrotropic, berry, fermi, &
                             cell_volume)
    !============================================================!
    !                                                            !
    !! Computes the following quantities:
    !!   (i) D tensor
    !!  (ii) K tensor
    !! (iii) C tensor
    !!  (iv) current-induced optical activity
    !!   (v) natural optical activity
    !                                                            !
    !============================================================!

    use pw90_parameters, only: gyrotropic_type, berry_type, postw90_oper_type, postw90_ham_type, &
      postw90_common_type
    use w90_comms, only: comms_reduce, w90commtype, mpirank, mpisize
    use w90_constants, only: dp, cmplx_0, twopi, pw90_physical_constants
    use w90_get_oper, only: get_HH_R, get_AA_R, get_BB_R, get_CC_R, get_SS_R
    use w90_io, only: io_error, io_file_unit, io_stopwatch
    use w90_parameters, only: fermi_data_type
    use w90_param_types, only: disentangle_type, k_point_type, parameter_input_type, kmesh_info_type
    use w90_utility, only: utility_det3

    implicit none

    ! passed variables
    integer, intent(in) :: num_bands, num_kpts, num_wann, nrpts, stdout
    integer, intent(inout) :: irvec(:, :), ndegen(:), rpt_origin

    real(kind=dp), intent(in) :: eigval(:, :), real_lattice(3, 3)
    real(kind=dp), intent(inout) :: crvec(:, :)
    real(kind=dp), intent(in) :: cell_volume

    complex(kind=dp), intent(in) :: u_matrix(:, :, :), v_matrix(:, :, :)

    type(disentangle_type), intent(in) :: dis_data
    type(kmesh_info_type), intent(in) :: kmesh_info
    type(k_point_type), intent(in) :: k_points
    type(parameter_input_type), intent(in) :: param_input
    type(postw90_common_type), intent(in) :: pw90_common
    type(postw90_ham_type), intent(in) :: pw90_ham
    type(postw90_oper_type), intent(in) :: postw90_oper
    type(pw90_physical_constants), intent(in) :: physics
    type(w90commtype), intent(in) :: comm
    type(gyrotropic_type), intent(in) :: gyrotropic
    type(berry_type), intent(in) :: berry
    type(fermi_data_type), intent(in) :: fermi

    character(len=50), intent(in) :: seedname

    ! local variables
    real(kind=dp), allocatable    :: gyro_K_spn(:, :, :)
    real(kind=dp), allocatable    :: gyro_DOS(:)
    real(kind=dp), allocatable    :: gyro_K_orb(:, :, :)
    real(kind=dp), allocatable    :: gyro_C(:, :, :)
    real(kind=dp), allocatable    :: gyro_D(:, :, :)
    real(kind=dp), allocatable    :: gyro_Dw(:, :, :, :)
    real(kind=dp), allocatable    :: gyro_NOA_spn(:, :, :, :)
    real(kind=dp), allocatable    :: gyro_NOA_orb(:, :, :, :)

    character(len=30) :: f_out_name_tmp
    character(len=30) :: units_tmp
    character(len=120) :: comment_tmp

    real(kind=dp)     :: kweight, kpt(3), &
                         db1, db2, db3, fac
    integer           :: loop_x, loop_y, loop_z, loop_xyz
    logical           :: eval_K, eval_C, eval_D, eval_Dw, eval_NOA, eval_spn, eval_DOS

    integer :: my_node_id, num_nodes
    my_node_id = mpirank(comm)
    num_nodes = mpisize(comm)

    if (fermi%n == 0) call io_error( &
      'Must specify one or more Fermi levels when gyrotropic=true', stdout, seedname)

    if (param_input%timing_level > 1 .and. param_input%iprint > 0) &
      call io_stopwatch('gyrotropic: prelims', 1, stdout, seedname)

    ! Mesh spacing in reduced coordinates
    !
    db1 = 1.0_dp/real(gyrotropic%kmesh(1), dp)
    db2 = 1.0_dp/real(gyrotropic%kmesh(2), dp)
    db3 = 1.0_dp/real(gyrotropic%kmesh(3), dp)

    eval_K = .false.
    eval_C = .false.
    eval_D = .false.
    eval_Dw = .false.
    eval_spn = .false.
    eval_NOA = .false.
    eval_DOS = .false.

    if (index(gyrotropic%task, '-k') > 0) eval_K = .true.
    if (index(gyrotropic%task, '-c') > 0) eval_C = .true.
    if (index(gyrotropic%task, '-d0') > 0) eval_D = .true.
    if (index(gyrotropic%task, '-dw') > 0) eval_Dw = .true.
    if (index(gyrotropic%task, '-spin') > 0) eval_spn = .true.
    if (index(gyrotropic%task, '-noa') > 0) eval_NOA = .true.
    if (index(gyrotropic%task, '-dos') > 0) eval_DOS = .true.
    if (index(gyrotropic%task, 'all') > 0) then
      eval_K = .true.
      eval_C = .true.
      eval_D = .true.
      eval_Dw = .true.
      if (param_input%spinors) eval_spn = .true.
      eval_NOA = .true.
      eval_DOS = .true.
    endif

    if (.not. (eval_K .or. eval_noa)) eval_spn = .false.

    if ((.not. param_input%spinors) .and. eval_spn) call io_error( &
      "spin contribution requested for gyrotropic, but the wavefunctions are not spinors", &
      stdout, seedname)

    ! Wannier matrix elements, allocations and initializations

    call get_HH_R(num_bands, num_kpts, num_wann, nrpts, ndegen, irvec, crvec, real_lattice, &
                  rpt_origin, eigval, u_matrix, v_matrix, dis_data, k_points, param_input, &
                  pw90_common, stdout, seedname, comm)
    if (eval_D .or. eval_Dw .or. eval_K .or. eval_NOA) then

      call get_AA_R(num_bands, num_kpts, num_wann, nrpts, irvec, eigval, v_matrix, berry, &
                    dis_data, kmesh_info, k_points, param_input, pw90_common, stdout, seedname, &
                    comm)
    endif

    if (eval_spn) then

      call get_SS_R(num_bands, num_kpts, num_wann, nrpts, irvec, eigval, v_matrix, dis_data, &
                    k_points, param_input, postw90_oper, stdout, seedname, comm)
    endif

    if (eval_K) then
      call get_BB_R(num_bands, num_kpts, num_wann, nrpts, irvec, eigval, v_matrix, dis_data, &
                    kmesh_info, k_points, param_input, pw90_common, stdout, seedname, comm)
      call get_CC_R(num_bands, num_kpts, num_wann, nrpts, irvec, eigval, v_matrix, dis_data, &
                    kmesh_info, k_points, param_input, postw90_oper, pw90_common, stdout, &
                    seedname, comm)
      allocate (gyro_K_orb(3, 3, fermi%n))
      gyro_K_orb = 0.0_dp
      if (eval_spn) then
        allocate (gyro_K_spn(3, 3, fermi%n))
        gyro_K_spn = 0.0_dp
      endif
    endif

    if (eval_D) then
      allocate (gyro_D(3, 3, fermi%n))
      gyro_D = 0.0_dp
    endif

    if (eval_DOS) then
      allocate (gyro_DOS(fermi%n))
      gyro_DOS = 0.0_dp
    endif

    if (eval_C) then
      allocate (gyro_C(3, 3, fermi%n))
      gyro_C = 0.0_dp
    endif

    if (eval_Dw) then
      allocate (gyro_Dw(3, 3, fermi%n, gyrotropic%nfreq))
      gyro_Dw = 0.0_dp
    endif

    if (eval_NOA) then
      allocate (gyro_NOA_orb(3, 3, fermi%n, gyrotropic%nfreq))
      gyro_NOA_orb = 0.0_dp
      if (eval_spn) then
        allocate (gyro_NOA_spn(3, 3, fermi%n, gyrotropic%nfreq))
        gyro_NOA_spn = 0.0_dp
      endif
    endif

    if (param_input%iprint > 0) then
      flush(stdout)
      write (stdout, '(/,/,1x,a)') 'Properties calculated in module  g y r o t r o p i c'
      write (stdout, '(1x,a)') '------------------------------------------'

      if (eval_D) write (stdout, '(/,3x,a)') '* D-tensor  --- Eq.2 of TAS17 '

      if (eval_dos) write (stdout, '(/,3x,a)') '* density of states '

      if (eval_K) then
        write (stdout, '(/,3x,a)') '* K-tensor  --- Eq.3 of TAS17 '
        if (eval_spn) then
          write (stdout, '(3x,a)') '    * including spin component '
        else
          write (stdout, '(3x,a)') '    * excluding spin component '
        endif
      endif

      if (eval_Dw) write (stdout, '(/,3x,a)') '* Dw-tensor  --- Eq.12 of TAS17 '

      if (eval_C) write (stdout, '(/,3x,a)') '* C-tensor  --- Eq.B6 of TAS17 '

      if (eval_NOA) then
        write (stdout, '(/,3x,a)') '* gamma-tensor of NOA --- Eq.C12 of TAS17 '
        if (eval_spn) then
          write (stdout, '(3x,a)') '    * including spin component '
        else
          write (stdout, '(3x,a)') '    * excluding spin component '
        endif
      endif

      if (berry%transl_inv) then
        if (eval_K) &
          call io_error('transl_inv=T disabled for K-tensor', stdout, seedname)
        write (stdout, '(/,1x,a)') &
          'Using a translationally-invariant discretization for the'
        write (stdout, '(1x,a)') &
          'band-diagonal Wannier matrix elements of r, etc.'
      endif

      if (param_input%timing_level > 1) then
        call io_stopwatch('gyrotropic: prelims', 2, stdout, seedname)
        call io_stopwatch('gyrotropic: k-interpolation', 1, stdout, seedname)
      endif

      write (stdout, '(1x,a20,3(i0,1x))') 'Interpolation grid: ', gyrotropic%kmesh(1:3)

      flush(stdout)

    end if ! param_input%iprint >0, aka "on_root"

    ! Do not read 'kpoint.dat'. Loop over a regular grid in the full BZ

    kweight = db1*db2*db3*utility_det3(gyrotropic%box)

    do loop_xyz = my_node_id, PRODUCT(gyrotropic%kmesh) - 1, num_nodes
      loop_x = loop_xyz/(gyrotropic%kmesh(2)*gyrotropic%kmesh(3))
      loop_y = (loop_xyz - loop_x*(gyrotropic%kmesh(2) &
                                   *gyrotropic%kmesh(3)))/gyrotropic%kmesh(3)
      loop_z = loop_xyz - loop_x*(gyrotropic%kmesh(2)*gyrotropic%kmesh(3)) &
               - loop_y*gyrotropic%kmesh(3)
      kpt(1) = loop_x*db1
      kpt(2) = loop_y*db2
      kpt(3) = loop_z*db3
      kpt(:) = gyrotropic%box_corner(:) + matmul(kpt, gyrotropic%box)

      call gyrotropic_get_k_list(num_bands, num_kpts, num_wann, nrpts, irvec, ndegen, rpt_origin, &
                                 eigval, real_lattice, crvec, u_matrix, v_matrix, dis_data, &
                                 k_points, param_input, pw90_common, comm, pw90_ham, kpt, stdout, &
                                 seedname, kweight, gyro_K_spn, gyro_K_orb, gyro_D, gyro_Dw, &
                                 gyro_C, gyro_DOS, gyro_NOA_orb, gyro_NOA_spn, eval_K, eval_D, &
                                 eval_Dw, eval_NOA, eval_spn, eval_C, eval_dos, gyrotropic, fermi)

    end do !loop_xyz

    ! Collect contributions from all nodes
    !
    if (eval_K) then
      call comms_reduce(gyro_K_orb(1, 1, 1), 3*3*fermi%n, 'SUM', stdout, seedname, comm)
      if (eval_spn) call comms_reduce(gyro_K_spn(1, 1, 1), 3*3*fermi%n, 'SUM', stdout, &
                                      seedname, comm)
    endif

    if (eval_D) &
      call comms_reduce(gyro_D(1, 1, 1), 3*3*fermi%n, 'SUM', stdout, seedname, comm)

    if (eval_C) &
      call comms_reduce(gyro_C(1, 1, 1), 3*3*fermi%n, 'SUM', stdout, seedname, comm)

    if (eval_Dw) &
      call comms_reduce(gyro_Dw(1, 1, 1, 1), 3*3*fermi%n*gyrotropic%nfreq, 'SUM', stdout, &
                        seedname, comm)

    if (eval_dos) &
      call comms_reduce(gyro_DOS(1), fermi%n, 'SUM', stdout, seedname, comm)

    if (eval_NOA) then
      call comms_reduce(gyro_NOA_orb(1, 1, 1, 1), 3*3*fermi%n*gyrotropic%nfreq, 'SUM', stdout, &
                        seedname, comm)
      if (eval_spn) call comms_reduce(gyro_NOA_spn(1, 1, 1, 1), 3*3*fermi%n*gyrotropic%nfreq, &
                                      'SUM', stdout, seedname, comm)
    endif

    if (param_input%iprint > 0) then

      if (param_input%timing_level > 1) call io_stopwatch('gyrotropic: k-interpolation', 2, &
                                                          stdout, seedname)
      write (stdout, '(1x,a)') ' '
      write (stdout, *) 'Calculation finished, writing results'
      flush(stdout)

      if (eval_K) then
        if (eval_spn) then
          ! At this point gme_spn_list contains
          ! (1/N) sum_k delta(E_kn-E_f).(d E_{kn}/d k_i).sigma_{kn,j}
          ! (units of length) in Angstroms.
          !        ====================================
          ! To get K in units of Ampere do the following:
          !        ====================================
          !   * Divide by V_c in Ang^3 to get a quantity with units of [L]^{-2}
          !   * Multiply by 10^20 to convert to SI
          !   * Multiply by -g_s.e.hbar/(4m_e) \simeq e.hbar/(2.m_e) in SI units
          ! ==============================
          ! fac = 10^20*e*hbar/(2.m_e.V_c)
          ! ==============================
          fac = -1.0e20_dp*physics%elem_charge_SI*physics%hbar_SI/(2.*physics%elec_mass_SI &
                                                                   *cell_volume)
          gyro_K_spn(:, :, :) = gyro_K_spn(:, :, :)*fac
          f_out_name_tmp = 'K_spin'
          units_tmp = "Ampere"
          comment_tmp = "spin part of the K tensor -- Eq. 3 of TAS17"
          call gyrotropic_outprint_tensor(stdout, seedname, gyrotropic, fermi, f_out_name_tmp, &
                                          arrEf=gyro_K_spn, units=units_tmp, comment=comment_tmp)
        endif  ! eval_K && eval_spin

        ! At this point gme_orb_list contains
        ! (1/N)sum_{k,n} delta(E_kn-E_f).(d E_{kn}/d k_i)
        !                           .Im[<del_k u_kn| x (H_k-E_kn)|del_k u_kn>]
        ! (units of energy times length^3) in eV.Ang^3.
        !        ====================================
        ! To get K  in units of Ampere do the following:
        !        ====================================
        !   * Divide by V_c in Ang^3 to get a quantity with units of eV
        !   * Multiply by 'e' in SI to convert to SI (Joules)
        !   * Multiply by e/(2.hbar) to get K in Ampere
        ! ====================================
        ! fac = e^2/(2.hbar.V_c)
        ! ====================================
        fac = physics%elem_charge_SI**2/(2.*physics%hbar_SI*cell_volume)
        gyro_K_orb(:, :, :) = gyro_K_orb(:, :, :)*fac

        f_out_name_tmp = 'K_orb'
        units_tmp = "Ampere"
        comment_tmp = "orbital part of the K tensor -- Eq. 3 of TAS17"
        call gyrotropic_outprint_tensor(stdout, seedname, gyrotropic, fermi, f_out_name_tmp, &
                                        arrEf=gyro_K_orb, units=units_tmp, comment=comment_tmp)
      endif ! eval_K

      if (eval_D) then
        fac = 1./cell_volume
        gyro_D(:, :, :) = gyro_D(:, :, :)*fac

        f_out_name_tmp = 'D'
        units_tmp = "dimensionless"
        comment_tmp = "the D tensor -- Eq. 2 of TAS17"
        call gyrotropic_outprint_tensor(stdout, seedname, gyrotropic, fermi, f_out_name_tmp, &
                                        arrEf=gyro_D, units=units_tmp, comment=comment_tmp)
      endif

      if (eval_Dw) then
        fac = 1./cell_volume
        gyro_Dw(:, :, :, :) = gyro_Dw(:, :, :, :)*fac

        f_out_name_tmp = 'tildeD'
        units_tmp = "dimensionless"
        comment_tmp = "the tildeD tensor -- Eq. 12 of TAS17"
        call gyrotropic_outprint_tensor(stdout, seedname, gyrotropic, fermi, f_out_name_tmp, &
                                        arrEfW=gyro_Dw, units=units_tmp, comment=comment_tmp)
      endif

      if (eval_C) then
        ! At this point gyro_C contains
        ! (1/N)sum_{k,n} delta(E_kn-E_f).(d E_{kn}/d k_i).(d E_{kn}/d k_j)
        ! (units of energy*length^2) in eV*Ang^2
        !
        ! To get it in Cab = e/h * (1/N*V_cell)sum_{k,n} delta(E_kn-E_f).(d E_{kn}/d k_i).(d E_{kn}/d k_j)
        ! in units Ampere/cm
        !
        ! divide by V_c in Ang^3 to get  eV/Ang
        ! multiply by 10^8*e in SI to get J/cm
        ! multiply by e/h in SI
        !
        fac = 1.0e+8_dp*physics%elem_charge_SI**2/(twopi*physics%hbar_SI*cell_volume)
        gyro_C(:, :, :) = gyro_C(:, :, :)*fac

        f_out_name_tmp = 'C'
        units_tmp = "Ampere/cm"
        comment_tmp = "the C tensor -- Eq. B6 of TAS17"
        call gyrotropic_outprint_tensor(stdout, seedname, gyrotropic, fermi, f_out_name_tmp, &
                                        arrEf=gyro_C, units=units_tmp, comment=comment_tmp)
      endif

      if (eval_noa) then
        ! at this point gyro_NOA_orb  is in eV^-1.Ang^3   !
        !  We want the result in angstrems  !
        !   * Divide by V_c in Ang^3 to make it eV^{-1}
        !   * Divide by e in SI to get J^{-1}
        !   * multiply by e^2/eps_0 to get meters
        !   *multiply dy 1e10 to get Ang
        fac = 1e+10_dp*physics%elem_charge_SI/(cell_volume*physics%eps0_SI)
        gyro_NOA_orb = gyro_NOA_orb*fac
        f_out_name_tmp = 'NOA_orb'
        units_tmp = "Ang"
        comment_tmp = "the tensor $gamma_{abc}^{orb}$ (Eq. C12,C14 of TAS17)"
        call gyrotropic_outprint_tensor(stdout, seedname, gyrotropic, fermi, f_out_name_tmp, &
                                        arrEfW=gyro_NOA_orb, units=units_tmp, comment=comment_tmp, &
                                        symmetrize=.false.)

        if (eval_spn) then
          ! at this point gyro_NOA_spn  is in eV^-2.Ang   !
          !  We want the result in angstrems  !
          !   * Divide by V_c in Ang^3 to make it (eV.Ang)^{-2}
          !   * multiply by 1e20/e^2 in SI to get (J.m)^{-2}
          !   * multiply by e^2/eps_0 to get (J.m)^{-1}
          !   *multiply dy hbar^2/m_e to get m
          !   *multiply by 1e10 to get Ang
          fac = 1e+30_dp*physics%hbar_SI**2/(cell_volume*physics%eps0_SI*physics%elec_mass_SI)
          gyro_NOA_spn = gyro_NOA_spn*fac
          f_out_name_tmp = 'NOA_spin'
          units_tmp = "Ang"
          comment_tmp = "the tensor $gamma_{abc}^{spin}$ (Eq. C12,C15 of TAS17)"
          call gyrotropic_outprint_tensor(stdout, seedname, gyrotropic, fermi, f_out_name_tmp, &
                                          arrEfW=gyro_NOA_spn, units=units_tmp, &
                                          comment=comment_tmp, symmetrize=.false.)
        endif
      endif  !eval_NOA

      if (eval_DOS) then
        ! At this point gyro_C contains
        ! (1/N)sum_{k,n} delta(E_kn-E_f)
        ! in units of eV^{-1}
        ! divide by V_c in Ang^3 to get units 1./(eV*Ang^3)
        gyro_DOS(:) = gyro_DOS(:)/cell_volume
        f_out_name_tmp = 'DOS'
        units_tmp = "eV^{-1}.Ang^{-3}"
        comment_tmp = "density of states"
        call gyrotropic_outprint_tensor(stdout, seedname, gyrotropic, fermi, f_out_name_tmp, &
                                        arrEf1d=gyro_DOS, units=units_tmp, comment=comment_tmp)
      endif

    end if !param_input%iprint >0, aka "on_root"

  end subroutine gyrotropic_main

  subroutine gyrotropic_get_k_list(num_bands, num_kpts, num_wann, nrpts, irvec, ndegen, &
                                   rpt_origin, eigval, real_lattice, crvec, u_matrix, v_matrix, &
                                   dis_data, k_points, param_input, pw90_common, comm, pw90_ham, &
                                   kpt, stdout, seedname, kweight, gyro_K_spn, gyro_K_orb, gyro_D, &
                                   gyro_Dw, gyro_C, gyro_DOS, gyro_NOA_orb, gyro_NOA_spn, eval_K, &
                                   eval_D, eval_Dw, eval_NOA, eval_spn, eval_C, eval_dos, &
                                   gyrotropic, fermi)
    !======================================================================!
    !                                                                      !
    ! Contribution from point k to the GME tensor, Eq.(9) of ZMS16,        !
    ! evaluated in the clean limit of omega.tau >> 1 where  it is real.    !
    ! The following two quantities are calculated (sigma = Pauli matrix):  !
    !                                                                      !
    ! gyro_K_spn_k = delta(E_kn-E_f).(d E_{kn}/d k_i).sigma_{kn,j}         !
    ! [units of length]                                                    !
    !                                                                      !
    ! gyro_K_orb_k = delta(E_kn-E_f).(d E_{kn}/d k_i).(2.hbar/e).m^orb_{kn,j} !
    ! [units of (length^3)*energy]                                         !
    !                                                                      !
    ! gyro_D_k = delta(E_kn-E_f).(d E_{kn}/d k_i).Omega_{kn,j}             !
    ! [units of length^3]                                                  !
    !                                                                      !
    ! gyro_Dw_k = delta(E_kn-E_f).(d E_{kn}/d k_i).tildeOmega_{kn,j}       !
    ! [units of length^3]                                                  !
    !                                                                      !
    ! gyro_C_k = delta(E_kn-E_f).(d E_{kn}/d k_i).(d E_{kn}/d k_j)         !
    ! [units of energy*length^3]                                           !
    !                                                                      !
    ! gyro_DOS_k = delta(E_kn-E_f)                                          !
    ! [units of 1/Energy]                                                  !
    !                                                                      !
    ! gme_NOA_orb_k = ?????                                           !
    !                                                                      !
    ! gme_NOA_spn_k = ??????                                           !
    !                                                                      !
    !======================================================================!

    use pw90_parameters, only: gyrotropic_type, postw90_common_type, postw90_ham_type
    use w90_comms, only: w90commtype, mpirank
    use w90_constants, only: dp, cmplx_0, cmplx_i
    use w90_get_oper, only: AA_R
    use w90_io, only: io_error, io_stopwatch, io_file_unit
    use w90_parameters, only: fermi_data_type
    use w90_param_types, only: disentangle_type, k_point_type, parameter_input_type
    use w90_postw90_common, only: pw90common_fourier_R_to_k, &
      pw90common_fourier_R_to_k_new_second_d, pw90common_get_occ, pw90common_fourier_R_to_k_vec
    use w90_spin, only: spin_get_S
    use w90_utility, only: utility_diagonalize, utility_rotate, utility_rotate_diag, utility_w0gauss
    use w90_wan_ham, only: wham_get_eig_deleig, wham_get_D_h

    implicit none

    ! Arguments
    !
    real(kind=dp), intent(in)                      :: kpt(3), kweight
    real(kind=dp), dimension(:, :, :), intent(inout)   :: gyro_K_spn, &
                                                          gyro_K_orb, &
                                                          gyro_D, &
                                                          gyro_C
    real(kind=dp), dimension(:, :, :, :), intent(inout) :: gyro_Dw, gyro_NOA_spn, gyro_NOA_orb
    real(kind=dp), dimension(:), intent(inout) :: gyro_DOS

    logical, intent(in) :: eval_K, eval_D, eval_Dw, &
                           eval_C, eval_NOA, eval_spn, eval_dos

    integer, intent(in) :: num_bands, num_kpts, num_wann, nrpts, stdout
    integer, intent(inout) :: irvec(:, :), ndegen(:), rpt_origin
    real(kind=dp), intent(in) :: eigval(:, :), real_lattice(3, 3)
    real(kind=dp), intent(inout) :: crvec(:, :)
    complex(kind=dp), intent(in) :: u_matrix(:, :, :), v_matrix(:, :, :)
    type(disentangle_type), intent(in) :: dis_data
    type(k_point_type), intent(in) :: k_points
    type(parameter_input_type), intent(in) :: param_input
    type(postw90_common_type), intent(in) :: pw90_common
    type(w90commtype), intent(in) :: comm
    type(postw90_ham_type), intent(in) :: pw90_ham
    type(gyrotropic_type), intent(in) :: gyrotropic
    type(fermi_data_type), intent(in) :: fermi
    character(len=50), intent(in) :: seedname

    complex(kind=dp), allocatable :: UU(:, :)
    complex(kind=dp), allocatable :: HH(:, :)
    complex(kind=dp), allocatable :: delHH(:, :, :)
    complex(kind=dp), allocatable :: SS(:, :, :)
    complex(kind=dp), allocatable :: AA(:, :, :)
    complex(kind=dp), allocatable :: D_h(:, :, :)

    real(kind=dp), allocatable :: curv_w_nk(:, :, :)

    integer          :: i, j, n, n1, ifermi
    real(kind=dp)    :: delta, occ(num_wann), &
                        eig(num_wann), del_eig(num_wann, 3), &
                        S(num_wann, 3), eta_smr, arg, &
                        orb_nk(3), curv_nk(3), &
                        imf_k(3, 3, 1), img_k(3, 3, 1), imh_k(3, 3, 1)
    logical          :: got_spin, got_orb_n

    allocate (UU(num_wann, num_wann))
    allocate (HH(num_wann, num_wann))
    allocate (delHH(num_wann, num_wann, 3))
    allocate (D_h(num_wann, num_wann, 3))

    if (eval_spn) allocate (SS(num_wann, num_wann, 3))

    call wham_get_eig_deleig(num_bands, num_kpts, nrpts, irvec, ndegen, rpt_origin, eigval, &
                             real_lattice, crvec, u_matrix, v_matrix, dis_data, k_points, &
                             param_input, pw90_common, comm, kpt, eig, del_eig, HH, delHH, UU, &
                             num_wann, pw90_ham, stdout, seedname)

    if (eval_Dw .or. eval_NOA) then
      allocate (AA(num_wann, num_wann, 3))
      call wham_get_D_h(delHH, UU, eig, D_h, num_wann)
      call pw90common_fourier_R_to_k_vec(stdout, seedname, kpt, AA_R, OO_true=AA)
      do i = 1, 3
        AA(:, :, i) = utility_rotate(AA(:, :, i), UU, num_wann)
      enddo
      AA = AA + cmplx_i*D_h ! Eq.(25) WYSV06
    endif

    if (eval_Dw) allocate (curv_w_nk(num_wann, gyrotropic%nfreq, 3))

    eta_smr = gyrotropic%smr_fixed_en_width

    got_spin = .false.

    do n1 = 1, gyrotropic%num_bands
      n = gyrotropic%band_list(n1)
      !
      ! ***ADJUSTABLE PARAMETER***
      ! avoid degeneracies
      !---------------------------------------------------
      if (n > 1) then
        if (eig(n) - eig(n - 1) <= gyrotropic%degen_thresh) cycle
      endif
      if (n < num_wann) then
        if (eig(n + 1) - eig(n) <= gyrotropic%degen_thresh) cycle
      endif
      !---------------------------------------------------
      got_orb_n = .false.
      do ifermi = 1, fermi%n
        arg = (eig(n) - fermi%energy_list(ifermi))/eta_smr
        !
        ! To save time: far from the Fermi surface, negligible contribution
        !
        !-------------------------
        if (abs(arg) > gyrotropic%smr_max_arg) cycle
        !-------------------------
        !
        ! Spin is computed for all bands simultaneously
        !
        if (eval_spn .and. .not. got_spin) then
          call spin_get_S(kpt, S, num_wann, stdout, seedname)
          got_spin = .true. ! Do it for only one value of ifermi and n
        endif
        ! Orbital quantities are computed for each band separately
        if (.not. got_orb_n) then
          if (eval_K) then
            ! Fake occupations: band n occupied, others empty
            occ = 0.0_dp
            occ(n) = 1.0_dp

            call berry_get_imfgh_klist(num_bands, num_kpts, nrpts, irvec, ndegen, rpt_origin, &
                                       eigval, real_lattice, u_matrix, v_matrix, dis_data, &
                                       k_points, comm, crvec, param_input, pw90_common, kpt, &
                                       num_wann, fermi, stdout, seedname, imf_k, img_k, imh_k, occ)
            do i = 1, 3
              orb_nk(i) = sum(imh_k(:, i, 1)) - sum(img_k(:, i, 1))
              curv_nk(i) = sum(imf_k(:, i, 1))
            enddo
          else if (eval_D) then
            occ = 0.0_dp
            occ(n) = 1.0_dp

            call berry_get_imf_klist(num_bands, num_kpts, nrpts, irvec, ndegen, rpt_origin, &
                                     eigval, real_lattice, u_matrix, v_matrix, dis_data, k_points, &
                                     comm, crvec, param_input, pw90_common, kpt, num_wann, fermi, &
                                     stdout, seedname, imf_k, occ)
            do i = 1, 3
              curv_nk(i) = sum(imf_k(:, i, 1))
            enddo
            got_orb_n = .true. ! Do it for only one value of ifermi
          endif

          if (eval_Dw) call gyrotropic_get_curv_w_k(eig, AA, curv_w_nk, gyrotropic)

          got_orb_n = .true. ! Do it for only one value of ifermi
        endif
        !
        delta = utility_w0gauss(arg, gyrotropic%smr_index, stdout, seedname)/eta_smr*kweight ! Broadened delta(E_nk-E_f)
        !
        ! Loop over Cartesian tensor components
        !
        do j = 1, 3
          if (eval_K .and. eval_spn) gyro_K_spn(:, j, ifermi) = &
            gyro_K_spn(:, j, ifermi) + del_eig(n, :)*S(n, j)*delta
          if (eval_K) gyro_K_orb(:, j, ifermi) = &
            gyro_K_orb(:, j, ifermi) + del_eig(n, :)*orb_nk(j)*delta
          if (eval_D) gyro_D(:, j, ifermi) = &
            gyro_D(:, j, ifermi) + del_eig(n, :)*curv_nk(j)*delta
          if (eval_Dw) then
            do i = 1, 3
              gyro_Dw(i, j, ifermi, :) = &
                gyro_Dw(i, j, ifermi, :) + del_eig(n, i)*delta*curv_w_nk(n, :, j)
            enddo
          endif
          if (eval_C) gyro_C(:, j, ifermi) = &
            gyro_C(:, j, ifermi) + del_eig(n, :)*del_eig(n, j)*delta
        enddo !j
        if (eval_dos) gyro_DOS(ifermi) = gyro_DOS(ifermi) + delta

      enddo !ifermi
    enddo !n

    if (eval_NOA) then
      if (eval_spn) then
        call gyrotropic_get_NOA_k(kpt, stdout, seedname, gyrotropic, num_wann, fermi, param_input, &
                                  kweight, eig, del_eig, AA, UU, gyro_NOA_orb, gyro_NOA_spn)
      else
        call gyrotropic_get_NOA_k(kpt, stdout, seedname, gyrotropic, num_wann, fermi, param_input, &
                                  kweight, eig, del_eig, AA, UU, gyro_NOA_orb)
      endif
    endif

  end subroutine gyrotropic_get_k_list

  subroutine gyrotropic_get_curv_w_k(eig, AA, curv_w_k, gyrotropic)
    !======================================================================!
    !                                                                      !
    ! calculation of the band-resolved                                     !
    ! frequency-dependent   berry curvature                                !
    !                                                                      !
    ! tildeOmega(w)=                                                        !
    !     -eps_{bcd}sum_m ( w_mn^2/(wmn^2-w^2)) *Im[A_{nm,c}A_{mn,d}       !
    !                                        !
    !======================================================================!

    use pw90_parameters, only: gyrotropic_type
    use w90_constants, only: dp

    implicit none

    ! Arguments
    !
    type(gyrotropic_type), intent(in) :: gyrotropic

    real(kind=dp), intent(in)                    :: eig(:)
    complex(kind=dp), intent(in)                 :: AA(:, :, :)
    real(kind=dp), dimension(:, :, :), intent(out) :: curv_w_k  ! (num_wann,n_freq,3)

    real(kind=dp), allocatable    :: multWre(:)
    integer          :: i, n, m, n1, m1
    real(kind=dp)    :: wmn

    allocate (multWre(gyrotropic%nfreq))
    curv_w_k(:, :, :) = 0_dp

    do n1 = 1, gyrotropic%num_bands
      n = gyrotropic%band_list(n1)
      do m1 = 1, gyrotropic%num_bands
        m = gyrotropic%band_list(m1)
        if (n == m) cycle
        wmn = eig(m) - eig(n)
        multWre(:) = real(wmn**2/(wmn**2 - gyrotropic%freq_list(:)**2))
        do i = 1, 3
          curv_w_k(n, :, i) = curv_w_k(n, :, i) - &
                              2_dp*aimag(AA(n, m, alpha_A(i))*AA(m, n, beta_A(i)))*multWre
        enddo
      enddo !m
    enddo !n

  end subroutine gyrotropic_get_curv_w_k

  subroutine gyrotropic_get_NOA_k(kpt, stdout, seedname, gyrotropic, num_wann, fermi, param_input, &
                                  kweight, eig, del_eig, AA, UU, gyro_NOA_orb, gyro_NOA_spn)
    !====================================================================!
    !                                                                    !
    ! Contribution from point k to the real (antisymmetric) part         !
    ! of the natural  complex interband optical conductivity             !
    !                                                                    !
    ! Re gyro_NOA_orb  =  SUM_{n,l}^{oe}  hbar^{-2}/(w_nl^2-w^2) *
    !   Re (  A_lnb Bnlac -Alna Bnlbc)
    !     -SUM_{n,l}^{oe}  hbar^{-2}(3*w_ln^2-w^2)/(w_nl^2-w^2)^2 *
    ! Im (  A_lnb Bnlac -Alna Bnlac)nm_a A_mn_b )
    ! [units of Ang^3/eV]                                                !

    ! [units of Ang^3]                                                !
    ! Re gyro_NOA_spn_{ab,c}  =  SUM_{n,l}^{oe}  hbar^{-2}/(w_nl^2-w^2) *
    !   Re (  A_lnb Bnlac -Alna Bnlbc)
    ! [units of Ang/eV^2]                                                !
    !
    !   here a,b  defined as epsilon_{abd}=1  (and NOA_dc tensor is saved)  !
    !====================================================================!

    use pw90_parameters, only: gyrotropic_type
    use w90_constants, only: dp, cmplx_1
    use w90_get_oper, only: SS_R
    use w90_io, only: io_time, io_error
    use w90_parameters, only: fermi_data_type, parameter_input_type
    use w90_postw90_common, only: pw90common_fourier_R_to_k_new
    use w90_spin, only: spin_get_S
    use w90_utility, only: utility_rotate

    implicit none

    ! Arguments
    !
    type(gyrotropic_type), intent(in) :: gyrotropic
    type(fermi_data_type), intent(in) :: fermi
    type(parameter_input_type), intent(in) :: param_input

    integer, intent(in) :: stdout
    integer, intent(in) :: num_wann
    real(kind=dp), intent(in)       :: kpt(3), kweight
    real(kind=dp), dimension(:, :, :, :), optional, intent(inout)    :: gyro_NOA_spn
    real(kind=dp), dimension(:, :, :, :), intent(inout)    :: gyro_NOA_orb
    complex(kind=dp), dimension(:, :, :), intent(in)       :: AA
    complex(kind=dp), dimension(:, :), intent(in)       :: UU
    real(kind=dp), dimension(:), intent(in)       :: eig
    real(kind=dp), dimension(:, :), intent(in)       :: del_eig
    character(len=50), intent(in)  :: seedname

    complex(kind=dp), allocatable :: Bnl_orb(:, :, :, :)
    complex(kind=dp), allocatable :: Bnl_spin(:, :, :, :)
    complex(kind=dp), allocatable :: SS(:, :, :)
    complex(kind=dp), allocatable :: S_h(:, :, :)

    complex(kind=dp)              :: multW1(gyrotropic%nfreq)
    real(kind=dp) :: multWe(gyrotropic%nfreq), multWm(gyrotropic%nfreq)

    integer ::  num_occ, num_unocc, occ_list(num_wann), unocc_list(num_wann)

    !real(kind=dp)    ::  wmn

    integer          :: j, n, l, n1, l1, a, b, c, ab, ifermi
    real(kind=dp)    :: wln

    if (present(gyro_NOA_spn)) then
      allocate (SS(num_wann, num_wann, 3))
      allocate (S_h(num_wann, num_wann, 3))
      do j = 1, 3 ! spin direction
        call pw90common_fourier_R_to_k_new(stdout, seedname, kpt, SS_R(:, :, :, j), OO=SS(:, :, j))
        S_h(:, :, j) = utility_rotate(SS(:, :, j), UU, num_wann)
      enddo
    endif

    do ifermi = 1, fermi%n

      num_occ = 0
      num_unocc = 0
      do n1 = 1, gyrotropic%num_bands
        n = gyrotropic%band_list(n1)
        if (eig(n) < fermi%energy_list(ifermi)) then
          num_occ = num_occ + 1
          occ_list(num_occ) = n
        elseif (eig(n) < gyrotropic%eigval_max) then
          num_unocc = num_unocc + 1
          unocc_list(num_unocc) = n
        endif
      enddo

      if (num_occ == 0) then
        if (param_input%iprint .ge. 2) &
          write (stdout, *) "WARNING no occupied bands included in the calculation for kpt=", &
          kpt, ", EF[", ifermi, "]=", fermi%energy_list(ifermi), "eV"
        cycle
      endif

      if (num_unocc == 0) then
        if (param_input%iprint .ge. 2) &
          write (stdout, *) "WARNING no unoccupied bands included in the calculation for kpt=", &
          kpt, ", EF[", ifermi, "]=", fermi%energy_list(ifermi), "eV"
        cycle
      endif

      allocate (Bnl_orb(num_occ, num_unocc, 3, 3))
      call gyrotropic_get_NOA_Bnl_orb(eig, del_eig, AA, num_occ, occ_list, num_unocc, unocc_list, &
                                      Bnl_orb, gyrotropic)

      if (present(gyro_NOA_spn)) then
        allocate (Bnl_spin(num_occ, num_unocc, 3, 3))
        call gyrotropic_get_NOA_Bnl_spin(S_h, num_occ, occ_list, num_unocc, unocc_list, Bnl_spin)
      endif

      do n1 = 1, num_occ
        n = occ_list(n1)
        do l1 = 1, num_unocc
          l = unocc_list(l1)

          wln = eig(l) - eig(n)
          multW1(:) = cmplx_1/(wln*wln - gyrotropic%freq_list(:)**2)
          multWm(:) = real(multW1)*kweight
          multWe(:) = real(-multW1(:)*(2*wln**2*multW1(:) + cmplx_1))*kweight
          do ab = 1, 3
            a = alpha_A(ab)
            b = beta_A(ab)
            do c = 1, 3
              gyro_NOA_orb(ab, c, ifermi, :) = &
                gyro_NOA_orb(ab, c, ifermi, :) + &
                multWm(:)*real(AA(l, n, b)*Bnl_orb(n1, l1, a, c) - &
                               AA(l, n, a)*Bnl_orb(n1, l1, b, c)) + &
                multWe(:)*(del_eig(n, c) + del_eig(l, c))*aimag(AA(n, l, a)*AA(l, n, b))

              if (present(gyro_NOA_spn)) &
                gyro_NOA_spn(ab, c, ifermi, :) = &
                gyro_NOA_spn(ab, c, ifermi, :) + &
                multWm(:)*real(AA(l, n, b)*Bnl_spin(n1, l1, a, c) - &
                               AA(l, n, a)*Bnl_spin(n1, l1, b, c))

            enddo ! c
          enddo ! ab
        enddo  ! l1
      enddo ! n1
      deallocate (Bnl_orb)
      if (present(gyro_NOA_spn)) deallocate (Bnl_spin)
    enddo !ifermi

  end subroutine gyrotropic_get_NOA_k

  subroutine gyrotropic_get_NOA_Bnl_orb(eig, del_eig, AA, num_occ, occ_list, num_unocc, &
                                        unocc_list, Bnl, gyrotropic)
    !====================================================================!
    !                                                                    !
    ! Calculating the matrix                                             !
    ! B_{nl,ac}(num_occ,num_unocc,3,3)=                    !
    !      -sum_m(  (E_m-E_n)A_nma*Amlc +(E_l-E_m)A_nmc*A_mla -          !
    !      -i( del_a (E_n+E_l) A_nlc                                     !
    !   in units eV*Ang^2                                                !
    !====================================================================!

    use pw90_parameters, only: gyrotropic_type
    use w90_constants, only: dp, cmplx_i, cmplx_0

    implicit none

    ! Arguments
    !
    type(gyrotropic_type), intent(in) :: gyrotropic
    integer, intent(in) ::num_occ, num_unocc
    integer, dimension(:), intent(in) ::occ_list, unocc_list
    real(kind=dp), dimension(:), intent(in) ::eig     !  n
    real(kind=dp), dimension(:, :), intent(in) ::del_eig !  n
    complex(kind=dp), dimension(:, :, :), intent(in) ::AA      !  n,l,a
    complex(kind=dp), dimension(:, :, :, :), intent(out)::Bnl     !   n,l,a,c
    integer n, m, l, a, c, n1, m1, l1

    Bnl(:, :, :, :) = cmplx_0

    do a = 1, 3
      do c = 1, 3
        do n1 = 1, num_occ
          n = occ_list(n1)
          do l1 = 1, num_unocc
            l = unocc_list(l1)
            Bnl(n1, l1, a, c) = -cmplx_i*(del_eig(n, a) + del_eig(l, a))*AA(n, l, c)
            do m1 = 1, gyrotropic%num_bands
              m = gyrotropic%band_list(m1)
              Bnl(n1, l1, a, c) = Bnl(n1, l1, a, c) + &
                                  (eig(n) - eig(m))*AA(n, m, a)*AA(m, l, c) - &
                                  (eig(l) - eig(m))*AA(n, m, c)*AA(m, l, a)
            enddo ! m1
          enddo !l1
        enddo !n1
      enddo !c
    enddo !a

  end subroutine gyrotropic_get_NOA_Bnl_orb

  subroutine gyrotropic_get_NOA_Bnl_spin(S_h, &
                                         num_occ, occ_list, num_unocc, unocc_list, Bnl)
    !====================================================================!
    !                                                                    !
    ! Calculating the matrix                                             !
    ! B_{nl,ac}^spin(num_occ,num_unocc,3,3)=                             !
    !      -i  eps_{abc}  < u_n | sigma_b | u_l >                        !
    !   ( dimensionless )                                                !
    !====================================================================!

    use w90_constants, only: dp, cmplx_i, cmplx_0

    implicit none

    ! Arguments
    !
    integer, intent(in) ::num_occ, num_unocc
    integer, dimension(:), intent(in) ::occ_list, unocc_list
    complex(kind=dp), dimension(:, :, :), intent(in) ::S_h     !  n,l,a
    complex(kind=dp), dimension(:, :, :, :), intent(out)::Bnl     !   n,l,a,c
    integer n, l, a, b, c, n1, l1

    Bnl(:, :, :, :) = cmplx_0

    do b = 1, 3
      c = alpha_A(b)
      a = beta_A(b)
      do n1 = 1, num_occ
        n = occ_list(n1)
        do l1 = 1, num_unocc
          l = unocc_list(l1)
          Bnl(n1, l1, a, c) = S_h(n, l, b)
        enddo !l1
      enddo !n1
    enddo !b

    Bnl = Bnl*(-cmplx_i)

  end subroutine gyrotropic_get_NOA_Bnl_spin

  subroutine gyrotropic_outprint_tensor(stdout, seedname, gyrotropic, fermi, f_out_name, arrEf, arrEF1D, arrEfW, &
                                        units, comment, symmetrize)
    use pw90_parameters, only: gyrotropic_type
    use w90_parameters, only: fermi_data_type
    use w90_io, only: io_file_unit

    implicit none

    type(gyrotropic_type), intent(in) :: gyrotropic
    type(fermi_data_type), intent(in) :: fermi

    integer, intent(in) :: stdout
    character(len=30), intent(in) :: f_out_name
    character(len=30), intent(in), optional :: units
    character(len=120), intent(in), optional :: comment
    real(kind=dp), dimension(:, :, :), intent(in), optional :: arrEf
    real(kind=dp), dimension(:, :, :, :), intent(in), optional :: arrEfW
    real(kind=dp), dimension(:), intent(in), optional :: arrEf1D
    logical, optional, intent(in) :: symmetrize
    character(len=50), intent(in)  :: seedname

    character(len=120)  :: file_name
    integer             :: i, file_unit
    logical             :: lsym

    lsym = .true.
    if (present(symmetrize)) then
      if (.not. symmetrize) lsym = .false.
    endif

    file_name = trim(seedname)//"-gyrotropic-"//trim(f_out_name)//".dat"
    file_name = trim(file_name)
    file_unit = io_file_unit()
    write (stdout, '(/,3x,a)') '* '//file_name
    open (file_unit, FILE=file_name, STATUS='UNKNOWN', FORM='FORMATTED')

    if (present(comment)) write (file_unit, *) "#"//trim(comment)
    if (present(units)) write (file_unit, *) "# in units of [ "//trim(units)//" ] "

    if (present(arrEf)) then
      call gyrotropic_outprint_tensor_w(fermi, file_unit, 0.0_dp, arr33N=arrEf, symmetrize=lsym)
    elseif (present(arrEfW)) then
      do i = 1, gyrotropic%nfreq
      call gyrotropic_outprint_tensor_w(fermi, file_unit, real(gyrotropic%freq_list(i)), arr33N=arrEfW(:, :, :, i), symmetrize=lsym)
      enddo
    elseif (present(arrEf1D)) then
      call gyrotropic_outprint_tensor_w(fermi, file_unit, 0.0_dp, arrN=arrEf1D)
    endif

    close (file_unit)

  end subroutine gyrotropic_outprint_tensor

  subroutine gyrotropic_outprint_tensor_w(fermi, file_unit, omega, arr33N, arrN, symmetrize)
    use w90_parameters, only: fermi_data_type

    implicit none

    type(fermi_data_type), intent(in) :: fermi

    integer, intent(in) :: file_unit
    real(kind=dp), intent(in) :: omega
    real(kind=dp), dimension(:, :, :), optional, intent(in) :: arr33N
    real(kind=dp), dimension(:), optional, intent(in) :: arrN
    ! symmetrize= True  - get symmetric and assimetric parts
    ! symmetrize= False - write the asymmetric (in xy) tensor gamma_{xyz}
    ! symmetrize not present - write as is
    logical, optional, intent(in) :: symmetrize

    real(kind=dp) ::  xx(fermi%n), yy(fermi%n), zz(fermi%n), &
                     xy(fermi%n), xz(fermi%n), yz(fermi%n), &
                     x(fermi%n), y(fermi%n), z(fermi%n)
    integer       ::  i
    logical lsym

    if (present(arr33N)) then
      lsym = .false.
      if (present(symmetrize)) lsym = symmetrize

      if (lsym) then
        ! Symmetric part
        xx = arr33N(1, 1, :)
        yy = arr33N(2, 2, :)
        zz = arr33N(3, 3, :)
        xy = (arr33N(1, 2, :) + arr33N(2, 1, :))/2.0_dp
        xz = (arr33N(1, 3, :) + arr33N(3, 1, :))/2.0_dp
        yz = (arr33N(2, 3, :) + arr33N(3, 2, :))/2.0_dp
        ! Antisymmetric part, in polar-vector form
        x = (arr33N(2, 3, :) - arr33N(3, 2, :))/2.0_dp
        y = (arr33N(3, 1, :) - arr33N(1, 3, :))/2.0_dp
        z = (arr33N(1, 2, :) - arr33N(2, 1, :))/2.0_dp
      else
        xx = arr33N(1, 1, :)
        yy = arr33N(2, 2, :)
        zz = arr33N(3, 3, :)
        xy = arr33N(1, 2, :)
        xz = arr33N(1, 3, :)
        yz = arr33N(2, 3, :)
        x = arr33N(3, 2, :)
        y = arr33N(3, 1, :)
        z = arr33N(2, 1, :)
      endif

      if (present(symmetrize)) then
        if (symmetrize) then
          write (file_unit, '(a1,29x,a1,38x,a14,37x,a2,14x,a15,14x,a1)') '#', "|", "symmetric part", "||", "asymmetric part", "|"
          write (file_unit, '(11a15)') '# EFERMI(eV)', "omega(eV)", 'xx', 'yy', 'zz', 'xy', 'xz', 'yz', 'x', 'y', 'z'
        else
          write (file_unit, '(11a15)') '# EFERMI(eV)', "omega(eV)", 'yzx', 'zxy', 'xyz', 'yzy', 'yzz', 'zxz', 'xyy', 'yzz', 'zxx'
        endif
      else
        write (file_unit, '(11a15)') '# EFERMI(eV)', "omega(eV)", 'xx', 'yy', 'zz', 'xy', 'xz', 'yz', 'zy', 'xz', 'yx'
      endif

      do i = 1, fermi%n
        write (file_unit, '(11E15.6)') fermi%energy_list(i), omega, xx(i), yy(i), zz(i), xy(i), xz(i), yz(i), x(i), y(i), z(i)
      enddo
    endif

    if (present(arrN)) then
      write (file_unit, '(2a15)') '# EFERMI(eV) '
      do i = 1, fermi%n
        write (file_unit, '(11E15.6)') fermi%energy_list(i), arrN(i)
      enddo

    endif
    write (file_unit, *)
    write (file_unit, *)
  end subroutine gyrotropic_outprint_tensor_w

end module w90_gyrotropic
