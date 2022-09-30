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
!  w90_wan_ham: Hamiltonian operations in Wannier basis      !
!                                                            !
!------------------------------------------------------------!

module w90_wan_ham

  !! This module contain operations on the Hamiltonian in the WF basis

  use w90_constants, only: dp
  use w90_error, only: w90_error_type, set_error_alloc, set_error_dealloc, set_error_fatal, &
    set_error_input, set_error_fatal, set_error_file

  implicit none

  private

  public :: wham_get_D_h
  public :: wham_get_D_h_P_value
  public :: wham_get_eig_deleig
  public :: wham_get_eig_deleig_TB_conv
  public :: wham_get_eig_UU_HH_AA_sc
  public :: wham_get_eig_UU_HH_AA_sc_TB_conv
  public :: wham_get_eig_UU_HH_JJlist
  public :: wham_get_occ_mat_list

contains

  !================================================!

  subroutine wham_get_D_h_a(delHH_a, UU, eig, ef, D_h_a, num_wann)
    !================================================!
    !
    !! Compute D^H_a=UU^dag.del_a UU (a=alpha,beta),
    !! using Eq.(24) of WYSV06
    !
    !================================================!

    use w90_constants, only: dp, cmplx_0
    use w90_utility, only: utility_rotate
    use w90_postw90_common, only: pw90common_get_occ

    ! arguments
    complex(kind=dp), intent(in)  :: delHH_a(:, :)
    complex(kind=dp), intent(in)  :: UU(:, :)
    complex(kind=dp), intent(out) :: D_h_a(:, :)

    real(kind=dp), intent(in) :: eig(:)
    real(kind=dp), intent(in) :: ef

    integer, intent(in) :: num_wann

    ! local variables
    complex(kind=dp), allocatable :: delHH_a_bar(:, :)
    real(kind=dp)                 :: occ(num_wann)
    integer                       :: n, m

    call pw90common_get_occ(ef, eig, occ, num_wann)

    allocate (delHH_a_bar(num_wann, num_wann))
    delHH_a_bar = utility_rotate(delHH_a, UU, num_wann)
    do m = 1, num_wann
      do n = 1, num_wann
        if (occ(n) > 0.999_dp .and. occ(m) < 0.001_dp) then
          D_h_a(n, m) = delHH_a_bar(n, m)/(eig(m) - eig(n))
        else
          D_h_a(n, m) = cmplx_0
        end if
      end do
    end do
    D_h_a = D_h_a - conjg(transpose(D_h_a))

  end subroutine wham_get_D_h_a

  !================================================!
  subroutine wham_get_D_h(delHH, D_h, UU, eig, num_wann)
    !================================================!
    !
    !! Compute D^H_a=UU^dag.del_a UU (a=x,y,z)
    !! using Eq.(24) of WYSV06
    !
    !================================================!

    ! TO DO: Implement version where energy denominators only connect
    !        occupied and empty states. In this case probably do not need
    !        to worry about avoiding small energy denominators

    use w90_constants, only: dp, cmplx_0
    use w90_utility, only: utility_rotate

    ! arguments
    complex(kind=dp), intent(in) :: delHH(:, :, :)
    complex(kind=dp), intent(in) :: UU(:, :)
    complex(kind=dp), intent(out) :: D_h(:, :, :)

    real(kind=dp), intent(in) :: eig(:)

    integer, intent(in) :: num_wann

    ! local variables
    complex(kind=dp), allocatable :: delHH_bar_i(:, :)
    integer                       :: n, m, i

    allocate (delHH_bar_i(num_wann, num_wann))
    D_h = cmplx_0
    do i = 1, 3
      delHH_bar_i(:, :) = utility_rotate(delHH(:, :, i), UU, num_wann)
      do m = 1, num_wann
        do n = 1, num_wann
          if (n == m .or. abs(eig(m) - eig(n)) < 1.0e-7_dp) cycle
          D_h(n, m, i) = delHH_bar_i(n, m)/(eig(m) - eig(n))
        end do
      end do
    enddo

  end subroutine wham_get_D_h

  !================================================!
  subroutine wham_get_D_h_P_value(pw90_berry, delHH, D_h, UU, eig, num_wann)
    !================================================!
    !
    !! Compute D^H_a=UU^dag.del_a UU (a=x,y,z)
    !! using Eq.(24) of WYSV06
    !  and prescription for energy denominator
    !  from BK81
    !
    !================================================!

    ! TO DO: Implement version where energy denominators only connect
    !        occupied and empty states. In this case probably do not need
    !        to worry about avoiding small energy denominators

    use w90_constants, only: dp, cmplx_0
    use w90_postw90_types, only: pw90_berry_mod_type !sc_eta
    use w90_utility, only: utility_rotate

    ! arguments
    type(pw90_berry_mod_type), intent(in) :: pw90_berry

    integer, intent(in) :: num_wann

    real(kind=dp), intent(in) :: eig(:)

    complex(kind=dp), intent(in)  :: delHH(:, :, :)
    complex(kind=dp), intent(in)    :: UU(:, :)
    complex(kind=dp), intent(out) :: D_h(:, :, :)

    ! local variables
    complex(kind=dp), allocatable :: delHH_bar_i(:, :)
    integer                       :: n, m, i
    real(kind=dp)                 :: deltaE

    allocate (delHH_bar_i(num_wann, num_wann))
    D_h = cmplx_0
    deltaE = 0.d0
    do i = 1, 3
      delHH_bar_i(:, :) = utility_rotate(delHH(:, :, i), UU, num_wann)
      do m = 1, num_wann
        do n = 1, num_wann
          if (n == m) cycle
          deltaE = eig(m) - eig(n)
          D_h(n, m, i) = delHH_bar_i(n, m)*(deltaE/(deltaE**(2) + pw90_berry%sc_eta**(2)))
        end do
      end do
    enddo

  end subroutine wham_get_D_h_P_value

  !================================================!
  subroutine wham_get_JJp_JJm_list(delHH, UU, eig, JJp_list, JJm_list, num_wann, &
                                   fermi_energy_list, occ)
    !================================================!
    !                                                !
    ! Compute JJ^+_a and JJ^-_a (a=Cartesian index)  !
    ! for a list of Fermi energies                   !
    !                                                !
    ! This routine is a replacement for              !
    ! wham_get_JJp_list and wham_getJJm_list.        !
    ! It computes both lists at once in a more       !
    ! efficient manner.                              !
    !                                                !
    !  Tsirkin:   added the optional occ parameter   !
    !                                                !
    !================================================!

    use w90_constants, only: dp, cmplx_0, cmplx_i
    use w90_utility, only: utility_rotate_new

    ! arguments
    real(kind=dp), allocatable, intent(in) :: fermi_energy_list(:)
    integer, intent(in) :: num_wann
    real(kind=dp), intent(in) :: eig(:)
    real(kind=dp), intent(in), optional, dimension(:) :: occ
    complex(kind=dp), intent(inout) :: delHH(:, :)
    complex(kind=dp), intent(in) :: UU(:, :)
    complex(kind=dp), intent(out) :: JJm_list(:, :, :)
    complex(kind=dp), intent(out) :: JJp_list(:, :, :)

    ! local variables
    integer :: n, m, ife, nfermi_loc
    real(kind=dp) :: fe

    if (present(occ)) then
      nfermi_loc = 1
    else
      nfermi_loc = 0
      if (allocated(fermi_energy_list)) nfermi_loc = size(fermi_energy_list)
    endif

    call utility_rotate_new(delHH, UU, num_wann)
    do ife = 1, nfermi_loc
      fe = fermi_energy_list(ife)
      do m = 1, num_wann
        do n = 1, num_wann
          if (present(occ)) then
            if (occ(m) < 0.5_dp .and. occ(n) > 0.5_dp) then
              JJm_list(n, m, ife) = cmplx_i*delHH(n, m)/(eig(m) - eig(n))
              JJp_list(m, n, ife) = cmplx_i*delHH(m, n)/(eig(n) - eig(m))
            else
              JJm_list(n, m, ife) = cmplx_0
              JJp_list(m, n, ife) = cmplx_0
            end if
          else
            if (eig(n) > fe .and. eig(m) < fe) then
              JJp_list(n, m, ife) = cmplx_i*delHH(n, m)/(eig(m) - eig(n))
              JJm_list(m, n, ife) = cmplx_i*delHH(m, n)/(eig(n) - eig(m))
            else
              JJp_list(n, m, ife) = cmplx_0
              JJm_list(m, n, ife) = cmplx_0
            endif
          endif
        enddo
      enddo
      call utility_rotate_new(JJp_list(:, :, ife), UU, num_wann, reverse=.true.)
      call utility_rotate_new(JJm_list(:, :, ife), UU, num_wann, reverse=.true.)
    end do

  end subroutine wham_get_JJp_JJm_list

  !================================================!
  subroutine wham_get_occ_mat_list(fermi_energy_list, f_list, g_list, UU, num_wann, error, comm, &
                                   eig, occ)
    !================================================!
    !
    !! Occupation matrix f, and g=1-f
    !! for a list of Fermi energies
    ! Tsirkin: !now optionally either eig or occ parameters may be supplied
    !    (Changed consistently the calls from the Berry module)
    !================================================!

    use w90_constants, only: dp, cmplx_0, cmplx_1
    use w90_postw90_common, only: pw90common_get_occ
    use w90_comms, only: w90_comm_type

    ! arguments
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    real(kind=dp), allocatable, intent(in) :: fermi_energy_list(:)

    integer, intent(in) :: num_wann

    real(kind=dp), intent(in), optional :: eig(:)
    real(kind=dp), intent(in), optional :: occ(:)

    complex(kind=dp), intent(in)  :: UU(:, :)
    complex(kind=dp), intent(out) :: f_list(:, :, :)
    complex(kind=dp), intent(out) :: g_list(:, :, :)

    ! local variables
    integer       :: n, m, i, if, nfermi_loc
    real(kind=dp), allocatable :: occ_list(:, :)

    if (present(occ)) then
      nfermi_loc = 1
    else
      nfermi_loc = 0
      if (allocated(fermi_energy_list)) nfermi_loc = size(fermi_energy_list)
    endif
    allocate (occ_list(num_wann, nfermi_loc))

    if (present(occ) .and. present(eig)) then
      call set_error_input(error, 'occ_list and eig cannot be both arguments in get_occ_mat_list', &
                           comm)
      return
    elseif (.not. present(occ) .and. .not. present(eig)) then
      call set_error_input(error, 'either occ_list or eig must be passed as arguments to get_occ_mat_list', comm)
      return
    endif

    if (present(occ)) then
      occ_list(:, 1) = occ(:)
    else
      do if = 1, nfermi_loc
        call pw90common_get_occ(fermi_energy_list(if), eig, occ_list(:, if), num_wann)
      enddo
    endif

    f_list = cmplx_0
    do if = 1, nfermi_loc
      do n = 1, num_wann
        do m = 1, num_wann
          do i = 1, num_wann
            f_list(n, m, if) = f_list(n, m, if) &
                               + UU(n, i)*occ_list(i, if)*conjg(UU(m, i))
          enddo
          g_list(n, m, if) = -f_list(n, m, if)
          if (m == n) g_list(n, n, if) = g_list(n, n, if) + cmplx_1
        enddo
      enddo
    enddo

  end subroutine wham_get_occ_mat_list

  !================================================!
  subroutine wham_get_deleig_a(deleig_a, eig, delHH_a, UU, num_wann, pw90_band_deriv_degen, &
                               error, comm)
    !================================================!
    !
    !! Band derivatives dE/dk_a
    !
    !================================================!

    use w90_constants, only: dp !, cmplx_0, cmplx_i
    use w90_utility, only: utility_diagonalize, utility_rotate, utility_rotate_diag
    use w90_postw90_types, only: pw90_band_deriv_degen_type
    use w90_comms, only: w90_comm_type

    ! arguments
    type(pw90_band_deriv_degen_type), intent(in) :: pw90_band_deriv_degen
    integer, intent(in) :: num_wann
    real(kind=dp), intent(in) :: eig(num_wann)
    real(kind=dp), intent(out) :: deleig_a(num_wann)
    complex(kind=dp), intent(in) :: delHH_a(:, :)
    complex(kind=dp), intent(in) :: UU(:, :)
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    ! local variables
    integer                       :: i, degen_min, degen_max, dim
    real(kind=dp)                 :: diff
    complex(kind=dp), allocatable :: delHH_bar_a(:, :), U_deg(:, :)

    allocate (delHH_bar_a(num_wann, num_wann))
    allocate (U_deg(num_wann, num_wann))

    if (pw90_band_deriv_degen%use_degen_pert) then

      delHH_bar_a = utility_rotate(delHH_a, UU, num_wann)

      ! Assuming that the energy eigenvalues are stored in eig(:) in
      ! increasing order (diff >= 0)

      i = 0
      do
        i = i + 1
        if (i > num_wann) exit
        if (i + 1 <= num_wann) then
          diff = eig(i + 1) - eig(i)
        else

          ! i-th is the highest band, and it is non-degenerate

          diff = pw90_band_deriv_degen%degen_thr + 1.0_dp
        end if
        if (diff < pw90_band_deriv_degen%degen_thr) then

          ! Bands i and i+1 are degenerate

          degen_min = i
          degen_max = degen_min + 1

          ! See if any higher bands are in the same degenerate group

          do
            if (degen_max + 1 > num_wann) exit
            diff = eig(degen_max + 1) - eig(degen_max)
            if (diff < pw90_band_deriv_degen%degen_thr) then
              degen_max = degen_max + 1
            else
              exit
            end if
          end do

          ! Bands from degen_min to degen_max are degenerate. Diagonalize
          ! the submatrix in Eq.(31) YWVS07 over this degenerate subspace.
          ! The eigenvalues are the band gradients

          dim = degen_max - degen_min + 1
          call utility_diagonalize(delHH_bar_a(degen_min:degen_max, degen_min:degen_max), dim, &
                                   deleig_a(degen_min:degen_max), U_deg(1:dim, 1:dim), error, comm)
          if (allocated(error)) return

          ! Scanned bands up to degen_max

          i = degen_max
        else

          ! Use non-degenerate form [Eq.(27) YWVS07] for current (i-th) band

          deleig_a(i) = real(delHH_bar_a(i, i), dp)
        end if
      end do

    else

      ! Use non-degenerate form for all bands

      deleig_a(:) = real(utility_rotate_diag(delHH_a(:, :), UU, num_wann), dp)

    end if

  end subroutine wham_get_deleig_a

  !================================================!
  subroutine wham_get_eig_deleig(dis_manifold, kpt_latt, pw90_band_deriv_degen, ws_region, &
                                 print_output, wannier_data, ws_distance, wigner_seitz, delHH, HH, &
                                 HH_R, u_matrix, UU, v_matrix, del_eig, eig, eigval, kpt, &
                                 real_lattice, scissors_shift, mp_grid, num_bands, num_kpts, &
                                 num_wann, num_valence_bands, effective_model, have_disentangled, &
                                 seedname, stdout, timer, error, comm)
    !================================================!
    !
    !! Given a k point, this function returns eigenvalues E and
    !! derivatives of the eigenvalues dE/dk_a, using wham_get_deleig_a
    !
    !================================================!

    use w90_constants, only: dp
    use w90_postw90_types, only: pw90_band_deriv_degen_type, wigner_seitz_type
    use w90_comms, only: w90_comm_type, mpirank
    use w90_constants, only: dp, cmplx_0
    use w90_get_oper, only: get_HH_R
    use w90_types, only: dis_manifold_type, print_output_type, wannier_data_type, &
      ws_region_type, ws_distance_type, timer_list_type
    use w90_postw90_common, only: pw90common_fourier_R_to_k_new_second_d, &
      pw90common_fourier_R_to_k
    use w90_utility, only: utility_diagonalize

    implicit none

    ! arguments
    type(dis_manifold_type), intent(in) :: dis_manifold
    real(kind=dp), intent(in) :: kpt_latt(:, :)
    type(pw90_band_deriv_degen_type), intent(in) :: pw90_band_deriv_degen
    type(print_output_type), intent(in) :: print_output
    type(ws_region_type), intent(in) :: ws_region
    type(w90_comm_type), intent(in) :: comm
    type(wannier_data_type), intent(in) :: wannier_data
    type(wigner_seitz_type), intent(inout) :: wigner_seitz
    type(ws_distance_type), intent(inout) :: ws_distance
    type(timer_list_type), intent(inout) :: timer
    type(w90_error_type), allocatable, intent(out) :: error

    integer, intent(in) :: num_wann, num_kpts, num_bands, num_valence_bands
    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: stdout

    real(kind=dp), intent(in) :: kpt(3)!! the three coordinates of the k point vector (in relative coordinates)
    real(kind=dp), intent(out) :: eig(num_wann)!! the calculated eigenvalues at kpt
    real(kind=dp), intent(out) :: del_eig(num_wann, 3)
    !! the calculated derivatives of the eigenvalues at kpt [first component: band; second component: 1,2,3
    !! for the derivatives along the three k directions]
    real(kind=dp), intent(in) :: eigval(:, :)
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    real(kind=dp), intent(in) :: scissors_shift

    complex(kind=dp), intent(out)   :: HH(:, :)
    !! the Hamiltonian matrix at kpt
    complex(kind=dp), intent(out) :: delHH(:, :, :)
    !! the delHH matrix (derivative of H) at kpt
    complex(kind=dp), intent(out)   :: UU(:, :)
    !! the rotation matrix that gives the eigenvectors of HH
    complex(kind=dp), intent(in) :: u_matrix(:, :, :), v_matrix(:, :, :)
    complex(kind=dp), allocatable, intent(inout) :: HH_R(:, :, :) !  <0n|r|Rm>

    character(len=50), intent(in)  :: seedname
    logical, intent(in) :: have_disentangled
    logical, intent(in) :: effective_model

    ! I call it to be sure that it has been called already once,
    ! and that HH_R contains the actual matrix.
    ! Further calls should return very fast.
    call get_HH_R(dis_manifold, kpt_latt, print_output, wigner_seitz, HH_R, u_matrix, v_matrix, &
                  eigval, real_lattice, scissors_shift, num_bands, num_kpts, num_wann, &
                  num_valence_bands, effective_model, have_disentangled, seedname, ws_distance, ws_region, &
                  stdout, timer, error, comm)
    if (allocated(error)) return

    call pw90common_fourier_R_to_k(ws_region, wannier_data, ws_distance, wigner_seitz, HH, HH_R, &
                                   kpt, real_lattice, mp_grid, 0, num_wann, error, comm)
    if (allocated(error)) return
    call utility_diagonalize(HH, num_wann, eig, UU, error, comm)
    if (allocated(error)) return
    call pw90common_fourier_R_to_k(ws_region, wannier_data, ws_distance, wigner_seitz, &
                                   delHH(:, :, 1), HH_R, kpt, real_lattice, mp_grid, 1, num_wann, &
                                   error, comm)
    if (allocated(error)) return
    call pw90common_fourier_R_to_k(ws_region, wannier_data, ws_distance, wigner_seitz, &
                                   delHH(:, :, 2), HH_R, kpt, real_lattice, mp_grid, 2, num_wann, &
                                   error, comm)
    if (allocated(error)) return
    call pw90common_fourier_R_to_k(ws_region, wannier_data, ws_distance, wigner_seitz, &
                                   delHH(:, :, 3), HH_R, kpt, real_lattice, mp_grid, 3, num_wann, &
                                   error, comm)
    if (allocated(error)) return
    call wham_get_deleig_a(del_eig(:, 1), eig, delHH(:, :, 1), UU, num_wann, &
                           pw90_band_deriv_degen, error, comm)
    if (allocated(error)) return
    call wham_get_deleig_a(del_eig(:, 2), eig, delHH(:, :, 2), UU, num_wann, &
                           pw90_band_deriv_degen, error, comm)
    if (allocated(error)) return
    call wham_get_deleig_a(del_eig(:, 3), eig, delHH(:, :, 3), UU, num_wann, &
                           pw90_band_deriv_degen, error, comm)
    if (allocated(error)) return

  end subroutine wham_get_eig_deleig

  !================================================!
  subroutine wham_get_eig_deleig_TB_conv(pw90_band_deriv_degen, delHH, UU, eig, del_eig, num_wann, &
                                         error, comm)
    !================================================!
    ! modified version of wham_get_eig_deleig for the TB convention
    ! avoids recalculating delHH and UU, works with input values
    !
    !! Given a k point, this function returns eigenvalues E and
    !! derivatives of the eigenvalues dE/dk_a, using wham_get_deleig_a
    !
    !================================================!

    use w90_postw90_types, only: pw90_band_deriv_degen_type
    use w90_comms, only: w90_comm_type

    ! arguments
    type(pw90_band_deriv_degen_type), intent(in) :: pw90_band_deriv_degen
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer, intent(in) :: num_wann

    real(kind=dp), intent(out) :: del_eig(num_wann, 3)
    real(kind=dp), intent(in) :: eig(num_wann)

    complex(kind=dp), intent(in) :: delHH(:, :, :)
    !! the delHH matrix (derivative of H) at kpt
    complex(kind=dp), intent(in) :: UU(:, :)
    !! the rotation matrix that gives the eigenvectors of HH

    call wham_get_deleig_a(del_eig(:, 1), eig, delHH(:, :, 1), UU, num_wann, &
                           pw90_band_deriv_degen, error, comm)
    if (allocated(error)) return
    call wham_get_deleig_a(del_eig(:, 2), eig, delHH(:, :, 2), UU, num_wann, &
                           pw90_band_deriv_degen, error, comm)
    if (allocated(error)) return
    call wham_get_deleig_a(del_eig(:, 3), eig, delHH(:, :, 3), UU, num_wann, &
                           pw90_band_deriv_degen, error, comm)
    if (allocated(error)) return

  end subroutine wham_get_eig_deleig_TB_conv

  !================================================!
  subroutine wham_get_eig_UU_HH_JJlist(dis_manifold, fermi_energy_list, kpt_latt, ws_region, &
                                       print_output, wannier_data, ws_distance, wigner_seitz, HH, &
                                       HH_R, JJm_list, JJp_list, u_matrix, UU, v_matrix, eig, &
                                       eigval, kpt, real_lattice, scissors_shift, mp_grid, &
                                       num_bands, num_kpts, num_wann, num_valence_bands, &
                                       effective_model, have_disentangled, seedname, stdout, &
                                       timer, error, comm, occ)
    !================================================!
    !
    !! Wrapper routine used to reduce number of Fourier calls
    !    Added the optional occ parameter
    !
    !================================================!

    use w90_constants, only: dp
    use w90_postw90_common, only: pw90common_fourier_R_to_k_new_second_d, &
      pw90common_fourier_R_to_k_new
    use w90_get_oper, only: get_HH_R
    use w90_utility, only: utility_diagonalize
    use w90_types, only: print_output_type, wannier_data_type, dis_manifold_type, &
      ws_region_type, ws_distance_type, timer_list_type
    use w90_comms, only: w90_comm_type, mpirank
    use w90_postw90_types, only: wigner_seitz_type

    implicit none

    ! arguments
    type(dis_manifold_type), intent(in) :: dis_manifold
    real(kind=dp), allocatable, intent(in) :: fermi_energy_list(:)
    real(kind=dp), intent(in) :: kpt_latt(:, :)
    type(print_output_type), intent(in) :: print_output
    type(ws_region_type), intent(in) :: ws_region
    type(w90_comm_type), intent(in) :: comm
    type(wannier_data_type), intent(in) :: wannier_data
    type(wigner_seitz_type), intent(inout) :: wigner_seitz
    type(ws_distance_type), intent(inout) :: ws_distance
    type(timer_list_type), intent(inout) :: timer
    type(w90_error_type), allocatable, intent(out) :: error

    integer, intent(in) :: num_wann, num_kpts, num_bands, num_valence_bands
    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: stdout

    real(kind=dp), intent(in) :: kpt(3), real_lattice(3, 3)
    real(kind=dp), intent(out) :: eig(:)
    real(kind=dp), intent(in) :: eigval(:, :)
    real(kind=dp), intent(in) :: scissors_shift
    real(kind=dp), intent(in), optional :: occ(:)

    complex(kind=dp), intent(out) :: UU(:, :)
    complex(kind=dp), intent(out) :: HH(:, :)
    complex(kind=dp), intent(out) :: JJp_list(:, :, :, :)
    complex(kind=dp), intent(out) :: JJm_list(:, :, :, :)
    complex(kind=dp), intent(in) :: u_matrix(:, :, :), v_matrix(:, :, :)
    complex(kind=dp), allocatable, intent(inout) :: HH_R(:, :, :) !  <0n|r|Rm>

    character(len=50), intent(in) :: seedname
    logical, intent(in) :: have_disentangled
    logical, intent(in) :: effective_model

    ! local variables
    integer                       :: i
    complex(kind=dp), allocatable :: delHH(:, :, :)

    call get_HH_R(dis_manifold, kpt_latt, print_output, wigner_seitz, HH_R, u_matrix, v_matrix, &
                  eigval, real_lattice, scissors_shift, num_bands, num_kpts, num_wann, &
                  num_valence_bands, effective_model, have_disentangled, seedname, ws_distance, ws_region, &
                  stdout, timer, error, comm)
    if (allocated(error)) return

    allocate (delHH(num_wann, num_wann, 3))
    call pw90common_fourier_R_to_k_new(ws_region, wannier_data, ws_distance, wigner_seitz, HH_R, &
                                       kpt, real_lattice, mp_grid, num_wann, error, comm, OO=HH, &
                                       OO_dx=delHH(:, :, 1), OO_dy=delHH(:, :, 2), &
                                       OO_dz=delHH(:, :, 3))
    if (allocated(error)) return

    call utility_diagonalize(HH, num_wann, eig, UU, error, comm)
    if (allocated(error)) return

    do i = 1, 3
      if (present(occ)) then
        call wham_get_JJp_JJm_list(delHH(:, :, i), UU, eig, JJp_list(:, :, :, i), &
                                   JJm_list(:, :, :, i), num_wann, fermi_energy_list, occ=occ)
      else
        call wham_get_JJp_JJm_list(delHH(:, :, i), UU, eig, JJp_list(:, :, :, i), &
                                   JJm_list(:, :, :, i), num_wann, fermi_energy_list)
      endif
    enddo

  end subroutine wham_get_eig_UU_HH_JJlist

  !================================================!
  subroutine wham_get_eig_UU_HH_AA_sc_TB_conv(pw90_berry, dis_manifold, kmesh_info, kpt_latt, &
                                              ws_region, print_output, wannier_data, ws_distance, &
                                              wigner_seitz, AA_R, HH, HH_da, HH_dadb, HH_R, &
                                              u_matrix, UU, v_matrix, eig, eigval, kpt, &
                                              real_lattice, scissors_shift, mp_grid, &
                                              num_bands, num_kpts, num_wann, num_valence_bands, &
                                              effective_model, have_disentangled, seedname, &
                                              stdout, timer, error, comm)
    !================================================!
    !
    ! modified version of wham_get_eig_UU_HH_AA_sc, calls routines
    ! satisfying the TB phase convention
    !
    !================================================!

    use w90_constants, only: dp
    use w90_get_oper, only: get_HH_R, get_AA_R_effective, get_AA_R
    use w90_postw90_common, only: pw90common_fourier_R_to_k_new_second_d_TB_conv
    use w90_types, only: print_output_type, wannier_data_type, dis_manifold_type, &
      kmesh_info_type, ws_region_type, ws_distance_type, timer_list_type
    use w90_utility, only: utility_diagonalize
    use w90_postw90_types, only: pw90_berry_mod_type, wigner_seitz_type
    use w90_comms, only: w90_comm_type, mpirank

    implicit none

    ! arguments
    type(pw90_berry_mod_type), intent(in) :: pw90_berry
    type(dis_manifold_type), intent(in) :: dis_manifold
    type(kmesh_info_type), intent(in) :: kmesh_info
    real(kind=dp), intent(in) :: kpt_latt(:, :)
    type(print_output_type), intent(in) :: print_output
    type(ws_region_type), intent(in) :: ws_region
    type(w90_comm_type), intent(in) :: comm
    type(wannier_data_type), intent(in) :: wannier_data
    type(wigner_seitz_type), intent(inout) :: wigner_seitz
    type(ws_distance_type), intent(inout) :: ws_distance
    type(timer_list_type), intent(inout) :: timer
    type(w90_error_type), allocatable, intent(out) :: error

    integer, intent(in) :: num_wann, num_kpts, num_bands, num_valence_bands
    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: stdout

    real(kind=dp), intent(out) :: eig(num_wann)
    real(kind=dp), intent(in) :: eigval(:, :)
    real(kind=dp), intent(in) :: kpt(3), real_lattice(3, 3)
    real(kind=dp), intent(in) :: scissors_shift

    complex(kind=dp), intent(out) :: UU(:, :)
    complex(kind=dp), intent(out) :: HH(:, :)
    complex(kind=dp), intent(out) :: HH_da(:, :, :)
    complex(kind=dp), intent(out) :: HH_dadb(:, :, :, :)
    complex(kind=dp), intent(in) :: u_matrix(:, :, :), v_matrix(:, :, :)
    complex(kind=dp), allocatable, intent(inout) :: HH_R(:, :, :) !  <0n|r|Rm>
    complex(kind=dp), allocatable, intent(inout) :: AA_R(:, :, :, :) ! <0n|r|Rm>

    character(len=50), intent(in)  :: seedname
    logical, intent(in) :: have_disentangled
    logical, intent(in) :: effective_model

    call get_HH_R(dis_manifold, kpt_latt, print_output, wigner_seitz, HH_R, u_matrix, v_matrix, &
                  eigval, real_lattice, scissors_shift, num_bands, num_kpts, num_wann, &
                  num_valence_bands, effective_model, have_disentangled, seedname, ws_distance, ws_region, &
                  stdout, timer, error, comm)
    if (allocated(error)) return

    if (effective_model) then
      call get_AA_R_effective(print_output, AA_R, HH_R, wigner_seitz%nrpts, num_wann, seedname, &
                              stdout, timer, error, comm)
    else
      call get_AA_R(pw90_berry, dis_manifold, kmesh_info, kpt_latt, print_output, wannier_data, AA_R, &
                    v_matrix, eigval, wigner_seitz, ws_distance, ws_region, num_bands, num_kpts, &
                    num_wann, have_disentangled, seedname, stdout, timer, error, comm)
    endif

    if (allocated(error)) return

    call pw90common_fourier_R_to_k_new_second_d_TB_conv(kpt, HH_R, AA_R, num_wann, ws_region, &
                                                        wannier_data, real_lattice, mp_grid, &
                                                        ws_distance, wigner_seitz, error, comm, &
                                                        OO=HH, OO_da=HH_da(:, :, :), &
                                                        OO_dadb=HH_dadb(:, :, :, :))
    if (allocated(error)) return
    call utility_diagonalize(HH, num_wann, eig, UU, error, comm)
    if (allocated(error)) return

  end subroutine wham_get_eig_UU_HH_AA_sc_TB_conv

  subroutine wham_get_eig_UU_HH_AA_sc(dis_manifold, kpt_latt, ws_region, print_output, &
                                      wannier_data, ws_distance, wigner_seitz, HH, HH_da, HH_dadb, &
                                      HH_R, u_matrix, UU, v_matrix, eig, eigval, kpt, &
                                      real_lattice, scissors_shift, mp_grid, num_bands, num_kpts, &
                                      num_wann, num_valence_bands, effective_model, &
                                      have_disentangled, seedname, stdout, timer, error, comm)
    !================================================!
    !
    !! Wrapper routine used to reduce number of Fourier calls
    !
    !================================================!

    use w90_constants, only: dp
    use w90_get_oper, only: get_HH_R
    use w90_postw90_common, only: pw90common_fourier_R_to_k_new_second_d
    use w90_utility, only: utility_diagonalize
    use w90_comms, only: w90_comm_type, mpirank
    use w90_types, only: print_output_type, wannier_data_type, dis_manifold_type, &
      ws_region_type, ws_distance_type, timer_list_type
    use w90_postw90_types, only: wigner_seitz_type

    implicit none

    ! arguments
    type(dis_manifold_type), intent(in) :: dis_manifold
    real(kind=dp), intent(in) :: kpt_latt(:, :)
    type(print_output_type), intent(in) :: print_output
    type(ws_region_type), intent(in) :: ws_region
    type(w90_comm_type), intent(in) :: comm
    type(wannier_data_type), intent(in) :: wannier_data
    type(wigner_seitz_type), intent(inout) :: wigner_seitz
    type(ws_distance_type), intent(inout) :: ws_distance
    type(timer_list_type), intent(inout) :: timer
    type(w90_error_type), allocatable, intent(out) :: error

    integer, intent(in) :: num_wann, num_kpts, num_bands, num_valence_bands
    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: stdout

    real(kind=dp), intent(out) :: eig(num_wann)
    real(kind=dp), intent(in) :: eigval(:, :)
    real(kind=dp), intent(in) :: kpt(3), real_lattice(3, 3)
    real(kind=dp), intent(in) :: scissors_shift

    complex(kind=dp), intent(out) :: UU(:, :)
    complex(kind=dp), intent(out) :: HH(:, :)
    complex(kind=dp), intent(out) :: HH_da(:, :, :)
    complex(kind=dp), intent(out) :: HH_dadb(:, :, :, :)
    complex(kind=dp), intent(in) :: u_matrix(:, :, :), v_matrix(:, :, :)
    complex(kind=dp), allocatable, intent(inout) :: HH_R(:, :, :) !  <0n|r|Rm>

    character(len=50), intent(in)  :: seedname
    logical, intent(in) :: have_disentangled
    logical, intent(in) :: effective_model

    call get_HH_R(dis_manifold, kpt_latt, print_output, wigner_seitz, HH_R, u_matrix, v_matrix, &
                  eigval, real_lattice, scissors_shift, num_bands, num_kpts, num_wann, &
                  num_valence_bands, effective_model, have_disentangled, seedname, ws_distance, ws_region, &
                  stdout, timer, error, comm)
    if (allocated(error)) return

    call pw90common_fourier_R_to_k_new_second_d(kpt, HH_R, num_wann, ws_region, wannier_data, &
                                                real_lattice, mp_grid, ws_distance, wigner_seitz, &
                                                error, comm, OO=HH, OO_da=HH_da(:, :, :), &
                                                OO_dadb=HH_dadb(:, :, :, :))
    if (allocated(error)) return
    call utility_diagonalize(HH, num_wann, eig, UU, error, comm)
    if (allocated(error)) return

  end subroutine wham_get_eig_UU_HH_AA_sc

end module w90_wan_ham
