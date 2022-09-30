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
!  w90_boltzwann: Boltzman transport                         !
!                                                            !
!------------------------------------------------------------!

module w90_boltzwann
  !================================================!
  !! Compute Boltzman tranport properties
  !!
  !! BoltzWann routines by
  !! G. Pizzi, D. Volja, B. Kozinsky, M. Fornari and N. Marzari
  !! August, 2012
  !!
  !! Affiliations:
  !! THEOS, EPFL, Station 12, 1015 Lausanne (Switzerland)
  !! DMSE, MIT, 77 Massachusetts Ave, Cambridge, MA, 02139
  !! Central Michigan University, Mount Pleasant, MI 48859
  !! Robert Bosch LLC, Cambridge, MA, 02139
  !!
  !!
  !! Please cite the following paper when publishing results
  !! obtained using the BoltzWann module:
  !!
  !![1] G. Pizzi, D. Volja, B. Kozinsky, M. Fornari, N. Marzari
  !!    Comp. Phys. Comm. 185, 422 (2014)
  !!    DOI: 10.1016/j.cpc.2013.09.015    (arXiv:1305.1587)
  !================================================!

  use w90_comms, only: mpisize, mpirank, comms_gatherv, comms_array_split, comms_reduce, &
    comms_allreduce, w90_comm_type
  use w90_constants, only: dp, pw90_physical_constants_type, min_smearing_binwidth_ratio
  use w90_dos, only: dos_get_k, dos_get_levelspacing
  use w90_utility, only: utility_inv3, utility_inv2
  use w90_error, only: w90_error_type, set_error_alloc, set_error_dealloc, set_error_fatal, &
    set_error_input, set_error_fatal, set_error_file

  implicit none

  private

  public :: boltzwann_main

  ! Constants to identify the six components of a tensor when it is stored in packed form
  integer, parameter :: XX = 1
  integer, parameter :: XY = 2
  integer, parameter :: YY = 3
  integer, parameter :: XZ = 4
  integer, parameter :: YZ = 5
  integer, parameter :: ZZ = 6

  character(len=74), parameter :: pub_string_1 = &
                                  "Please cite the following paper when publishing results obtained using    "
  character(len=74), parameter :: pub_string_2 = &
                                  "the BoltzWann module:                                                     "
  character(len=74), parameter :: pub_string_3 = &
                                  "G. Pizzi, D. Volja, B. Kozinsky, M. Fornari, and N. Marzari,              "
  character(len=74), parameter :: pub_string_4 = &
                                  "Comp. Phys. Comm. 185, 422 (2014); DOI:10.1016/j.cpc.2013.09.015          "

contains

  !================================================!

  subroutine boltzwann_main(pw90_boltzwann, dis_manifold, pw90_dos, kpt_latt, &
                            pw90_band_deriv_degen, postw90_oper, pw90_spin, physics, ws_region, &
                            w90_system, wannier_data, ws_distance, wigner_seitz, print_output, &
                            HH_R, SS_R, v_matrix, u_matrix, eigval, real_lattice, scissors_shift, &
                            mp_grid, num_wann, num_bands, num_kpts, effective_model, &
                            have_disentangled, spin_decomp, seedname, stdout, timer, error, comm)
    !================================================!
    !! This is the main routine of the BoltzWann module.
    !! It calculates the transport coefficients using the Boltzmann transport equation.
    !!
    !! It produces six files that contain:
    !!
    !!  1. the Transport Distribution function (TDF) in units of 1/hbar^2 * eV*fs/angstrom
    !!  2. the electrical conductivity in SI units (1/Ohm/m)
    !!  3. the tensor sigma*S (sigma=el.cond., S=seebeck) in SI units (Ampere/meter/K)
    !!  4. the Seebeck coefficient in SI units (V/K)
    !!  5. the thermal conductivity in SI units (W/meter/K)
    !!  6. if requested, the density of states
    !!
    !! Files from 2 to 4 are output on a grid of (mu,T) points, where mu is the chemical potential in eV and
    !! T is the temperature in Kelvin. The grid is defined in the input.
    !================================================!

    use w90_constants, only: dp
    use w90_io, only: io_stopwatch_start, io_stopwatch_stop
    use w90_comms, only: comms_bcast, w90_comm_type, mpirank
    use w90_types, only: dis_manifold_type, print_output_type, wannier_data_type, &
      ws_region_type, w90_system_type, ws_distance_type, timer_list_type
    use w90_postw90_types, only: pw90_boltzwann_type, pw90_spin_mod_type, &
      pw90_band_deriv_degen_type, pw90_dos_mod_type, pw90_oper_read_type, wigner_seitz_type

    implicit none

    ! arguments
    type(pw90_boltzwann_type), intent(in) :: pw90_boltzwann
    type(dis_manifold_type), intent(in) :: dis_manifold
    type(pw90_dos_mod_type), intent(in) :: pw90_dos
    type(pw90_band_deriv_degen_type), intent(in) :: pw90_band_deriv_degen
    type(pw90_oper_read_type), intent(in) :: postw90_oper
    type(pw90_spin_mod_type), intent(in) :: pw90_spin
    type(print_output_type), intent(in) :: print_output
    type(pw90_physical_constants_type), intent(in) :: physics
    type(ws_region_type), intent(in) :: ws_region
    type(w90_comm_type), intent(in) :: comm
    type(w90_system_type), intent(in) :: w90_system
    type(wannier_data_type), intent(in) :: wannier_data
    type(wigner_seitz_type), intent(inout) :: wigner_seitz
    type(ws_distance_type), intent(inout) :: ws_distance
    type(timer_list_type), intent(inout) :: timer
    type(w90_error_type), allocatable, intent(out) :: error

    complex(kind=dp), allocatable, intent(inout) :: HH_R(:, :, :) !  <0n|r|Rm>
    complex(kind=dp), allocatable, intent(inout) :: SS_R(:, :, :, :) ! <0n|sigma_x,y,z|Rm>
    complex(kind=dp), intent(in) :: v_matrix(:, :, :), u_matrix(:, :, :)

    real(kind=dp), intent(in) :: eigval(:, :)
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    real(kind=dp), intent(in) :: scissors_shift
    real(kind=dp), intent(in) :: kpt_latt(:, :)

    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: num_wann, num_bands, num_kpts
    integer, intent(in) :: stdout

    character(len=50), intent(in) :: seedname
    logical, intent(in) :: have_disentangled
    logical, intent(in) :: spin_decomp
    logical, intent(in) :: effective_model

    ! local variables
    integer :: TempNumPoints, MuNumPoints, TDFEnergyNumPoints
    integer :: i, j, ierr, EnIdx, TempIdx, MuIdx
    real(kind=dp), allocatable :: TempArray(:), MuArray(:), KTArray(:)
    real(kind=dp), allocatable :: TDF(:, :, :) ! (coordinate,Energy)
    real(kind=dp), allocatable :: TDFEnergyArray(:)
    real(kind=dp), allocatable :: IntegrandArray(:, :) ! (coordinate, Energy) at a given T and mu
    real(kind=dp) :: SigmaS_FP(3, 3), ThisElCond(3, 3), ElCondInverse(3, 3), ThisSeebeck(3, 3)
    real(kind=dp) :: ThisElCond2d(2, 2), ElCondInverse2d(2, 2)
    !real(kind=dp), dimension(6) :: ElCondTimesSeebeck
    real(kind=dp), allocatable :: ElCond(:, :, :) ! (coordinate,Temp, mu)
    real(kind=dp), allocatable :: SigmaS(:, :, :) ! (coordinate,Temp, mu)
    real(kind=dp), allocatable :: Seebeck(:, :, :) ! (coordinate,Temp, mu)
    real(kind=dp), allocatable :: Kappa(:, :, :) ! (coordinate,Temp, mu)
    real(kind=dp), allocatable :: LocalElCond(:, :) ! (coordinate,Temp+mu combined index)
    real(kind=dp), allocatable :: LocalSigmaS(:, :) ! (coordinate,Temp+mu combined index)
    real(kind=dp), allocatable :: LocalSeebeck(:, :) ! (coordinate,Temp+mu combined index)
    real(kind=dp), allocatable :: LocalKappa(:, :) ! (coordinate,Temp+mu combined index)
    real(kind=dp) :: Determinant
    integer :: tdf_unit, elcond_unit, sigmas_unit, seebeck_unit, kappa_unit, ndim

    integer :: LocalIdx, GlobalIdx
    ! I also add 3 times the smearing on each side of the TDF energy array to take into account also possible smearing effects
    real(kind=dp), parameter :: TDF_exceeding_energy_times_smr = 3._dp
    real(kind=dp) :: TDF_exceeding_energy
    real(kind=dp) :: cell_volume
    integer :: NumberZeroDet

    integer, allocatable :: counts(:), displs(:)
    integer :: my_node_id, num_nodes
    logical :: on_root = .false.

    my_node_id = mpirank(comm)
    num_nodes = mpisize(comm)
    if (my_node_id == 0) on_root = .true.
    allocate (counts(0:num_nodes - 1))
    allocate (displs(0:num_nodes - 1))

    cell_volume = real_lattice(1, 1)*(real_lattice(2, 2)*real_lattice(3, 3) - real_lattice(3, 2)*real_lattice(2, 3)) + &
                  real_lattice(1, 2)*(real_lattice(2, 3)*real_lattice(3, 1) - real_lattice(3, 3)*real_lattice(2, 1)) + &
                  real_lattice(1, 3)*(real_lattice(2, 1)*real_lattice(3, 2) - real_lattice(3, 1)*real_lattice(2, 2))

    if (print_output%iprint > 0 .and. print_output%timing_level > 0) call io_stopwatch_start('boltzwann_main', timer)

    if (print_output%iprint > 0) then
      write (stdout, *)
      write (stdout, '(1x,a)') '*---------------------------------------------------------------------------*'
      write (stdout, '(1x,a)') '|                   Boltzmann Transport (BoltzWann module)                  |'
      write (stdout, '(1x,a)') '*---------------------------------------------------------------------------*'
      write (stdout, '(1x,a)') '| '//pub_string_1//'|'
      write (stdout, '(1x,a)') '| '//pub_string_2//'|'
      write (stdout, '(1x,a)') '| '//pub_string_3//'|'
      write (stdout, '(1x,a)') '| '//pub_string_4//'|'
      write (stdout, '(1x,a)') '*---------------------------------------------------------------------------*'
      write (stdout, *)
    end if

    if (print_output%iprint > 0) then
      if (pw90_boltzwann%dir_num_2d /= 0) then
        write (stdout, '(1x,a)') '>                                                                           <'
        write (stdout, '(1x,a)') '> NOTE! Using the 2D version for the calculation of the Seebeck             <'
        if (pw90_boltzwann%dir_num_2d == 1) then
          write (stdout, '(1x,a)') '>       coefficient, where the non-periodic direction is x.                 <'
        elseif (pw90_boltzwann%dir_num_2d == 2) then
          write (stdout, '(1x,a)') '>       coefficient, where the non-periodic direction is y.                 <'
        elseif (pw90_boltzwann%dir_num_2d == 3) then
          write (stdout, '(1x,a)') '>       coefficient, where the non-periodic direction is z.                 <'
        end if
        write (stdout, '(1x,a)') '>                                                                           <'
        write (stdout, '(1x,a)') ''
      end if
    end if

    ! separate error condition from info printout above (which may be avoided entirely in lib mode?)
    if (pw90_boltzwann%dir_num_2d < 0 .or. pw90_boltzwann%dir_num_2d > 3) then
      call set_error_input(error, 'Unrecognized value of pw90_boltzwann_2d_dir_num', comm)
      return
    endif

    ! I precalculate the TempArray and the MuArray
    TempNumPoints = int(floor((pw90_boltzwann%temp_max - pw90_boltzwann%temp_min)/pw90_boltzwann%temp_step)) + 1
    allocate (TempArray(TempNumPoints), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating TempArray in boltzwann_main', comm)
      return
    endif
    do i = 1, TempNumPoints
      TempArray(i) = pw90_boltzwann%temp_min + real(i - 1, dp)*pw90_boltzwann%temp_step
    end do

    ! This array contains the same temperatures of the TempArray, but multiplied by k_boltzmann, in units of eV
    allocate (KTArray(TempNumPoints), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating KTArray in boltzwann_main', comm)
      return
    endif
    ! (k_B in eV/kelvin is equal to k_B_SI / elem_charge_SI)
    KTArray = TempArray*physics%k_B_SI/physics%elem_charge_SI

    MuNumPoints = int(floor((pw90_boltzwann%mu_max - pw90_boltzwann%mu_min)/pw90_boltzwann%mu_step)) + 1
    allocate (MuArray(MuNumPoints), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating MuArray in boltzwann_main', comm)
      return
    endif
    do i = 1, MuNumPoints
      MuArray(i) = pw90_boltzwann%mu_min + real(i - 1, dp)*pw90_boltzwann%mu_step
    end do

    if (pw90_boltzwann%tdf_smearing%use_adaptive) then
      call set_error_input(error, 'Adaptive smearing not allowed in Boltzwann TDF', comm)
      return
    endif
    ! I precalculate the TDFEnergyArray
    ! I assume that dis_win_min and dis_win_max are set to sensible values, related to the max and min energy
    ! This is true if the .eig file is present. I can assume its presence since we need it to interpolate the
    ! bands.
    ! I also add 3 times the smearing on each side of the TDF energy array to take into account also possible smearing effects,
    ! or at least 0.2 eV
    TDF_exceeding_energy = max(TDF_exceeding_energy_times_smr*pw90_boltzwann%tdf_smearing%fixed_width, 0.2_dp)
    TDFEnergyNumPoints = int(floor((dis_manifold%win_max - dis_manifold%win_min &
                                    + 2._dp*TDF_exceeding_energy)/pw90_boltzwann%tdf_energy_step)) + 1
    if (TDFEnergyNumPoints .eq. 1) TDFEnergyNumPoints = 2
    allocate (TDFEnergyArray(TDFEnergyNumPoints), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating TDFEnergyArray in boltzwann_main', comm)
      return
    endif
    do i = 1, TDFEnergyNumPoints
      TDFEnergyArray(i) = dis_manifold%win_min - TDF_exceeding_energy &
                          + real(i - 1, dp)*pw90_boltzwann%tdf_energy_step
    end do

    if (spin_decomp) then
      ndim = 3
    else
      ndim = 1
    end if

    ! I allocate the array for the TDF
    allocate (TDF(6, TDFEnergyNumPoints, ndim), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating TDF in boltzwann_main', comm)
      return
    endif

    ! I call the subroutine that calculates the Transport Distribution Function
    call calcTDFandDOS(pw90_boltzwann, dis_manifold, pw90_dos, kpt_latt, postw90_oper, &
                       pw90_band_deriv_degen, pw90_spin, ws_region, print_output, wannier_data, &
                       ws_distance, wigner_seitz, HH_R, SS_R, u_matrix, v_matrix, eigval, &
                       real_lattice, TDF, TDFEnergyArray, cell_volume, scissors_shift, mp_grid, &
                       num_bands, num_kpts, num_wann, w90_system%num_valence_bands, &
                       w90_system%num_elec_per_state, effective_model, have_disentangled, &
                       spin_decomp, seedname, stdout, timer, error, comm)
    ! The TDF array contains now the TDF, or more precisely
    ! hbar^2 * TDF in units of eV * fs / angstrom

    ! I print on file the TDF
    if (on_root) then
      open (newunit=tdf_unit, file=trim(seedname)//'_tdf.dat')
      write (tdf_unit, '(A)') "# Written by the BoltzWann module of the Wannier90 code."
      write (tdf_unit, '(A)') "# Transport distribution function (in units of 1/hbar^2 * eV * fs / angstrom)"// &
        " vs energy in eV"
      write (tdf_unit, '(A)') "# Content of the columns:"
      write (tdf_unit, '(A)') "# Energy TDF_xx TDF_xy TDF_yy TDF_xz TDF_yz TDF_zz"
      write (tdf_unit, '(A)') '#   (if spin decomposition is required, 12 further columns are provided, with the 6'
      write (tdf_unit, '(A)') '#    components of the TDF for the spin up, followed by those for the spin down)'
      do i = 1, size(TDFEnergyArray)
        if (ndim .eq. 1) then
          write (tdf_unit, 101) TDFEnergyArray(i), TDF(:, i, 1)
        else
          write (tdf_unit, 102) TDFEnergyArray(i), TDF(:, i, :)
        end if
      end do
      close (tdf_unit)
      if (print_output%iprint > 1) &
        write (stdout, '(3X,A)') "Transport distribution function written on the "//trim(seedname)//"_tdf.dat file."
    end if

    ! *********************************************************************************
    ! I got the TDF and I printed it. Now I use it to calculate the transport properties.

    if (on_root .and. (print_output%timing_level > 0)) call io_stopwatch_start('boltzwann_main: calc_props', timer)

    ! I obtain the counts and displs arrays, which tell how I should partition a big array
    ! on the different nodes.
    call comms_array_split(TempNumPoints*MuNumPoints, counts, displs, comm)

    ! I allocate the arrays for the spectra
    ! Allocate at least 1 entry
    allocate (LocalElCond(6, max(1, counts(my_node_id))), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating LocalElCond in boltzwann_main', comm)
      return
    endif
    allocate (LocalSigmaS(6, max(1, counts(my_node_id))), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating LocalSigmaS in boltzwann_main', comm)
      return
    endif
    allocate (LocalSeebeck(9, max(1, counts(my_node_id))), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating LocalSeebeck in boltzwann_main', comm)
      return
    endif
    allocate (LocalKappa(6, max(1, counts(my_node_id))), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating LocalKappa in boltzwann_main', comm)
      return
    endif
    LocalSigmaS = 0._dp
    LocalElCond = 0._dp
    LocalSeebeck = 0._dp
    LocalKappa = 0._dp

    ! I allocate the array that I will use to store the functions to be integrated
    allocate (IntegrandArray(6, TDFEnergyNumPoints), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating FermiDerivArray in boltzwann_main', comm)
      return
    endif

    NumberZeroDet = 0
    ! Now, I calculate the various spectra for all mu and T values
    do LocalIdx = 1, counts(my_node_id)
      ! GlobalIdx is an index from 0 to TempNumPoints*MuNumPoints-1
      GlobalIdx = displs(my_node_id) + (LocalIdx - 1)
      ! MuIdx goes from 1 to MuNumPoints
      MuIdx = GlobalIdx/TempNumPoints + 1
      ! TempIdx goes from 1 to TempNumPoints
      TempIdx = GlobalIdx - TempNumPoints*(MuIdx - 1) + 1

      ! For the calculation of the properties, I only use the total TDF (i.e., unresolved for spin)
      IntegrandArray = TDF(:, :, 1)
      do EnIdx = 1, TDFEnergyNumPoints
        IntegrandArray(:, EnIdx) = IntegrandArray(:, EnIdx)* &
                                   MinusFermiDerivative(E=TDFEnergyArray(EnIdx), mu=MuArray(MuIdx), KT=KTArray(TempIdx))
      end do
      ! Now, IntegrandArray contains (-dn/dE) * TDF_ij(E), where n is the Fermi distribution function
      ! Its integral is ElCond_ij/e^2
      LocalElCond(:, LocalIdx) = sum(IntegrandArray, DIM=2)*pw90_boltzwann%tdf_energy_step
      ! ElCond contains now (hbar^2/e^2) * sigma in eV*fs/angstrom, where sigma is the conductivity tensor
      ! (note that MinusFermiDerivative is in units of 1/eV, so that when I perform the integration
      ! ElCond has the same units of TDF)

      ! I store in ThisElCond the conductivity tensor in standard format
      ! Store both the upper and lower triangle
      do j = 1, 3
        do i = 1, j
          ThisElCond(i, j) = LocalElCond(i + ((j - 1)*j)/2, LocalIdx)
          ThisElCond(j, i) = LocalElCond(i + ((j - 1)*j)/2, LocalIdx)
        end do
      end do

      ! I calculate the inverse matrix of the conductivity
      if (pw90_boltzwann%dir_num_2d /= 0) then
        ! Invert only the appropriate 2x2 submatrix
        if (pw90_boltzwann%dir_num_2d == 1) then
          ThisElCond2d(1, 1) = ThisElCond(2, 2)
          ThisElCond2d(1, 2) = ThisElCond(2, 3)
          ThisElCond2d(2, 1) = ThisElCond(3, 2)
          ThisElCond2d(2, 2) = ThisElCond(3, 3)
        elseif (pw90_boltzwann%dir_num_2d == 2) then
          ThisElCond2d(1, 1) = ThisElCond(1, 1)
          ThisElCond2d(1, 2) = ThisElCond(1, 3)
          ThisElCond2d(2, 1) = ThisElCond(3, 1)
          ThisElCond2d(2, 2) = ThisElCond(3, 3)
        elseif (pw90_boltzwann%dir_num_2d == 3) then
          ThisElCond2d(1, 1) = ThisElCond(1, 1)
          ThisElCond2d(1, 2) = ThisElCond(1, 2)
          ThisElCond2d(2, 1) = ThisElCond(2, 1)
          ThisElCond2d(2, 2) = ThisElCond(2, 2)
          ! I do not do the else case because this was already checked at the beginning
          ! of this routine
        end if
        call utility_inv2(ThisElCond2d, ElCondInverse2d, Determinant)
        ElCondInverse = 0._dp ! Other elements must be set to zero
        if (pw90_boltzwann%dir_num_2d == 1) then
          ElCondInverse(2, 2) = ElCondInverse2d(1, 1)
          ElCondInverse(2, 3) = ElCondInverse2d(1, 2)
          ElCondInverse(3, 2) = ElCondInverse2d(2, 1)
          ElCondInverse(3, 3) = ElCondInverse2d(2, 2)
        elseif (pw90_boltzwann%dir_num_2d == 2) then
          ElCondInverse(1, 1) = ElCondInverse2d(1, 1)
          ElCondInverse(1, 3) = ElCondInverse2d(1, 2)
          ElCondInverse(3, 1) = ElCondInverse2d(2, 1)
          ElCondInverse(3, 3) = ElCondInverse2d(2, 2)
        elseif (pw90_boltzwann%dir_num_2d == 3) then
          ElCondInverse(1, 1) = ElCondInverse2d(1, 1)
          ElCondInverse(1, 2) = ElCondInverse2d(1, 2)
          ElCondInverse(2, 1) = ElCondInverse2d(2, 1)
          ElCondInverse(2, 2) = ElCondInverse2d(2, 2)
          ! I do not do the else case because this was already checked at the beginning
          ! of this routine
        end if
      else
        call utility_inv3(ThisElCond, ElCondInverse, Determinant)
      end if
      if (Determinant .eq. 0._dp) then
        NumberZeroDet = NumberZeroDet + 1
        ! If the determinant is zero (i.e., zero conductivity along a given direction)
        ! I set the Inverse to zero, so that the Seebeck is zero. I will also issue
        ! a warning.
        ElCondInverse = 0._dp
      else
        ! The routine returns the adjoint of ThisElCond; in order to get
        ! the inverse, I have to calculate ElCondInverse / Determinant
        ElCondInverse = ElCondInverse/Determinant
      end if

      ! Now, I multiply IntegrandArray by (E-mu): then, IntegrandArray contains
      ! (-dn/dE) * TDF_ij(E) * (E-mu) and its integral is (ElCond*Seebeck)_ij * T / e, where
      ! T is the temperature
      do EnIdx = 1, TDFEnergyNumPoints
        IntegrandArray(:, EnIdx) = IntegrandArray(:, EnIdx)*(TDFEnergyArray(EnIdx) - MuArray(MuIdx))
      end do

      ! I store in SigmaS_FP the product of the two tensors in full-packed format
      LocalSigmaS(:, LocalIdx) = sum(IntegrandArray, DIM=2)*pw90_boltzwann%tdf_energy_step/TempArray(TempIdx)
      do j = 1, 3
        do i = 1, j
          ! Both upper and lower diagonal
          SigmaS_FP(i, j) = LocalSigmaS(i + ((j - 1)*j)/2, LocalIdx)
          SigmaS_FP(j, i) = LocalSigmaS(i + ((j - 1)*j)/2, LocalIdx)
        end do
      end do
      ! Now, LocalSigmaS (and SigmaS_FP) contain
      ! [ElCond*Seebeck] * hbar^2/e in units of eV^2*fs/angstrom/kelvin

      ! I calculate ElCond^(-1) . (LocalSigmaS) = Seebeck in fully-packed format and then
      ! store it in the LocalSeebeck array (not in packed format because, unless S and sigma
      ! commute, there is no a-priori reason for which S should be symmetric - even if
      ! probably one can find physical reasons).

      ! I invert the sign because the electron charge is < 0
      ThisSeebeck = -matmul(ElCondInverse, SigmaS_FP)
      ! Reshuffle in 1D; order: xx, xy, xz, yx, yy, yz, zx, zy, zz
      ! Note that this is different
      do j = 1, 3
        do i = 1, 3
          LocalSeebeck((i - 1)*3 + j, LocalIdx) = ThisSeebeck(i, j)
        end do
      end do

      ! Now, Seebeck contains the Seebeck coefficient in volt / Kelvin. In fact:
      ! - ElCond contains (hbar^2/e^2) * sigma in eV*fs/angstrom, where sigma is the value of the
      !   conductivity, i.e. it is the conductivity in units of (e^2/hbar^2) * eV * fs / angstrom
      ! - ElCondInverse is in thus in units of (hbar^2/e^2) / eV / fs * angstrom
      ! - ElCondTimesSeebeck is in units of  e / hbar^2 * eV^2 * fs / angstrom / Kelvin
      ! therefore ThisSeebeck, which has the units of ElCondInverse * ElCondTimesSeebeck, i.e.
      ! [(hbar^2/e^2) / eV / fs * angstrom] * [ e / hbar^2 * eV^2 * fs / angstrom / kelvin] =
      ! eV/e/kelvin = volt/kelvin

      ! Now, I multiply IntegrandArray by (E-mu): then, IntegrandArray contains
      ! (-dn/dE) * TDF_ij(E) * (E-mu)^2 and its integral is (Kappa)_ij * T
      do EnIdx = 1, TDFEnergyNumPoints
        IntegrandArray(:, EnIdx) = IntegrandArray(:, EnIdx)*(TDFEnergyArray(EnIdx) - MuArray(MuIdx))
      end do
      LocalKappa(:, LocalIdx) = sum(IntegrandArray, DIM=2)*pw90_boltzwann%tdf_energy_step/TempArray(TempIdx)
      ! Kappa contains now the thermal conductivity in units of
      ! 1/hbar^2 * eV^3*fs/angstrom/kelvin
    end do
    ! I check if there were (mu,T) pairs for which we got sigma = 0
    call comms_reduce(NumberZeroDet, 1, 'SUM', error, comm)
    if (allocated(error)) return
    if (on_root) then
      if ((NumberZeroDet .gt. 0)) then
        write (stdout, '(1X,A,I0,A)') "> WARNING! There are ", NumberZeroDet, " (mu,T) pairs for which the electrical"
        write (stdout, '(1X,A)') ">          conductivity has zero determinant."
        write (stdout, '(1X,A)') ">          Seebeck coefficient set to zero for those pairs."
        write (stdout, '(1X,A)') ">          Check if this is physical or not."
        write (stdout, '(1X,A)') ">          (If you are dealing with a 2D system, set the pw90_boltzwann_2d_dir flag.)"
        write (stdout, '(1X,A)') ""
      end if
    end if

    ! Now, I multiply by the correct factors to obtain the tensors in SI units

    ! **** Electrical conductity ****
    ! Now, ElCond is in units of (e^2/hbar^2) * eV * fs / angstrom
    ! I want the conductivity in units of 1/Ohm/meter (i.e., SI units). The conversion factor is then
    ! e^2/hbar^2 * eV*fs/angstrom * Ohm * meter; we use the constants in constants.F90 by noting that:
    ! Ohm = V/A = V*s/C, where A=ampere, C=coulomb
    ! Then we have: (e/C)^3/hbar^2 * V*s^2 * V * C^2 * (meter / angstrom)  * (fs / s)
    ! Now: e/C = elem_charge_SI; CV=Joule, CV/s=Watt, hbar/Watt = hbar_SI;
    ! moreover meter / angstrom = 1e10, fs / s = 1e-15 so that we finally get the
    ! CONVERSION FACTOR: elem_charge_SI**3 / (hbar_SI**2) * 1.e-5_dp
    LocalElCond = LocalElCond*physics%elem_charge_SI**3/(physics%hbar_SI**2)*1.e-5_dp
    ! THIS IS NOW THE ELECTRICAL CONDUCTIVITY IN SI UNITS, i.e. in 1/Ohm/meter

    ! *** Sigma * S ****
    ! Again, as above or below for Kappa, the conversion factor is
    ! * elem_charge_SI**3 / (hbar_SI**2) * 1.e-5_dp
    ! and brings the result to Ampere/m/K
    LocalSigmaS = LocalSigmaS*physics%elem_charge_SI**3/(physics%hbar_SI**2)*1.e-5_dp

    ! **** Seebeck coefficient ****
    ! THE SEEECK COEFFICIENTS IS ALREADY IN volt/kelvin, so nothing has to be done

    ! **** K coefficient (approx. the Thermal conductivity) ****
    ! Now, Kappa is in units of 1/hbar^2 * eV^3*fs/angstrom/K
    ! I want it to be in units of W/m/K (i.e., SI units). Then conversion factor is then
    ! 1/hbar^2 * eV^3 * fs/angstrom/K / W * m * K =
    ! 1/hbar^2 * eV^3 * fs / W * (m/angstrom) =
    ! 1/hbar^2 * C^3 * V^3 * s / W * [e/C]^3 * (m/angstrom) * (fs / s) =
    ! 1/hbar^2 * J^2 * [e/C]^3 * (m/angstrom) * (fs / s) =
    ! elem_charge_SI**3 / (hbar_SI**2) * 1.e-5_dp, i.e. the same conversion factor as above
    LocalKappa = LocalKappa*physics%elem_charge_SI**3/(physics%hbar_SI**2)*1.e-5_dp
    ! THIS IS NOW THE THERMAL CONDUCTIVITY IN SI UNITS, i.e. in W/meter/K

    ! Now I send the different pieces to the local node
    if (on_root) then
      allocate (ElCond(6, TempNumPoints, MuNumPoints), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating ElCond in boltzwann_main', comm)
        return
      endif
      allocate (SigmaS(6, TempNumPoints, MuNumPoints), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating SigmaS in boltzwann_main', comm)
        return
      endif
      allocate (Seebeck(9, TempNumPoints, MuNumPoints), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating Seebeck in boltzwann_main', comm)
        return
      endif
      allocate (Kappa(6, TempNumPoints, MuNumPoints), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating Kappa in boltzwann_main', comm)
        return
      endif
    else
      ! In principle, this should not be needed, because we use ElCond,
      ! Seebeck and Kappa only on the root node. However, since all
      ! processors call comms_gatherv a few lines below, and one argument
      ! is ElCond(1,1,1), some compilers complain.
      allocate (ElCond(1, 1, 1), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating ElCond in boltzwann_main (2)', comm)
        return
      endif
      allocate (SigmaS(1, 1, 1), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating SigmaS in boltzwann_main (2)', comm)
        return
      endif
      allocate (Seebeck(1, 1, 1), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating Seebeck in boltzwann_main (2)', comm)
        return
      endif
      allocate (Kappa(1, 1, 1), stat=ierr)
      if (ierr /= 0) then
        call set_error_alloc(error, 'Error in allocating Kappa in boltzwann_main (2)', comm)
        return
      endif
    end if

    ! The 6* factors are due to the fact that for each (T,mu) pair we have 6 components (xx,xy,yy,xz,yz,zz)
    ! NOTE THAT INSTEAD SEEBECK IS A FULL MATRIX AND HAS 9 COMPONENTS!
    call comms_gatherv(LocalElCond, 6*counts(my_node_id), ElCond, 6*counts, 6*displs, error, comm)
    if (allocated(error)) return
    call comms_gatherv(LocalSigmaS, 6*counts(my_node_id), SigmaS, 6*counts, 6*displs, error, comm)
    if (allocated(error)) return
    call comms_gatherv(LocalSeebeck, 9*counts(my_node_id), Seebeck, 9*counts, 9*displs, error, comm)
    if (allocated(error)) return
    call comms_gatherv(LocalKappa, 6*counts(my_node_id), Kappa, 6*counts, 6*displs, error, comm)
    if (allocated(error)) return

    if (on_root .and. (print_output%timing_level > 0)) call io_stopwatch_stop('boltzwann_main: calc_props', timer)

    ! Open files and print
    if (on_root) then
      open (newunit=elcond_unit, file=trim(seedname)//'_elcond.dat')
      write (elcond_unit, '(A)') "# Written by the BoltzWann module of the Wannier90 code."
      write (elcond_unit, '(A)') "# [Electrical conductivity in SI units, i.e. in 1/Ohm/m]"
      write (elcond_unit, '(A)') "# Mu(eV) Temp(K) ElCond_xx ElCond_xy ElCond_yy ElCond_xz ElCond_yz ElCond_zz"
      do MuIdx = 1, MuNumPoints
        do TempIdx = 1, TempNumPoints
          write (elcond_unit, 103) MuArray(MuIdx), TempArray(TempIdx), ElCond(:, TempIdx, MuIdx)
        end do
      end do
      close (elcond_unit)
      if (print_output%iprint > 1) &
        write (stdout, '(3X,A)') "Electrical conductivity written on the "//trim(seedname)//"_elcond.dat file."

      open (newunit=sigmas_unit, file=trim(seedname)//'_sigmas.dat')
      write (sigmas_unit, '(A)') "# Written by the BoltzWann module of the Wannier90 code."
      write (sigmas_unit, '(A)') "# [(Electrical conductivity * Seebeck coefficient) in SI units, i.e. in Ampere/m/K]"
      write (sigmas_unit, '(A)') "# Mu(eV) Temp(K) (Sigma*S)_xx (Sigma*S)_xy (Sigma*S)_yy (Sigma*S)_xz (Sigma*S)_yz (Sigma*S)_zz"
      do MuIdx = 1, MuNumPoints
        do TempIdx = 1, TempNumPoints
          write (sigmas_unit, 103) MuArray(MuIdx), TempArray(TempIdx), SigmaS(:, TempIdx, MuIdx)
        end do
      end do
      close (sigmas_unit)
      if (print_output%iprint > 1) write (stdout, '(3X,A)') &
        "sigma*S (sigma=el. conductivity, S=Seebeck coeff.) written on the "//trim(seedname)//"_sigmas.dat file."

      open (newunit=seebeck_unit, file=trim(seedname)//'_seebeck.dat')
      write (seebeck_unit, '(A)') "# Written by the BoltzWann module of the Wannier90 code."
      write (seebeck_unit, '(A)') "# [Seebeck coefficient in SI units, i.e. in V/K]"
      write (seebeck_unit, '(A)') &
        "# Mu(eV) Temp(K) Seebeck_xx Seebeck_xy Seebeck_xz Seebeck_yx Seebeck_yy Seebeck_yz Seebeck_zx Seebeck_zy Seebeck_zz"
      do MuIdx = 1, MuNumPoints
        do TempIdx = 1, TempNumPoints
          write (seebeck_unit, 104) MuArray(MuIdx), TempArray(TempIdx), Seebeck(:, TempIdx, MuIdx)
        end do
      end do
      close (seebeck_unit)
      if (print_output%iprint > 1) &
        write (stdout, '(3X,A)') "Seebeck coefficient written on the "//trim(seedname)//"_seebeck.dat file."

      open (newunit=kappa_unit, file=trim(seedname)//'_kappa.dat')
      write (kappa_unit, '(A)') "# Written by the BoltzWann module of the Wannier90 code."
      write (kappa_unit, '(A)') "# [K coefficient in SI units, i.e. in W/m/K]"
      write (kappa_unit, '(A)') "# [the K coefficient is defined in the documentation, and is an ingredient of"
      write (kappa_unit, '(A)') "#  the thermal conductivity. See the docs for further information.]"
      write (kappa_unit, '(A)') "# Mu(eV) Temp(K) Kappa_xx Kappa_xy Kappa_yy Kappa_xz Kappa_yz Kappa_zz"
      do MuIdx = 1, MuNumPoints
        do TempIdx = 1, TempNumPoints
          write (kappa_unit, 103) MuArray(MuIdx), TempArray(TempIdx), Kappa(:, TempIdx, MuIdx)
        end do
      end do
      close (kappa_unit)
      if (print_output%iprint > 1) &
        write (stdout, '(3X,A)') "K coefficient written on the "//trim(seedname)//"_kappa.dat file."
    end if

    if (on_root) then
      write (stdout, '(3X,A)') "Transport properties calculated."
      write (stdout, *)
      write (stdout, '(1x,a)') '*---------------------------------------------------------------------------*'
      write (stdout, '(1x,a)') '|                        End of the BoltzWann module                        |'
      write (stdout, '(1x,a)') '*---------------------------------------------------------------------------*'
      write (stdout, *)
    end if

    ! Before ending, I deallocate memory
    deallocate (TempArray, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating TempArray in boltzwann_main', comm)
      return
    endif
    deallocate (KTArray, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating KTArray in boltzwann_main', comm)
      return
    endif
    deallocate (MuArray, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating MuArray in boltzwann_main', comm)
      return
    endif
    deallocate (TDFEnergyArray, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating TDFEnergyArray in boltzwann_main', comm)
      return
    endif
    deallocate (TDF, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating TDF in boltzwann_main', comm)
      return
    endif
    deallocate (LocalElCond, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating LocalElCond in boltzwann_main', comm)
      return
    endif
    deallocate (LocalSigmaS, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating LocalSigmaS in boltzwann_main', comm)
      return
    endif
    deallocate (LocalSeebeck, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating LocalSeebeck in boltzwann_main', comm)
      return
    endif
    deallocate (LocalKappa, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating LocalKappa in boltzwann_main', comm)
      return
    endif

    deallocate (ElCond, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating ElCond in boltzwann_main', comm)
      return
    endif
    deallocate (SigmaS, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating SigmaS in boltzwann_main', comm)
      return
    endif
    deallocate (Seebeck, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating Seebeck in boltzwann_main', comm)
      return
    endif
    deallocate (Kappa, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating Kappa in boltzwann_main', comm)
      return
    endif

    deallocate (IntegrandArray, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating IntegrandArray in boltzwann_main', comm)
      return
    endif

    if (on_root .and. (print_output%timing_level > 0)) call io_stopwatch_stop('boltzwann_main', timer)

101 FORMAT(7G18.10)
102 FORMAT(19G18.10)
103 FORMAT(8G18.10)
104 FORMAT(11G18.10)

  end subroutine boltzwann_main

  !================================================!
  subroutine calcTDFandDOS(pw90_boltzwann, dis_manifold, pw90_dos, kpt_latt, postw90_oper, &
                           pw90_band_deriv_degen, pw90_spin, ws_region, print_output, &
                           wannier_data, ws_distance, wigner_seitz, HH_R, SS_R, u_matrix, &
                           v_matrix, eigval, real_lattice, TDF, TDFEnergyArray, cell_volume, &
                           scissors_shift, mp_grid, num_bands, num_kpts, num_wann, &
                           num_valence_bands, num_elec_per_state, effective_model, &
                           have_disentangled, spin_decomp, seedname, stdout, timer, error, comm)
    !================================================!
    !! This routine calculates the Transport Distribution Function $$\sigma_{ij}(\epsilon)$$ (TDF)
    !! in units of 1/hbar^2 * eV*fs/angstrom, and possibly the DOS.
    !!
    !! The TDFEnergyArray must be already allocated and initialized with the
    !! energies in eV before calling this routine.
    !!
    !!  The TDF array must be already allocated with dimensions 6 * size(TDFEnergyArray) * ndim, before calling
    !! this routine, where ndim=1 if spin_decomp==false, or ndim=3 if spin_decomp==true. This is not checked.
    !!
    !!  If run in parallel, at the end each processor will have a copy of the full TDF array
    !!
    !!  We assume that the TDFEnergyArray is uniformely spaced and that it has at least
    !!       two elements (no checks are performed on this).
    !!
    !! Note that the order of indices of TDF is different w.r.t. the DOS array (the energy is not the first but
    !!       the second index)
    !!
    !!  If the input flag pw90_boltzwann_bandshift is set to .true., the code will also shift the
    !!       conduction bands by a given amount, as defined by the pw90_boltzwann_bandshift_energyshift
    !!       and pw90_boltzwann_bandshift_firstband input flags.
    !!
    !================================================!

    use w90_constants, only: dp
    use w90_comms, only: comms_bcast, w90_comm_type, mpirank
    use w90_io, only: io_stopwatch_start, io_stopwatch_stop
    use w90_utility, only: utility_recip_lattice_base
    use w90_get_oper, only: get_HH_R, get_SS_R
    use w90_types, only: print_output_type, wannier_data_type, dis_manifold_type, &
      ws_region_type, ws_distance_type, timer_list_type
    use w90_postw90_types, only: pw90_boltzwann_type, pw90_spin_mod_type, &
      pw90_band_deriv_degen_type, pw90_dos_mod_type, pw90_oper_read_type, wigner_seitz_type
    use w90_readwrite, only: w90_readwrite_get_smearing_type
    use w90_wan_ham, only: wham_get_eig_deleig

    implicit none

    ! arguments
    type(pw90_boltzwann_type), intent(in) :: pw90_boltzwann
    type(dis_manifold_type), intent(in) :: dis_manifold
    type(pw90_dos_mod_type), intent(in) :: pw90_dos
    real(kind=dp), intent(in) :: kpt_latt(:, :)
    type(pw90_band_deriv_degen_type), intent(in) :: pw90_band_deriv_degen
    type(pw90_oper_read_type), intent(in) :: postw90_oper
    type(pw90_spin_mod_type), intent(in) :: pw90_spin
    type(print_output_type), intent(in) :: print_output
    type(ws_region_type), intent(in) :: ws_region
    type(w90_comm_type), intent(in) :: comm
    type(wannier_data_type), intent(in) :: wannier_data
    type(wigner_seitz_type), intent(inout) :: wigner_seitz
    type(ws_distance_type), intent(inout) :: ws_distance
    type(timer_list_type), intent(inout) :: timer
    type(w90_error_type), allocatable, intent(out) :: error

    integer, intent(in) :: num_wann, num_bands, num_kpts, num_valence_bands, num_elec_per_state
    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: stdout

    real(kind=dp), intent(out) :: TDF(:, :, :) ! (coordinate,Energy,spin)
    !! The TDF(i,EnIdx,spin) output array, where:
    !!        - i is an index from 1 to 6 giving the component of the symmetric tensor
    !!          $$ Sigma_{ij}(\eps) $$,
    !!          where 1=xx, 2=xy, 3=yy, 4=xz, 5=yz, 6=zz
    !!          as defined by the module constants XX, XY, ...
    !!          (in this way the mapping (i,j) -> i+((j-1)*j)/2 is satisfied for the packed storage of the
    !!          upper triangle [i<=j]).
    !!        - EnIdx is the index of the energies; the corresponding energy is given by
    !!          TDFEnergyArray(EndIdx) array (in eV).
    !!        - Spin may be only 1 if spin_decomp=.false. If it is instead true, 1 contains the total TDF,
    !!          2 the spin-up component and 3 the spin-up component
    real(kind=dp), intent(in) :: TDFEnergyArray(:)
    !! TDFEnergyArray The array with the energies for which the TDF is calculated, in eV
    ! Comments:
    ! issue warnings if going outside of the energy window
    ! check that we actually get hbar*velocity in eV*angstrom
    real(kind=dp), intent(in) :: eigval(:, :)
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    real(kind=dp), intent(in) :: cell_volume
    real(kind=dp), intent(in) :: scissors_shift

    complex(kind=dp), intent(in) :: v_matrix(:, :, :), u_matrix(:, :, :)
    complex(kind=dp), allocatable, intent(inout) :: HH_R(:, :, :) !  <0n|r|Rm>
    complex(kind=dp), allocatable, intent(inout) :: SS_R(:, :, :, :) ! <0n|sigma_x,y,z|Rm>

    character(len=50), intent(in) :: seedname
    logical, intent(in) :: have_disentangled
    logical, intent(in) :: spin_decomp
    logical, intent(in) :: effective_model

    ! local variables
    real(kind=dp) :: recip_lattice(3, 3), volume
    real(kind=dp) :: kpt(3), orig_kpt(3)
    integer :: loop_tot, loop_x, loop_y, loop_z, ierr

    complex(kind=dp), allocatable :: HH(:, :)
    complex(kind=dp), allocatable :: delHH(:, :, :)
    complex(kind=dp), allocatable :: UU(:, :)
    real(kind=dp) :: del_eig(num_wann, 3)
    real(kind=dp) :: eig(num_wann), levelspacing_k(num_wann)

    real(kind=dp), allocatable :: DOS_EnergyArray(:)
    real(kind=dp), allocatable :: DOS_k(:, :), TDF_k(:, :, :)
    real(kind=dp), allocatable :: DOS_all(:, :)
    real(kind=dp) :: kweight
    integer :: ndim, DOS_NumPoints, i, j, k, EnIdx

    character(len=20) :: numfieldsstr
    integer :: boltzdos_unit, NumPtsRefined

    real(kind=dp), parameter :: SPACING_THRESHOLD = 1.e-3
    real(kind=dp) :: min_spacing, max_spacing

    integer :: my_node_id, num_nodes
    logical :: on_root = .false.

    my_node_id = mpirank(comm)
    num_nodes = mpisize(comm)
    if (my_node_id == 0) on_root = .true.

    if (print_output%iprint > 0 .and. (print_output%timing_level > 0)) call io_stopwatch_start('calcTDF', timer)
    if (print_output%iprint > 0) then
      if (pw90_boltzwann%calc_also_dos) then
        write (stdout, '(3X,A)') "Calculating Transport Distribution function (TDF) and DOS..."
      else
        write (stdout, '(3X,A)') "Calculating Transport Distribution function (TDF)..."
      end if
    end if

    ! I call once the routine to calculate the Hamiltonian in real-space <0n|H|Rm>

    call get_HH_R(dis_manifold, kpt_latt, print_output, wigner_seitz, HH_R, u_matrix, v_matrix, &
                  eigval, real_lattice, scissors_shift, num_bands, num_kpts, num_wann, &
                  num_valence_bands, effective_model, have_disentangled, seedname, ws_distance, ws_region, &
                  stdout, timer, error, comm)
    if (allocated(error)) return

    if (spin_decomp) then
      ndim = 3

      call get_SS_R(dis_manifold, kpt_latt, print_output, postw90_oper, SS_R, v_matrix, eigval, &
                    wigner_seitz, ws_distance, ws_region, num_bands, num_kpts, num_wann, &
                    have_disentangled, seedname, stdout, timer, error, comm)
      if (allocated(error)) return

    else
      ndim = 1
    end if

    ! Some initial checks
    if (size(TDF, 1) /= 6 .or. size(TDF, 2) /= size(TDFEnergyArray) .or. size(TDF, 3) /= ndim) then
      call set_error_input(error, 'Wrong size for the TDF array in calcTDF', comm)
      return
    end if

    ! I zero the TDF array before starting
    TDF = 0._dp

    allocate (TDF_k(6, size(TDFEnergyArray), ndim), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating TDF_k in calcTDF', comm)
      return
    endif

    allocate (HH(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating HH in calcTDF', comm)
      return
    endif
    allocate (delHH(num_wann, num_wann, 3), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating delHH in calcTDF', comm)
      return
    endif
    allocate (UU(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating UU in calcTDF', comm)
      return
    endif

    DOS_NumPoints = int(floor((pw90_boltzwann%dos_energy_max - pw90_boltzwann%dos_energy_min)/pw90_boltzwann%dos_energy_step)) + 1
    if (DOS_NumPoints .eq. 1) DOS_NumPoints = 2
    allocate (DOS_EnergyArray(DOS_NumPoints), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating DOS_EnergyArray in calcTDF', comm)
      return
    endif
    do i = 1, DOS_NumPoints
      DOS_EnergyArray(i) = pw90_boltzwann%dos_energy_min + real(i - 1, dp)*pw90_boltzwann%dos_energy_step
    end do

    allocate (DOS_k(size(DOS_EnergyArray), ndim), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating DOS_k in calcTDF', comm)
      return
    endif
    allocate (DOS_all(size(DOS_EnergyArray), ndim), stat=ierr)
    if (ierr /= 0) then
      call set_error_alloc(error, 'Error in allocating DOS_all in calcTDF', comm)
      return
    endif
    dos_all = 0.0_dp

    ! I open the output files
    if (pw90_boltzwann%calc_also_dos .and. on_root) then
      open (newunit=boltzdos_unit, file=trim(seedname)//'_boltzdos.dat')
    end if

    if (pw90_boltzwann%calc_also_dos .and. on_root .and. (print_output%iprint > 1)) then
      write (stdout, '(5X,A)') "Smearing for DOS: "
      if (pw90_boltzwann%dos_smearing%use_adaptive) then
        write (stdout, '(7X,A)') trim(w90_readwrite_get_smearing_type(pw90_boltzwann%dos_smearing%type_index))//", adaptive"
      else
        if (pw90_boltzwann%dos_smearing%fixed_width/(DOS_EnergyArray(2) - DOS_EnergyArray(1)) < &
            min_smearing_binwidth_ratio) then
          write (stdout, '(7X,A)') "Unsmeared (use smearing width larger than bin width to smear)"
        else
          write (stdout, '(7X,A,G18.10)') trim(w90_readwrite_get_smearing_type(pw90_boltzwann%dos_smearing%type_index))// &
            ", non-adaptive, width (eV) =", pw90_boltzwann%dos_smearing%fixed_width
        end if
      end if
    end if
    if (pw90_boltzwann%calc_also_dos .and. pw90_boltzwann%dos_smearing%use_adaptive .and. &
        (pw90_boltzwann%dos_smearing%fixed_width .ne. 0._dp) .and. on_root) then
      write (stdout, '(5X,A)') "*** WARNING! pw90_boltzwann_dos_smr_fixed_en_width ignored since you chose"
      write (stdout, '(5X,A)') "             an adaptive smearing."
    end if

    if (on_root .and. (print_output%iprint > 1)) then
      if (pw90_boltzwann%tdf_smearing%fixed_width/(TDFEnergyArray(2) - TDFEnergyArray(1)) &
          < min_smearing_binwidth_ratio) then
        write (stdout, '(5X,A)') "Smearing for TDF: "
        write (stdout, '(7X,A)') "Unsmeared (use smearing width larger than bin width to smear)"
      else
        write (stdout, '(5X,A)') "Smearing for TDF: "
        write (stdout, '(7X,A,G18.10)') &
          trim(w90_readwrite_get_smearing_type(pw90_boltzwann%tdf_smearing%type_index))//", non-adaptive, width (eV) =", &
          pw90_boltzwann%tdf_smearing%fixed_width
      end if
    end if

    if (on_root) then
      write (stdout, '(5X,A,I0,A,I0,A,I0)') "k-grid used for band interpolation in BoltzWann: ", &
        pw90_boltzwann%kmesh%mesh(1), 'x', pw90_boltzwann%kmesh%mesh(2), 'x', pw90_boltzwann%kmesh%mesh(3)
      write (stdout, '(5X,A,I1)') "Number of electrons per state: ", num_elec_per_state
      write (stdout, '(5X,A,G18.10)') "Relaxation time (fs): ", pw90_boltzwann%relax_time
      if (print_output%iprint > 1) then
        write (stdout, '(5X,A,G18.10)') "Energy step for TDF (eV): ", pw90_boltzwann%tdf_energy_step
      end if
    end if
    kweight = 1.0_dp/real(PRODUCT(pw90_boltzwann%kmesh%mesh), kind=dp)

    if (pw90_boltzwann%bandshift .and. on_root) then
      write (stdout, '(5X,A,I0,A,G18.10,A)') "Shifting energy bands with index >= ", pw90_boltzwann%bandshift_firstband, " by ", &
        pw90_boltzwann%bandshift_energyshift, " eV."
    end if

    call utility_recip_lattice_base(real_lattice, recip_lattice, volume)
    NumPtsRefined = 0
    min_spacing = 1.e10_dp ! very large initial value
    max_spacing = 0.e0_dp
    ! I loop over all kpoints
    do loop_tot = my_node_id, PRODUCT(pw90_boltzwann%kmesh%mesh) - 1, num_nodes

      ! I get the coordinates for the x,y,z components starting from a single loop variable
      ! (which is better for parallelization purposes)
      ! Important! This works only if loop_tot starts from ZERO and ends with
      !            PRODUCT(pw90_boltzwann_kmesh)-1, so be careful when parallelizing
      loop_x = loop_tot/(pw90_boltzwann%kmesh%mesh(2)*pw90_boltzwann%kmesh%mesh(3))
      loop_y = (loop_tot - loop_x*(pw90_boltzwann%kmesh%mesh(2)*pw90_boltzwann%kmesh%mesh(3)))/pw90_boltzwann%kmesh%mesh(3)
      loop_z = loop_tot - loop_x*(pw90_boltzwann%kmesh%mesh(2)*pw90_boltzwann%kmesh%mesh(3)) - loop_y*pw90_boltzwann%kmesh%mesh(3)

      ! kpt(i) is in in the [0,d-1]/d range, with d=pw90_boltzwann_kmesh(i)
      kpt(1) = (real(loop_x, dp)/real(pw90_boltzwann%kmesh%mesh(1), dp))
      kpt(2) = (real(loop_y, dp)/real(pw90_boltzwann%kmesh%mesh(2), dp))
      kpt(3) = (real(loop_z, dp)/real(pw90_boltzwann%kmesh%mesh(3), dp))

      ! Here I get the band energies and the velocities
      call wham_get_eig_deleig(dis_manifold, kpt_latt, pw90_band_deriv_degen, ws_region, print_output, wannier_data, &
                               ws_distance, wigner_seitz, delHH, HH, HH_R, u_matrix, UU, v_matrix, &
                               del_eig, eig, eigval, kpt, real_lattice, scissors_shift, mp_grid, &
                               num_bands, num_kpts, num_wann, num_valence_bands, effective_model, &
                               have_disentangled, seedname, stdout, timer, error, comm)
      if (allocated(error)) return

      call dos_get_levelspacing(del_eig, pw90_boltzwann%kmesh%mesh, levelspacing_k, num_wann, recip_lattice)

      ! Here I apply a scissor operator to the conduction bands, if required in the input
      if (pw90_boltzwann%bandshift) then
        eig(pw90_boltzwann%bandshift_firstband:) = eig(pw90_boltzwann%bandshift_firstband:) + pw90_boltzwann%bandshift_energyshift
      end if

      call TDF_kpt(pw90_boltzwann, ws_region, pw90_spin, wannier_data, ws_distance, wigner_seitz, &
                   HH_R, SS_R, del_eig, eig, TDFEnergyArray, kpt, real_lattice, TDF_k, mp_grid, &
                   num_wann, num_elec_per_state, spin_decomp, error, comm)
      if (allocated(error)) return

      ! As above, the sum of TDF_k * kweight amounts to calculate
      ! spin_degeneracy * V_cell/(2*pi)^3 * \int_BZ d^3k
      ! so that we divide by the cell_volume (in Angstrom^3) to have
      ! the correct integral
      TDF = TDF + TDF_k*kweight/cell_volume

      ! DOS part !

      if (pw90_boltzwann%calc_also_dos) then
        if (pw90_boltzwann%dos_smearing%use_adaptive) then

          ! This may happen if at least one band has zero derivative (along all three directions)
          ! Then I substitute this point with its 8 neighbors (+/- 1/4 of the spacing with the next point on the grid
          ! on each of the three directions)
          min_spacing = min(min_spacing, minval(abs(levelspacing_k)))
          max_spacing = max(max_spacing, maxval(abs(levelspacing_k)))
          if (any(abs(levelspacing_k) < SPACING_THRESHOLD)) then
            orig_kpt = kpt
            NumPtsRefined = NumPtsRefined + 1
            do i = -1, 1, 2
              do j = -1, 1, 2
                do k = -1, 1, 2
                  kpt = orig_kpt + &
                        (/real(i, kind=dp)/real(pw90_boltzwann%kmesh%mesh(1), dp)/4._dp, &
                          real(j, kind=dp)/real(pw90_boltzwann%kmesh%mesh(2), dp)/4._dp, &
                          real(k, kind=dp)/real(pw90_boltzwann%kmesh%mesh(3), dp)/4._dp/)
                  call wham_get_eig_deleig(dis_manifold, kpt_latt, pw90_band_deriv_degen, &
                                           ws_region, print_output, wannier_data, ws_distance, &
                                           wigner_seitz, delHH, HH, HH_R, u_matrix, UU, v_matrix, &
                                           del_eig, eig, eigval, kpt, real_lattice, &
                                           scissors_shift, mp_grid, num_bands, num_kpts, num_wann, &
                                           num_valence_bands, effective_model, have_disentangled, &
                                           seedname, stdout, timer, error, comm)
                  if (allocated(error)) return

                  call dos_get_levelspacing(del_eig, pw90_boltzwann%kmesh%mesh, levelspacing_k, &
                                            num_wann, recip_lattice)
                  call dos_get_k(num_elec_per_state, ws_region, kpt, DOS_EnergyArray, eig, dos_k, &
                                 num_wann, wannier_data, real_lattice, mp_grid, pw90_dos, &
                                 spin_decomp, pw90_spin, ws_distance, wigner_seitz, HH_R, SS_R, &
                                 pw90_boltzwann%dos_smearing, error, comm, &
                                 levelspacing_k=levelspacing_k)
                  if (allocated(error)) return

                  ! I divide by 8 because I'm substituting a point with its 8 neighbors
                  dos_all = dos_all + dos_k*kweight/8.
                end do
              end do
            end do
          else
            call dos_get_k(num_elec_per_state, ws_region, kpt, DOS_EnergyArray, eig, dos_k, &
                           num_wann, wannier_data, real_lattice, mp_grid, pw90_dos, spin_decomp, &
                           pw90_spin, ws_distance, wigner_seitz, HH_R, SS_R, &
                           pw90_boltzwann%dos_smearing, error, comm, levelspacing_k=levelspacing_k)
            if (allocated(error)) return

            dos_all = dos_all + dos_k*kweight
          end if
        else
          call dos_get_k(num_elec_per_state, ws_region, kpt, DOS_EnergyArray, eig, dos_k, &
                         num_wann, wannier_data, real_lattice, mp_grid, pw90_dos, spin_decomp, &
                         pw90_spin, ws_distance, wigner_seitz, HH_R, SS_R, &
                         pw90_boltzwann%dos_smearing, error, comm)
          if (allocated(error)) return

          ! This sum multiplied by kweight amounts to calculate
          ! spin_degeneracy * V_cell/(2*pi)^3 * \int_BZ d^3k
          ! So that the DOS will be in units of 1/eV, normalized so that
          ! \int_{-\infty}^{\infty} DOS(E) dE = Num.Electrons
          dos_all = dos_all + dos_k*kweight
        end if
      end if

    end do

    ! I sum the results of the calculation for the DOS on the root node only
    ! (I have to print the results only)
    if (pw90_boltzwann%calc_also_dos) then
      call comms_reduce(DOS_all(1, 1), size(DOS_all), 'SUM', error, comm)
      if (allocated(error)) return
      call comms_reduce(NumPtsRefined, 1, 'SUM', error, comm)
      if (allocated(error)) return
      call comms_reduce(min_spacing, 1, 'MIN', error, comm)
      if (allocated(error)) return
      call comms_reduce(max_spacing, 1, 'MAX', error, comm)
      if (allocated(error)) return
    end if
    ! I sum the results of the calculation on all nodes, and I store them on all
    ! nodes (because for the following, each node will do a different calculation,
    ! each of which will require the whole knowledge of the TDF array)
    call comms_allreduce(TDF(1, 1, 1), size(TDF), 'SUM', error, comm)
    if (allocated(error)) return

    if (pw90_boltzwann%calc_also_dos .and. on_root) then
      write (boltzdos_unit, '(A)') "# Written by the BoltzWann module of the Wannier90 code."
      write (boltzdos_unit, '(A)') "# The first column."
      if (pw90_boltzwann%dos_smearing%use_adaptive) then
        write (boltzdos_unit, '(A)') '# The second column is the adaptively-smeared DOS'
        write (boltzdos_unit, '(A)') '# (see Yates et al., PRB 75, 195121 (2007)'
        if (spin_decomp) then
          write (boltzdos_unit, '(A)') '# The third column is the spin-up projection of the DOS'
          write (boltzdos_unit, '(A)') '# The fourth column is the spin-down projection of the DOS'
        end if
        write (boltzdos_unit, '(A,1X,G14.6)') '# Smearing coefficient: ', pw90_boltzwann%dos_smearing%adaptive_prefactor
        write (boltzdos_unit, '(A,I0,A,I0)') '# Number of points refined: ', NumPtsRefined, &
          ' out of ', product(pw90_boltzwann%kmesh%mesh)
        write (boltzdos_unit, '(A,G18.10,A,G18.10,A)') '# (Min spacing: ', min_spacing, &
          ', max spacing: ', max_spacing, ')'
      else
        if (pw90_boltzwann%dos_smearing%fixed_width/(DOS_EnergyArray(2) - DOS_EnergyArray(1)) < min_smearing_binwidth_ratio) then
          write (boltzdos_unit, '(A)') '# The second column is the unsmeared DOS.'
        else
          write (boltzdos_unit, '(A,G14.6,A)') '# The second column is the DOS for a fixed smearing of ', &
            pw90_boltzwann%dos_smearing%fixed_width, ' eV.'
        end if
      end if
      write (boltzdos_unit, '(A,1X,G14.6)') '# Cell volume (ang^3): ', cell_volume

      write (boltzdos_unit, '(A)') '# Energy(eV) DOS [DOS DOS ...]'

      ! I save a string with the number of fields to print
      write (numfieldsstr, '(I0)') 1 + ndim
      do EnIdx = 1, size(DOS_EnergyArray)
        write (boltzdos_unit, '(1X,'//trim(numfieldsstr)//'G18.10)') &
          DOS_EnergyArray(EnIdx), dos_all(EnIdx, :)
      end do
    end if

    if (on_root .and. (print_output%timing_level > 0)) call io_stopwatch_stop('calcTDF', timer)
    if (on_root) then
      if (pw90_boltzwann%calc_also_dos) then
        write (stdout, '(3X,A)') "TDF and DOS calculated."
      else
        write (stdout, '(3X,A)') "TDF calculated."
      end if
    end if
    if (on_root) write (stdout, *)

    if (on_root .and. pw90_boltzwann%calc_also_dos) then
      close (boltzdos_unit)
      if (print_output%iprint > 1) write (stdout, '(3X,A)') "DOS written on the "//trim(seedname)//"_boltzdos.dat file."
    end if

    deallocate (HH, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating HH in calcTDF', comm)
      return
    endif
    deallocate (delHH, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating delHH in calcTDF', comm)
      return
    endif
    deallocate (UU, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating UU in calcTDF', comm)
      return
    endif
    deallocate (DOS_EnergyArray, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating DOS_EnergyArray in calcTDF', comm)
      return
    endif
    deallocate (DOS_k, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating DOS_k in calcTDF', comm)
      return
    endif
    deallocate (DOS_all, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating DOS_all in calcTDF', comm)
      return
    endif
    deallocate (TDF_k, stat=ierr)
    if (ierr /= 0) then
      call set_error_dealloc(error, 'Error in deallocating TDF_k in calcTDF', comm)
      return
    endif

  end subroutine calcTDFandDOS

  !================================================!
  function MinusFermiDerivative(E, mu, KT)
    !================================================!
    !> This function calculates -dn(E)/dE, where n(E) is the Fermi distribution function.
    !>
    !> \param E  Energy at which we want to calculate -dn(E)/dE, in 1/eV
    !> \param mu Chemical potential in eV
    !> \param KT k_Boltzmann * Temperature, in eV
    !================================================!

    real(kind=dp), intent(in) :: E
    real(kind=dp), intent(in) :: mu
    real(kind=dp), intent(in) :: KT
    real(kind=dp) :: MinusFermiDerivative

    ! I do not put stopwatches here because it would slow down the calculation by orders of magnitude

    ! MaxExp is the maximum value to be used for the exp function.
    ! The function is taken to be zero for x>MaxExp. This value is chosen so that
    ! the function is truncated when its value is smaller than about 1.e-16
    real(kind=dp), parameter :: MaxExp = 36._dp
    real(kind=dp) :: MyExp

    MyExp = (E - mu)/KT
    if (abs(MyExp) > MaxExp) then
      MinusFermiDerivative = 0._dp
    else
      MinusFermiDerivative = 1._dp/KT*exp(MyExp)/((exp(MyExp) + 1._dp)**2)
    end if

  end function MinusFermiDerivative

  !================================================!
  subroutine TDF_kpt(pw90_boltzwann, ws_region, pw90_spin, wannier_data, ws_distance, &
                     wigner_seitz, HH_R, SS_R, deleig_k, eig_k, EnergyArray, kpt, real_lattice, &
                     TDF_k, mp_grid, num_wann, num_elec_per_state, spin_decomp, error, comm)
    !================================================!
    !! This subroutine calculates the contribution to the TDF of a single k point
    !!
    !!  This routine does not use the adaptive smearing; in fact, for non-zero temperatures
    !!       one often doesn't even need to smear. It simply uses a standard smearing as defined by
    !!       the variables pw90_boltzwann_TDF_smr_fixed_en_width and pw90_boltzwann_TDF_smr_index
    !!
    !! still to do: adapt spin_get_nk to read in input the UU rotation matrix
    !!
    !! This routine simply provides the dos contribution of a given
    !!       point. This must be externally summed after proper weighting.
    !!       The weight factor (for a full BZ sampling with N^3 points) is 1/N^3/cell_volume
    !!       if we want to calculate 1/(2*pi)^3 * \int_{BZ} d^3 k
    !! The only factor that is included INSIDE this routine is the spin degeneracy
    !!       factor (=2 if spinors is .false., =1 if spinors is .true.)
    !! The EnergyArray is assumed to be evenly spaced (and the energy spacing
    !!       is taken from EnergyArray(2)-EnergyArray(1))
    !! The routine is assuming that EnergyArray has at least two elements.
    !!  The meaning of the three indices of the TDF_k array is different with respect to
    !!       those of the dos_k array returned by the dos_get_k routine
    !! The TDF_k array must have dimensions 6 * size(EnergyArray) * ndim, where
    !!       ndim=1 if spin_decomp==false, or ndim=3 if spin_decomp==true. This is not checked.
    !!
    !================================================!

    use w90_constants, only: dp, smearing_cutoff, min_smearing_binwidth_ratio
    use w90_utility, only: utility_w0gauss
    use w90_types, only: print_output_type, wannier_data_type, ws_region_type, ws_distance_type
    use w90_postw90_types, only: pw90_boltzwann_type, pw90_spin_mod_type, wigner_seitz_type
    use w90_spin, only: spin_get_nk
    use w90_utility, only: utility_w0gauss
    use w90_comms, only: w90_comm_type

    implicit none

    ! arguments
    type(pw90_boltzwann_type), intent(in) :: pw90_boltzwann
    type(ws_region_type), intent(in) :: ws_region
    type(pw90_spin_mod_type), intent(in) :: pw90_spin
    type(wannier_data_type), intent(in) :: wannier_data
    type(ws_distance_type), intent(inout) :: ws_distance
    type(wigner_seitz_type), intent(in) :: wigner_seitz
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm

    integer, intent(in) :: num_wann
    integer, intent(in) :: mp_grid(3)

    real(kind=dp), intent(in) :: kpt(3)
    !! the three coordinates of the k point vector whose DOS contribution we
    !! want to calculate (in relative coordinates)
    real(kind=dp), intent(in) :: EnergyArray(:)
    !! array with the energy grid on which to calculate the DOS (in eV)
    !! It must have at least two elements
    real(kind=dp), intent(in) :: eig_k(:)
    !! array with the eigenvalues at the given k point (in eV)
    real(kind=dp), intent(in) :: deleig_k(:, :)
    !! array with the band derivatives at the given k point
    !! (in eV * angstrom / (2pi) as internally given by the code)
    !! already corrected in case of degeneracies, as returned by the
    !! wham_get_deleig_a routine
    real(kind=dp), intent(out) :: TDF_k(:, :, :)
    !! TDF_k array in which the contribution is stored. Three dimensions:
    !! TDF_k(ij, energyidx, spinidx), where:
    !! - ij indexes the components of the TDF (symmetric) tensor (1=XX, 2=XY, ...);
    !!   see the global constants defined in the module
    !! - energyidx is the index of the energies, corresponding to the one
    !!   of the EnergyArray array;
    !!  - spinidx=1 contains the total dos; if if spin_decomp==.true., then
    !!  spinidx=2 and spinidx=3 contain the spin-up and spin-down contributions to the DOS
    real(kind=dp), intent(in) :: real_lattice(3, 3)

    complex(kind=dp), allocatable, intent(inout) :: HH_R(:, :, :) !  <0n|r|Rm>
    complex(kind=dp), allocatable, intent(inout) :: SS_R(:, :, :, :) ! <0n|sigma_x,y,z|Rm>

    logical, intent(in) :: spin_decomp
    integer, intent(in) :: num_elec_per_state

    ! local variables
    real(kind=dp) :: smear, arg ! Adaptive smearing
    real(kind=dp) :: rdum, spn_nk(num_wann), alpha_sq, beta_sq
    real(kind=dp) :: binwidth, r_num_elec_per_state
    integer :: BandIdx, loop_f, min_f, max_f
    logical :: DoSmearing

    r_num_elec_per_state = real(num_elec_per_state, kind=dp)

    ! Get spin projections for every band
    !
    if (spin_decomp) then
      call spin_get_nk(ws_region, pw90_spin, wannier_data, ws_distance, wigner_seitz, HH_R, SS_R, &
                       kpt, real_lattice, spn_nk, mp_grid, num_wann, error, comm)
      if (allocated(error)) return
    endif

    binwidth = EnergyArray(2) - EnergyArray(1)

    TDF_k = 0.0_dp
    do BandIdx = 1, num_wann
      if (spin_decomp) then
        ! Contribution to spin-up DOS of Bloch spinor with component
        ! (alpha,beta) with respect to the chosen quantization axis
        alpha_sq = (1.0_dp + spn_nk(BandIdx))/2.0_dp ! |alpha|^2
        ! Contribution to spin-down DOS
        beta_sq = 1.0_dp - alpha_sq ! |beta|^2 = 1 - |alpha|^2
      end if

      ! Do not use an adaptive smearing here, it would require the knowledge of second derivatives
      ! And typically, when working at not too small temperatures, smearing is not needed

      ! Faster optimization: I precalculate the indices
      ! Value of the smearing in eV; default = 0 eV, i.e. no smearing
      smear = pw90_boltzwann%tdf_smearing%fixed_width
      if (smear/binwidth < min_smearing_binwidth_ratio) then
        min_f = max(nint((eig_k(BandIdx) - EnergyArray(1))/ &
                         (EnergyArray(size(EnergyArray)) - EnergyArray(1)) &
                         *real(size(EnergyArray) - 1, kind=dp)) + 1, 1)
        max_f = min(nint((eig_k(BandIdx) - EnergyArray(1))/ &
                         (EnergyArray(size(EnergyArray)) - EnergyArray(1)) &
                         *real(size(EnergyArray) - 1, kind=dp)) + 1, size(EnergyArray))
        DoSmearing = .false.
      else
        min_f = max(nint((eig_k(BandIdx) - smearing_cutoff*smear - EnergyArray(1))/ &
                         (EnergyArray(size(EnergyArray)) - EnergyArray(1)) &
                         *real(size(EnergyArray) - 1, kind=dp)) + 1, 1)
        max_f = min(nint((eig_k(BandIdx) + smearing_cutoff*smear - EnergyArray(1))/ &
                         (EnergyArray(size(EnergyArray)) - EnergyArray(1)) &
                         *real(size(EnergyArray) - 1, kind=dp)) + 1, size(EnergyArray))
        DoSmearing = .true.
      end if

      do loop_f = min_f, max_f
        if (DoSmearing) then
          arg = (EnergyArray(loop_f) - eig_k(BandIdx))/smear
          rdum = utility_w0gauss(arg, pw90_boltzwann%tdf_smearing%type_index, error, comm)/smear
          if (allocated(error)) return
        else
          rdum = 1._dp/(EnergyArray(2) - EnergyArray(1))
        end if
        !
        ! Contribution to total DOS
        !
        TDF_k(XX, loop_f, 1) = TDF_k(XX, loop_f, 1) + rdum* &
                               r_num_elec_per_state*deleig_k(BandIdx, 1)*deleig_k(BandIdx, 1)
        TDF_k(XY, loop_f, 1) = TDF_k(XY, loop_f, 1) + rdum* &
                               r_num_elec_per_state*deleig_k(BandIdx, 1)*deleig_k(BandIdx, 2)
        TDF_k(YY, loop_f, 1) = TDF_k(YY, loop_f, 1) + rdum* &
                               r_num_elec_per_state*deleig_k(BandIdx, 2)*deleig_k(BandIdx, 2)
        TDF_k(XZ, loop_f, 1) = TDF_k(XZ, loop_f, 1) + rdum* &
                               r_num_elec_per_state*deleig_k(BandIdx, 1)*deleig_k(BandIdx, 3)
        TDF_k(YZ, loop_f, 1) = TDF_k(YZ, loop_f, 1) + rdum* &
                               r_num_elec_per_state*deleig_k(BandIdx, 2)*deleig_k(BandIdx, 3)
        TDF_k(ZZ, loop_f, 1) = TDF_k(ZZ, loop_f, 1) + rdum* &
                               r_num_elec_per_state*deleig_k(BandIdx, 3)*deleig_k(BandIdx, 3)

        ! I don't put num_elec_per_state here below: if we are calculating the spin decomposition,
        ! we should be doing a calcultation with spin-orbit, and thus num_elec_per_state=1!
        if (spin_decomp) then
          ! Spin-up contribution
          TDF_k(XX, loop_f, 2) = TDF_k(XX, loop_f, 2) + rdum* &
                                 alpha_sq*deleig_k(BandIdx, 1)*deleig_k(BandIdx, 1)
          TDF_k(XY, loop_f, 2) = TDF_k(XY, loop_f, 2) + rdum* &
                                 alpha_sq*deleig_k(BandIdx, 1)*deleig_k(BandIdx, 2)
          TDF_k(YY, loop_f, 2) = TDF_k(YY, loop_f, 2) + rdum* &
                                 alpha_sq*deleig_k(BandIdx, 2)*deleig_k(BandIdx, 2)
          TDF_k(XZ, loop_f, 2) = TDF_k(XZ, loop_f, 2) + rdum* &
                                 alpha_sq*deleig_k(BandIdx, 1)*deleig_k(BandIdx, 3)
          TDF_k(YZ, loop_f, 2) = TDF_k(YZ, loop_f, 2) + rdum* &
                                 alpha_sq*deleig_k(BandIdx, 2)*deleig_k(BandIdx, 3)
          TDF_k(ZZ, loop_f, 2) = TDF_k(ZZ, loop_f, 2) + rdum* &
                                 alpha_sq*deleig_k(BandIdx, 3)*deleig_k(BandIdx, 3)
          ! Spin-down contribution
          TDF_k(XX, loop_f, 3) = TDF_k(XX, loop_f, 3) + rdum* &
                                 beta_sq*deleig_k(BandIdx, 1)*deleig_k(BandIdx, 1)
          TDF_k(XY, loop_f, 3) = TDF_k(XY, loop_f, 3) + rdum* &
                                 beta_sq*deleig_k(BandIdx, 1)*deleig_k(BandIdx, 2)
          TDF_k(YY, loop_f, 3) = TDF_k(YY, loop_f, 3) + rdum* &
                                 beta_sq*deleig_k(BandIdx, 2)*deleig_k(BandIdx, 2)
          TDF_k(XZ, loop_f, 3) = TDF_k(XZ, loop_f, 3) + rdum* &
                                 beta_sq*deleig_k(BandIdx, 1)*deleig_k(BandIdx, 3)
          TDF_k(YZ, loop_f, 3) = TDF_k(YZ, loop_f, 3) + rdum* &
                                 beta_sq*deleig_k(BandIdx, 2)*deleig_k(BandIdx, 3)
          TDF_k(ZZ, loop_f, 3) = TDF_k(ZZ, loop_f, 3) + rdum* &
                                 beta_sq*deleig_k(BandIdx, 3)*deleig_k(BandIdx, 3)
        end if
      end do
    end do !loop over bands

    ! I multiply it here, since I am assuming a constant relaxation time, independent of the band index
    ! (actually, it is also independent of k)
    TDF_k = TDF_k*pw90_boltzwann%relax_time

  end subroutine TDF_kpt

end module w90_boltzwann
