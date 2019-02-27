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

module w90_boltzwann
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
  !============================================================!
  use w90_constants
  use w90_parameters, only: &
    boltz_mu_min, boltz_mu_max, boltz_mu_step, boltz_temp_min, boltz_temp_max, boltz_temp_step, &
    boltz_kmesh_spacing, boltz_kmesh, boltz_tdf_energy_step, boltz_relax_time, &
    boltz_bandshift, boltz_bandshift_firstband, boltz_bandshift_energyshift, &
    timing_level, dis_win_min, dis_win_max, spin_decomp, boltz_dos_adpt_smr, &
    boltz_dos_adpt_smr_fac, boltz_dos_smr_fixed_en_width, boltz_2d_dir_num, &
    boltz_tdf_smr_fixed_en_width, cell_volume, num_elec_per_state, iprint
  use w90_io, only: io_error, stdout, io_stopwatch, io_file_unit, seedname
  use w90_utility, only: utility_inv3, utility_inv2
  use w90_postw90_common
  use w90_comms
  use w90_dos, only: dos_get_k, dos_get_levelspacing
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

  subroutine boltzwann_main()
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
    integer :: TempNumPoints, MuNumPoints, TDFEnergyNumPoints
    integer :: i, j, ierr, EnIdx, TempIdx, MuIdx
    real(kind=dp), dimension(:), allocatable :: TempArray, MuArray, KTArray
    real(kind=dp), dimension(:, :, :), allocatable :: TDF ! (coordinate,Energy)
    real(kind=dp), dimension(:), allocatable :: TDFEnergyArray
    real(kind=dp), dimension(:, :), allocatable :: IntegrandArray ! (coordinate, Energy) at a given T and mu
    real(kind=dp), dimension(3, 3)                  :: SigmaS_FP, ThisElCond, ElCondInverse, ThisSeebeck
    real(kind=dp), dimension(2, 2)                  :: ThisElCond2d, ElCondInverse2d
    real(kind=dp), dimension(6)                    :: ElCondTimesSeebeck
    real(kind=dp), dimension(:, :, :), allocatable :: ElCond ! (coordinate,Temp, mu)
    real(kind=dp), dimension(:, :, :), allocatable :: SigmaS ! (coordinate,Temp, mu)
    real(kind=dp), dimension(:, :, :), allocatable :: Seebeck ! (coordinate,Temp, mu)
    real(kind=dp), dimension(:, :, :), allocatable :: Kappa ! (coordinate,Temp, mu)
    real(kind=dp), dimension(:, :), allocatable :: LocalElCond ! (coordinate,Temp+mu combined index)
    real(kind=dp), dimension(:, :), allocatable :: LocalSigmaS ! (coordinate,Temp+mu combined index)
    real(kind=dp), dimension(:, :), allocatable :: LocalSeebeck ! (coordinate,Temp+mu combined index)
    real(kind=dp), dimension(:, :), allocatable :: LocalKappa ! (coordinate,Temp+mu combined index)
    real(kind=dp)                                  :: Determinant
    integer :: tdf_unit, elcond_unit, sigmas_unit, seebeck_unit, kappa_unit, ndim

    ! Needed to split an array on different nodes
    integer, dimension(0:num_nodes - 1) :: counts
    integer, dimension(0:num_nodes - 1) :: displs
    integer :: LocalIdx, GlobalIdx
    ! I also add 3 times the smearing on each side of the TDF energy array to take into account also possible smearing effects
    real(kind=dp), parameter :: TDF_exceeding_energy_times_smr = 3._dp
    real(kind=dp) :: TDF_exceeding_energy
    integer :: NumberZeroDet

    if (on_root .and. (timing_level > 0)) call io_stopwatch('boltzwann_main', 1)

    if (on_root) then
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

    if (on_root) then
      if (boltz_2d_dir_num /= 0) then
        write (stdout, '(1x,a)') '>                                                                           <'
        write (stdout, '(1x,a)') '> NOTE! Using the 2D version for the calculation of the Seebeck             <'
        if (boltz_2d_dir_num == 1) then
          write (stdout, '(1x,a)') '>       coefficient, where the non-periodic direction is x.                 <'
        elseif (boltz_2d_dir_num == 2) then
          write (stdout, '(1x,a)') '>       coefficient, where the non-periodic direction is y.                 <'
        elseif (boltz_2d_dir_num == 3) then
          write (stdout, '(1x,a)') '>       coefficient, where the non-periodic direction is z.                 <'
        else
          call io_error('Unrecognized value of boltz_2d_dir_num')
        end if
        write (stdout, '(1x,a)') '>                                                                           <'
        write (stdout, '(1x,a)') ''
      end if
    end if

    ! I precalculate the TempArray and the MuArray
    TempNumPoints = int(floor((boltz_temp_max - boltz_temp_min)/boltz_temp_step)) + 1
    allocate (TempArray(TempNumPoints), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating TempArray in boltzwann_main')
    do i = 1, TempNumPoints
      TempArray(i) = boltz_temp_min + real(i - 1, dp)*boltz_temp_step
    end do

    ! This array contains the same temperatures of the TempArray, but multiplied by k_boltzmann, in units of eV
    allocate (KTArray(TempNumPoints), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating KTArray in boltzwann_main')
    ! (k_B in eV/kelvin is equal to k_B_SI / elem_charge_SI)
    KTArray = TempArray*k_B_SI/elem_charge_SI

    MuNumPoints = int(floor((boltz_mu_max - boltz_mu_min)/boltz_mu_step)) + 1
    allocate (MuArray(MuNumPoints), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating MuArray in boltzwann_main')
    do i = 1, MuNumPoints
      MuArray(i) = boltz_mu_min + real(i - 1, dp)*boltz_mu_step
    end do

    ! I precalculate the TDFEnergyArray
    ! I assume that dis_win_min and dis_win_max are set to sensible values, related to the max and min energy
    ! This is true if the .eig file is present. I can assume its presence since we need it to interpolate the
    ! bands.
    ! I also add 3 times the smearing on each side of the TDF energy array to take into account also possible smearing effects,
    ! or at least 0.2 eV
    TDF_exceeding_energy = max(TDF_exceeding_energy_times_smr*boltz_TDF_smr_fixed_en_width, 0.2_dp)
    TDFEnergyNumPoints = int(floor((dis_win_max - dis_win_min + 2._dp*TDF_exceeding_energy)/boltz_tdf_energy_step)) + 1
    if (TDFEnergyNumPoints .eq. 1) TDFEnergyNumPoints = 2
    allocate (TDFEnergyArray(TDFEnergyNumPoints), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating TDFEnergyArray in boltzwann_main')
    do i = 1, TDFEnergyNumPoints
      TDFEnergyArray(i) = dis_win_min - TDF_exceeding_energy + real(i - 1, dp)*boltz_tdf_energy_step
    end do

    if (spin_decomp) then
      ndim = 3
    else
      ndim = 1
    end if

    ! I allocate the array for the TDF
    allocate (TDF(6, TDFEnergyNumPoints, ndim), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating TDF in boltzwann_main')

    ! I call the subroutine that calculates the Transport Distribution Function
    call calcTDFandDOS(TDF, TDFEnergyArray)
    ! The TDF array contains now the TDF, or more precisely
    ! hbar^2 * TDF in units of eV * fs / angstrom

    ! I print on file the TDF
    if (on_root) then
      tdf_unit = io_file_unit()
      open (unit=tdf_unit, file=trim(seedname)//'_tdf.dat')
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
      if (iprint > 1) write (stdout, '(3X,A)') "Transport distribution function written on the "//trim(seedname)//"_tdf.dat file."
    end if

    ! *********************************************************************************
    ! I got the TDF and I printed it. Now I use it to calculate the transport properties.

    if (on_root .and. (timing_level > 0)) call io_stopwatch('boltzwann_main: calc_props', 1)

    ! I obtain the counts and displs arrays, which tell how I should partition a big array
    ! on the different nodes.
    call comms_array_split(TempNumPoints*MuNumPoints, counts, displs)

    ! I allocate the arrays for the spectra
    ! Allocate at least 1 entry
    allocate (LocalElCond(6, max(1, counts(my_node_id))), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating LocalElCond in boltzwann_main')
    allocate (LocalSigmaS(6, max(1, counts(my_node_id))), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating LocalSigmaS in boltzwann_main')
    allocate (LocalSeebeck(9, max(1, counts(my_node_id))), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating LocalSeebeck in boltzwann_main')
    allocate (LocalKappa(6, max(1, counts(my_node_id))), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating LocalKappa in boltzwann_main')
    LocalElCond = 0._dp
    LocalSeebeck = 0._dp
    LocalKappa = 0._dp

    ! I allocate the array that I will use to store the functions to be integrated
    allocate (IntegrandArray(6, TDFEnergyNumPoints), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating FermiDerivArray in boltzwann_main')

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
      LocalElCond(:, LocalIdx) = sum(IntegrandArray, DIM=2)*boltz_tdf_energy_step
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
      if (boltz_2d_dir_num /= 0) then
        ! Invert only the appropriate 2x2 submatrix
        if (boltz_2d_dir_num == 1) then
          ThisElCond2d(1, 1) = ThisElCond(2, 2)
          ThisElCond2d(1, 2) = ThisElCond(2, 3)
          ThisElCond2d(2, 1) = ThisElCond(3, 2)
          ThisElCond2d(2, 2) = ThisElCond(3, 3)
        elseif (boltz_2d_dir_num == 2) then
          ThisElCond2d(1, 1) = ThisElCond(1, 1)
          ThisElCond2d(1, 2) = ThisElCond(1, 3)
          ThisElCond2d(2, 1) = ThisElCond(3, 1)
          ThisElCond2d(2, 2) = ThisElCond(3, 3)
        elseif (boltz_2d_dir_num == 3) then
          ThisElCond2d(1, 1) = ThisElCond(1, 1)
          ThisElCond2d(1, 2) = ThisElCond(1, 2)
          ThisElCond2d(2, 1) = ThisElCond(2, 1)
          ThisElCond2d(2, 2) = ThisElCond(2, 2)
          ! I do not do the else case because this was already checked at the beginning
          ! of this routine
        end if
        call utility_inv2(ThisElCond2d, ElCondInverse2d, Determinant)
        ElCondInverse = 0._dp ! Other elements must be set to zero
        if (boltz_2d_dir_num == 1) then
          ElCondInverse(2, 2) = ElCondInverse2d(1, 1)
          ElCondInverse(2, 3) = ElCondInverse2d(1, 2)
          ElCondInverse(3, 2) = ElCondInverse2d(2, 1)
          ElCondInverse(3, 3) = ElCondInverse2d(2, 2)
        elseif (boltz_2d_dir_num == 2) then
          ElCondInverse(1, 1) = ElCondInverse2d(1, 1)
          ElCondInverse(1, 3) = ElCondInverse2d(1, 2)
          ElCondInverse(3, 1) = ElCondInverse2d(2, 1)
          ElCondInverse(3, 3) = ElCondInverse2d(2, 2)
        elseif (boltz_2d_dir_num == 3) then
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
      LocalSigmaS(:, LocalIdx) = sum(IntegrandArray, DIM=2)*boltz_tdf_energy_step/TempArray(TempIdx)
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
      LocalKappa(:, LocalIdx) = sum(IntegrandArray, DIM=2)*boltz_tdf_energy_step/TempArray(TempIdx)
      ! Kappa contains now the thermal conductivity in units of
      ! 1/hbar^2 * eV^3*fs/angstrom/kelvin
    end do
    ! I check if there were (mu,T) pairs for which we got sigma = 0
    call comms_reduce(NumberZeroDet, 1, 'SUM')
    if (on_root) then
      if ((NumberZeroDet .gt. 0)) then
        write (stdout, '(1X,A,I0,A)') "> WARNING! There are ", NumberZeroDet, " (mu,T) pairs for which the electrical"
        write (stdout, '(1X,A)') ">          conductivity has zero determinant."
        write (stdout, '(1X,A)') ">          Seebeck coefficient set to zero for those pairs."
        write (stdout, '(1X,A)') ">          Check if this is physical or not."
        write (stdout, '(1X,A)') ">          (If you are dealing with a 2D system, set the boltz_2d_dir flag.)"
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
    LocalElCond = LocalElCond*elem_charge_SI**3/(hbar_SI**2)*1.e-5_dp
    ! THIS IS NOW THE ELECTRICAL CONDUCTIVITY IN SI UNITS, i.e. in 1/Ohm/meter

    ! *** Sigma * S ****
    ! Again, as above or below for Kappa, the conversion factor is
    ! * elem_charge_SI**3 / (hbar_SI**2) * 1.e-5_dp
    ! and brings the result to Ampere/m/K
    LocalSigmaS = LocalSigmaS*elem_charge_SI**3/(hbar_SI**2)*1.e-5_dp

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
    LocalKappa = LocalKappa*elem_charge_SI**3/(hbar_SI**2)*1.e-5_dp
    ! THIS IS NOW THE THERMAL CONDUCTIVITY IN SI UNITS, i.e. in W/meter/K

    ! Now I send the different pieces to the local node
    if (on_root) then
      allocate (ElCond(6, TempNumPoints, MuNumPoints), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating ElCond in boltzwann_main')
      allocate (SigmaS(6, TempNumPoints, MuNumPoints), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating SigmaS in boltzwann_main')
      allocate (Seebeck(9, TempNumPoints, MuNumPoints), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating Seebeck in boltzwann_main')
      allocate (Kappa(6, TempNumPoints, MuNumPoints), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating Kappa in boltzwann_main')
    else
      ! In principle, this should not be needed, because we use ElCond,
      ! Seebeck and Kappa only on the root node. However, since all
      ! processors call comms_gatherv a few lines below, and one argument
      ! is ElCond(1,1,1), some compilers complain.
      allocate (ElCond(1, 1, 1), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating ElCond in boltzwann_main (2)')
      allocate (SigmaS(1, 1, 1), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating SigmaS in boltzwann_main (2)')
      allocate (Seebeck(1, 1, 1), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating Seebeck in boltzwann_main (2)')
      allocate (Kappa(1, 1, 1), stat=ierr)
      if (ierr /= 0) call io_error('Error in allocating Kappa in boltzwann_main (2)')
    end if

    ! The 6* factors are due to the fact that for each (T,mu) pair we have 6 components (xx,xy,yy,xz,yz,zz)
    ! NOTE THAT INSTEAD SEEBECK IS A FULL MATRIX AND HAS 9 COMPONENTS!
    call comms_gatherv(LocalElCond, 6*counts(my_node_id), ElCond, 6*counts, 6*displs)
    call comms_gatherv(LocalSigmaS, 6*counts(my_node_id), SigmaS, 6*counts, 6*displs)
    call comms_gatherv(LocalSeebeck, 9*counts(my_node_id), Seebeck, 9*counts, 9*displs)
    call comms_gatherv(LocalKappa, 6*counts(my_node_id), Kappa, 6*counts, 6*displs)

    if (on_root .and. (timing_level > 0)) call io_stopwatch('boltzwann_main: calc_props', 2)

    ! Open files and print
    if (on_root) then
      elcond_unit = io_file_unit()
      open (unit=elcond_unit, file=trim(seedname)//'_elcond.dat')
      write (elcond_unit, '(A)') "# Written by the BoltzWann module of the Wannier90 code."
      write (elcond_unit, '(A)') "# [Electrical conductivity in SI units, i.e. in 1/Ohm/m]"
      write (elcond_unit, '(A)') "# Mu(eV) Temp(K) ElCond_xx ElCond_xy ElCond_yy ElCond_xz ElCond_yz ElCond_zz"
      do MuIdx = 1, MuNumPoints
        do TempIdx = 1, TempNumPoints
          write (elcond_unit, 103) MuArray(MuIdx), TempArray(TempIdx), ElCond(:, TempIdx, MuIdx)
        end do
      end do
      close (elcond_unit)
      if (iprint > 1) write (stdout, '(3X,A)') "Electrical conductivity written on the "//trim(seedname)//"_elcond.dat file."

      sigmas_unit = io_file_unit()
      open (unit=sigmas_unit, file=trim(seedname)//'_sigmas.dat')
      write (sigmas_unit, '(A)') "# Written by the BoltzWann module of the Wannier90 code."
      write (sigmas_unit, '(A)') "# [(Electrical conductivity * Seebeck coefficient) in SI units, i.e. in Ampere/m/K]"
      write (sigmas_unit, '(A)') "# Mu(eV) Temp(K) (Sigma*S)_xx (Sigma*S)_xy (Sigma*S)_yy (Sigma*S)_xz (Sigma*S)_yz (Sigma*S)_zz"
      do MuIdx = 1, MuNumPoints
        do TempIdx = 1, TempNumPoints
          write (sigmas_unit, 103) MuArray(MuIdx), TempArray(TempIdx), SigmaS(:, TempIdx, MuIdx)
        end do
      end do
      close (sigmas_unit)
      if (iprint > 1) write (stdout, '(3X,A)') &
        "sigma*S (sigma=el. conductivity, S=Seebeck coeff.) written on the "//trim(seedname)//"_sigmas.dat file."

      seebeck_unit = io_file_unit()
      open (unit=seebeck_unit, file=trim(seedname)//'_seebeck.dat')
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
      if (iprint > 1) write (stdout, '(3X,A)') "Seebeck coefficient written on the "//trim(seedname)//"_seebeck.dat file."

      kappa_unit = io_file_unit()
      open (unit=kappa_unit, file=trim(seedname)//'_kappa.dat')
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
      if (iprint > 1) write (stdout, '(3X,A)') "K coefficient written on the "//trim(seedname)//"_kappa.dat file."
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
    if (ierr /= 0) call io_error('Error in deallocating TempArray in boltzwann_main')
    deallocate (KTArray, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating KTArray in boltzwann_main')
    deallocate (MuArray, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating MuArray in boltzwann_main')
    deallocate (TDFEnergyArray, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating TDFEnergyArray in boltzwann_main')
    deallocate (TDF, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating TDF in boltzwann_main')
    deallocate (LocalElCond, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating LocalElCond in boltzwann_main')
    deallocate (LocalSigmaS, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating LocalSigmaS in boltzwann_main')
    deallocate (LocalSeebeck, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating LocalSeebeck in boltzwann_main')
    deallocate (LocalKappa, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating LocalKappa in boltzwann_main')

    deallocate (ElCond, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating ElCond in boltzwann_main')
    deallocate (SigmaS, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating SigmaS in boltzwann_main')
    deallocate (Seebeck, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating Seebeck in boltzwann_main')
    deallocate (Kappa, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating Kappa in boltzwann_main')

    deallocate (IntegrandArray, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating IntegrandArray in boltzwann_main')

    if (on_root .and. (timing_level > 0)) call io_stopwatch('boltzwann_main', 2)

101 FORMAT(7G18.10)
102 FORMAT(19G18.10)
103 FORMAT(8G18.10)
104 FORMAT(11G18.10)

  end subroutine boltzwann_main

  subroutine calcTDFandDOS(TDF, TDFEnergyArray)
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
    !!  If the input flag boltz_bandshift is set to .true., the code will also shift the
    !!       conduction bands by a given amount, as defined by the boltz_bandshift_energyshift
    !!       and boltz_bandshift_firstband input flags.
    !!
    use w90_get_oper, only: get_HH_R, get_SS_R, HH_R
    use w90_parameters, only: num_wann, boltz_calc_also_dos, &
      boltz_dos_energy_step, boltz_dos_energy_min, boltz_dos_energy_max, &
      boltz_dos_adpt_smr, boltz_dos_smr_fixed_en_width, boltz_dos_adpt_smr_fac, &
      boltz_dos_adpt_smr_max, &
      param_get_smearing_type, boltz_dos_smr_index, boltz_tdf_smr_index
    use w90_utility, only: utility_diagonalize
    use w90_wan_ham, only: wham_get_eig_deleig

    real(kind=dp), dimension(:, :, :), intent(out)   :: TDF ! (coordinate,Energy,spin)
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
    real(kind=dp), dimension(:), intent(in)      :: TDFEnergyArray
    !! TDFEnergyArray The array with the energies for which the TDF is calculated, in eV
    ! Comments:
    ! issue warnings if going outside of the energy window
    ! check that we actually get hbar*velocity in eV*angstrom

    real(kind=dp), dimension(3) :: kpt, orig_kpt
    integer :: loop_tot, loop_x, loop_y, loop_z, ierr

    complex(kind=dp), allocatable :: HH(:, :)
    complex(kind=dp), allocatable :: delHH(:, :, :)
    complex(kind=dp), allocatable :: UU(:, :)
    real(kind=dp) :: del_eig(num_wann, 3)
    real(kind=dp) :: eig(num_wann), levelspacing_k(num_wann)

    real(kind=dp), allocatable :: DOS_EnergyArray(:)
    real(kind=dp), allocatable :: DOS_k(:, :), TDF_k(:, :, :)
    real(kind=dp), allocatable :: DOS_all(:, :)
    real(kind=dp)              :: kweight
    integer                    :: ndim, DOS_NumPoints, i, j, k, EnIdx

    character(len=20)          :: numfieldsstr
    integer :: boltzdos_unit, num_s_steps, NumPtsRefined

    real(kind=dp), parameter :: SPACING_THRESHOLD = 1.e-3
    real(kind=dp) :: min_spacing, max_spacing

    if (on_root .and. (timing_level > 0)) call io_stopwatch('calcTDF', 1)
    if (on_root) then
      if (boltz_calc_also_dos) then
        write (stdout, '(3X,A)') "Calculating Transport Distribution function (TDF) and DOS..."
      else
        write (stdout, '(3X,A)') "Calculating Transport Distribution function (TDF)..."
      end if
    end if

    ! I call once the routine to calculate the Hamiltonian in real-space <0n|H|Rm>
    call get_HH_R
    if (spin_decomp) then
      ndim = 3
      call get_SS_R
    else
      ndim = 1
    end if

    ! Some initial checks
    if (size(TDF, 1) /= 6 .or. size(TDF, 2) /= size(TDFEnergyArray) .or. size(TDF, 3) /= ndim) then
      call io_error('Wrong size for the TDF array in calcTDF')
    end if

    ! I zero the TDF array before starting
    TDF = 0._dp

    allocate (TDF_k(6, size(TDFEnergyArray), ndim), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating TDF_k in calcTDF')

    allocate (HH(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating HH in calcTDF')
    allocate (delHH(num_wann, num_wann, 3), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating delHH in calcTDF')
    allocate (UU(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating UU in calcTDF')

    DOS_NumPoints = int(floor((boltz_dos_energy_max - boltz_dos_energy_min)/boltz_dos_energy_step)) + 1
    if (DOS_NumPoints .eq. 1) DOS_NumPoints = 2
    allocate (DOS_EnergyArray(DOS_NumPoints), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating DOS_EnergyArray in calcTDF')
    do i = 1, DOS_NumPoints
      DOS_EnergyArray(i) = boltz_dos_energy_min + real(i - 1, dp)*boltz_dos_energy_step
    end do

    allocate (DOS_k(size(DOS_EnergyArray), ndim), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating DOS_k in calcTDF')
    allocate (DOS_all(size(DOS_EnergyArray), ndim), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating DOS_all in calcTDF')
    dos_all = 0.0_dp

    ! I open the output files
    if (boltz_calc_also_dos .and. on_root) then
      boltzdos_unit = io_file_unit()
      open (unit=boltzdos_unit, file=trim(seedname)//'_boltzdos.dat')
    end if

    if (boltz_calc_also_dos .and. on_root .and. (iprint > 1)) then
      write (stdout, '(5X,A)') "Smearing for DOS: "
      if (boltz_dos_adpt_smr) then
        write (stdout, '(7X,A)') trim(param_get_smearing_type(boltz_dos_smr_index))//", adaptive"
      else
        if (boltz_dos_smr_fixed_en_width/(DOS_EnergyArray(2) - DOS_EnergyArray(1)) < &
            min_smearing_binwidth_ratio) then
          write (stdout, '(7X,A)') "Unsmeared (use smearing width larger than bin width to smear)"
        else
          write (stdout, '(7X,A,G18.10)') trim(param_get_smearing_type(boltz_dos_smr_index))// &
            ", non-adaptive, width (eV) =", boltz_dos_smr_fixed_en_width
        end if
      end if
    end if
    if (boltz_calc_also_dos .and. boltz_dos_adpt_smr .and. (boltz_dos_smr_fixed_en_width .ne. 0._dp) .and. on_root) then
      write (stdout, '(5X,A)') "*** WARNING! boltz_dos_smr_fixed_en_width ignored since you chose"
      write (stdout, '(5X,A)') "             an adaptive smearing."
    end if

    if (on_root .and. (iprint > 1)) then
      if (boltz_TDF_smr_fixed_en_width/(TDFEnergyArray(2) - TDFEnergyArray(1)) < min_smearing_binwidth_ratio) then
        write (stdout, '(5X,A)') "Smearing for TDF: "
        write (stdout, '(7X,A)') "Unsmeared (use smearing width larger than bin width to smear)"
      else
        write (stdout, '(5X,A)') "Smearing for TDF: "
        write (stdout, '(7X,A,G18.10)') &
          trim(param_get_smearing_type(boltz_TDF_smr_index))//", non-adaptive, width (eV) =", &
          boltz_TDF_smr_fixed_en_width
      end if
    end if

    if (on_root) then
      write (stdout, '(5X,A,I0,A,I0,A,I0)') "k-grid used for band interpolation in BoltzWann: ", &
        boltz_kmesh(1), 'x', boltz_kmesh(2), 'x', boltz_kmesh(3)
      write (stdout, '(5X,A,I1)') "Number of electrons per state: ", num_elec_per_state
      write (stdout, '(5X,A,G18.10)') "Relaxation time (fs): ", boltz_relax_time
      if (iprint > 1) then
        write (stdout, '(5X,A,G18.10)') "Energy step for TDF (eV): ", boltz_tdf_energy_step
      end if
    end if
    kweight = 1.0_dp/real(PRODUCT(boltz_kmesh), kind=dp)

    if (boltz_bandshift .and. on_root) then
      write (stdout, '(5X,A,I0,A,G18.10,A)') "Shifting energy bands with index >= ", boltz_bandshift_firstband, " by ", &
        boltz_bandshift_energyshift, " eV."
    end if

    NumPtsRefined = 0
    min_spacing = 1.e10_dp ! very large initial value
    max_spacing = 0.e0_dp
    ! I loop over all kpoints
    do loop_tot = my_node_id, PRODUCT(boltz_kmesh) - 1, num_nodes

      ! I get the coordinates for the x,y,z components starting from a single loop variable
      ! (which is better for parallelization purposes)
      ! Important! This works only if loop_tot starts from ZERO and ends with
      !            PRODUCT(boltz_kmesh)-1, so be careful when parallelizing
      loop_x = loop_tot/(boltz_kmesh(2)*boltz_kmesh(3))
      loop_y = (loop_tot - loop_x*(boltz_kmesh(2)*boltz_kmesh(3)))/boltz_kmesh(3)
      loop_z = loop_tot - loop_x*(boltz_kmesh(2)*boltz_kmesh(3)) - loop_y*boltz_kmesh(3)

      ! kpt(i) is in in the [0,d-1]/d range, with d=boltz_kmesh(i)
      kpt(1) = (real(loop_x, dp)/real(boltz_kmesh(1), dp))
      kpt(2) = (real(loop_y, dp)/real(boltz_kmesh(2), dp))
      kpt(3) = (real(loop_z, dp)/real(boltz_kmesh(3), dp))

      ! Here I get the band energies and the velocities
      call wham_get_eig_deleig(kpt, eig, del_eig, HH, delHH, UU)
      call dos_get_levelspacing(del_eig, boltz_kmesh, levelspacing_k)

      ! Here I apply a scissor operator to the conduction bands, if required in the input
      if (boltz_bandshift) then
        eig(boltz_bandshift_firstband:) = eig(boltz_bandshift_firstband:) + boltz_bandshift_energyshift
      end if

      call TDF_kpt(kpt, TDFEnergyArray, eig, del_eig, TDF_k)
      ! As above, the sum of TDF_k * kweight amounts to calculate
      ! spin_degeneracy * V_cell/(2*pi)^3 * \int_BZ d^3k
      ! so that we divide by the cell_volume (in Angstrom^3) to have
      ! the correct integral
      TDF = TDF + TDF_k*kweight/cell_volume

      ! DOS part !

      if (boltz_calc_also_dos) then
        if (boltz_dos_adpt_smr) then

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
                        (/real(i, kind=dp)/real(boltz_kmesh(1), dp)/4._dp, &
                          real(j, kind=dp)/real(boltz_kmesh(2), dp)/4._dp, &
                          real(k, kind=dp)/real(boltz_kmesh(3), dp)/4._dp/)
                  call wham_get_eig_deleig(kpt, eig, del_eig, HH, delHH, UU)
                  call dos_get_levelspacing(del_eig, boltz_kmesh, levelspacing_k)
                  call dos_get_k(kpt, DOS_EnergyArray, eig, dos_k, &
                                 smr_index=boltz_dos_smr_index, &
                                 adpt_smr_fac=boltz_dos_adpt_smr_fac, &
                                 adpt_smr_max=boltz_dos_adpt_smr_max, &
                                 levelspacing_k=levelspacing_k)
                  ! I divide by 8 because I'm substituting a point with its 8 neighbors
                  dos_all = dos_all + dos_k*kweight/8.
                end do
              end do
            end do
          else
            call dos_get_k(kpt, DOS_EnergyArray, eig, dos_k, &
                           smr_index=boltz_dos_smr_index, &
                           adpt_smr_fac=boltz_dos_adpt_smr_fac, &
                           adpt_smr_max=boltz_dos_adpt_smr_max, &
                           levelspacing_k=levelspacing_k)
            dos_all = dos_all + dos_k*kweight
          end if
        else
          call dos_get_k(kpt, DOS_EnergyArray, eig, dos_k, &
                         smr_index=boltz_dos_smr_index, &
                         smr_fixed_en_width=boltz_dos_smr_fixed_en_width)
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
    if (boltz_calc_also_dos) then
      call comms_reduce(DOS_all(1, 1), size(DOS_all), 'SUM')
      call comms_reduce(NumPtsRefined, 1, 'SUM')
      call comms_reduce(min_spacing, 1, 'MIN')
      call comms_reduce(max_spacing, 1, 'MAX')
    end if
    ! I sum the results of the calculation on all nodes, and I store them on all
    ! nodes (because for the following, each node will do a different calculation,
    ! each of which will require the whole knowledge of the TDF array)
    call comms_allreduce(TDF(1, 1, 1), size(TDF), 'SUM')

    if (boltz_calc_also_dos .and. on_root) then
      write (boltzdos_unit, '(A)') "# Written by the BoltzWann module of the Wannier90 code."
      write (boltzdos_unit, '(A)') "# The first column."
      if (boltz_dos_adpt_smr) then
        write (boltzdos_unit, '(A)') '# The second column is the adaptively-smeared DOS'
        write (boltzdos_unit, '(A)') '# (see Yates et al., PRB 75, 195121 (2007)'
        if (spin_decomp) then
          write (boltzdos_unit, '(A)') '# The third column is the spin-up projection of the DOS'
          write (boltzdos_unit, '(A)') '# The fourth column is the spin-down projection of the DOS'
        end if
        write (boltzdos_unit, '(A,1X,G14.6)') '# Smearing coefficient: ', boltz_dos_adpt_smr_fac
        write (boltzdos_unit, '(A,I0,A,I0)') '# Number of points refined: ', NumPtsRefined, &
          ' out of ', product(boltz_kmesh)
        write (boltzdos_unit, '(A,G18.10,A,G18.10,A)') '# (Min spacing: ', min_spacing, &
          ', max spacing: ', max_spacing, ')'
      else
        if (boltz_dos_smr_fixed_en_width/(DOS_EnergyArray(2) - DOS_EnergyArray(1)) < min_smearing_binwidth_ratio) then
          write (boltzdos_unit, '(A)') '# The second column is the unsmeared DOS.'
        else
          write (boltzdos_unit, '(A,G14.6,A)') '# The second column is the DOS for a fixed smearing of ', &
            boltz_dos_smr_fixed_en_width, ' eV.'
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

    if (on_root .and. (timing_level > 0)) call io_stopwatch('calcTDF', 2)
    if (on_root) then
      if (boltz_calc_also_dos) then
        write (stdout, '(3X,A)') "TDF and DOS calculated."
      else
        write (stdout, '(3X,A)') "TDF calculated."
      end if
    end if
    if (on_root) write (stdout, *)

    if (on_root .and. boltz_calc_also_dos) then
      close (boltzdos_unit)
      if (iprint > 1) write (stdout, '(3X,A)') "DOS written on the "//trim(seedname)//"_boltzdos.dat file."
    end if

    deallocate (HH, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating HH in calcTDF')
    deallocate (delHH, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating delHH in calcTDF')
    deallocate (UU, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating UU in calcTDF')
    deallocate (DOS_EnergyArray, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating DOS_EnergyArray in calcTDF')
    deallocate (DOS_k, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating DOS_k in calcTDF')
    deallocate (DOS_all, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating DOS_all in calcTDF')
    deallocate (TDF_k, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating TDF_k in calcTDF')

  end subroutine calcTDFandDOS

  !> This function calculates -dn(E)/dE, where n(E) is the Fermi distribution function.
  !>
  !> \param E  Energy at which we want to calculate -dn(E)/dE, in 1/eV
  !> \param mu Chemical potential in eV
  !> \param KT k_Boltzmann * Temperature, in eV
  function MinusFermiDerivative(E, mu, KT)
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

  subroutine TDF_kpt(kpt, EnergyArray, eig_k, deleig_k, TDF_k)
    !! This subroutine calculates the contribution to the TDF of a single k point
    !!
    !!  This routine does not use the adaptive smearing; in fact, for non-zero temperatures
    !!       one often doesn't even need to smear. It simply uses a standard smearing as defined by
    !!       the variables boltz_TDF_smr_fixed_en_width and boltz_TDF_smr_index
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
    use w90_constants, only: dp, smearing_cutoff, min_smearing_binwidth_ratio
    use w90_utility, only: utility_w0gauss
    use w90_parameters, only: num_wann, spin_decomp, num_elec_per_state, &
      boltz_TDF_smr_fixed_en_width, boltz_TDF_smr_index, boltz_relax_time
    use w90_spin, only: spin_get_nk

    ! Arguments
    !
    real(kind=dp), dimension(3), intent(in)      :: kpt
    !! the three coordinates of the k point vector whose DOS contribution we
    !! want to calculate (in relative coordinates)
    real(kind=dp), dimension(:), intent(in)      :: EnergyArray
    !! array with the energy grid on which to calculate the DOS (in eV)
    !! It must have at least two elements
    real(kind=dp), dimension(:), intent(in)      :: eig_k
    !! array with the eigenvalues at the given k point (in eV)
    real(kind=dp), dimension(:, :), intent(in)    :: deleig_k
    !! array with the band derivatives at the given k point
    !! (in eV * angstrom / (2pi) as internally given by the code)
    !! already corrected in case of degeneracies, as returned by the
    !! wham_get_deleig_a routine
    real(kind=dp), dimension(:, :, :), intent(out) :: TDF_k
    !! TDF_k array in which the contribution is stored. Three dimensions:
    !! TDF_k(ij, energyidx, spinidx), where:
    !! - ij indexes the components of the TDF (symmetric) tensor (1=XX, 2=XY, ...);
    !!   see the global constants defined in the module
    !! - energyidx is the index of the energies, corresponding to the one
    !!   of the EnergyArray array;
    !!  - spinidx=1 contains the total dos; if if spin_decomp==.true., then
    !!  spinidx=2 and spinidx=3 contain the spin-up and spin-down contributions to the DOS

    ! Adaptive smearing
    !
    real(kind=dp) :: smear, arg

    ! Misc/Dummy
    !
    integer          :: BandIdx, loop_f, min_f, max_f
    real(kind=dp)    :: rdum, spn_nk(num_wann), alpha_sq, beta_sq
    real(kind=dp)    :: binwidth, r_num_elec_per_state
    logical          :: DoSmearing

    r_num_elec_per_state = real(num_elec_per_state, kind=dp)

    ! Get spin projections for every band
    !
    if (spin_decomp) call spin_get_nk(kpt, spn_nk)

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
      smear = boltz_TDF_smr_fixed_en_width
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
          rdum = utility_w0gauss(arg, boltz_TDF_smr_index)/smear
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
    TDF_k = TDF_k*boltz_relax_time

  end subroutine TDF_kpt

end module w90_boltzwann
