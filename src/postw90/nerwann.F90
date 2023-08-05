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
!  w90_nerwann: Boltzman transport under a magnetic field    !
!                                                            !
!------------------------------------------------------------!
module w90_nerwann
  !================================================!
  !! Compute Thermomagnetic tranport properties
  !!
  !! NerWann routines by
  !! S.E. Rezaei, M. Zebajadi, K. Esfarjani 
  !! Jan, 2021
  !!
  !! Affiliations:
  !! University of Virginia, Charlottesville, VA (USA)
  !!
  !!
  !! Please cite the following paper when publishing results
  !! obtained using the NerWann module:
  !!
  !! S.E. Rezaei, M. Zebajadi, K. Esfarjani 
  !!    Comput. Mater. Sci. (2022)
  !!    DOI: 10.1016/j.commatsci.2022.111412
  !================================================!

  use w90_comms, only: mpisize, mpirank, comms_gatherv, comms_array_split,comms_reduce, &
    comms_allreduce, w90comm_type
  use w90_constants, only: dp, pw90_physical_constants_type,min_smearing_binwidth_ratio,dkpt
  use w90_io, only: io_error, io_stopwatch, io_file_unit
  use w90_utility, only: utility_inv3, utility_inv2

  implicit none

  private

  public :: nerwann_main

  ! Constants to identify the nine components of a tensor 
  integer, parameter :: XX = 1
  integer, parameter :: XY = 2
  integer, parameter :: YY = 3
  integer, parameter :: XZ = 4
  integer, parameter :: YZ = 5
  integer, parameter :: ZZ = 6
  integer, parameter :: YX = 7 
  integer, parameter :: ZX = 8 
  integer, parameter :: ZY = 9 

  character(len=74), parameter :: pub_string_1 = &
                                  "Please cite the following paper when publishing results obtained using    "
  character(len=74), parameter :: pub_string_2 = &
                                  "the Nerwann module "
  character(len=74), parameter :: pub_string_3 = &
                                  "S.E. Rezaei, M. Zebajadi, K. Esfarjani,              "
  character(len=74), parameter :: pub_string_4 = &
                                  "Comput. Mater. Sci. (2022)10.1016/j.commatsci.2022.111412         "

contains

  subroutine nerwann_main(pw90_nerwann, dis_manifold, kpt_latt, &
                            pw90_band_deriv_degen, postw90_oper, pw90_spin,physics, ws_region, &
                            w90_system, wannier_data, ws_distance, wigner_seitz,print_output, &
                            HH_R, SS_R, v_matrix, u_matrix, eigval,real_lattice, scissors_shift, &
                            mp_grid, num_wann, num_bands, num_kpts,effective_model, &
                            have_disentangled, spin_decomp, seedname, stdout,comm)
    !! This is the main routine of the NerWan module.
    !! It computes the transport coefficients under a magnetic field using the Boltzmann transport equation.
    !!
    !!  Four files will be generated:
    !!
    !!  1. The total Transport Distribution function (TDF) in units of 1/hbar^2 * eV*fs/angstrom
    !!  2. The Hall coefficient in SI units (m^3/C)
    !!  3. The Nernst Coefficient in SI units (V/K)
    !!  4. The Ettingshausen Coefficient in SI units (m.K/Amp)
    !!
  !================================================!

    use w90_constants, only: dp,dkpt
    use w90_io, only: io_file_unit, io_error, io_stopwatch
    use w90_comms, only: comms_bcast, w90comm_type, mpirank
    use w90_types, only: dis_manifold_type, print_output_type,wannier_data_type, &
      ws_region_type, w90_system_type, ws_distance_type
    use w90_postw90_types, only: pw90_nerwann_type, pw90_spin_mod_type, &
      pw90_band_deriv_degen_type, pw90_oper_read_type, wigner_seitz_type

    implicit none

    ! arguments
    type(pw90_nerwann_type), intent(in) :: pw90_nerwann
    type(dis_manifold_type), intent(in) :: dis_manifold
    type(pw90_band_deriv_degen_type), intent(in) :: pw90_band_deriv_degen
    type(pw90_oper_read_type), intent(in) :: postw90_oper
    type(pw90_spin_mod_type), intent(in) :: pw90_spin
    type(print_output_type), intent(in) :: print_output
    type(pw90_physical_constants_type), intent(in) :: physics
    type(ws_region_type), intent(in) :: ws_region
    type(w90comm_type), intent(in) :: comm
    type(w90_system_type), intent(in) :: w90_system
    type(wannier_data_type), intent(in) :: wannier_data
    type(wigner_seitz_type), intent(inout) :: wigner_seitz
    type(ws_distance_type), intent(inout) :: ws_distance

    complex(kind=dp), allocatable, intent(inout) :: HH_R(:, :, :) !  <0n|r|Rm>
    complex(kind=dp), allocatable, intent(inout) :: SS_R(:, :, :, :) !<0n|sigma_x,y,z|Rm>
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

    integer :: TDFEnergyNumPoints
    integer :: i, j, ierr, EnIdx
    real(kind=dp), allocatable :: TDFtotEnergyArr(:)
    integer :: tdf1_unit,ndim
    integer :: tdf2_unit
!Total TDF 1st+2nd along Bz 
    real(kind=dp), allocatable :: TDF1totz(:, :, :),TDF2totz(:, :, :)! (coordinate,Energy) 

    integer :: LocalIdx, GlobalIdx
    ! I also add 3 times the smearing on each side of the TDF energy array to take into account also possible smearing effects
    real(kind=dp), parameter :: TDF_exceeding_energy_times_smr = 3._dp
    real(kind=dp) :: TDF_exceeding_energy
    real(kind=dp) :: cell_volume

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

    if (print_output%iprint > 0 .and. print_output%timing_level > 0) call io_stopwatch('nerwann_main', 1, stdout, seedname)

    if (print_output%iprint > 0) then
      write (stdout, *)
      write (stdout, '(1x,a)') '*---------------------------------------------------------------------------*'
      write (stdout, '(1x,a)') '|             Boltzmann Transport under a magnetic field (NerWan module)    |'
      write (stdout, '(1x,a)') '*---------------------------------------------------------------------------*'
      write (stdout, '(1x,a)') '| '//pub_string_1//'|'
      write (stdout, '(1x,a)') '| '//pub_string_2//'|'
      write (stdout, '(1x,a)') '| '//pub_string_3//'|'
      write (stdout, '(1x,a)') '| '//pub_string_4//'|'
      write (stdout, '(1x,a)') '*---------------------------------------------------------------------------*'
      write (stdout, *)
    end if 


    if (pw90_nerwann%tdf_smearing%use_adaptive) then
      call io_error('Adaptive smearing not allowed in Nerwann TDF', stdout, seedname)
    endif
    ! TDFtotEnergyArr
    TDF_exceeding_energy = max(TDF_exceeding_energy_times_smr*pw90_nerwann%tdf_smearing%fixed_width, 0.2_dp)
    TDFEnergyNumPoints = int(floor((dis_manifold%win_max - dis_manifold%win_min &
                                    + 2._dp*TDF_exceeding_energy)/pw90_nerwann%tdf_energy_step)) + 1
    if (TDFEnergyNumPoints .eq. 1) TDFEnergyNumPoints = 2
    allocate (TDFtotEnergyArr(TDFEnergyNumPoints), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating TDFtotEnergyArr in nerwann_main',stdout, seedname)
    do i = 1, TDFEnergyNumPoints
      TDFtotEnergyArr(i) = dis_manifold%win_min - TDF_exceeding_energy &
                          + real(i - 1, dp)*pw90_nerwann%tdf_energy_step
    end do

    if (spin_decomp) then
      ndim = 3
    else
      ndim = 1
    end if

    !Allocation for  TDF tensor with 9 components
    allocate (TDF1totz(9, TDFEnergyNumPoints, ndim), stat=ierr) !total TDF 1st+2nd along Bz
    if (ierr /= 0) call io_error('Error in allocating TDF1totz in nerwann_main', stdout, seedname)
    allocate (TDF2totz(9, TDFEnergyNumPoints, ndim), stat=ierr) !total TDF 1st+2nd along Bz
    if (ierr /= 0) call io_error('Error in allocating TDF2totz in nerwann_main', stdout, seedname)
	call calcTDFtot(pw90_nerwann, dis_manifold, kpt_latt, postw90_oper, pw90_band_deriv_degen, pw90_spin, &
                       ws_region, print_output, wannier_data, ws_distance, wigner_seitz, HH_R, SS_R, u_matrix, &
                       v_matrix, eigval, real_lattice, TDF1totz,TDF2totz, TDFtotEnergyArr, &
                       cell_volume, scissors_shift, mp_grid, num_bands, num_kpts, num_wann, &
                       w90_system%num_valence_bands, w90_system%num_elec_per_state, effective_model, &
                       have_disentangled, spin_decomp, seedname, stdout, comm)

!================================!
!Total TDF 1st term
!================================!
    if (on_root) then
      tdf1_unit = io_file_unit()
      open (unit=tdf1_unit, file=trim(seedname)//'_tdf1.dat')
      write (tdf1_unit, '(A)') "# Developed by  NerWann module of the Wannier90 code."
     write (tdf1_unit, '(A)') "# Total Transport distribution function (in SI units of C^2/kg/m^3)"// &
       " vs energy in eV"
      write (tdf1_unit, '(A)') "# Content of the columns:"
      write (tdf1_unit, '(A)') "# Energy xxz xyz yyz xzz yzz zzz yxz zxz zyz"
      write (tdf1_unit, '(A)') '#   (if spin decomposition is required, 18 further columns are provided, with the 9'
      write (tdf1_unit, '(A)') '#    components of the TDFtot for the spin up, followed by those for the spin down)'
      do i = 1, size(TDFtotEnergyArr)
        if (ndim .eq. 1) then
          write (tdf1_unit, 104) TDFtotEnergyArr(i), TDF1totz(:, i, 1)
        else
          write (tdf1_unit, 104) TDFtotEnergyArr(i), TDF1totz(:, i, :)
        end if
      end do
      close (tdf1_unit)
	    if (print_output%iprint > 1) &
        write (stdout, '(3X,A)') "Total Transport distribution function along Bz written on the " &
&              //trim(seedname)//"_tdf1.dat file."
    end if

!================================!
!Total TDF 1st term
!================================!
    if (on_root) then
      tdf2_unit = io_file_unit()
      open (unit=tdf2_unit, file=trim(seedname)//'_tdf2.dat')
      write (tdf2_unit, '(A)') "# Developed by  NerWann module of the Wannier90 code."
     write (tdf2_unit, '(A)') "# Total Transport distribution function (in SI units of C^2/kg/m^3/S)"// &
       " vs energy in eV"
      write (tdf2_unit, '(A)') "# Content of the columns:"
      write (tdf2_unit, '(A)') "# Energy xxz xyz yyz xzz yzz zzz yxz zxz zyz"
      write (tdf2_unit, '(A)') '#   (if spin decomposition is required, 18 further columns are provided, with the 9'
      write (tdf2_unit, '(A)') '#    components of the TDFtot for the spin up, followed by those for the spin down)'
      do i = 1, size(TDFtotEnergyArr)
        if (ndim .eq. 1) then
          write (tdf2_unit, 104) TDFtotEnergyArr(i), TDF2totz(:, i, 1)
        else
          write (tdf2_unit, 104) TDFtotEnergyArr(i), TDF2totz(:, i, :)
        end if
      end do
      close (tdf2_unit)
	    if (print_output%iprint > 1) &
        write (stdout, '(3X,A)') "Total Transport distribution function along Bz written on the " &
&              //trim(seedname)//"_tdf2.dat file."
    end if


    ! *********************************************************************************
    !The total TDF is computed and will be used to calculate the thermomagnetic responses

    if (on_root .and. (print_output%timing_level > 0)) call io_stopwatch('nerwann_main: calc_props', 1, stdout, seedname)

    ! I obtain the counts and displs arrays, which tell how I should partition a big array
    ! on the different nodes.



!!!Allocation for  response functions along Bz 
!!!END  


    do LocalIdx = 1, counts(my_node_id)
      GlobalIdx = displs(my_node_id) + (LocalIdx - 1)


    end do


    if (on_root .and. (print_output%timing_level > 0)) call io_stopwatch('nerwann_main: calc_props', 2, stdout, seedname)

    ! Open files and print
    if (on_root) then
      write (stdout, '(3X,A)') "Thermomagnetic coefficients calculated."
      write (stdout, *)
      write (stdout, '(1x,a)') '*---------------------------------------------------------------------------*'
      write (stdout, '(1x,a)') '|                        End of the NerWann module                        |'
      write (stdout, '(1x,a)') '*---------------------------------------------------------------------------*'
      write (stdout, *)
    end if

    ! Before ending, I deallocate memory
    deallocate (TDFtotEnergyArr, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating TDFtotEnergyArr in nerwann_main', stdout, seedname)
    deallocate (TDF1totz, stat=ierr) 
    if (ierr /= 0) call io_error('Error in deallocating TDF1totz in nerwann_main', stdout, seedname)
    deallocate (TDF2totz, stat=ierr) 
    if (ierr /= 0) call io_error('Error in deallocating TDF2totz in nerwann_main', stdout, seedname)


    if (on_root .and. (print_output%timing_level > 0)) call io_stopwatch('nerwann_main', 2, stdout, seedname)

101 FORMAT(7G18.10)
102 FORMAT(19G18.10)
103 FORMAT(8G18.10)
104 FORMAT(11G18.10)
  end subroutine nerwann_main

  subroutine calcTDFtot(pw90_nerwann, dis_manifold, kpt_latt, postw90_oper, pw90_band_deriv_degen, &
                           pw90_spin, ws_region, print_output, wannier_data, ws_distance, wigner_seitz, HH_R, &
                           SS_R, u_matrix, v_matrix, eigval, real_lattice, TDF1totz,TDF2totz,TDFtotEnergyArr, &
                           cell_volume, scissors_shift, mp_grid, num_bands, &
                           num_kpts, num_wann, num_valence_bands, num_elec_per_state, &
                           effective_model, have_disentangled, spin_decomp, seedname, stdout, comm)
!============================!
! This subroutine is aimed to calculate the Toal Transport distribution Function
! In SI units of m^2*C^3/S^3
    !!
    !! The TDFtotEnergyArr must be already allocated and initialized with the
    !! energies in eV before calling this routine.
    !!
    !!  The TDFtot array is allocated with dimensions 9 *size(TDFtotEnergyArr) * ndim, before calling
    !! this routine, where ndim=1 if spin_decomp==false, or ndim=3 if
    !spin_decomp==true. This is not checked.
    !!
    !!  If run in parallel, at the end each processor will have a copy of the full TDFtot array
    !!
    !!  If the input flag pw90_nerwann_bandshift is set to .true., the code will move
    !!  all bands above the valence band (pw90_nerwann_bandshift_firstband) as much as the pw90_nerwann_bandshift_energyshift

    use w90_constants, only: dp, dkpt
    use w90_comms, only: comms_bcast, w90comm_type, mpirank
    use w90_io, only: io_file_unit, io_error, io_stopwatch
    use w90_utility, only: utility_recip_lattice_base
    use w90_get_oper, only: get_HH_R, get_SS_R
    use w90_types, only: print_output_type, wannier_data_type, dis_manifold_type, &
      ws_region_type, ws_distance_type
    use w90_postw90_types, only: pw90_nerwann_type, pw90_spin_mod_type, &
      pw90_band_deriv_degen_type, pw90_oper_read_type, wigner_seitz_type
    use w90_readwrite, only: w90_readwrite_get_smearing_type
    use w90_wan_ham, only: wham_get_eig_deleig, wham_get_eig_deltwoeig, & 
      Omega_operator

!Variables
    implicit none

    ! arguments
    type(pw90_nerwann_type), intent(in) :: pw90_nerwann
    type(dis_manifold_type), intent(in) :: dis_manifold
    real(kind=dp), intent(in) :: kpt_latt(:, :)
    type(pw90_band_deriv_degen_type), intent(in) :: pw90_band_deriv_degen
    type(pw90_oper_read_type), intent(in) :: postw90_oper
    type(pw90_spin_mod_type), intent(in) :: pw90_spin
    type(print_output_type), intent(in) :: print_output
    type(ws_region_type), intent(in) :: ws_region
    type(w90comm_type), intent(in) :: comm
    type(wannier_data_type), intent(in) :: wannier_data
    type(wigner_seitz_type), intent(inout) :: wigner_seitz
    type(ws_distance_type), intent(inout) :: ws_distance

    integer, intent(in) :: num_wann, num_bands, num_kpts, num_valence_bands, num_elec_per_state
    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: stdout

    real(kind=dp), intent(out)   :: TDF1totz(:,:,:),TDF2totz(:,:,:)
    real(kind=dp), intent(in)      :: TDFtotEnergyArr(:)

    !! TDFtotEnergyArr The array with the energies for which the TDF is calculated, in eV

! Dummy Variables
    
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
    real(kind=dp) :: kpt(3), orig_kpt(3), kredpp(3), krednn(3), & 
                     kcart(3), nkdk(3)
    integer :: loop_tot, loop_x, loop_y, loop_z, ierr

    complex(kind=dp), allocatable :: HH(:, :), HHpp(:, :), HHnn(:, :) 
    complex(kind=dp), allocatable :: delHH(:, :,:),delHHnn(:,:,:),delHHpp(:,:,:)
    complex(kind=dp), allocatable :: UU(:, :), UUnn(:,:), UUpp(:,:)
    real(kind=dp) :: vel(num_wann, 3), velnn(num_wann, 3), velpp(num_wann, 3) 
    real(kind=dp) ::dvdk(num_wann,3,3) 
    real(kind=dp) :: eig(num_wann), eigpp(num_wann),eignn(num_wann), levelspacing_k(num_wann) 
    real(kind=dp)::omga_bndx(num_wann),omga_bndy(num_wann),omga_bndz(num_wann) 

    real(kind=dp), allocatable :: TDF1_kz(:,:,:),TDF2_kz(:,:,:)
    real(kind=dp)              :: kweight
    integer                    :: ndim, i, j, k, EnIdx

    character(len=20)          :: numfieldsstr
    integer ::  NumPtsRefined

    real(kind=dp), parameter :: SPACING_THRESHOLD = 1.e-3
    real(kind=dp) :: min_spacing, max_spacing
    type(pw90_physical_constants_type) :: physics

    integer :: my_node_id, num_nodes
    logical :: on_root = .false.
	
    my_node_id = mpirank(comm)
    num_nodes = mpisize(comm)
	if (my_node_id == 0) on_root = .true.

    if (print_output%iprint > 0 .and. (print_output%timing_level > 0)) call io_stopwatch('calcTDFtot', 1, stdout, seedname)
    if (print_output%iprint > 0) then
        write (stdout, '(3X,A)') "Calculating Total Transport Distribution function (TDFtot)..."
    end if

    ! I call once the routine to calculate the Hamiltonian in real-space <0n|H|Rm>
    call get_HH_R(dis_manifold, kpt_latt, print_output, wigner_seitz, HH_R,u_matrix, v_matrix, eigval, &
                  real_lattice, scissors_shift, num_bands, num_kpts, num_wann,num_valence_bands, &
                  effective_model, have_disentangled, seedname, stdout, comm)

    if (spin_decomp) then
      ndim = 3
      call get_SS_R(dis_manifold, kpt_latt, print_output, postw90_oper, SS_R, v_matrix, eigval, &
                    wigner_seitz%irvec, wigner_seitz%nrpts, num_bands, num_kpts, num_wann, have_disentangled, seedname, stdout, &
                    comm)
    else
      ndim = 1
    end if

    ! Initial check of sizes
    if (size(TDF1totz, 1) /= 9 .or. size(TDF1totz, 2) /= size(TDFtotEnergyArr) .or. size(TDF1totz, 3) /= ndim) then
      call io_error('Wrong size for the TDF1totz array in calcTDFtot', stdout, seedname)
    end if
    if (size(TDF2totz, 1) /= 9 .or. size(TDF2totz, 2) /= size(TDFtotEnergyArr) .or. size(TDF2totz, 3) /= ndim) then
      call io_error('Wrong size for the TDF2totz array in calcTDFtot', stdout, seedname)
    end if

    ! I zero the TDF array before starting
    TDF1totz = 0.0_dp !Setting TDFtotz to Zero 
    TDF2totz = 0.0_dp !Setting TDFtotz to Zero 


    allocate (TDF1_kz(9, size(TDFtotEnergyArr), ndim), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating TDF1_kz in calcTDFtot', stdout, seedname) 
    allocate (TDF2_kz(9, size(TDFtotEnergyArr), ndim), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating TDF2_kz in calcTDFtot', stdout, seedname)
    allocate (HH(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating HH in calcTDFtot', stdout, seedname)
    allocate (HHpp(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating HHpp in calcTDFtot', stdout, seedname)
    allocate (HHnn(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating HHnn in calcTDFtot', stdout, seedname)
    allocate (delHH(num_wann, num_wann, 3), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating delHH in calcTDFtot', stdout, seedname)
    allocate (delHHpp(num_wann, num_wann, 3), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating delHHpp in calcTDFtot', stdout, seedname)
    allocate (delHHnn(num_wann, num_wann, 3), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating delHHnn in calcTDFtot', stdout, seedname)
    allocate (UU(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating UU in calcTDFtot', stdout, seedname)
    allocate (UUpp(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating UUpp in calcTDFtot', stdout, seedname)
    allocate (UUnn(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating UUnn in calcTDFtot', stdout, seedname)


    if (on_root .and. (print_output%iprint > 1)) then
      if (pw90_nerwann%tdf_smearing%fixed_width/(TDFtotEnergyArr(2) - TDFtotEnergyArr(1)) &
		  < min_smearing_binwidth_ratio) then
        write (stdout, '(5X,A)') "Smearing for TDF: "
        write (stdout, '(7X,A)') "Unsmeared (use smearing width larger than bin width to smear)"
      else
        write (stdout, '(5X,A)') "Smearing for TDF: "
        write (stdout, '(7X,A,G18.10)') &
          trim(w90_readwrite_get_smearing_type(pw90_nerwann%tdf_smearing%type_index))//", non-adaptive, width (eV) =", &
          pw90_nerwann%tdf_smearing%fixed_width
      end if
    end if

    if (on_root) then
      write (stdout, '(5X,A,I0,A,I0,A,I0)') "k-grid used for band interpolation in NerWann: ", &
        pw90_nerwann%kmesh%mesh(1), 'x', pw90_nerwann%kmesh%mesh(2), 'x', pw90_nerwann%kmesh%mesh(3)
      write (stdout, '(5X,A,I1)') "Number of electrons per state: ", num_elec_per_state
      write (stdout, '(5X,A,G18.10,G18.10,G18.10)') "External magnetic field in NerWann: ", &
        pw90_nerwann%bext(1), pw90_nerwann%bext(2), pw90_nerwann%bext(3)
      if (print_output%iprint > 1) then
        write (stdout, '(5X,A,G18.10)') "Energy step for TDF (eV): ", pw90_nerwann%tdf_energy_step
      end if
    end if
    kweight = 1.0_dp/real(PRODUCT(pw90_nerwann%kmesh%mesh), kind=dp)

    if (pw90_nerwann%bandshift .and. on_root) then
      write (stdout, '(5X,A,I0,A,G18.10,A)') "Shifting energy bands with index >= ", pw90_nerwann%bandshift_firstband, " by ", &
        pw90_nerwann%bandshift_energyshift, " eV."
    end if

call utility_recip_lattice_base(real_lattice, recip_lattice, volume)

    NumPtsRefined = 0
    min_spacing = 1.e10_dp ! very large initial value
    max_spacing = 0.e0_dp

    ! loop over kpoints

    do loop_tot = my_node_id, PRODUCT(pw90_nerwann%kmesh%mesh) - 1, num_nodes

      ! coordinates for the x,y,z components starting from a single loop variable
      ! This works only if loop_tot starts from ZERO and ends with
      !            PRODUCT(pw90_nerwann_kmesh)-1, so be careful when parallelizing
      loop_x = loop_tot/(pw90_nerwann%kmesh%mesh(2)*pw90_nerwann%kmesh%mesh(3))
      loop_y = (loop_tot - loop_x*(pw90_nerwann%kmesh%mesh(2)*pw90_nerwann%kmesh%mesh(3)))/pw90_nerwann%kmesh%mesh(3)
      loop_z = loop_tot - loop_x*(pw90_nerwann%kmesh%mesh(2)*pw90_nerwann%kmesh%mesh(3)) - loop_y*pw90_nerwann%kmesh%mesh(3)

      ! kpt(i) is in in the [0,d-1]/d range, with d=pw90_nerwann_kmesh(i)
      kpt(1) = (real(loop_x, dp)/real(pw90_nerwann%kmesh%mesh(1), dp))
      kpt(2) = (real(loop_y, dp)/real(pw90_nerwann%kmesh%mesh(2), dp))
      kpt(3) = (real(loop_z, dp)/real(pw90_nerwann%kmesh%mesh(3), dp))

      ! Here I get the band energies and the velocities
      call wham_get_eig_deleig(dis_manifold, kpt_latt, pw90_band_deriv_degen, ws_region, print_output, wannier_data, &
                               ws_distance, wigner_seitz, delHH, HH, HH_R, u_matrix, UU, v_matrix, &
                               vel, eig, eigval, kpt, real_lattice, scissors_shift, mp_grid, &
                               num_bands, num_kpts, num_wann, num_valence_bands, effective_model, &
                               have_disentangled, seedname, stdout, comm)

      call  wham_get_eig_deltwoeig(dis_manifold, kpt_latt, pw90_band_deriv_degen,print_output, &
                                  ws_region,wannier_data,ws_distance,wigner_seitz,delHHnn, delHHpp, &
                                  HHnn, HHpp, HH_R, u_matrix, UUnn, UUpp,v_matrix, velnn, &
                                  velpp, eignn, eigpp,eigval, kpt,krednn,kredpp, kcart,nkdk, &
                                  real_lattice,scissors_shift,mp_grid,num_bands, num_kpts, &
                                  num_wann, num_valence_bands, effective_model,have_disentangled, &
                                  seedname, stdout, comm, recip_lattice,dvdk)
      call Omega_operator(dis_manifold, kpt_latt, pw90_band_deriv_degen,print_output, &
                                  ws_region,wannier_data,ws_distance,wigner_seitz,delHH, delHHnn, &
                                  delHHpp, HH, HHnn, HHpp, HH_R, u_matrix, UU, UUnn,UUpp,v_matrix, &
                                  vel, velnn, velpp, eig, eignn, eigpp, eigval, kpt, krednn, kredpp, & 
                                  kcart, nkdk, real_lattice, scissors_shift, mp_grid, num_bands, &
                                  num_kpts, num_wann, num_valence_bands,effective_model, & 
                                  have_disentangled, seedname, stdout, comm, recip_lattice, dvdk, &
                                  omga_bndx, omga_bndy, omga_bndz, pw90_nerwann, physics)


      ! The scissor operator is applied to shift all bands above the valence band, if required in the input  

      if (pw90_nerwann%bandshift) then
        eig(pw90_nerwann%bandshift_firstband:) = eig(pw90_nerwann%bandshift_firstband:) + pw90_nerwann%bandshift_energyshift
      end if
!!Calling Transport Functions
      call TDFtot_kpt(pw90_nerwann, ws_region, pw90_spin, wannier_data,ws_distance, wigner_seitz, HH_R, SS_R, &
                   eig, vel, omga_bndx, omga_bndy,omga_bndz, TDFtotEnergyArr,kpt,&
                   real_lattice, TDF1_kz,TDF2_kz,mp_grid,num_wann, num_elec_per_state, physics, spin_decomp, seedname, stdout)


      TDF1totz = TDF1totz + TDF1_kz*kweight/cell_volume 
      TDF2totz = TDF2totz + TDF2_kz*kweight/cell_volume 
    end do

    ! I sum the results of the calculation on all nodes, and I store them on all
    ! nodes (because for the following, each node will do a different calculation,
    ! each of which will require the whole knowledge of the TDF array)
    call comms_allreduce(TDF1totz(1, 1, 1), size(TDF1totz), 'SUM', stdout, seedname, comm)
    call comms_allreduce(TDF2totz(1, 1, 1), size(TDF2totz), 'SUM', stdout, seedname, comm)

    if (on_root .and. (print_output%timing_level > 0)) call io_stopwatch('calcTDFtot', 2, stdout, seedname)
    if (on_root) then
        write (stdout, '(3X,A)') "TDF calculated."
    end if
    if (on_root) write (stdout, *)


    deallocate (HH, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating HH in calcTDFtot', stdout, seedname)
    deallocate (HHpp, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating HHpp in calcTDFtot', stdout, seedname)
    deallocate (HHnn, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating HHnn in calcTDFtot', stdout, seedname)
    deallocate (delHH, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating delHH in calcTDFtot', stdout, seedname)
    deallocate (delHHpp, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating delHHpp in calcTDFtot', stdout, seedname)
    deallocate (delHHnn, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating delHHnn in calcTDFtot', stdout, seedname)
    deallocate (UU, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating UU in calcTDFtot', stdout, seedname)
    deallocate (UUpp, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating UUpp in calcTDFtot', stdout, seedname)
    deallocate (UUnn, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating UUnn in calcTDFtot', stdout, seedname)
    deallocate (TDF1_kz, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating TDF1_kz in calcTDFtot', stdout, seedname) 
    deallocate (TDF2_kz, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating TDF2_kz in calcTDFtot', stdout, seedname) 
  end subroutine calcTDFtot





  function MinusFermiDerivative(E, mu, KT)
    !================================================!
    !> This function calculates -dn(E)/dE, where n(E) is the Fermi distribution
    !function.
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
  

subroutine TDFtot_kpt(pw90_nerwann, ws_region, pw90_spin, wannier_data,ws_distance, wigner_seitz, HH_R, SS_R, &
                      eig_k,vel_k, omg_bnd1, omg_bnd2, omg_bnd3,EnergyArray,kpt,&
                      real_lattice,TDF1_kz,TDF2_kz,mp_grid, num_wann, num_elec_per_state, physics, spin_decomp,seedname, stdout)
    !================================================!
    !! This subroutine calculates the contribution to the TDF of a single k point
    !!
    use w90_constants, only: dp, smearing_cutoff, min_smearing_binwidth_ratio
    use w90_utility, only: utility_w0gauss
    use w90_spin, only: spin_get_nk
    use w90_types, only: print_output_type, wannier_data_type, ws_region_type, &
      ws_distance_type
    use w90_postw90_types, only: pw90_nerwann_type, pw90_spin_mod_type, wigner_seitz_type


    ! Arguments
    !
    implicit none

    ! arguments
    type(pw90_nerwann_type), intent(in) :: pw90_nerwann
    type(ws_region_type), intent(in) :: ws_region
    type(pw90_spin_mod_type), intent(in) :: pw90_spin
    type(wannier_data_type), intent(in) :: wannier_data
    type(ws_distance_type), intent(inout) :: ws_distance
    type(wigner_seitz_type), intent(in) :: wigner_seitz

    integer, intent(in) :: num_wann
    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: stdout

    real(kind=dp), intent(in) :: kpt(3)
    real(kind=dp), intent(in)      :: EnergyArray(:)
    real(kind=dp), intent(in)      :: eig_k(:)
    real(kind=dp), intent(in)    :: vel_k(:, :)
    real(kind=dp), intent(in)::omg_bnd1(:),omg_bnd2(:),omg_bnd3(:)
    logical, intent(in) :: spin_decomp
    integer, intent(in) :: num_elec_per_state
    real(kind=dp), intent(out):: TDF1_kz(:, :, :),TDF2_kz(:, :, :)
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    type(pw90_physical_constants_type), intent(in) :: physics

    complex(kind=dp), allocatable, intent(inout) :: HH_R(:, :, :) !  <0n|r|Rm>
    complex(kind=dp), allocatable, intent(inout) :: SS_R(:, :, :, :) ! <0n|sigma_x,y,z|Rm>

    character(len=50), intent(in) :: seedname

    ! local variables
    real(kind=dp) :: smear, arg ! Adaptive smearing
    real(kind=dp) :: rdum, spn_nk(num_wann), alpha_sq, beta_sq
    real(kind=dp) :: binwidth, r_num_elec_per_state
    integer :: BandIdx, loop_f, min_f, max_f
    logical :: DoSmearing    

    r_num_elec_per_state = real(num_elec_per_state, kind=dp)

    ! Get spin projections for every band
    !
    if (spin_decomp) call spin_get_nk(ws_region, pw90_spin, wannier_data, ws_distance, wigner_seitz, &
                                      HH_R, SS_R, kpt, real_lattice, spn_nk, &
                                      mp_grid, num_wann, seedname, stdout)

    binwidth = EnergyArray(2) - EnergyArray(1)

    TDF1_kz=0.0_dp !Putting initial value zero
    TDF2_kz=0.0_dp !Putting initial value zero

    do BandIdx = 1, num_wann
      if (spin_decomp) then
        ! Contribution to spin-up DOS of Bloch spinor with component
        ! (alpha,beta) with respect to the chosen quantization axis
        alpha_sq = (1.0_dp + spn_nk(BandIdx))/2.0_dp ! |alpha|^2
        ! Contribution to spin-down DOS
        beta_sq = 1.0_dp - alpha_sq ! |beta|^2 = 1 - |alpha|^2
      end if

      ! Do not use an adaptive smearing here, it would require the knowledge of second derivatives

      ! Faster optimization: I precalculate the indices
      ! Value of the smearing in eV; default = 0 eV, i.e. no smearing
      smear = pw90_nerwann%tdf_smearing%fixed_width
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
          rdum = utility_w0gauss(arg, pw90_nerwann%tdf_smearing%type_index, stdout, seedname)/smear
        else
          rdum = 1._dp/(EnergyArray(2) - EnergyArray(1))
        end if

!  write(*,6)'loop,rdum,occ,vel=',loop_f,rdum,r_num_elec_per_state,vel(BandIdx,1:3) !Added by Dr.Keivan main write
!Total TDF 1st+2nd along Bz
        TDF1_kz(XX, loop_f, 1) = TDF1_kz(XX, loop_f,1)+rdum*r_num_elec_per_state*& 
                        vel_k(BandIdx,1)*vel_k(BandIdx,1)*physics%elem_charge_SI**3/physics%hbar_SI**2*1.e10_dp 
        TDF2_kz(XX, loop_f, 1) = TDF2_kz(XX, loop_f,1)+rdum*r_num_elec_per_state*vel_k(BandIdx,1)*omg_bnd1(BandIdx)*& 
                                physics%elem_charge_SI**4/physics%hbar_SI**3*1.e-10_dp
        TDF1_kz(XY, loop_f, 1) = TDF1_kz(XY, loop_f,1)+rdum*r_num_elec_per_state*& 
                        vel_k(BandIdx,1)*vel_k(BandIdx,2)*physics%elem_charge_SI**3/physics%hbar_SI**2*1.e10_dp 
        TDF2_kz(XY, loop_f, 1) = TDF2_kz(XY, loop_f,1)+rdum*r_num_elec_per_state*vel_k(BandIdx,1)*omg_bnd2(BandIdx)*& 
                                physics%elem_charge_SI**4/physics%hbar_SI**3*1.e-10_dp
        TDF1_kz(YX, loop_f, 1) = TDF1_kz(YX, loop_f,1)+rdum*r_num_elec_per_state*& 
                        vel_k(BandIdx,2)*vel_k(BandIdx,1)*physics%elem_charge_SI**3/physics%hbar_SI**2*1.e10_dp 
        TDF2_kz(YX, loop_f, 1) = TDF2_kz(YX, loop_f,1)+rdum*r_num_elec_per_state*vel_k(BandIdx,2)*omg_bnd1(BandIdx)*& 
                                physics%elem_charge_SI**4/physics%hbar_SI**3*1.e-10_dp
        TDF1_kz(YY, loop_f, 1) = TDF1_kz(YY, loop_f,1)+rdum*r_num_elec_per_state*& 
                        vel_k(BandIdx,2)*vel_k(BandIdx,2)*physics%elem_charge_SI**3/physics%hbar_SI**2*1.e10_dp 
        TDF2_kz(YY, loop_f, 1) = TDF2_kz(YY, loop_f,1)+rdum*r_num_elec_per_state*vel_k(BandIdx,2)*omg_bnd2(BandIdx)*& 
                                physics%elem_charge_SI**4/physics%hbar_SI**3*1.e-10_dp
        TDF1_kz(XZ, loop_f, 1) = TDF1_kz(XZ, loop_f,1)+rdum*r_num_elec_per_state*& 
                        vel_k(BandIdx,1)*vel_k(BandIdx,3)*physics%elem_charge_SI**3/physics%hbar_SI**2*1.e10_dp 
        TDF2_kz(XZ, loop_f, 1) = TDF2_kz(XZ, loop_f,1)+rdum*r_num_elec_per_state*vel_k(BandIdx,1)*omg_bnd3(BandIdx)*& 
                                physics%elem_charge_SI**4/physics%hbar_SI**3*1.e-10_dp
        TDF1_kz(YZ, loop_f, 1) = TDF1_kz(YZ, loop_f,1)+rdum*r_num_elec_per_state*& 
                        vel_k(BandIdx,2)*vel_k(BandIdx,3)*physics%elem_charge_SI**3/physics%hbar_SI**2*1.e10_dp 
        TDF2_kz(YZ, loop_f, 1) = TDF2_kz(YZ, loop_f,1)+rdum*r_num_elec_per_state*vel_k(BandIdx,2)*omg_bnd3(BandIdx)*& 
                                physics%elem_charge_SI**4/physics%hbar_SI**3*1.e-10_dp
        TDF1_kz(ZZ, loop_f, 1) = TDF1_kz(ZZ, loop_f,1)+rdum*r_num_elec_per_state*& 
                        vel_k(BandIdx,3)*vel_k(BandIdx,3)*physics%elem_charge_SI**3/physics%hbar_SI**2*1.e10_dp 
        TDF2_kz(ZZ, loop_f, 1) = TDF2_kz(ZZ, loop_f,1)+rdum*r_num_elec_per_state*vel_k(BandIdx,3)*omg_bnd3(BandIdx)*& 
                                physics%elem_charge_SI**4/physics%hbar_SI**3*1.e-10_dp
        TDF1_kz(ZX, loop_f, 1) = TDF1_kz(ZX, loop_f,1)+rdum*r_num_elec_per_state*& 
                        vel_k(BandIdx,3)*vel_k(BandIdx,1)*physics%elem_charge_SI**3/physics%hbar_SI**2*1.e10_dp 
        TDF2_kz(ZX, loop_f, 1) = TDF2_kz(ZX, loop_f,1)+rdum*r_num_elec_per_state*vel_k(BandIdx,3)*omg_bnd1(BandIdx)*& 
                                physics%elem_charge_SI**4/physics%hbar_SI**3*1.e-10_dp
        TDF1_kz(ZY, loop_f, 1) = TDF1_kz(ZY, loop_f,1)+rdum*r_num_elec_per_state*& 
                        vel_k(BandIdx,3)*vel_k(BandIdx,2)*physics%elem_charge_SI**3/physics%hbar_SI**2*1.e10_dp 
        TDF2_kz(ZY, loop_f, 1) = TDF2_kz(ZY, loop_f,1)+rdum*r_num_elec_per_state*vel_k(BandIdx,3)*omg_bnd2(BandIdx)*& 
                                physics%elem_charge_SI**4/physics%hbar_SI**3*1.e-10_dp

        if (spin_decomp) then

          ! Spin-up contribution
        TDF1_kz(XX, loop_f,2) = TDF1_kz(XX, loop_f,2)+rdum*alpha_sq*&
                        vel_k(BandIdx,1)*vel_k(BandIdx,1)*physics%elem_charge_SI**3/physics%hbar_SI**2*1.e10_dp
        TDF2_kz(XX, loop_f,2) = TDF2_kz(XX,loop_f,2)+rdum*alpha_sq*vel_k(BandIdx,1)*omg_bnd1(BandIdx)*&
                                physics%elem_charge_SI**4/physics%hbar_SI**3*1.e-10_dp
        TDF1_kz(XY, loop_f,2) = TDF1_kz(XY, loop_f,2)+rdum*alpha_sq*&
                        vel_k(BandIdx,1)*vel_k(BandIdx,2)*physics%elem_charge_SI**3/physics%hbar_SI**2*1.e10_dp
        TDF2_kz(XY, loop_f,1) = TDF2_kz(XY,loop_f,2)+rdum*alpha_sq*vel_k(BandIdx,1)*omg_bnd2(BandIdx)*&
                                physics%elem_charge_SI**4/physics%hbar_SI**3*1.e-10_dp
        TDF1_kz(YX, loop_f,2) = TDF1_kz(YX, loop_f,2)+rdum*alpha_sq*&
                        vel_k(BandIdx,2)*vel_k(BandIdx,1)*physics%elem_charge_SI**3/physics%hbar_SI**2*1.e10_dp
        TDF2_kz(YX, loop_f,2) = TDF2_kz(YX,loop_f,2)+rdum*alpha_sq*vel_k(BandIdx,2)*omg_bnd1(BandIdx)*&
                                physics%elem_charge_SI**4/physics%hbar_SI**3*1.e-10_dp
        TDF1_kz(YY, loop_f,2) = TDF1_kz(YY, loop_f,2)+rdum*alpha_sq*&
                        vel_k(BandIdx,2)*vel_k(BandIdx,2)*physics%elem_charge_SI**3/physics%hbar_SI**2*1.e10_dp
        TDF2_kz(YY, loop_f,2) = TDF2_kz(YY,loop_f,2)+rdum*alpha_sq*vel_k(BandIdx,2)*omg_bnd2(BandIdx)*&
                                physics%elem_charge_SI**4/physics%hbar_SI**3*1.e-10_dp
        TDF1_kz(XZ, loop_f,2) = TDF1_kz(XZ, loop_f,2)+rdum*alpha_sq*&
                        vel_k(BandIdx,1)*vel_k(BandIdx,3)*physics%elem_charge_SI**3/physics%hbar_SI**2*1.e10_dp
        TDF2_kz(XZ, loop_f,2) = TDF2_kz(XZ,loop_f,2)+rdum*alpha_sq*vel_k(BandIdx,1)*omg_bnd3(BandIdx)*&
                                physics%elem_charge_SI**4/physics%hbar_SI**3*1.e-10_dp
        TDF1_kz(YZ, loop_f,2) = TDF1_kz(YZ, loop_f,2)+rdum*alpha_sq*&
                        vel_k(BandIdx,2)*vel_k(BandIdx,3)*physics%elem_charge_SI**3/physics%hbar_SI**2*1.e10_dp
        TDF2_kz(YZ, loop_f,2) = TDF2_kz(YZ,loop_f,2)+rdum*alpha_sq*vel_k(BandIdx,2)*omg_bnd3(BandIdx)*&
                                physics%elem_charge_SI**4/physics%hbar_SI**3*1.e-10_dp
        TDF1_kz(ZZ, loop_f,2) = TDF1_kz(ZZ, loop_f,2)+rdum*alpha_sq*&
                        vel_k(BandIdx,3)*vel_k(BandIdx,3)*physics%elem_charge_SI**3/physics%hbar_SI**2*1.e10_dp
        TDF2_kz(ZZ, loop_f,2) = TDF2_kz(ZZ,loop_f,2)+rdum*alpha_sq*vel_k(BandIdx,3)*omg_bnd3(BandIdx)*&
                                physics%elem_charge_SI**4/physics%hbar_SI**3*1.e-10_dp
        TDF1_kz(ZX, loop_f,2) = TDF1_kz(ZX, loop_f,2)+rdum*alpha_sq*&
                        vel_k(BandIdx,3)*vel_k(BandIdx,1)*physics%elem_charge_SI**3/physics%hbar_SI**2*1.e10_dp
        TDF2_kz(ZX, loop_f,2) = TDF2_kz(ZX,loop_f,2)+rdum*alpha_sq*vel_k(BandIdx,3)*omg_bnd1(BandIdx)*&
                                physics%elem_charge_SI**4/physics%hbar_SI**3*1.e-10_dp
        TDF1_kz(ZY, loop_f,2) = TDF1_kz(ZY, loop_f,2)+rdum*alpha_sq*&
                        vel_k(BandIdx,3)*vel_k(BandIdx,2)*physics%elem_charge_SI**3/physics%hbar_SI**2*1.e10_dp
        TDF2_kz(ZY, loop_f,2) = TDF2_kz(ZY,loop_f,2)+rdum*alpha_sq*vel_k(BandIdx,3)*omg_bnd2(BandIdx)*&
                                physics%elem_charge_SI**4/physics%hbar_SI**3*1.e-10_dp
       	!Spin-down contribution
        TDF1_kz(XX, loop_f,3) = TDF1_kz(XX, loop_f,3)+rdum*beta_sq*&
                        vel_k(BandIdx,1)*vel_k(BandIdx,1)*physics%elem_charge_SI**3/physics%hbar_SI**2*1.e10_dp
        TDF2_kz(XX, loop_f,3) = TDF2_kz(XX,loop_f,3)+rdum*beta_sq*vel_k(BandIdx,1)*omg_bnd1(BandIdx)*&
                                physics%elem_charge_SI**4/physics%hbar_SI**3*1.e-10_dp
        TDF1_kz(XY, loop_f,3) = TDF1_kz(XY, loop_f,3)+rdum*beta_sq*&
                        vel_k(BandIdx,1)*vel_k(BandIdx,2)*physics%elem_charge_SI**3/physics%hbar_SI**2*1.e10_dp
        TDF2_kz(XY, loop_f,3) = TDF2_kz(XY,loop_f,3)+rdum*beta_sq*vel_k(BandIdx,1)*omg_bnd2(BandIdx)*&
                                physics%elem_charge_SI**4/physics%hbar_SI**3*1.e-10_dp
        TDF1_kz(YX, loop_f,3) = TDF1_kz(YX, loop_f,3)+rdum*beta_sq*&
                        vel_k(BandIdx,2)*vel_k(BandIdx,1)*physics%elem_charge_SI**3/physics%hbar_SI**2*1.e10_dp
        TDF2_kz(YX, loop_f,3) = TDF2_kz(YX,loop_f,3)+rdum*beta_sq*vel_k(BandIdx,2)*omg_bnd1(BandIdx)*&
                                physics%elem_charge_SI**4/physics%hbar_SI**3*1.e-10_dp
        TDF1_kz(YY, loop_f,3) = TDF1_kz(YY, loop_f,3)+rdum*beta_sq*&
                        vel_k(BandIdx,2)*vel_k(BandIdx,2)*physics%elem_charge_SI**3/physics%hbar_SI**2*1.e10_dp
        TDF2_kz(YY, loop_f,3) = TDF2_kz(YY,loop_f,3)+rdum*beta_sq*vel_k(BandIdx,2)*omg_bnd2(BandIdx)*&
                                physics%elem_charge_SI**4/physics%hbar_SI**3*1.e-10_dp
        TDF1_kz(XZ, loop_f,3) = TDF1_kz(XZ, loop_f,3)+rdum*beta_sq*&
                        vel_k(BandIdx,1)*vel_k(BandIdx,3)*physics%elem_charge_SI**3/physics%hbar_SI**2*1.e10_dp
        TDF2_kz(XZ, loop_f,3) = TDF2_kz(XZ,loop_f,3)+rdum*beta_sq*vel_k(BandIdx,1)*omg_bnd3(BandIdx)*&
                                physics%elem_charge_SI**4/physics%hbar_SI**3*1.e-10_dp
        TDF1_kz(YZ, loop_f,3) = TDF1_kz(YZ, loop_f,3)+rdum*beta_sq*&
                        vel_k(BandIdx,2)*vel_k(BandIdx,3)*physics%elem_charge_SI**3/physics%hbar_SI**2*1.e10_dp
        TDF2_kz(YZ, loop_f,3) = TDF2_kz(YZ,loop_f,3)+rdum*beta_sq*vel_k(BandIdx,2)*omg_bnd3(BandIdx)*&
                                physics%elem_charge_SI**4/physics%hbar_SI**3*1.e-10_dp
        TDF1_kz(ZZ, loop_f,3) = TDF1_kz(ZZ, loop_f,3)+rdum*beta_sq*&
                        vel_k(BandIdx,3)*vel_k(BandIdx,3)*physics%elem_charge_SI**3/physics%hbar_SI**2*1.e10_dp
        TDF2_kz(ZZ, loop_f,3) = TDF2_kz(ZZ,loop_f,3)+rdum*beta_sq*vel_k(BandIdx,3)*omg_bnd3(BandIdx)*&
                                physics%elem_charge_SI**4/physics%hbar_SI**3*1.e-10_dp
        TDF1_kz(ZX, loop_f,3) = TDF1_kz(ZX, loop_f,3)+rdum*beta_sq*&
                        vel_k(BandIdx,3)*vel_k(BandIdx,1)*physics%elem_charge_SI**3/physics%hbar_SI**2*1.e10_dp
        TDF2_kz(ZX, loop_f,3) = TDF2_kz(ZX,loop_f,3)+rdum*beta_sq*vel_k(BandIdx,3)*omg_bnd1(BandIdx)*&
                                physics%elem_charge_SI**4/physics%hbar_SI**3*1.e-10_dp
        TDF1_kz(ZY, loop_f,3) = TDF1_kz(ZY, loop_f,3)+rdum*beta_sq*&
                        vel_k(BandIdx,3)*vel_k(BandIdx,2)*physics%elem_charge_SI**3/physics%hbar_SI**2*1.e10_dp
        TDF2_kz(ZY, loop_f,3) = TDF2_kz(ZY,loop_f,3)+rdum*beta_sq*vel_k(BandIdx,3)*omg_bnd2(BandIdx)*&
                                physics%elem_charge_SI**4/physics%hbar_SI**3*1.e-10_dp
        end if
      end do
    end do !loop over bands

  end subroutine TDFtot_kpt

 
end module w90_nerwann
