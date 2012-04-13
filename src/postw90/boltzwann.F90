!******
! HEADER HERE:
!******
! TODOs:
! * check results for thermal conductivity
! * Check input variables name consistency with DOS module
! * put timings and memory allocations
! * undestand behavior if the number of steps is zero! (as it is: division by zero)
! * add check of inputs somewhere
! * read for TODOs
! * check comments
! * implement rigid bandshift
! * parameters related to smearing (2x, DOS and TDF)
! * in the documentation, say that it also uses the parameters of the ham. interpolation, as
!   for instance those to choose if we want to take into account degeneracies when calculating the
!   velocities
! * check that we get the same results with and without spin-orbit! 
! * document that for spin-degeneracy one has to use spinors=.true. in input of one gets a factor of 2 error

module w90_boltzwann
  use w90_constants
  use w90_parameters, only : &
       boltz_calc_also_dos, boltz_dos_energy_step, boltz_dos_min_energy, boltz_dos_max_energy, &
       boltz_mu_min, boltz_mu_max, boltz_mu_step, boltz_temp_min, boltz_temp_max, boltz_temp_step, &
       boltz_interp_mesh_spacing, boltz_interp_mesh, boltz_tdf_energy_step, &
       boltz_bandshift, boltz_bandshift_firstband, boltz_bandshift_energyshift, &
       timing_level, dis_win_min, dis_win_max 
  use w90_io, only : io_error,stdout,io_stopwatch,io_file_unit,seedname  
  use w90_utility,    only : utility_inv3
  use w90_wanint_common
  implicit none

  private 
  public :: boltzwann_main

contains 

  !> This is the main routine of the BoltzWann module.
  !> It calculates the transport coefficients using the Boltzmann transport equation.
  !>
  !> It produces five files that contain:
  !> 1. the Transport Distribution function (TDF) in units of 1/hbar^2 * eV*ps/angstrom
  !> 2. the electrical conductivity in SI units (1/Ohm/m)
  !> 3. the Seebeck coefficient in SI units (V/K)
  !> 4. the thermal conductivity in SI units (W/meter/K)
  !> 5. if requested, the density of states
  !> Files from 2 to 4 are output on a grid of (mu,T) points, where mu is the chemical potential in eV and
  !> T is the temperature in Kelvin. The grid is defined in the input.
  subroutine boltzwann_main()
    integer :: TempNumPoints, MuNumPoints, TDFEnergyNumPoints
    integer :: i, j, ierr, EnIdx, TempIdx, MuIdx
    real(kind=dp), dimension(:),     allocatable :: TempArray, MuArray, KTArray
    real(kind=dp), dimension(:,:),   allocatable :: TDF ! (coordinate,Energy)
    real(kind=dp), dimension(:),     allocatable :: TDFEnergyArray
    real(kind=dp), dimension(:,:),   allocatable :: IntegrandArray ! (coordinate, Energy) at a given T and mu
    real(kind=dp), dimension(3,3)                :: ElCondTimesSeebeckFP, ThisElCond, ElCondInverse, ThisSeebeck
    real(kind=dp), dimension(6)                  :: ElCondTimesSeebeck
    real(kind=dp), dimension(:,:,:), allocatable :: ElCond ! (coordinate,Temp, mu)
    real(kind=dp), dimension(:,:,:), allocatable :: Seebeck ! (coordinate,Temp, mu)
    real(kind=dp), dimension(:,:,:), allocatable :: ThermCond ! (coordinate,Temp, mu)
    real(kind=dp)                                :: Determinant
    integer :: tdf_unit, elcond_unit, seebeck_unit, boltzdos_unit, thermcond_unit

    if(on_root .and. (timing_level>0)) call io_stopwatch('boltzwann_main',1)
    ! TODO here: (only on node 0): print initial banner

    ! I open the output files
    if (on_root) then
       tdf_unit=io_file_unit()
       open(unit=tdf_unit,file=trim(seedname)//'.tdf')  
       
       elcond_unit = io_file_unit()
       open(unit=elcond_unit,file=trim(seedname)//'.elcond')  
       
       seebeck_unit =io_file_unit()
       open(unit=seebeck_unit,file=trim(seedname)//'.seebeck')  
       
       thermcond_unit =io_file_unit()
       open(unit=thermcond_unit,file=trim(seedname)//'.thermcond')  
       
       if (boltz_calc_also_dos) then
          boltzdos_unit=io_file_unit()
          open(unit=boltzdos_unit,file=trim(seedname)//'.boltzdos')  
       end if
    end if

    ! I precalculate the TempArray and the MuArray
    TempNumPoints = int(floor((boltz_temp_max-boltz_temp_min)/boltz_temp_step))+1
    allocate(TempArray(TempNumPoints),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating TempArray in boltzwann_main')
    do i=1,TempNumPoints
       TempArray(i) = boltz_temp_min + real(i-1,dp)*boltz_temp_step
    end do
    
    ! This array contains the same temperatures of the TempArray, but multiplied by k_boltzmann, in units of eV
    allocate(KTArray(TempNumPoints),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating KTArray in boltzwann_main')
    ! (k_B in eV/kelvin is equal to k_B_SI / elem_charge_SI)
    KTArray = TempArray * k_B_SI / elem_charge_SI

    MuNumPoints = int(floor((boltz_mu_max-boltz_mu_min)/boltz_mu_step))+1
    allocate(MuArray(MuNumPoints),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating MuArray in boltzwann_main')
    do i=1,MuNumPoints
       MuArray(i) = boltz_mu_min + real(i-1,dp)*boltz_mu_step
    end do

    ! I precalculate the TDFEnergyArray
    ! I assume that dis_win_min and dis_win_max are set to sensible values, related to the max and min energy
    ! This is true if the .eig file is present. I can assume its presence since we need it to interpolate the
    ! bands.
    ! TODO: when we put smearing, understand if we want to enlarge the window by some multiple of the smearing
    TDFEnergyNumPoints = int(floor((dis_win_max-dis_win_min)/boltz_tdf_energy_step))+1
    allocate(TDFEnergyArray(TDFEnergyNumPoints),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating TDFEnergyArray in boltzwann_main')
    do i=1,TDFEnergyNumPoints
       TDFEnergyArray(i) = dis_win_min + real(i-1,dp)*boltz_tdf_energy_step
    end do

    ! I allocate the array for the TDF
    allocate(TDF(6,TDFEnergyNumPoints),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating TDF in boltzwann_main')

    ! I allocate the array that I will use to store the functions to be integrated
    allocate(IntegrandArray(6,TDFEnergyNumPoints),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating FermiDerivArray in boltzwann_main')

    ! I call the subroutine that calculates the Transport Distribution Function
    call calcTDF(TDF,TDFEnergyArray)
    ! The TDF array contains now the TDF, or more precisely
    ! hbar^2 * TDF in units of eV * ps / angstrom

    ! I print on file the TDF
    if (on_root) then
       ! TODO: Write more complete header with authors, ...
       write(tdf_unit,'(A)') "# Transport distribution function (in units of 1/hbar^2 * eV * ps / angstrom)" // &
            " vs energy in eV"
       write(tdf_unit,'(A)') "# Energy TDF_xx TDF_xy TDF_yy TDF_xz TDF_yz TDF_zz"
       do i=1,size(TDFEnergyArray)
          write(tdf_unit,101) TDFEnergyArray(i), TDF(:,i)
       end do
    end if

    ! I allocate the arrays for the spectra
    allocate(ElCond(6,TempNumPoints,MuNumPoints),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating ElCond in boltzwann_main')
    allocate(Seebeck(6,TempNumPoints,MuNumPoints),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating Seebeck in boltzwann_main')
    allocate(ThermCond(6,TempNumPoints,MuNumPoints),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating Seebeck in boltzwann_main')
    ElCond = 0._dp
    Seebeck = 0._dp
    ThermCond = 0._dp

    ! Now, I calculate the various spectra for all mu and T values
    do MuIdx=1,MuNumPoints
!       write(0,*) MuIdx, '/', MuNumPoints, TDFEnergyNumPoints
       do TempIdx = 1,TempNumPoints
          IntegrandArray = TDF
          do EnIdx=1,TDFEnergyNumPoints
             IntegrandArray(:,EnIdx) = IntegrandArray(:,EnIdx) * &
                  MinusFermiDerivative(E=TDFEnergyArray(EnIdx),mu=MuArray(MuIdx),KT=KTArray(TempIdx))
          end do
          ! Now, IntegrandArray contains (-dn/dE) * TDF_ij(E), where n is the Fermi distribution function
          ! Its integral is ElCond_ij/e^2
          ElCond(:,TempIdx,MuIdx) = sum(IntegrandArray,DIM=2)*boltz_tdf_energy_step
          ! ElCond contains now (hbar^2/e^2) * sigma in eV*ps/angstrom, where sigma is the conductivity tensor 
          ! (note that MinusFermiDerivative is in units of 1/eV, so that when I perform the integration
          ! ElCond has the same units of TDF)
          
          ! I store in ThisElCond the conductivity tensor in standard format
          do j=1,3
             do i=1,j
                ThisElCond(i,j)=ElCond(i+((j-1)*j)/2,TempIdx,MuIdx)
             end do
          end do
          ! I calculate the inverse matrix of the conductivity 
          call utility_inv3(ThisElCond,ElCondInverse,Determinant)
          ! TODO: check if Determinant == 0!
          ! The routine returns the adjoint of ThisElCond; in order to get 
          ! the inverse, I have to calculate ElCondInverse / Determinant
          ElCondInverse = ElCondInverse / Determinant

          ! Now, I multiply IntegrandArray by (E-mu): then, IntegrandArray contains
          ! (-dn/dE) * TDF_ij(E) * (E-mu) and its integral is (ElCond*Seebeck)_ij * T / e, where
          ! T is the temperature
          do EnIdx=1,TDFEnergyNumPoints
             IntegrandArray(:,EnIdx) = IntegrandArray(:,EnIdx) * (TDFEnergyArray(EnIdx) -MuArray(MuIdx))
          end do

          ! I store in ElCondTimesSeebeckFP the product of the two tensors in full-packed format
          ElCondTimesSeebeck = sum(IntegrandArray,DIM=2)*boltz_tdf_energy_step / TempArray(TempIdx)
          do j=1,3
             do i=1,j
                ElCondTimesSeebeckFP(i,j)=ElCondTimesSeebeck(i+((j-1)*j)/2)
             end do
          end do
          ! Now, ElCondTimesSeebeck (and ElCondTimesSeebeckFP) contain
          ! [ElCond*Seebeck] * hbar^2/e in units of eV^2*ps/angstrom/kelvin

          ! I calculate ElCond^(-1) . (ElCondTimesSeebeck) = Seebeck in standard format, and then
          ! store it in the Seebeck array (in packed format)
          ThisSeebeck = matmul(ElCondInverse,ElCondTimesSeebeckFP) 
          do j=1,3
             do i=1,j
                Seebeck(i+((j-1)*j)/2,TempIdx,MuIdx) = ThisSeebeck(i,j)
             end do
          end do
          ! Now, Seebeck contains the Seebeck coefficient in volt / Kelvin. In fact:
          ! - ElCond contains (hbar^2/e^2) * sigma in eV*ps/angstrom, where sigma is the value of the
          !   conductivity, i.e. it is the conductivity in units of (e^2/hbar^2) * eV * ps / angstrom
          ! - ElCondInverse is in thus in units of (hbar^2/e^2) / eV / ps * angstrom
          ! - ElCondTimesSeebeck is in units of  e / hbar^2 * eV^2 * ps / angstrom / Kelvin
          ! therefore ThisSeebeck, which has the units of ElCondInverse * ElCondTimesSeebeck, i.e.
          ! [(hbar^2/e^2) / eV / ps * angstrom] * [ e / hbar^2 * eV^2 * ps / angstrom / kelvin] = 
          ! eV/e/kelvin = volt/kelvin
          
          ! Now, I multiply IntegrandArray by (E-mu): then, IntegrandArray contains
          ! (-dn/dE) * TDF_ij(E) * (E-mu)^2 and its integral is (ThermCond)_ij * T
          do EnIdx=1,TDFEnergyNumPoints
             IntegrandArray(:,EnIdx) = IntegrandArray(:,EnIdx) * (TDFEnergyArray(EnIdx) -MuArray(MuIdx))
          end do
          ThermCond(:,TempIdx,MuIdx) = sum(IntegrandArray,DIM=2)*boltz_tdf_energy_step / TempArray(TempIdx)
          ! ThermCond contains now the thermal conductivity in units of
          ! 1/hbar^2 * eV^3*ps/angstrom/kelvin
       end do
    end do
    ! Now, I multiply by the correct factors to obtain the tensors in SI units
    
    ! **** Electrical conductity ****
    ! Now, ElCond is in units of (e^2/hbar^2) * eV * ps / angstrom
    ! I want the conductivity in units of 1/Ohm/meter (i.e., SI units). The conversion factor is then
    ! e^2/hbar^2 * eV*ps/angstrom * Ohm * meter; we use the constants in constants.F90 by noting that:
    ! Ohm = V/A = V*s/C, where A=ampere, C=coulomb
    ! Then we have: (e/C)^3/hbar^2 * V*s^2 * V * C^2 * (meter / angstrom)  * (ps / s)
    ! Now: e/C = elem_charge_SI; CV=Joule, CV/s=Watt, hbar/Watt = hbar_SI;
    ! moreover meter / angstrom = 1e10, ps / s = 1e-12 so that we finally get the
    ! CONVERSION FACTOR: elem_charge_SI**3 / (hbar_SI**2) * 1.e-2_dp
    ElCond = ElCond * elem_charge_SI**3 / (hbar_SI**2) * 1.e-2_dp
    ! THIS IS NOW THE ELECTRICAL CONDUCTIVITY IN SI UNITS, i.e. in 1/Ohm/meter

    ! **** Seebeck coefficient ****
    ! THE SEEECK COEFFICIENTS IS ALREADY IN volt/kelvin, so nothing has to be done

    ! **** Thermal conductivity ****
    ! Now, ThermCond is in units of 1/hbar^2 * eV^3*ps/angstrom/K
    ! I want it to be in units of W/m/K (i.e., SI units). Then conversion factor is then
    ! 1/hbar^2 * eV^3*ps/angstrom/K / W * m * K =
    ! 1/hbar^2 * eV^3 * ps / W * (m/angstrom) =
    ! 1/hbar^2 * C^3 * V^3 * s / W * [e/C]^3 * (m/angstrom) * (ps / s) =
    ! 1/hbar^2 * J^2 * [e/C]^3 * (m/angstrom) * (ps / s) =
    ! elem_charge_SI**3 / (hbar_SI**2) * 1.e-2_dp, i.e. the same conversion factor as above
    ThermCond = ThermCond * elem_charge_SI**3 / (hbar_SI**2) * 1.e-2_dp
    ! THIS IS NOW THE THERMAL CONDUCTIVITY IN SI UNITS, i.e. in W/meter/K
    
    ! TODO: printing here!
      
    ! Before ending, I close files and deallocate memory
    if (on_root) then
       close(tdf_unit)
       close(elcond_unit)
       close(seebeck_unit)
       if (boltz_calc_also_dos) close(boltzdos_unit)
    end if

    deallocate(TempArray,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating TempArray in boltzwann_main')
    deallocate(KTArray,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating KTArray in boltzwann_main')
    deallocate(MuArray,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating MuArray in boltzwann_main')
    deallocate(TDFEnergyArray,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating TDFEnergyArray in boltzwann_main')
    deallocate(TDF,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating TDF in boltzwann_main')
    deallocate(ElCond,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating ElCond in boltzwann_main')
    deallocate(Seebeck,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating Seebeck in boltzwann_main')
    deallocate(ThermCond,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating ThermCond in boltzwann_main')
    deallocate(IntegrandArray,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating IntegrandArray in boltzwann_main')

    if(on_root .and. (timing_level>0)) call io_stopwatch('boltzwann_main',2)

101 FORMAT(7G18.10)

  end subroutine boltzwann_main

  !> This routine calculates the Transport Distribution Function \f$ \Sigma_{ij}(\eps) \f$ (TDF)
  !> in units of 1/hbar^2 * eV*ps/angstrom.
  !> 
  !> \note The TDFEnergyArray must be already allocated and initialized with the energies in eV before calling
  !>       this routine
  !> 
  !> \note The TDF array must be already allocated with dimensions 6*N before calling this routine, where
  !>       N is the dimension of the TDFEnergyArray array. 
  !>
  !> \note We assume that the TDFEnergyArray is uniformely spaced (no checks are performed on this).
  !> 
  !> \param TDF The TDF(i,EnIdx) output array, where:
  !>  - i is an index from 1 to 6 giving the component of the symmetric tensor \f$ Sigma_{ij}(\eps) \f$,
  !>    where 1=xx, 2=xy, 3=yy, 4=xz, 5=yz, 6=zz
  !>    (in this way the mapping (i,j) -> i+((j-1)*j)/2 is satisfied for the packed storage of the
  !>    upper triangle [i<=j]).
  !>  - EnIdx is the index of the energies; the corresponding energy is given by
  !>    TDFEnergyArray(EndIdx) array (in eV).
  !> \param TDFEnergyArray The array with the energies for which the TDF is calculated, in eV
  subroutine calcTDF(TDF,TDFEnergyArray)
    use w90_get_oper, only: get_HH_R, HH_R
    use w90_parameters, only    : num_wann
    use w90_utility, only : utility_diagonalize
    use w90_wanint_common, only: fourier_R_to_k
    use w90_wan_ham, only: get_deleig_a

    real(kind=dp), dimension(:,:), intent(out)   :: TDF ! (coordinate,Energy)
    real(kind=dp), dimension(:), intent(in)      :: TDFEnergyArray
    ! Comments:
    ! issue warnings if going outside of the energy window
    ! check that we actually get hbar*velocity in eV*angstrom
    ! TODO: should I also calc the DOS?
    ! comments on smearing
    ! PARTODO: comments on which processor the TDF is actually stored at the end

    real(kind=dp), dimension(3) :: kpt
    integer :: loop_tot, loop_x, loop_y, loop_z, ierr

    complex(kind=dp), allocatable :: HH(:,:)
    complex(kind=dp), allocatable :: delHH(:,:,:)
    complex(kind=dp), allocatable :: UU(:,:)
    real(kind=dp) :: del_eig(num_wann,3)
    real(kind=dp) :: eig(num_wann)

    if(on_root .and. (timing_level>0)) call io_stopwatch('calcTDF',1)

    ! Some initial checks
    if (size(TDF,1)/=6 .or. size(TDF,2)/=size(TDFEnergyArray)) then
       call io_error('Wrong size for the TDF array in calcTDF')
    end if

    ! I zero the TDF array before starting
    TDF = 0._dp
    
    ! I call once the routine to calculate the Hamiltonian in real-space <0n|H|Rm>
    ! TODO: to be checked: if I have to check for the presence of the .eig file
    call get_HH_R()

    allocate(HH(num_wann,num_wann),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating HH in calcTDF')
    allocate(delHH(num_wann,num_wann,3),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating delHH in calcTDF')
    allocate(UU(num_wann,num_wann),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating UU in calcTDF')    

    ! TODO: print on output the grid that I am using (and min_ksep value, if boltz_interp_mesh_spacing >0)
    
    print*, my_node_id, boltz_interp_mesh, boltz_interp_mesh_spacing

    ! I loop over all kpoints
    do loop_tot=0,PRODUCT(boltz_interp_mesh)-1

       ! I get the coordinates for the x,y,z components starting from a single loop variable
       ! (which is better for parallelization purposes)
       ! Important! This works only if loop_tot starts from ZERO and ends with 
       !            PRODUCT(boltz_interp_mesh)-1, so be careful when parallelizing
       loop_x=loop_tot/(boltz_interp_mesh(2)*boltz_interp_mesh(3))
       loop_y=(loop_tot-loop_x*(boltz_interp_mesh(2)*boltz_interp_mesh(3)))/boltz_interp_mesh(3)
       loop_z=loop_tot-loop_x*(boltz_interp_mesh(2)*boltz_interp_mesh(3))-loop_y*boltz_interp_mesh(3)

       ! kpt(i) is in in the [0,d-1]/d range, with d=boltz_interp_mesh(i)
       kpt(1)=(real(loop_x,dp)/real(boltz_interp_mesh(1),dp))
       kpt(2)=(real(loop_y,dp)/real(boltz_interp_mesh(2),dp))
       kpt(3)=(real(loop_z,dp)/real(boltz_interp_mesh(3),dp))
       
       if (boltz_calc_also_dos) then
          !TODO!!
       end if
       
       ! Here I have to get the band energies and the velocities
       call fourier_R_to_k(kpt,HH_R,HH,0) 
       call utility_diagonalize(HH,num_wann,eig,UU) 
       call fourier_R_to_k(kpt,HH_R,delHH(:,:,1),1) 
       call fourier_R_to_k(kpt,HH_R,delHH(:,:,2),2) 
       call fourier_R_to_k(kpt,HH_R,delHH(:,:,3),3) 
       call get_deleig_a(del_eig(:,1),eig,delHH(:,:,1),UU)
       call get_deleig_a(del_eig(:,2),eig,delHH(:,:,2),UU)
       call get_deleig_a(del_eig(:,3),eig,delHH(:,:,3),UU)
       

!!$       do BandIdx=1, num_wann
!!$          ! Maybe we should check once and for all if it makes sense to broaden or if the broadening
!!$          ! is too small!
!!$          MinEn=eig(BandIdx)-constant*broadening
!!$          MaxEn=eig(BandIdx)+constant*broadening
!!$          ! We assume here that the TDFEnergyArray is uniformely spaced
!!$          MinEnIdx = int((MinEn-TDFEnergyArray(1))/(TDFEnergyArray(size(TDFEnergyArray))-TDFEnergyArray(1)) * &
!!$               (size(TDFEnergyArray)-1)) + 1
!!$          MaxEnIdx = int((MaxEn-TDFEnergyArray(1))/(TDFEnergyArray(size(TDFEnergyArray))-TDFEnergyArray(1)) * &
!!$               (size(TDFEnergyArray)-1)) + 1
!!$          MinEnIdx = max(1,MinEnIdx)
!!$          MaxEnIdx = min(size(TDFEnergyArray),MaxEnIdx)
!!$          do EnIdx=MinEnIdx, MaxEnIdx
!!$             ! TODO: Factors must contain: tau; spin-degeneracy; dos in k space; ...
!!$             ! for spin-degeneracy: use spinors=.true. in input and elec_per_state variable
!!$             ! TDF_xx
!!$             TDF(EnIdx,1) = TDF(EndIdx,1) + Factors * & 
!!$                  ! normGauss gives the value of a normalized gaussian of proper width with center the origin
!!$                  normalizedGaussian(TDFEnergyArray(EnIdx)-eig(BandIdx)) * &
!!$                  del_eig(BandIdx, 1) * del_eig (BandIdx, 1)  
!!$          ! TODO: repeat the same for the other 5 components
!!$       end do
!!$

    end do

    if(on_root .and. (timing_level>0)) call io_stopwatch('calcTDF',2)

    deallocate(HH,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating HH in calcTDF')
    deallocate(delHH,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating delHH in calcTDF')
    deallocate(UU,stat=ierr)
    if (ierr/=0) call io_error('Error in deallocating UU in calcTDF')

  end subroutine calcTDF

  !> This function calculates -dn(E)/dE, where n(E) is the Fermi distribution function.
  !>
  !> \note I do not put stopwatches here because it would slow down the calculation by orders of magnitude
  !> 
  !> \param E  Energy at which we want to calculate -dn(E)/dE, in eV
  !> \param mu Chemical potential in eV
  !> \param KT k_Boltzmann * Temperature in eV
  function MinusFermiDerivative(E,mu,KT)
    real(kind=dp), intent(in) :: E
    real(kind=dp), intent(in) :: mu
    real(kind=dp), intent(in) :: KT
    real(kind=dp) :: MinusFermiDerivative

    ! This is used as the maximum value to be used for the exp function. 
    ! exp(-x) for x>MaxExp is assumed to be zero. This value si chose so that
    ! the function is truncated when its value is smaller than about 1.e-16
    real(kind=dp), parameter :: MaxExp = 36._dp
    real(kind=dp) :: MyExp

    MyExp = (E - mu) / KT
    if (abs(MyExp) > MaxExp) then
       MinusFermiDerivative = 0._dp
    else
       MinusFermiDerivative = 1/KT * exp(MyExp)/((exp(MyExp)+1)**2)
    end if

  end function MinusFermiDerivative

!!$  ! This happens to be slower
!!$  function MinusFermiDerivative2(EArray,mu,KT)
!!$    real(kind=dp), dimension(:), intent(in) :: EArray
!!$    real(kind=dp), intent(in) :: mu
!!$    real(kind=dp), intent(in) :: KT
!!$    real(kind=dp), dimension(size(EArray)) :: MinusFermiDerivative2
!!$
!!$    ! This is used as the maximum value to be used for the exp function. 
!!$    ! exp(-x) for x>MaxExp is assumed to be zero. This value si chose so that
!!$    ! the function is truncated when its value is smaller than about 1.e-16
!!$    real(kind=dp), parameter :: MaxExp = 36._dp
!!$    real(kind=dp), dimension(size(EArray)) :: MyExpArray
!!$    
!!$    MyExpArray = (EArray - mu) / KT
!!$
!!$    MinusFermiDerivative2 = 0._dp
!!$    where (abs(MyExpArray) < MaxExp) &
!!$         MinusFermiDerivative2 = 1/KT * exp(MyExpArray)/((exp(MyExpArray)+1)**2)
!!$
!!$  end function MinusFermiDerivative2

end module w90_boltzwann
