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

module w90_dos
  !! Compute Density of States
  use w90_constants, only: dp

  implicit none

  private

  public :: dos_main, dos_get_levelspacing, dos_get_k

  integer       :: num_freq
  !! Number of sampling points
  real(kind=dp) :: d_omega
  !! Step between energies

contains

  !=========================================================!
  !                   PUBLIC PROCEDURES                     !
  !=========================================================!

  subroutine dos_main
    !=======================================================!
    !                                                       !
    !! Computes the electronic density of states. Can
    !! resolve into up-spin and down-spin parts, project
    !! onto selected Wannier orbitals, and use adaptive
    !! broadening, as in PRB 75, 195121 (2007) [YWVS07].
    !                                                       !
    !=======================================================!

    use w90_io, only: io_error, io_file_unit, io_date, io_stopwatch, &
      seedname, stdout
    use w90_comms, only: on_root, num_nodes, my_node_id, comms_reduce
    use w90_postw90_common, only: num_int_kpts_on_node, int_kpts, weight, &
      pw90common_fourier_R_to_k
    use w90_parameters, only: num_wann, dos_energy_min, dos_energy_max, &
      dos_energy_step, timing_level, &
      wanint_kpoint_file, dos_kmesh, &
      dos_smr_index, dos_adpt_smr, dos_adpt_smr_fac, &
      dos_adpt_smr_max, spin_decomp, &
      dos_smr_fixed_en_width, &
      dos_project, num_dos_project
    use w90_get_oper, only: get_HH_R, get_SS_R, HH_R
    use w90_wan_ham, only: wham_get_eig_deleig
    use w90_utility, only: utility_diagonalize

    ! 'dos_k' contains contrib. from one k-point,
    ! 'dos_all' from all nodes/k-points (first summed on one node and
    ! then reduced (i.e. summed) over all nodes)
    !
    real(kind=dp), allocatable :: dos_k(:, :)
    real(kind=dp), allocatable :: dos_all(:, :)

    real(kind=dp)    :: kweight, kpt(3), omega
    integer          :: i, loop_x, loop_y, loop_z, loop_tot, ifreq
    integer          :: dos_unit, ndim, ierr
    real(kind=dp), dimension(:), allocatable :: dos_energyarray

    complex(kind=dp), allocatable :: HH(:, :)
    complex(kind=dp), allocatable :: delHH(:, :, :)
    complex(kind=dp), allocatable :: UU(:, :)
    real(kind=dp) :: del_eig(num_wann, 3)
    real(kind=dp) :: eig(num_wann), levelspacing_k(num_wann)

    num_freq = nint((dos_energy_max - dos_energy_min)/dos_energy_step) + 1
    if (num_freq == 1) num_freq = 2
    d_omega = (dos_energy_max - dos_energy_min)/(num_freq - 1)

    allocate (dos_energyarray(num_freq), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating dos_energyarray in ' &
                                 //'dos subroutine')
    do ifreq = 1, num_freq
      dos_energyarray(ifreq) = dos_energy_min + real(ifreq - 1, dp)*d_omega
    end do

    allocate (HH(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating HH in dos')
    allocate (delHH(num_wann, num_wann, 3), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating delHH in dos')
    allocate (UU(num_wann, num_wann), stat=ierr)
    if (ierr /= 0) call io_error('Error in allocating UU in dos')

    call get_HH_R
    if (spin_decomp) then
      ndim = 3
      call get_SS_R
    else
      ndim = 1
    end if

    allocate (dos_k(num_freq, ndim))
    allocate (dos_all(num_freq, ndim))

    if (on_root) then

      if (timing_level > 1) call io_stopwatch('dos', 1)

!       write(stdout,'(/,1x,a)') '============'
!       write(stdout,'(1x,a)')   'Calculating:'
!       write(stdout,'(1x,a)')   '============'

      write (stdout, '(/,/,1x,a)') &
        'Properties calculated in module  d o s'
      write (stdout, '(1x,a)') &
        '--------------------------------------'

      if (num_dos_project == num_wann) then
        write (stdout, '(/,3x,a)') '* Total density of states (_dos)'
      else
        write (stdout, '(/,3x,a)') &
          '* Density of states projected onto selected WFs (_dos)'
        write (stdout, '(3x,a)') 'Selected WFs |Rn> are:'
        do i = 1, num_dos_project
          write (stdout, '(5x,a,2x,i3)') 'n =', dos_project(i)
        enddo
      endif

      write (stdout, '(/,5x,a,f9.4,a,f9.4,a)') &
        'Energy range: [', dos_energy_min, ',', dos_energy_max, '] eV'

      write (stdout, '(/,5x,a,(f6.3,1x))') &
        'Adaptive smearing width prefactor: ', &
        dos_adpt_smr_fac

      write (stdout, '(/,/,1x,a20,3(i0,1x))') 'Interpolation grid: ', &
        dos_kmesh(1:3)

    end if

    dos_all = 0.0_dp

    if (wanint_kpoint_file) then
      !
      ! Unlike for optical properties, this should always work for the DOS
      !
      if (on_root) write (stdout, '(/,1x,a)') 'Sampling the irreducible BZ only'

      ! Loop over k-points on the irreducible wedge of the Brillouin zone,
      ! read from file 'kpoint.dat'
      !
      do loop_tot = 1, num_int_kpts_on_node(my_node_id)
        kpt(:) = int_kpts(:, loop_tot)
        if (dos_adpt_smr) then
          call wham_get_eig_deleig(kpt, eig, del_eig, HH, delHH, UU)
          call dos_get_levelspacing(del_eig, dos_kmesh, levelspacing_k)
          call dos_get_k(kpt, dos_energyarray, eig, dos_k, &
                         smr_index=dos_smr_index, &
                         adpt_smr_fac=dos_adpt_smr_fac, &
                         adpt_smr_max=dos_adpt_smr_max, &
                         levelspacing_k=levelspacing_k, &
                         UU=UU)
        else
          call pw90common_fourier_R_to_k(kpt, HH_R, HH, 0)
          call utility_diagonalize(HH, num_wann, eig, UU)
          call dos_get_k(kpt, dos_energyarray, eig, dos_k, &
                         smr_index=dos_smr_index, &
                         smr_fixed_en_width=dos_smr_fixed_en_width, &
                         UU=UU)
        end if
        dos_all = dos_all + dos_k*weight(loop_tot)
      end do

    else

      if (on_root) write (stdout, '(/,1x,a)') 'Sampling the full BZ'

      kweight = 1.0_dp/real(PRODUCT(dos_kmesh), kind=dp)
      do loop_tot = my_node_id, PRODUCT(dos_kmesh) - 1, num_nodes
        loop_x = loop_tot/(dos_kmesh(2)*dos_kmesh(3))
        loop_y = (loop_tot - loop_x*(dos_kmesh(2) &
                                     *dos_kmesh(3)))/dos_kmesh(3)
        loop_z = loop_tot - loop_x*(dos_kmesh(2)*dos_kmesh(3)) &
                 - loop_y*dos_kmesh(3)
        kpt(1) = real(loop_x, dp)/real(dos_kmesh(1), dp)
        kpt(2) = real(loop_y, dp)/real(dos_kmesh(2), dp)
        kpt(3) = real(loop_z, dp)/real(dos_kmesh(3), dp)
        if (dos_adpt_smr) then
          call wham_get_eig_deleig(kpt, eig, del_eig, HH, delHH, UU)
          call dos_get_levelspacing(del_eig, dos_kmesh, levelspacing_k)
          call dos_get_k(kpt, dos_energyarray, eig, dos_k, &
                         smr_index=dos_smr_index, &
                         adpt_smr_fac=dos_adpt_smr_fac, &
                         adpt_smr_max=dos_adpt_smr_max, &
                         levelspacing_k=levelspacing_k, &
                         UU=UU)
        else
          call pw90common_fourier_R_to_k(kpt, HH_R, HH, 0)
          call utility_diagonalize(HH, num_wann, eig, UU)
          call dos_get_k(kpt, dos_energyarray, eig, dos_k, &
                         smr_index=dos_smr_index, &
                         smr_fixed_en_width=dos_smr_fixed_en_width, &
                         UU=UU)
        end if
        dos_all = dos_all + dos_k*kweight
      end do

    end if

    ! Collect contributions from all nodes
    !
    call comms_reduce(dos_all(1, 1), num_freq*ndim, 'SUM')

    if (on_root) then
      write (stdout, '(1x,a)') 'Output data files:'
      write (stdout, '(/,3x,a)') trim(seedname)//'-dos.dat'
      dos_unit = io_file_unit()
      open (dos_unit, FILE=trim(seedname)//'-dos.dat', STATUS='UNKNOWN', &
            FORM='FORMATTED')
      do ifreq = 1, num_freq
        omega = dos_energyarray(ifreq)
        write (dos_unit, '(4E16.8)') omega, dos_all(ifreq, :)
      enddo
      close (dos_unit)
      if (timing_level > 1) call io_stopwatch('dos', 2)
    end if

    deallocate (HH, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating HH in dos_main')
    deallocate (delHH, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating delHH in dos_main')
    deallocate (UU, stat=ierr)
    if (ierr /= 0) call io_error('Error in deallocating UU in dos_main')

  end subroutine dos_main

  ! =========================================================================

! The next routine is commented. It should be working (apart for a
! missing broadcast at the very end, see comments there).  However,
! it should be debugged, and probably the best thing is to avoid to
! resample the BZ, but rather use the calculated DOS (maybe it can be
! something that is done at the end of the DOS routine?)
!~  subroutine find_fermi_level
!~    !==============================================!
!~    !                                              !
!~    ! Finds the Fermi level by integrating the DOS !
!~    !                                              !
!~    !==============================================!
!~
!~    use w90_io, only            : stdout,io_error
!~    use w90_comms
!~    use w90_postw90_common, only : max_int_kpts_on_node,num_int_kpts_on_node,&
!~         int_kpts,weight
!~    use w90_parameters, only    : nfermi,fermi_energy_list,&
!~         num_valence_bands,&
!~         num_wann,dos_num_points,dos_min_energy,&
!~         dos_max_energy,dos_energy_step,&
!~         wanint_kpoint_file
!~
!~#ifdef MPI
!~    include 'mpif.h'
!~#endif
!~
!~    real(kind=dp) :: kpt(3),sum_max_node,sum_max_all,&
!~         sum_mid_node,sum_mid_all,emin,emax,emid,&
!~         emin_node(0:num_nodes-1),emax_node(0:num_nodes-1),&
!~         ef
!~    integer       :: loop_x,loop_y,loop_z,loop_kpt,loop_nodes,&
!~         loop_iter,ierr,num_int_kpts,ikp
!~
!~    real(kind=dp), allocatable :: eig_node(:,:)
!~    real(kind=dp), allocatable :: levelspacing_node(:,:)
!~
!~    if(on_root) write(stdout,'(/,a)') 'Finding the value of the Fermi level'
!~
!~    if(.not.wanint_kpoint_file) then
!~       !
!~       ! Already done in pw90common_wanint_get_kpoint_file if
!~       ! wanint_kpoint_file=.true.
!~       !
!~       allocate(num_int_kpts_on_node(0:num_nodes-1))
!~       num_int_kpts=dos_num_points**3
!~       !
!~       ! Local k-point counter on each node (lazy way of doing it, there is
!~       ! probably a smarter way)
!~       !
!~       ikp=0
!~       do loop_kpt=my_node_id,num_int_kpts-1,num_nodes
!~          ikp=ikp+1
!~       end do
!~       num_int_kpts_on_node(my_node_id)=ikp
!~#ifdef MPI
!~       call MPI_reduce(ikp,max_int_kpts_on_node,1,MPI_integer,&
!~            MPI_MAX,0,MPI_COMM_WORLD,ierr)
!~#else
!~       max_int_kpts_on_node=ikp
!~#endif
!~       call comms_bcast(max_int_kpts_on_node,1)
!~    end if
!~
!~    allocate(eig_node(num_wann,max_int_kpts_on_node),stat=ierr)
!~    if (ierr/=0)&
!~         call io_error('Error in allocating eig_node in find_fermi_level')
!~    eig_node=0.0_dp
!~    allocate(levelspacing_node(num_wann,max_int_kpts_on_node),stat=ierr)
!~    if (ierr/=0)&
!~         call io_error('Error in allocating levelspacing in find_fermi_level')
!~    levelspacing_node=0.0_dp
!~
!~    if(wanint_kpoint_file) then
!~       if(on_root) write(stdout,'(/,1x,a)') 'Sampling the irreducible BZ only'
!~       do loop_kpt=1,num_int_kpts_on_node(my_node_id)
!~          kpt(:)=int_kpts(:,loop_kpt)
!~          call dos_get_eig_levelspacing_k(kpt,eig_node(:,loop_kpt),&
!~               levelspacing_node(:,loop_kpt))
!~       end do
!~    else
!~       if (on_root)&
!~            write(stdout,'(/,1x,a)') 'Sampling the full BZ (not using symmetry)'
!~       allocate(weight(max_int_kpts_on_node),stat=ierr)
!~       if (ierr/=0)&
!~            call io_error('Error in allocating weight in find_fermi_level')
!~       weight=0.0_dp
!~       ikp=0
!~       do loop_kpt=my_node_id,num_int_kpts-1,num_nodes
!~          ikp=ikp+1
!~          loop_x=loop_kpt/dos_num_points**2
!~          loop_y=(loop_kpt-loop_x*dos_num_points**2)/dos_num_points
!~          loop_z=loop_kpt-loop_x*dos_num_points**2-loop_y*dos_num_points
!~          kpt(1)=real(loop_x,dp)/dos_num_points
!~          kpt(2)=real(loop_y,dp)/dos_num_points
!~          kpt(3)=real(loop_z,dp)/dos_num_points
!~          weight(ikp)=1.0_dp/dos_num_points**3
!~          call dos_get_eig_levelspacing_k(kpt,eig_node(:,ikp),&
!~               levelspacing_node(:,ikp))
!~       end do
!~    end if
!~
!~    ! Find minimum and maximum band energies within projected subspace
!~    !
!~    emin_node(my_node_id)=&
!~         minval(eig_node(1,1:num_int_kpts_on_node(my_node_id)))
!~    emax_node(my_node_id)=&
!~         maxval(eig_node(num_wann,1:num_int_kpts_on_node(my_node_id)))
!~    if(.not.on_root) then
!~       call comms_send(emin_node(my_node_id),1,root_id)
!~       call comms_send(emax_node(my_node_id),1,root_id)
!~    else
!~       do loop_nodes=1,num_nodes-1
!~          call comms_recv(emin_node(loop_nodes),1,loop_nodes)
!~          call comms_recv(emax_node(loop_nodes),1,loop_nodes)
!~       end do
!~       emin=minval(emin_node)
!~       emax=maxval(emax_node)
!~    end if
!~    call comms_bcast(emin,1)
!~    call comms_bcast(emax,1)
!~
!~    ! Check that the Fermi level lies within the projected subspace
!~    !
!~    sum_max_node=count_states(emax,eig_node,levelspacing_node,&
!~         num_int_kpts_on_node(my_node_id))
!~#ifdef MPI
!~    call MPI_reduce(sum_max_node,sum_max_all,1,MPI_DOUBLE_PRECISION,&
!~         MPI_SUM,0,MPI_COMM_WORLD,ierr)
!~#else
!~    sum_max_all=sum_max_node
!~#endif
!~    if(on_root) then
!~       if(num_valence_bands>sum_max_all) then
!~          write(stdout,*) 'Something wrong in find_fermi_level:'
!~          write(stdout,*)&
!~               '   Fermi level does not lie within projected subspace'
!~          write(stdout,*) 'num_valence_bands= ',num_valence_bands
!~          write(stdout,*) 'sum_max_all= ',sum_max_all
!~          stop 'Stopped: see output file'
!~       end if
!~    end if
!~
!~    ! Now interval search for the Fermi level
!~    !
!~    do loop_iter=1,1000
!~       emid=(emin+emax)/2.0_dp
!~       sum_mid_node=count_states(emid,eig_node,levelspacing_node,&
!~            num_int_kpts_on_node(my_node_id))
!~#ifdef MPI
!~       call MPI_reduce(sum_mid_node,sum_mid_all,1,MPI_DOUBLE_PRECISION,&
!~            MPI_SUM,0,MPI_COMM_WORLD,ierr)
!~#else
!~       sum_mid_all=sum_mid_node
!~#endif
!~       ! This is needed because MPI_reduce only returns sum_mid_all to the
!~       ! root (To understand: could we use MPI_Allreduce instead?)
!~       !
!~       call comms_bcast(sum_mid_all,1)
!~       if(abs(sum_mid_all-num_valence_bands) < 1.e-10_dp) then
!~          !
!~          ! NOTE: Here should assign a value to an entry in a fermi-level
!~          !       vector. Then at the end average over adaptive smearing
!~          !       widths and broadcast the result
!~          !
!~          ef=emid
!~          exit
!~       elseif((sum_mid_all-num_valence_bands) < -1.e-10_dp) then
!~          emin=emid
!~       else
!~          emax=emid
!~       end if
!~    end do
!~
!~    nfermi=1
!~    allocate(fermi_energy_list(1))
!~    fermi_energy_list(1)=ef
!~    !~~ PROBABLY HERE YOU MAY WANT TO BROADCAST THE ABOVE TWO VARIABLES~!
!~    if(on_root) then
!~       write(stdout,*) ' '
!~       write(stdout,'(1x,a,f10.6,a)')&
!~            'Fermi energy = ',ef, ' eV'
!~       write(stdout,'(1x,a)')&
!~            '---------------------------------------------------------'
!~    end if
!~
!~  end subroutine find_fermi_level

  !> This subroutine calculates the contribution to the DOS of a single k point
  !>
  !> \todo still to do: adapt spin_get_nk to read in input the UU rotation matrix
  !>
  !> \note This routine simply provides the dos contribution of a given
  !>       point. This must be externally summed after proper weighting.
  !>       The weight factor (for a full BZ sampling with N^3 points) is 1/N^3 if we want
  !>       the final DOS to be normalized to the total number of electrons.
  !> \note The only factor that is included INSIDE this routine is the spin degeneracy
  !>       factor (=num_elec_per_state variable)
  !> \note The EnergyArray is assumed to be evenly spaced (and the energy spacing
  !>       is taken from EnergyArray(2)-EnergyArray(1))
  !> \note The routine is assuming that EnergyArray has at least two elements.
  !> \note The dos_k array must have dimensions size(EnergyArray) * ndim, where
  !>       ndim=1 if spin_decomp==false, or ndim=3 if spin_decomp==true. This is not checked.
  !> \note If smearing/binwidth < min_smearing_binwidth_ratio,
  !>       no smearing is applied (for that k point)
  !>
  !> \param kpt         the three coordinates of the k point vector whose DOS contribution we
  !>                    want to calculate (in relative coordinates)
  !> \param EnergyArray array with the energy grid on which to calculate the DOS (in eV)
  !>                    It must have at least two elements
  !> \param eig_k       array with the eigenvalues at the given k point (in eV)
  !> \param dos_k       array in which the contribution is stored. Three dimensions:
  !>                    dos_k(energyidx, spinidx), where:
  !>                    - energyidx is the index of the energies, corresponding to the one
  !>                      of the EnergyArray array;
  !>                    - spinidx=1 contains the total dos; if if spin_decomp==.true., then
  !>                      spinidx=2 and spinidx=3 contain the spin-up and spin-down contributions to the DOS
  !> \param smr_index  index that tells the kind of smearing
  !> \param smr_fixed_en_width optional parameter with the fixed energy for smearing, in eV. Can be provided only if the
  !>                    levelspacing_k parameter is NOT given
  !> \param adpt_smr_fac optional parameter with the factor for the adaptive smearing. Can be provided only if the
  !>                    levelspacing_k parameter IS given
  !> \param levelspacing_k optional array with the level spacings, i.e. how much each level changes
  !>                    near a given point of the interpolation mesh, as given by the
  !>                    dos_get_levelspacing() routine
  !>                    If present: adaptive smearing
  !>                    If not present: fixed-energy-width smearing
  subroutine dos_get_k(kpt, EnergyArray, eig_k, dos_k, smr_index, &
                       smr_fixed_en_width, adpt_smr_fac, adpt_smr_max, levelspacing_k, UU)
    use w90_io, only: io_error
    use w90_constants, only: dp, smearing_cutoff, min_smearing_binwidth_ratio
    use w90_utility, only: utility_w0gauss
    use w90_parameters, only: num_wann, spin_decomp, num_elec_per_state, &
      num_dos_project, dos_project
    use w90_spin, only: spin_get_nk

    ! Arguments
    !
    real(kind=dp), dimension(3), intent(in)               :: kpt
    real(kind=dp), dimension(:), intent(in)               :: EnergyArray
    real(kind=dp), dimension(:), intent(in)               :: eig_k
    real(kind=dp), dimension(:, :), intent(out)            :: dos_k
    integer, intent(in)                                   :: smr_index
    real(kind=dp), intent(in), optional                    :: smr_fixed_en_width
    real(kind=dp), intent(in), optional                    :: adpt_smr_fac
    real(kind=dp), intent(in), optional                    :: adpt_smr_max
    real(kind=dp), dimension(:), intent(in), optional      :: levelspacing_k
    complex(kind=dp), dimension(:, :), intent(in), optional :: UU

    ! Adaptive smearing
    !
    real(kind=dp) :: eta_smr, arg

    ! Misc/Dummy
    !
    integer          :: i, j, loop_f, min_f, max_f, num_s_steps
    real(kind=dp)    :: rdum, spn_nk(num_wann), alpha_sq, beta_sq
    real(kind=dp)    :: binwidth, r_num_elec_per_state
    logical          :: DoSmearing

    if (present(levelspacing_k)) then
      if (present(smr_fixed_en_width)) &
        call io_error('Cannot call doskpt with levelspacing_k and ' &
                      //'with smr_fixed_en_width parameters together')
      if (.not. (present(adpt_smr_fac))) &
        call io_error('Cannot call doskpt with levelspacing_k and ' &
                      //'without adpt_smr_fac parameter')
      if (.not. (present(adpt_smr_max))) &
        call io_error('Cannot call doskpt with levelspacing_k and ' &
                      //'without adpt_smr_max parameter')
    else
      if (present(adpt_smr_fac)) &
        call io_error('Cannot call doskpt without levelspacing_k and ' &
                      //'with adpt_smr_fac parameter')
      if (present(adpt_smr_max)) &
        call io_error('Cannot call doskpt without levelspacing_k and ' &
                      //'with adpt_smr_max parameter')
      if (.not. (present(smr_fixed_en_width))) &
        call io_error('Cannot call doskpt without levelspacing_k and ' &
                      //'without smr_fixed_en_width parameter')
    end if

    r_num_elec_per_state = real(num_elec_per_state, kind=dp)

    ! Get spin projections for every band
    !
    if (spin_decomp) call spin_get_nk(kpt, spn_nk)

    binwidth = EnergyArray(2) - EnergyArray(1)

    dos_k = 0.0_dp
    do i = 1, num_wann
      if (spin_decomp) then
        ! Contribution to spin-up DOS of Bloch spinor with component
        ! (alpha,beta) with respect to the chosen quantization axis
        alpha_sq = (1.0_dp + spn_nk(i))/2.0_dp ! |alpha|^2
        ! Contribution to spin-down DOS
        beta_sq = 1.0_dp - alpha_sq ! |beta|^2 = 1 - |alpha|^2
      end if

      if (.not. present(levelspacing_k)) then
        eta_smr = smr_fixed_en_width
      else
        ! Eq.(35) YWVS07
        eta_smr = min(levelspacing_k(i)*adpt_smr_fac, adpt_smr_max)
!          eta_smr=max(eta_smr,min_smearing_binwidth_ratio) !! No: it would render the next if always false
      end if

      ! Faster optimization: I precalculate the indices
      if (eta_smr/binwidth < min_smearing_binwidth_ratio) then
        min_f = max(nint((eig_k(i) - EnergyArray(1))/ &
                         (EnergyArray(size(EnergyArray)) - EnergyArray(1)) &
                         *real(size(EnergyArray) - 1, kind=dp)) + 1, 1)
        max_f = min(nint((eig_k(i) - EnergyArray(1))/ &
                         (EnergyArray(size(EnergyArray)) - EnergyArray(1)) &
                         *real(size(EnergyArray) - 1, kind=dp)) + 1, size(EnergyArray))
        DoSmearing = .false.
      else
        min_f = max(nint((eig_k(i) - smearing_cutoff*eta_smr - EnergyArray(1))/ &
                         (EnergyArray(size(EnergyArray)) - EnergyArray(1)) &
                         *real(size(EnergyArray) - 1, kind=dp)) + 1, 1)
        max_f = min(nint((eig_k(i) + smearing_cutoff*eta_smr - EnergyArray(1))/ &
                         (EnergyArray(size(EnergyArray)) - EnergyArray(1)) &
                         *real(size(EnergyArray) - 1, kind=dp)) + 1, size(EnergyArray))
        DoSmearing = .true.
      end if

      do loop_f = min_f, max_f
        ! kind of smearing read from input (internal smearing_index variable)
        if (DoSmearing) then
          arg = (EnergyArray(loop_f) - eig_k(i))/eta_smr
          rdum = utility_w0gauss(arg, smr_index)/eta_smr
        else
          rdum = 1._dp/(EnergyArray(2) - EnergyArray(1))
        end if

        !
        ! Contribution to total DOS
        !
        if (num_dos_project == num_wann) then
          !
          ! Total DOS (default): do not loop over j, to save time
          !
          dos_k(loop_f, 1) = dos_k(loop_f, 1) + rdum*r_num_elec_per_state
          ! [GP] I don't put num_elec_per_state here below: if we are
          ! calculating the spin decomposition, we should be doing a
          ! calcultation with spin-orbit, and thus num_elec_per_state=1!
          if (spin_decomp) then
            ! Spin-up contribution
            dos_k(loop_f, 2) = dos_k(loop_f, 2) + rdum*alpha_sq
            ! Spin-down contribution
            dos_k(loop_f, 3) = dos_k(loop_f, 3) + rdum*beta_sq
          end if
        else ! 0<num_dos_project<num_wann
          !
          ! Partial DOS, projected onto the WFs with indices
          ! n=dos_project(1:num_dos_project)
          !
          do j = 1, num_dos_project
            dos_k(loop_f, 1) = dos_k(loop_f, 1) + rdum*r_num_elec_per_state &
                               *abs(UU(dos_project(j), i))**2
            if (spin_decomp) then
              ! Spin-up contribution
              dos_k(loop_f, 2) = dos_k(loop_f, 2) &
                                 + rdum*alpha_sq*abs(UU(dos_project(j), i))**2
              ! Spin-down contribution
              dos_k(loop_f, 3) = dos_k(loop_f, 3) &
                                 + rdum*beta_sq*abs(UU(dos_project(j), i))**2
            end if
          enddo
        endif
      enddo !loop_f
    end do !loop over bands

  end subroutine dos_get_k

  ! =========================================================================

  subroutine dos_get_levelspacing(del_eig, kmesh, levelspacing)
    !! This subroutine calculates the level spacing, i.e. how much the level changes
    !! near a given point of the interpolation mesh
    use w90_parameters, only: num_wann
    use w90_postw90_common, only: pw90common_kmesh_spacing

    real(kind=dp), dimension(num_wann, 3), intent(in) :: del_eig
    !! Band velocities, already corrected when degeneracies occur
    integer, dimension(3), intent(in)                :: kmesh
    !! array of three integers, giving the number of k points along
    !! each of the three directions defined by the reciprocal lattice vectors
    real(kind=dp), dimension(num_wann), intent(out)  :: levelspacing
    !! On output, the spacing for each of the bands (in eV)

    real(kind=dp) :: Delta_k
    integer :: band

    Delta_k = pw90common_kmesh_spacing(kmesh)
    do band = 1, num_wann
      levelspacing(band) = &
        sqrt(dot_product(del_eig(band, :), del_eig(band, :)))*Delta_k
    end do

  end subroutine dos_get_levelspacing

! Next routine is commented; it is the older version
!~  subroutine dos_get_eig_levelspacing_k(kpt,eig,levelspacing)
!~
!~    use w90_constants, only     : dp,cmplx_0,cmplx_i,twopi
!~    use w90_io, only            : io_error
!~    use w90_utility, only   : utility_diagonalize
!~    use w90_parameters, only    : num_wann,dos_num_points
!~    use w90_postw90_common, only : pw90common_fourier_R_to_k,pw90common_kmesh_spacing
!~    use w90_get_oper, only      : HH_R
!~    use w90_wan_ham, only   : wham_get_deleig_a
!~
!~    ! Arguments
!~    !
!~    real(kind=dp), intent(in)  :: kpt(3)
!~    real(kind=dp), intent(out) :: eig(num_wann)
!~    real(kind=dp), intent(out) :: levelspacing(num_wann)
!~
!~    complex(kind=dp), allocatable :: HH(:,:)
!~    complex(kind=dp), allocatable :: delHH(:,:,:)
!~    complex(kind=dp), allocatable :: UU(:,:)
!~
!~    ! Adaptive smearing
!~    !
!~    real(kind=dp) :: del_eig(num_wann,3),Delta_k
!~
!~    integer          :: i
!~
!~    allocate(HH(num_wann,num_wann))
!~    allocate(delHH(num_wann,num_wann,3))
!~    allocate(UU(num_wann,num_wann))
!~
!~    call pw90common_fourier_R_to_k(kpt,HH_R,HH,0)
!~    call utility_diagonalize(HH,num_wann,eig,UU)
!~    call pw90common_fourier_R_to_k(kpt,HH_R,delHH(:,:,1),1)
!~    call pw90common_fourier_R_to_k(kpt,HH_R,delHH(:,:,2),2)
!~    call pw90common_fourier_R_to_k(kpt,HH_R,delHH(:,:,3),3)
!~    call wham_get_deleig_a(del_eig(:,1),eig,delHH(:,:,1),UU)
!~    call wham_get_deleig_a(del_eig(:,2),eig,delHH(:,:,2),UU)
!~    call wham_get_deleig_a(del_eig(:,3),eig,delHH(:,:,3),UU)
!~
!~    Delta_k=pw90common_kmesh_spacing(dos_num_points)
!~    do i=1,num_wann
!~       levelspacing(i)=&
!~            sqrt(dot_product(del_eig(i,:),del_eig(i,:)))*Delta_k
!~    end do
!~
!~  end subroutine dos_get_eig_levelspacing_k

!~  !=========================================================!
!~  !                   PRIVATE PROCEDURES                    !
!~  !=========================================================!
!~
!~
!~  function count_states(energy,eig,levelspacing,npts)
!~
!~    use w90_constants, only     : dp,cmplx_0,cmplx_i,twopi
!~    use w90_utility, only       : utility_wgauss
!~    use w90_postw90_common, only : weight
!~    use w90_parameters, only    : num_wann,dos_adpt_smr_fac
!~
!~    real(kind=dp) :: count_states
!~
!~    ! Arguments
!~    !
!~    real(kind=dp)                  :: energy
!~    real(kind=dp), dimension (:,:) :: eig
!~    real(kind=dp), dimension (:,:) :: levelspacing
!~    integer                        :: npts
!~
!~    ! Misc/Dummy
!~    !
!~    integer       :: loop_k,i
!~    real(kind=dp) :: sum,eta_smr,arg
!~
!~    count_states=0.0_dp
!~    do loop_k=1,npts
!~       sum=0.0_dp
!~       do i=1,num_wann
!~          eta_smr=levelspacing(i,loop_k)*dos_adpt_smr_fac
!~          arg=(energy-eig(i,loop_k))/eta_smr
!~          !
!~          ! For Fe and a 125x125x125 interpolation mesh, E_f=12.6306 with M-P
!~          ! smearing, and E_f=12.6512 with F-D smearing
!~          !
!~          !          sum=sum+utility_wgauss(arg,-99) ! Fermi-Dirac
!~          sum=sum+utility_wgauss(arg,1)    ! Methfessel-Paxton case
!~       end do
!~       count_states=count_states+weight(loop_k)*sum
!~    end do
!~
!~  end function count_states

end module w90_dos

