!===================================================================
subroutine get_recip_lattice(real_lat, recip_lat)  !
!==================================================================!
!                                                                  !
!!  Calculates the reciprical lattice vectors and the cell volume
!                                                                  !
!===================================================================
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(kind=dp), parameter :: pi = 3.141592653589793238462643383279_dp
  real(kind=dp), parameter :: twopi = 2*pi
  real(kind=dp), parameter :: eps5 = 1.0e-5_dp
  real(kind=dp) :: real_lat(3, 3)
  real(kind=dp) :: recip_lat(3, 3)
  real(kind=dp) :: volume

  recip_lat(1, 1) = real_lat(2, 2)*real_lat(3, 3) - real_lat(3, 2)*real_lat(2, 3)
  recip_lat(1, 2) = real_lat(2, 3)*real_lat(3, 1) - real_lat(3, 3)*real_lat(2, 1)
  recip_lat(1, 3) = real_lat(2, 1)*real_lat(3, 2) - real_lat(3, 1)*real_lat(2, 2)
  recip_lat(2, 1) = real_lat(3, 2)*real_lat(1, 3) - real_lat(1, 2)*real_lat(3, 3)
  recip_lat(2, 2) = real_lat(3, 3)*real_lat(1, 1) - real_lat(1, 3)*real_lat(3, 1)
  recip_lat(2, 3) = real_lat(3, 1)*real_lat(1, 2) - real_lat(1, 1)*real_lat(3, 2)
  recip_lat(3, 1) = real_lat(1, 2)*real_lat(2, 3) - real_lat(2, 2)*real_lat(1, 3)
  recip_lat(3, 2) = real_lat(1, 3)*real_lat(2, 1) - real_lat(2, 3)*real_lat(1, 1)
  recip_lat(3, 3) = real_lat(1, 1)*real_lat(2, 2) - real_lat(2, 1)*real_lat(1, 2)

  volume = real_lat(1, 1)*recip_lat(1, 1) + &
           real_lat(1, 2)*recip_lat(1, 2) + &
           real_lat(1, 3)*recip_lat(1, 3)

  recip_lat = twopi*recip_lat/volume

  return

end subroutine get_recip_lattice

program test_library
  !! NOTE! THIS PROGRAM ONLY WORKS IN SERIAL FOR NOW
  !! (even if there are some stubs that could make you think
  !! it works in parallel...)
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: num_nnmax = 12

  !! Wannier-setup params
  character(len=100) :: seed__name
  integer, dimension(3) :: mp_grid_loc
  integer :: num_kpts_loc
  real(kind=dp), dimension(3, 3) :: real_lattice_loc
  real(kind=dp), dimension(3, 3) :: recip_lattice_loc
  real(kind=dp), allocatable :: kpt_latt_loc(:, :)
  integer :: num_bands_tot
  integer :: num_atoms_loc
  character(len=4), allocatable :: atom_symbols_loc(:)
  real(kind=dp), allocatable :: atoms_cart_loc(:, :)
  logical :: gamma_only_loc
  logical :: spinors_loc
  ! Intent out below here
  integer :: nntot_loc
  integer, allocatable :: nnlist_loc(:, :)
  integer, allocatable :: nncell_loc(:, :, :)
  integer :: num_bands_loc
  integer :: num_wann_loc
  real(kind=dp), allocatable :: proj_site_loc(:, :)
  integer, allocatable :: proj_l_loc(:)
  integer, allocatable :: proj_m_loc(:)
  integer, allocatable :: proj_radial_loc(:)
  real(kind=dp), allocatable :: proj_z_loc(:, :)
  real(kind=dp), allocatable :: proj_x_loc(:, :)
  real(kind=dp), allocatable :: proj_zona_loc(:)
  integer, allocatable :: exclude_bands_loc(:)
  integer, allocatable :: proj_s_loc(:)
  real(kind=dp), allocatable :: proj_s_qaxis_loc(:, :)

  !! Wannier-run params
  complex(kind=dp), allocatable :: M_matrix_loc(:, :, :, :)
  complex(kind=dp), allocatable :: A_matrix_loc(:, :, :)
  real(kind=dp), allocatable :: eigenvalues_loc(:, :)
  ! Intent out below here
  complex(kind=dp), allocatable :: U_matrix_loc(:, :, :)
  complex(kind=dp), allocatable :: U_matrix_opt_loc(:, :, :)
  logical, allocatable :: lwindow_loc(:, :)
  real(kind=dp), allocatable :: wann_centres_loc(:, :)
  real(kind=dp), allocatable :: wann_spreads_loc(:)
  real(kind=dp), dimension(3) :: spread_loc

  integer :: i, j, k, l, ierr, mmn_in, amn_in, nn, inn
  logical :: nn_found
  character(len=50) :: dummy, verbosity
  real(kind=dp) :: re_tmp, im_tmp

  integer, parameter :: stdout = 6
  integer, parameter :: root_id = 0
  integer :: nb_tmp, nkp_tmp, nntot_tmp, nw_tmp, num_mmn, num_amn, ncount
  integer :: nkp, nkp2, nnl, nnm, nnn, m, n
  real(kind=dp) :: m_real, m_imag, a_real, a_imag
  complex(kind=dp), allocatable :: mmn_tmp(:, :)

  integer :: num_nodes, my_node_id

  logical :: verbose

  NAMELIST /PARAMS/ seed__name, mp_grid_loc, num_bands_tot, gamma_only_loc, spinors_loc, verbosity

#ifdef MPI
  include 'mpif.h'

  call mpi_init(ierr)
  if (ierr .ne. 0) then
    write (0, '(/a/)') '# MPI initialisation error'
    stop 1
  end if

  call mpi_comm_rank(mpi_comm_world, my_node_id, ierr)
  call mpi_comm_size(mpi_comm_world, num_nodes, ierr)
  if (my_node_id == 0) then
    print *, "# COMPILED IN PARALLEL, RUNNING ON ", num_nodes, " NODES"
  end if
#else
  num_nodes = 1
  my_node_id = 0
  print *, "# COMPILED IN SERIAL"
#endif

  verbosity = 'low'
  OPEN (unit=100, file='PARAMS', status='old', action='read')
  READ (UNIT=100, NML=PARAMS)
  CLOSE (100)
  num_kpts_loc = PRODUCT(mp_grid_loc)

  if (TRIM(verbosity) == 'high') then
    verbose = .true.
  elseif (TRIM(verbosity) == 'low') then
    verbose = .false.
  else
    write (0, *) "# INVALID VERBOSITY VALUE, can be only 'low' or 'high'"
    stop 1
  end if

  ! First line: num_atoms
  ! Next num_atoms lines: positions
  ! next num_atoms lines: coordinates
  OPEN (unit=100, file='POSITIONS', status='old', action='read')
  READ (100, *) num_atoms_loc
  allocate (atom_symbols_loc(num_atoms_loc), stat=ierr)
  if (ierr /= 0) then
    write (0, *) "# ERROR DURING ALLOCATION"
    stop 1
  end if
  allocate (atoms_cart_loc(3, num_atoms_loc), stat=ierr)
  if (ierr /= 0) then
    write (0, *) "# ERROR DURING ALLOCATION"
    stop 1
  end if
  READ (100, *) atom_symbols_loc
  DO i = 1, num_atoms_loc
    READ (100, *) atoms_cart_loc(:, i)
  END DO
  CLOSE (100)

  ! 3x3 cell, each row is a vector (so transpose w.r.t. fortran)
  OPEN (unit=100, file='CELL', status='old', action='read')
  DO i = 1, 3
    ! Index inverted to read in directly the transpose
    READ (100, *) real_lattice_loc(i, :)
  END DO
  !real_lattice_loc = TRANSPOSE(real_lattice_loc)
  CLOSE (100)
  CALL get_recip_lattice(real_lattice_loc, recip_lattice_loc)

  ! each line: the coordinates of the kpoints
  allocate (kpt_latt_loc(3, num_kpts_loc))
  OPEN (unit=100, file='KPOINTS', status='old', action='read')
  DO i = 1, num_kpts_loc
    READ (100, *) kpt_latt_loc(:, i)
  END DO
  CLOSE (100)

  if (verbose) then
    print *, "# INPUTS READ."

    print *, '#seed_name:', trim(seed__name)
    print *, '#num_atoms_loc:', num_atoms_loc
    print *, '#num_kpts_loc:', num_kpts_loc
    print *, '#num_nnmax:', num_nnmax
    print *, '#num_bands_tot:', num_bands_tot
  end if

  allocate (nnlist_loc(num_kpts_loc, num_nnmax), stat=ierr)
  if (ierr /= 0) then
    write (0, *) "# ERROR DURING ALLOCATION"
    stop 1
  end if
  allocate (nncell_loc(3, num_kpts_loc, num_nnmax), stat=ierr)
  if (ierr /= 0) then
    write (0, *) "# ERROR DURING ALLOCATION"
    stop 1
  end if
  allocate (proj_site_loc(3, num_bands_tot), stat=ierr)
  if (ierr /= 0) then
    write (0, *) "# ERROR DURING ALLOCATION"
    stop 1
  end if
  allocate (proj_l_loc(num_bands_tot), stat=ierr)
  if (ierr /= 0) then
    write (0, *) "# ERROR DURING ALLOCATION"
    stop 1
  end if
  allocate (proj_m_loc(num_bands_tot), stat=ierr)
  if (ierr /= 0) then
    write (0, *) "# ERROR DURING ALLOCATION"
    stop 1
  end if
  allocate (proj_radial_loc(num_bands_tot), stat=ierr)
  if (ierr /= 0) then
    write (0, *) "# ERROR DURING ALLOCATION"
    stop 1
  end if
  allocate (proj_z_loc(3, num_bands_tot), stat=ierr)
  if (ierr /= 0) then
    write (0, *) "# ERROR DURING ALLOCATION"
    stop 1
  end if
  allocate (proj_x_loc(3, num_bands_tot), stat=ierr)
  if (ierr /= 0) then
    write (0, *) "# ERROR DURING ALLOCATION"
    stop 1
  end if
  allocate (proj_zona_loc(num_bands_tot), stat=ierr)
  if (ierr /= 0) then
    write (0, *) "# ERROR DURING ALLOCATION"
    stop 1
  end if
  allocate (exclude_bands_loc(num_bands_tot), stat=ierr)
  if (ierr /= 0) then
    write (0, *) "# ERROR DURING ALLOCATION"
    stop 1
  end if
  allocate (proj_s_loc(num_bands_tot), stat=ierr)
  if (ierr /= 0) then
    write (0, *) "# ERROR DURING ALLOCATION"
    stop 1
  end if
  allocate (proj_s_qaxis_loc(3, num_bands_tot), stat=ierr)
  if (ierr /= 0) then
    write (0, *) "# ERROR DURING ALLOCATION"
    stop 1
  end if

  if (verbose) then
    print *, '#seed_name:', trim(seed__name)
    print *, '#mp_grid_loc:', mp_grid_loc
    print *, '#num_kpts_loc:', num_kpts_loc
    print *, '#real_lattice_loc(vec1):', real_lattice_loc(1, :)
    print *, '#real_lattice_loc(vec2):', real_lattice_loc(2, :)
    print *, '#real_lattice_loc(vec3):', real_lattice_loc(3, :)
    !print*, '#recip_lattice_loc(vec1):', recip_lattice_loc(:,1)
    !print*, '#recip_lattice_loc(vec2):', recip_lattice_loc(:,2)
    !print*, '#recip_lattice_loc(vec3):', recip_lattice_loc(:,3)
    !print*, '#kpt_latt_loc:', kpt_latt_loc
    !print*, '#atom_symbols_loc:', atom_symbols_loc
    !print*, '#atoms_cart_loc:', atoms_cart_loc
    print *, '#gamma_only_loc:', gamma_only_loc
    print *, '#spinors_loc:', spinors_loc
  end if

  call wannier_setup(seed__name, mp_grid_loc, num_kpts_loc, &
                     real_lattice_loc, recip_lattice_loc, kpt_latt_loc, num_bands_tot, &
                     num_atoms_loc, atom_symbols_loc, atoms_cart_loc, gamma_only_loc, spinors_loc, &
                     nntot_loc, nnlist_loc, nncell_loc, num_bands_loc, num_wann_loc, &
                     proj_site_loc, proj_l_loc, proj_m_loc, proj_radial_loc, proj_z_loc, &
                     proj_x_loc, proj_zona_loc, exclude_bands_loc, proj_s_loc, proj_s_qaxis_loc)

  if (verbose) then
    print *, "# WANNIER_SETUP CALLED."

    print *, "# nntot_loc", nntot_loc
    !print*, "# nnlist_loc", nnlist_loc
    !print*, "# nncell_loc", nncell_loc
    print *, '#num_bands_loc:', num_bands_loc
    print *, '#num_wann_loc:', num_wann_loc
    !print*, '#proj_site_loc:', proj_site_loc
  end if

  allocate (M_matrix_loc(num_bands_loc, num_bands_loc, nntot_loc, num_kpts_loc), stat=ierr)
  if (ierr /= 0) then
    write (0, *) "# ERROR DURING ALLOCATION"
    stop 1
  end if
  M_matrix_loc = 0.

  allocate (A_matrix_loc(num_bands_loc, num_wann_loc, num_kpts_loc), stat=ierr)
  if (ierr /= 0) then
    write (0, *) "# ERROR DURING ALLOCATION"
    stop 1
  end if
  A_matrix_loc = 0.

  ! Use the usual .mmn and .amn files

  if (my_node_id == root_id) then
    ! Read M_matrix_orig from file
    mmn_in = 100 ! Unit number
    open (unit=mmn_in, file=trim(seed__name)//'.mmn', &
          form='formatted', status='old', action='read', err=101)

    if (verbose) write (stdout, '(/a)', advance='no') '# Reading overlaps from '//trim(seed__name)//'.mmn    : '

    ! Read the comment line
    read (mmn_in, '(a)', err=103, end=103) dummy
    if (verbose) write (stdout, '(a)') '# '//trim(dummy)

    ! Read the number of bands, k-points and nearest neighbours
    read (mmn_in, *, err=103, end=103) nb_tmp, nkp_tmp, nntot_tmp

    ! Checks
    if (nb_tmp .ne. num_bands_loc) then
      write (stdout, '(/a/)') '# '//trim(seed__name)//'.mmn has not the right number of bands'
      stop 1
    end if
    if (nkp_tmp .ne. num_kpts_loc) then
      write (stdout, '(/a/)') '# '//trim(seed__name)//'.mmn has not the right number of k-points'
      stop 1
    end if
    if (nntot_tmp .ne. nntot_loc) then
      write (stdout, '(/a/)') '# '//trim(seed__name)//'.mmn has not the right number of nearest neighbours'
      stop 1
    end if

    ! Read the overlaps
    num_mmn = num_kpts_loc*nntot_loc
    allocate (mmn_tmp(num_bands_loc, num_bands_loc), stat=ierr)
    if (ierr /= 0) then
      write (stdout, '(/a/)') '# Error in allocating mmn_tmp in overlap_read'
      stop 1
    end if
    do ncount = 1, num_mmn
      read (mmn_in, *, err=103, end=103) nkp, nkp2, nnl, nnm, nnn
      do n = 1, num_bands_loc
        do m = 1, num_bands_loc
          read (mmn_in, *, err=103, end=103) m_real, m_imag
          mmn_tmp(m, n) = cmplx(m_real, m_imag, kind=dp)
        enddo
      enddo
      nn = 0
      nn_found = .false.
      do inn = 1, nntot_loc
        if ((nkp2 .eq. nnlist_loc(nkp, inn)) .and. &
            (nnl .eq. nncell_loc(1, nkp, inn)) .and. &
            (nnm .eq. nncell_loc(2, nkp, inn)) .and. &
            (nnn .eq. nncell_loc(3, nkp, inn))) then
          if (.not. nn_found) then
            nn_found = .true.
            nn = inn
          else
            write (stdout, '(/a/)') '# Error reading '//trim(seed__name)// &
              '.mmn. More than one matching nearest neighbour found'
            stop 1
          endif
        endif
      end do
      if (nn .eq. 0) then
        write (stdout, '(/a,i8,2i5,i4,2x,3i3)') &
          '# Error reading '//trim(seed__name)//'.mmn, Neighbor not found:', &
          ncount, nkp, nkp2, nn, nnl, nnm, nnn

        !do i=1, num_kpts_loc
        !    write(stdout, *) nnlist_loc(i, :)
        !end do

        write (stdout, *) '# nkp:', nkp
        write (stdout, *) '# nkp2:', nkp2
        do inn = 1, nntot_loc
          write (stdout, '(a,I1,a,I10,a)', advance='no') '# nnlist(nkp,', inn, '):', nnlist_loc(nkp, inn), '  ---> '
          write (stdout, '(a,I1,a,3I5)') '# nncell(nkp,', inn, ',:):', nncell_loc(:, nkp, inn)
        enddo

        write (stdout, *) '# KPOINTS:'
        do i = 1, num_kpts_loc
          write (stdout, '(A,3G18.10)') '#', kpt_latt_loc(:, i)
        END DO
        !write(stdout, *) '_____'
        !write(stdout, *) nncell_loc
        stop 1
      end if
      m_matrix_loc(:, :, nn, nkp) = mmn_tmp(:, :)
    end do
    deallocate (mmn_tmp, stat=ierr)
    if (ierr /= 0) then
      write (stdout, '(/a/)') '# Error in deallocating mmn_tmp in overlap_read'
      stop 1
    end if
    close (mmn_in)
  end if

  ! Broadcast
#ifdef MPI
  call MPI_bcast(m_matrix_loc(1, 1, 1, 1), num_bands_loc*num_bands_loc*nntot_loc*num_kpts_loc, &
                 MPI_double_complex, root_id, mpi_comm_world, ierr)
  if (ierr .ne. MPI_success) then
    write (stdout, '(/a/)') '# Error in comms_bcast_cmplx'
    stop 1
  end if
#endif
  if (my_node_id == root_id) then
    ! Read A_matrix from file wannier.amn
    amn_in = 100 ! Unit number
    open (unit=amn_in, file=trim(seed__name)//'.amn', form='formatted', status='old', err=102)

    if (verbose) write (stdout, '(/a)', advance='no') '# Reading projections from '//trim(seed__name)//'.amn : '

    ! Read the comment line
    read (amn_in, '(a)', err=104, end=104) dummy
    if (verbose) write (stdout, '(a)') '# '//trim(dummy)

    ! Read the number of bands, k-points and wannier functions
    read (amn_in, *, err=104, end=104) nb_tmp, nkp_tmp, nw_tmp

    ! Checks
    if (nb_tmp .ne. num_bands_loc) then
      write (stdout, '(/a/)') '# '//trim(seed__name)//'.amn has not the right number of bands'
      stop 1
    end if
    if (nkp_tmp .ne. num_kpts_loc) then
      write (stdout, '(/a/)') '# '//trim(seed__name)//'.amn has not the right number of k-points'
      stop 1
    end if
    if (nw_tmp .ne. num_wann_loc) then
      write (stdout, '(/a/)') '# '//trim(seed__name)//'.amn has not the right number of Wannier functions'
      stop 1
    end if

    ! Read the projections
    num_amn = num_bands_loc*num_wann_loc*num_kpts_loc
    do ncount = 1, num_amn
      read (amn_in, *, err=104, end=104) m, n, nkp, a_real, a_imag
      a_matrix_loc(m, n, nkp) = cmplx(a_real, a_imag, kind=dp)
    end do
    close (amn_in)
  end if

#ifdef MPI
  call MPI_bcast(a_matrix_loc(1, 1, 1), num_bands_loc*num_wann_loc*num_kpts_loc, &
                 MPI_double_complex, root_id, mpi_comm_world, ierr)
  if (ierr .ne. MPI_success) then
    write (stdout, '(/a/)') '# Error in comms_bcast_cmplx'
    stop 1
  end if
#endif

  allocate (eigenvalues_loc(num_bands_loc, num_kpts_loc), stat=ierr)
  if (ierr /= 0) then
    write (0, *) "# ERROR DURING ALLOCATION"
    stop 1
  end if
  OPEN (unit=100, file='EIG', status='old', action='read')
  do j = 1, num_kpts_loc
    do i = 1, num_bands_loc
      read (100, '(2G25.17)') eigenvalues_loc(i, j)
    end do
  end do
  CLOSE (100)

  allocate (U_matrix_loc(num_wann_loc, num_wann_loc, num_kpts_loc), stat=ierr)
  if (ierr /= 0) then
    write (0, *) "# ERROR DURING ALLOCATION"
    stop 1
  end if
  allocate (U_matrix_opt_loc(num_bands_loc, num_wann_loc, num_kpts_loc), stat=ierr)
  if (ierr /= 0) then
    write (0, *) "# ERROR DURING ALLOCATION"
    stop 1
  end if
  allocate (lwindow_loc(num_bands_loc, num_kpts_loc), stat=ierr)
  if (ierr /= 0) then
    write (0, *) "# ERROR DURING ALLOCATION"
    stop 1
  end if
  allocate (wann_centres_loc(3, num_wann_loc), stat=ierr)
  if (ierr /= 0) then
    write (0, *) "# ERROR DURING ALLOCATION"
    stop 1
  end if
  allocate (wann_spreads_loc(num_wann_loc), stat=ierr)
  if (ierr /= 0) then
    write (0, *) "# ERROR DURING ALLOCATION"
    stop 1
  end if

  if (verbose) print *, "# SECOND ALLOCATION PHASE COMPLETED."

  call wannier_run(seed__name, mp_grid_loc, num_kpts_loc, &
                   real_lattice_loc, recip_lattice_loc, kpt_latt_loc, num_bands_loc, &
                   num_wann_loc, nntot_loc, num_atoms_loc, atom_symbols_loc, &
                   atoms_cart_loc, gamma_only_loc, M_matrix_loc, A_matrix_loc, eigenvalues_loc, &
                   U_matrix_loc, U_matrix_opt_loc, lwindow_loc, wann_centres_loc, &
                   wann_spreads_loc, spread_loc)

  if (verbose) print *, "# WANNIER_RUN CALLED."

  print *, "# Wannier run completed without errors."

  OPEN (unit=123, file='results.dat', action='write')
  do i = 1, num_wann_loc
    write (123, '(4G18.10)') wann_centres_loc(:, i), wann_spreads_loc(i)
  end do
  CLOSE (123)
  print *, "# MLWF centre_x centre_y centre_z spread writte in results.dat"

  stop

101 write (stdout, '(/a/)') '# Error: Problem opening input file .mmn'
  stop 1
102 write (stdout, '(/a/)') '# Error: Problem opening input file .amn'
  stop 1
103 write (stdout, '(/a/)') '# Error: Problem reading input file .mmn'
  stop 1
104 write (stdout, '(/a/)') '# Error: Problem reading input file .amn'
  stop 1

end program test_library
