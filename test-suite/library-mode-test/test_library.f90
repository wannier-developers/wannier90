!===================================================================
subroutine get_recip_lattice (real_lat,recip_lat)  !
!==================================================================!
!                                                                  !
!!  Calculates the reciprical lattice vectors and the cell volume 
!                                                                  !
!===================================================================
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    real(kind=dp), parameter :: pi=3.141592653589793238462643383279_dp
    real(kind=dp), parameter :: twopi = 2*pi
    real(kind=dp), parameter :: eps5 = 1.0e-5_dp
    real(kind=dp) :: real_lat (3, 3)
    real(kind=dp) :: recip_lat (3, 3)  
    real(kind=dp) :: volume

    recip_lat(1,1)=real_lat(2,2)*real_lat(3,3)-real_lat(3,2)*real_lat(2,3)
    recip_lat(1,2)=real_lat(2,3)*real_lat(3,1)-real_lat(3,3)*real_lat(2,1)
    recip_lat(1,3)=real_lat(2,1)*real_lat(3,2)-real_lat(3,1)*real_lat(2,2)
    recip_lat(2,1)=real_lat(3,2)*real_lat(1,3)-real_lat(1,2)*real_lat(3,3)
    recip_lat(2,2)=real_lat(3,3)*real_lat(1,1)-real_lat(1,3)*real_lat(3,1)
    recip_lat(2,3)=real_lat(3,1)*real_lat(1,2)-real_lat(1,1)*real_lat(3,2)
    recip_lat(3,1)=real_lat(1,2)*real_lat(2,3)-real_lat(2,2)*real_lat(1,3)
    recip_lat(3,2)=real_lat(1,3)*real_lat(2,1)-real_lat(2,3)*real_lat(1,1)
    recip_lat(3,3)=real_lat(1,1)*real_lat(2,2)-real_lat(2,1)*real_lat(1,2)

    volume=real_lat(1,1)*recip_lat(1,1) + &
            real_lat(1,2)*recip_lat(1,2) + &
            real_lat(1,3)*recip_lat(1,3)  

    recip_lat=twopi*recip_lat/volume

    return

end subroutine get_recip_lattice

program test_library
    implicit none

    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: num_nnmax=12

    !! Wannier-setup params
    character(len=100) :: seed__name
    integer, dimension(3) :: mp_grid_loc
    integer :: num_kpts_loc
    real(kind=dp), dimension(3,3) :: real_lattice_loc
    real(kind=dp), dimension(3,3) :: recip_lattice_loc
    real(kind=dp), allocatable :: kpt_latt_loc(:,:)
    integer :: num_bands_tot
    integer :: num_atoms_loc
    character(len=4), allocatable :: atom_symbols_loc(:)
    real(kind=dp), allocatable :: atoms_cart_loc(:,:)
    logical :: gamma_only_loc
    logical :: spinors_loc
    ! Intent out below here
    integer :: nntot_loc
    integer, allocatable :: nnlist_loc(:,:)
    integer, allocatable :: nncell_loc(:,:,:)
    integer :: num_bands_loc
    integer :: num_wann_loc
    real(kind=dp), allocatable :: proj_site_loc(:,:)
    integer, allocatable :: proj_l_loc(:)
    integer, allocatable :: proj_m_loc(:)
    integer, allocatable :: proj_radial_loc(:)
    real(kind=dp), allocatable :: proj_z_loc(:,:)
    real(kind=dp), allocatable :: proj_x_loc(:,:)
    real(kind=dp), allocatable :: proj_zona_loc(:)
    integer, allocatable :: exclude_bands_loc(:)
    integer, allocatable :: proj_s_loc(:)
    real(kind=dp), allocatable :: proj_s_qaxis_loc(:,:)

    !! Wannier-run params
    complex(kind=dp), allocatable :: M_matrix_loc(:,:,:,:)
    complex(kind=dp), allocatable :: A_matrix_loc(:,:,:)
    real(kind=dp), allocatable :: eigenvalues_loc(:,:)
    ! Intent out below here
    complex(kind=dp), allocatable :: U_matrix_loc(:,:,:)
    complex(kind=dp), allocatable :: U_matrix_opt_loc(:,:,:)
    logical, allocatable :: lwindow_loc(:,:)
    real(kind=dp), allocatable :: wann_centres_loc(:,:)
    real(kind=dp), allocatable :: wann_spreads_loc(:)
    real(kind=dp), dimension(3) :: spread_loc

    integer :: i, j, k, l, ierr
    character(len=50) :: dummy
    real(kind=dp) :: re_tmp, im_tmp

    NAMELIST / PARAMS / seed__name, mp_grid_loc, num_bands_tot, gamma_only_loc, spinors_loc

    print*, "STARTING..."

    OPEN(unit=100, file='PARAMS', status='old', action='read')
    READ (UNIT=100, NML=PARAMS)
    CLOSE(100)
    num_kpts_loc = PRODUCT(mp_grid_loc)

    ! First line: num_atoms
    ! Next num_atoms lines: positions
    ! next num_atoms lines: coordinates
    OPEN(unit=100, file='POSITIONS', status='old', action='read')
    READ(100,*) num_atoms_loc
    allocate(atom_symbols_loc(num_atoms_loc),stat=ierr)
    if (ierr/=0) then
        write(0,*) "ERROR DURING ALLOCATION"
        stop 1
    end if
    allocate(atoms_cart_loc(3,num_atoms_loc),stat=ierr)
    if (ierr/=0) then
        write(0,*) "ERROR DURING ALLOCATION"
        stop 1
    end if
    READ(100,*) atom_symbols_loc
    DO i=1, num_atoms_loc
        READ(100,*) atoms_cart_loc(:, i)
    END DO
    CLOSE(100)

    ! 3x3 cell, each row is a vector (so transpose w.r.t. fortran)
    OPEN(unit=100, file='CELL', status='old', action='read')
    DO i=1, 3
        READ(100,*) real_lattice_loc(:, i)
    END DO
    !real_lattice_loc = TRANSPOSE(real_lattice_loc)
    CLOSE(100)
    CALL get_recip_lattice (real_lattice_loc,recip_lattice_loc)

    ! each line: the coordinates of the kpoints
    allocate(kpt_latt_loc(3,num_kpts_loc))
    OPEN(unit=100, file='KPOINTS', status='old', action='read')
    DO i=1, num_kpts_loc
        READ(100,*) kpt_latt_loc(:, i)
    END DO
    CLOSE(100)
   
    print*, "INPUTS READ."

    print*, 'seed_name:', trim(seed__name)
    print*, 'num_atoms_loc:', num_atoms_loc
    print*, 'num_kpts_loc:', num_kpts_loc
    print*, 'num_nnmax:', num_nnmax
    print*, 'num_bands_tot:', num_bands_tot

    allocate(nnlist_loc(num_kpts_loc,num_nnmax),stat=ierr)
    if (ierr/=0) then
        write(0,*) "ERROR DURING ALLOCATION"
        stop 1
    end if
    allocate(nncell_loc(3,num_kpts_loc,num_nnmax),stat=ierr)
    if (ierr/=0) then
        write(0,*) "ERROR DURING ALLOCATION"
        stop 1
    end if
    allocate(proj_site_loc(3,num_bands_tot),stat=ierr)
    if (ierr/=0) then
        write(0,*) "ERROR DURING ALLOCATION"
        stop 1
    end if
    allocate(proj_l_loc(num_bands_tot),stat=ierr)
    if (ierr/=0) then
        write(0,*) "ERROR DURING ALLOCATION"
        stop 1
    end if
    allocate(proj_m_loc(num_bands_tot),stat=ierr)
    if (ierr/=0) then
        write(0,*) "ERROR DURING ALLOCATION"
        stop 1
    end if
    allocate(proj_radial_loc(num_bands_tot),stat=ierr)
    if (ierr/=0) then
        write(0,*) "ERROR DURING ALLOCATION"
        stop 1
    end if
    allocate(proj_z_loc(3,num_bands_tot),stat=ierr)
    if (ierr/=0) then
        write(0,*) "ERROR DURING ALLOCATION"
        stop 1
    end if
    allocate(proj_x_loc(3,num_bands_tot),stat=ierr)
    if (ierr/=0) then
        write(0,*) "ERROR DURING ALLOCATION"
        stop 1
    end if
    allocate(proj_zona_loc(num_bands_tot),stat=ierr)
    if (ierr/=0) then
        write(0,*) "ERROR DURING ALLOCATION"
        stop 1
    end if
    allocate(exclude_bands_loc(num_bands_tot),stat=ierr)
    if (ierr/=0) then
        write(0,*) "ERROR DURING ALLOCATION"
        stop 1
    end if
    allocate(proj_s_loc(num_bands_tot),stat=ierr)
    if (ierr/=0) then
        write(0,*) "ERROR DURING ALLOCATION"
        stop 1
    end if
    allocate(proj_s_qaxis_loc(3,num_bands_tot),stat=ierr)
    if (ierr/=0) then
        write(0,*) "ERROR DURING ALLOCATION"
        stop 1
    end if

    print*, 'seed_name:', trim(seed__name)
    print*, 'mp_grid_loc:', mp_grid_loc
    print*, 'num_kpts_loc:', num_kpts_loc
    print*, 'real_lattice_loc(vec1):', real_lattice_loc(1,:)
    print*, 'real_lattice_loc(vec2):', real_lattice_loc(2,:)
    print*, 'real_lattice_loc(vec3):', real_lattice_loc(3,:)
    !print*, 'recip_lattice_loc(vec1):', recip_lattice_loc(:,1)
    !print*, 'recip_lattice_loc(vec2):', recip_lattice_loc(:,2)
    !print*, 'recip_lattice_loc(vec3):', recip_lattice_loc(:,3)
    !print*, 'kpt_latt_loc:', kpt_latt_loc
    !print*, 'atom_symbols_loc:', atom_symbols_loc
    !print*, 'atoms_cart_loc:', atoms_cart_loc
    print*, 'gamma_only_loc:', gamma_only_loc
    print*, 'spinors_loc:', spinors_loc

    call wannier_setup(seed__name,mp_grid_loc,num_kpts_loc,&
        real_lattice_loc,recip_lattice_loc,kpt_latt_loc,num_bands_tot, &
        num_atoms_loc,atom_symbols_loc,atoms_cart_loc, gamma_only_loc,spinors_loc, &
        nntot_loc,nnlist_loc,nncell_loc,num_bands_loc,num_wann_loc, &
        proj_site_loc,proj_l_loc,proj_m_loc,proj_radial_loc,proj_z_loc, &
        proj_x_loc,proj_zona_loc,exclude_bands_loc,proj_s_loc,proj_s_qaxis_loc)

    print*, "WANNIER_SETUP CALLED."

    print*, "nntot_loc", nntot_loc
    !print*, "nnlist_loc", nnlist_loc
    !print*, "nncell_loc", nncell_loc
    print*, 'num_bands_loc:', num_bands_loc
    print*, 'num_wann_loc:', num_wann_loc
    !print*, 'proj_site_loc:', proj_site_loc
    
    ! The next three matrices should be generated from the .mmn, ... files
    ! running wannier90.x on gaas.win using 'devel_flag=print_raw_mmn_amn_eig_and_exit' in the input
    allocate(M_matrix_loc(num_bands_loc,num_bands_loc,nntot_loc,num_kpts_loc),stat=ierr)
    if (ierr/=0) then
        write(0,*) "ERROR DURING ALLOCATION"
        stop 1
    end if
    OPEN(unit=100, file='MMN', status='old', action='read')
    do l=1,num_kpts_loc
        do k=1,nntot_loc
            do j=1,num_bands_loc
                do i=1,num_bands_loc
                    read(100,'(2G25.17)') re_tmp, im_tmp
                    m_matrix_loc(i,j,k,l) = COMPLEX(re_tmp, im_tmp)
                end do
            end do
        end do
    end do
    CLOSE(100)

    allocate(A_matrix_loc(num_bands_loc,num_wann_loc,num_kpts_loc),stat=ierr)
    if (ierr/=0) then
        write(0,*) "ERROR DURING ALLOCATION"
        stop 1
    end if
    OPEN(unit=100, file='AMN', status='old', action='read')
    do k=1,num_kpts_loc
        do j=1,num_wann_loc
            do i=1,num_bands_loc
                read(100,'(2G25.17)') re_tmp, im_tmp
                a_matrix_loc(i,j,k) = COMPLEX(re_tmp, im_tmp)
            end do
        end do
    end do
    CLOSE(100)

    allocate(eigenvalues_loc(num_bands_loc,num_kpts_loc),stat=ierr)
    if (ierr/=0) then
        write(0,*) "ERROR DURING ALLOCATION"
        stop 1
    end if
    OPEN(unit=100, file='EIG', status='old', action='read')
    do j=1,num_kpts_loc
        do i=1,num_bands_loc
            read(100,'(2G25.17)') eigenvalues_loc(i,j)
        end do
    end do
    CLOSE(100)

    allocate(U_matrix_loc(num_wann_loc,num_wann_loc,num_kpts_loc),stat=ierr)
    if (ierr/=0) then
        write(0,*) "ERROR DURING ALLOCATION"
        stop 1
    end if
    allocate(U_matrix_opt_loc(num_bands_loc,num_wann_loc,num_kpts_loc),stat=ierr)
    if (ierr/=0) then
        write(0,*) "ERROR DURING ALLOCATION"
        stop 1
    end if
    allocate(lwindow_loc(num_bands_loc,num_kpts_loc),stat=ierr)
    if (ierr/=0) then
        write(0,*) "ERROR DURING ALLOCATION"
        stop 1
    end if
    allocate(wann_centres_loc(3,num_wann_loc),stat=ierr)
    if (ierr/=0) then
        write(0,*) "ERROR DURING ALLOCATION"
        stop 1
    end if
    allocate(wann_spreads_loc(num_wann_loc),stat=ierr)
    if (ierr/=0) then
        write(0,*) "ERROR DURING ALLOCATION"
        stop 1
    end if

    print*, "SECOND ALLOCATION PHASE COMPLETED."

    call wannier_run(seed__name,mp_grid_loc,num_kpts_loc, &
        real_lattice_loc,recip_lattice_loc,kpt_latt_loc,num_bands_loc, &
        num_wann_loc,nntot_loc,num_atoms_loc,atom_symbols_loc, &
        atoms_cart_loc,gamma_only_loc,M_matrix_loc,A_matrix_loc,eigenvalues_loc, &
        U_matrix_loc,U_matrix_opt_loc,lwindow_loc,wann_centres_loc, &
        wann_spreads_loc,spread_loc)

    print*, "WANNIER_RUN CALLED."

end program test_library
