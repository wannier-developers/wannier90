!-*- mode: F90; mode: font-lock; column-number-mode: true -*-!
! Copyright (C) 2004,2006 Jonathan Yates, Arash Mostofi,
!            Nicola Marzari, Ivo Souza, David Vanderbilt
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------

module parameters

  use constants, only : dp
  use io,        only : stdout,maxlen

  implicit none

  private

  !Input
  integer,           public, save :: iprint
  character(len=20), public, save :: energy_unit
  character(len=20), public, save :: length_unit
  logical,           public, save :: wvfn_formatted
  integer,           public, save :: spin
  integer,           public, save :: num_bands
  integer,           public, save :: num_dump_cycles
  integer,           public, save :: num_print_cycles
  character(len=50), public, save :: devel_flag
  integer, allocatable, public,save :: exclude_bands(:)  
  integer,           public, save :: num_wann
  integer,           public, save :: mp_grid(3)
  logical,           public, save :: automatic_mp_grid
  real(kind=dp),     public, save :: dis_win_min
  real(kind=dp),     public, save :: dis_win_max
  real(kind=dp),     public, save :: dis_froz_min
  real(kind=dp),     public, save :: dis_froz_max
  integer,           public, save :: dis_num_iter
  real(kind=dp),     public, save :: dis_mix_ratio
  real(kind=dp),     public, save :: dis_conv_tol
  integer,           public, save :: dis_conv_window
  integer,           public, save :: num_iter
  integer,           public, save :: num_cg_steps
  real(kind=dp),     public, save :: conv_tol
  integer,           public, save :: conv_window
  logical,           public, save :: wannier_plot
  integer, allocatable, public,save :: wannier_plot_list(:)
  integer,           public, save :: wannier_plot_supercell
  character(len=20), public, save :: wannier_plot_format
  logical,           public, save :: bands_plot
  integer,           public, save :: bands_num_points
  character(len=20), public, save :: bands_plot_format
  logical,           public, save :: fermi_surface_plot
  integer,           public, save :: fermi_surface_num_points
  character(len=20), public, save :: fermi_surface_plot_format
  real(kind=dp),     public, save :: fermi_energy
  logical,           public, save :: slice_plot
  integer,           public, save :: slice_num_points
  character(len=20), public, save :: slice_plot_format
  logical,           public, save :: dos_plot
  integer,           public, save :: dos_num_points
  real(kind=dp),     public, save :: dos_energy_step
  real(kind=dp),     public, save :: dos_gaussian_width
  character(len=20), public, save :: dos_plot_format
  integer,           public, save :: num_nnmax    !expert  - max nearest neighbours
  integer,           public, save :: num_shells
  integer, allocatable, public,save :: shell_list(:)
  integer,           public, save :: insign
  real(kind=dp), allocatable,    public, save :: kpt_latt(:,:) !kpoints in lattice vecs
  real(kind=dp),     public, save :: real_lattice(3,3)
  logical,           public, save :: postproc_setup
  logical,           public, save :: cp_pp ! Car-Parinello post-proc flag
  logical,           public, save :: calc_only_A
  character(len=20), public, save :: restart
  logical,           public, save :: write_r2mn
  logical,           public, save :: guiding_centres
  integer,           public, save :: num_guide_cycles
  integer,           public, save :: num_no_guide_iter

  ! Restarts
  real(kind=dp),     public, save :: omega_invariant
  character(len=20), public, save :: checkpoint
  logical,           public, save :: have_disentangled

  ! Atom sites
  real(kind=dp), allocatable,     public, save :: atoms_pos_frac(:,:,:)
  real(kind=dp), allocatable,     public, save :: atoms_pos_cart(:,:,:)
  integer, allocatable,           public, save :: atoms_species_num(:)  
  character(len=2), allocatable,  public, save :: atoms_label(:)
  integer,                        public, save :: num_atoms
  integer,                        public, save :: num_species

  ! Projections
  real(kind=dp), allocatable,     public, save :: proj_site(:,:)
  integer, allocatable,           public, save :: proj_l(:)  
  integer, allocatable,           public, save :: proj_m(:)  
  real(kind=dp), allocatable,     public, save :: proj_z(:,:)
  real(kind=dp), allocatable,     public, save :: proj_x(:,:)
  integer, allocatable,           public, save :: proj_radial(:)  
  real(kind=dp), allocatable,     public, save :: proj_zona(:)
  real(kind=dp), allocatable,     public, save :: proj_box(:)


  !parameters dervied from input
  integer,           public, save :: num_kpts
  real(kind=dp), allocatable,    public, save ::wtkpt(:)
  real(kind=dp),     public, save :: recip_lattice(3,3)
  real(kind=dp),     public, save :: cell_volume
  real(kind=dp),     public, save :: real_metric(3,3)
  real(kind=dp),     public, save :: recip_metric(3,3)
  integer,           public, save :: bands_num_spec_points  
  character(len=1), allocatable,    public, save ::bands_label(:)
  real(kind=dp), allocatable,    public, save ::bands_spec_points(:,:)
  real(kind=dp), allocatable,    public, save ::kpt_cart(:,:) !kpoints in cartesians
  logical,           public, save :: disentanglement
  real(kind=dp),     public, save :: lenconfac
  integer,           public, save :: num_wannier_plot
  integer,           public, save :: num_exclude_bands


  ! kmesh parameters (set in kmesh)

  integer,       public, save              :: nnh           ! the number of b-directions (bka)
  integer,       public, save              :: nntot         ! total number of neighbours for each k-point
  integer,       public, save, allocatable :: nnlist(:,:)   ! list of neighbours for each k-point
  integer,       public, save, allocatable :: neigh(:,:)    
  integer,       public, save, allocatable :: nncell(:,:,:) ! gives BZ of each neighbour of each k-point
  real(kind=dp), public, save              :: wbtot
  real(kind=dp), public, save, allocatable :: wb(:)       ! weights associated with neighbours of each k-point
  real(kind=dp), public, save, allocatable :: bk(:,:,:)     ! the b-vectors that go from each k-point to its neighbours
  real(kind=dp), public, save, allocatable :: bka(:,:)      ! the b-directions from 1st k-point to its neighbours


  ! disentangle parameters
  integer, public, save, allocatable :: ndimwin(:)
  logical, public, save, allocatable :: lwindow(:,:)
  logical, public, save :: frozen_states

  ! a_matrix and m_matrix_orig can be calculated internally from bloch states
  ! or read in from an ab-initio grid
  ! a_matrix      = projection of trial orbitals on bloch states
  ! m_matrix_orig = overlap of bloch states

  complex(kind=dp), allocatable, save, public :: a_matrix(:,:,:)
  complex(kind=dp), allocatable, save, public :: m_matrix_orig(:,:,:,:)
  real(kind=dp),    allocatable, save, public :: eigval(:,:)

  ! u_matrix_opt gives the num_wann dimension optimal subspace from the
  ! original bloch states

  complex(kind=dp), allocatable, save, public :: u_matrix_opt(:,:,:)

  ! u_matrix gives the unitary rotations from the optimal subspace to the
  ! optimally smooth states. 
  ! m_matrix we store here, becuase it is needed for restart of wannierise

  complex(kind=dp), allocatable, save, public :: u_matrix(:,:,:)
  complex(kind=dp), allocatable, save, public :: m_matrix(:,:,:,:)

  ! The maximum number of shells we need to satisfy B1 condition in kmesh
  integer, parameter, public :: max_shells=6

  !private data
  integer                            :: num_lines
  character(len=maxlen), allocatable :: in_data(:)


  public :: param_read
  public :: param_write
  public :: param_dealloc
  public :: param_write_header
  public :: param_write_um
  public :: param_read_um
  public :: param_write_chkpt
  public :: param_read_chkpt

contains



  !==================================================================!
  subroutine param_read ( )
  !==================================================================!
  !                                                                  !
  ! Read parameters and calculate derived values                     !
  !                                                                  !
  !===================================================================  
    use constants, only : bohr   
    use utility, only : utility_recip_lattice,utility_compute_metric
    use io,      only : io_error,io_file_unit,seedname
    implicit none

    !local variables
    integer :: nkp,i,j,n,k,i_temp,i_temp2,unit,loop
    logical :: found,found2,eig_found,lunits
    character(len=6) :: spin_str

    call param_in_file

    !%%%%%%%%%%%%%%%%
    !System variables
    !%%%%%%%%%%%%%%%%

    iprint          =  1             ! Verbosity
    call param_get_keyword('iprint',found,i_value=iprint)

    energy_unit     =  'ev'          !
    call param_get_keyword('energy_unit',found,c_value=energy_unit)

    length_unit     =  'ang'         !
    lenconfac=1.0_dp
    call param_get_keyword('length_unit',found,c_value=length_unit)
    if (length_unit.ne.'ang' .and. length_unit.ne.'bohr') &
         call io_error('Error: value of length_unit not recognised in param_read')
    if (length_unit.eq.'bohr') lenconfac=1.0_dp/bohr

    wvfn_formatted  =  .false.       ! formatted or "binary" file
    call param_get_keyword('wvfn_formatted',found,l_value=wvfn_formatted)

    spin=1
    call param_get_keyword('spin',found,c_value=spin_str)
    if(found) then
       if(index(spin_str,'up')>0) then
          spin=1
       elseif(index(spin_str,'down')>0) then
          spin=2
       else
          call io_error('Error: unrecognised value of spin found: '//trim(spin_str))
       end if
    end if

    num_wann      =   -99
    call param_get_keyword('num_wann',found,i_value=num_wann)
    if (.not. found) then
       call io_error('Error: You must specify num_wann')
    elseif(num_wann<=0) then
       call io_error('Error: num_wann must be greater than zero')
    endif

    num_bands       =   -1   
    call param_get_keyword('num_bands',found,i_value=num_bands)
    if(.not.found) num_bands=num_wann
    if(found .and. num_bands<num_wann) then
       call io_error('Error: num_bands must be greater than or equal to num_wann')
    endif

    num_dump_cycles =   100          ! frequency to write backups at
    call param_get_keyword('num_dump_cycles',found,i_value=num_dump_cycles)
    if (num_dump_cycles<0) call io_error('Error: num_dump_cycles must be positive')       

    num_print_cycles =   1          ! frequency to write at
    call param_get_keyword('num_print_cycles',found,i_value=num_print_cycles)
    if (num_print_cycles<0) call io_error('Error: num_print_cycles must be positive')       

    devel_flag      =   ' '          !       
    call param_get_keyword('devel_flag',found,c_value=devel_flag)


    num_exclude_bands=0
    call param_get_vector_length('exclude_bands',found,num_exclude_bands)
    if(found) then
       if(num_exclude_bands<1) call io_error('Error: problem reading exclude_bands')
       allocate(exclude_bands(num_exclude_bands))
       call param_get_keyword_vector('exclude_bands',found,num_exclude_bands,i_value=exclude_bands)
       if (any(wannier_plot_list<1)  ) &
            call io_error('Error: exclude_bands must contain positive numbers')
    end if

    mp_grid=-99
    call param_get_keyword_vector('mp_grid',found,3,i_value=mp_grid)
    if (any(mp_grid(:)==-99)) then
       call io_error('Error: You must specify dimensions of the Monkhorst-Pack grid by setting mp_grid')
    elseif (any(mp_grid<1)) then
       call io_error('Error: mp_grid must be greater than zero')
    end if
    num_kpts= mp_grid(1)*mp_grid(2)*mp_grid(3)

    automatic_mp_grid = .false.
    call param_get_keyword('automatic_mp_grid',found,l_value=automatic_mp_grid)

    postproc_setup = .false.            ! set to true to write .nnkp file and exit
    call param_get_keyword('postproc_setup',found,l_value=postproc_setup)

    cp_pp = .false.                  ! set to true if doing CP post-processing
    call param_get_keyword('cp_pp',found,l_value=cp_pp)

    calc_only_A = .false.
    call param_get_keyword('calc_only_A',found,l_value=calc_only_A)

    restart = ' '
    call param_get_keyword('restart',found,c_value=restart)

    write_r2mn = .false.
    call param_get_keyword('write_r2mn',found,l_value=write_r2mn)

    !%%%%%%%%%%%
    ! Wannierise
    !%%%%%%%%%%%

    num_iter          = 100    
    call param_get_keyword('num_iter',found,i_value=num_iter)
    if (num_iter<0) call io_error('Error: num_iter must be positive')       

    num_cg_steps      =   5
    call param_get_keyword('num_cg_steps',found,i_value=num_cg_steps)
    if (num_cg_steps<0) call io_error('Error: num_cg_steps must be positive')       

    conv_tol=0.0_dp
    call param_get_keyword('conv_tol',found,r_value=conv_tol)
    if (conv_tol<0.d0) call io_error('Error: conv_tol must be positive')

    conv_window=3
    call param_get_keyword('conv_window',found,i_value=conv_window)
    if (conv_window<0) call io_error('Error: conv_window must be positive')

    guiding_centres=.false.
    call param_get_keyword('guiding_centres',found,l_value=guiding_centres)

    num_guide_cycles=1
    call param_get_keyword('num_guide_cycles',found,i_value=num_guide_cycles)
    if (num_guide_cycles<0) call io_error('Error: num_guide_cycles must be >= 0')

    num_no_guide_iter=0
    call param_get_keyword('num_no_guide_iter',found,i_value=num_no_guide_iter)
    if (num_no_guide_iter<0) call io_error('Error: num_no_guide_iter must be >= 0')
    

    !%%%%%%%%%
    ! Plotting
    !%%%%%%%%%

    wannier_plot              = .false.
    call param_get_keyword('wannier_plot',found,l_value=wannier_plot)

    wannier_plot_supercell    = 2
    call param_get_keyword('wannier_plot_supercell',found,i_value=wannier_plot_supercell)
    if (wannier_plot_supercell<0) call io_error('Error: wannier_plot_supercell must be positive')       

    wannier_plot_format       = 'xcrysden'
    call param_get_keyword('wannier_plot_format',found,c_value=wannier_plot_format)

    call param_get_vector_length('wannier_plot_list',found,num_wannier_plot)
    if(found) then
       if(num_wannier_plot<1) call io_error('Error: problem reading wannier_plot_list')
       allocate(wannier_plot_list(num_wannier_plot))
       call param_get_keyword_vector('wannier_plot_list',found,num_wannier_plot,i_value=wannier_plot_list)
       if (any(wannier_plot_list<1) .or. any(wannier_plot_list>num_wann) ) &
            call io_error('Error: wannier_plot_list asks for a non-valid wannier function to be plotted')
    else
       ! we plot all wannier functions
       num_wannier_plot=num_wann
       allocate(wannier_plot_list(num_wannier_plot))
       do loop=1,num_wann
          wannier_plot_list(loop)=loop
       end do
    end if


    bands_plot                = .false.
    call param_get_keyword('bands_plot',found,l_value=bands_plot)

    bands_num_points          = 100
    call param_get_keyword('bands_num_points',found,i_value=bands_num_points)
    if (bands_num_points<0) call io_error('Error: bands_num_points must be positive')       

    bands_plot_format         = 'gnuplot'
    call param_get_keyword('bands_plot_format',found,c_value=bands_plot_format)

    call param_get_block_length('kpoint_path',found,i_temp)
    if (found) then
       bands_num_spec_points=i_temp*2
       allocate(bands_label(bands_num_spec_points))
       allocate(bands_spec_points(3,bands_num_spec_points))
       call param_get_keyword_kpath
    end if
    if(.not.found .and. bands_plot) &
         call io_error('A bandstructure plot has been requested but there is no kpoint_path block') 

    fermi_surface_plot        =  .false.
    call param_get_keyword('fermi_surface_plot',found,l_value=fermi_surface_plot)

    fermi_surface_num_points  = 50
    call param_get_keyword('fermi_surface_num_points',found,i_value=fermi_surface_num_points)
    if (fermi_surface_num_points<0) call io_error('Error: fermi_surface_num_points must be positive')       
    fermi_surface_plot_format = 'xcrysden'
    call param_get_keyword('fermi_surface_plot_format',found,c_value=fermi_surface_plot_format)

    fermi_energy=0.0_dp
    call param_get_keyword('fermi_energy',found,r_value=fermi_energy)

    slice_plot                = .false.
    call param_get_keyword('slice_plot',found,l_value=slice_plot)

    slice_num_points          = 50
    call param_get_keyword('slice_num_points',found,i_value=slice_num_points)
    if (slice_num_points<0) call io_error('Error: slice_num_points must be positive')       

    slice_plot_format         = 'plotmv'
    call param_get_keyword('slice_plot_format',found,c_value=slice_plot_format)

    dos_plot                  = .false.
    call param_get_keyword('dos_plot',found,l_value=dos_plot)

    dos_num_points            = 50
    call param_get_keyword('dos_num_points',found,i_value=dos_num_points)
    if (dos_num_points<0) call io_error('Error: dos_num_points must be positive')       

    dos_energy_step           = 0.01_dp
    call param_get_keyword('dos_energy_step',found,r_value=dos_energy_step)

    dos_gaussian_width        = 0.1_dp
    call param_get_keyword('dos_gaussian_width',found,r_value=dos_gaussian_width)

    dos_plot_format           = 'gnuplot'
    call param_get_keyword('dos_plot_format',found,c_value=dos_plot_format)

    !%%%%%%%%%%%%%%%%
    ! Disentanglement
    !%%%%%%%%%%%%%%%%

    disentanglement=.false.
    if(num_bands>num_wann) disentanglement=.true.

    ! Read the eigenvalues from wannier.eig
    allocate(eigval(num_bands,num_kpts))

    inquire(file=trim(seedname)//'.eig',exist=eig_found)
    if(.not. eig_found) then
       if ( disentanglement.and.(.not.postproc_setup) ) then
          call io_error('No '//trim(seedname)//'.eig file found. Needed for disentanglement')
       else if ((bands_plot .or. dos_plot .or. fermi_surface_plot) .and.(.not.postproc_setup) ) then
          call io_error('No '//trim(seedname)//'.eig file found. Needed for interpolation')
       end if
    else
       unit=io_file_unit()
       open(unit=unit,file=trim(seedname)//'.eig',form='formatted',status='old',err=105)
       do k=1,num_kpts
          do n=1,num_bands
             read(unit,*) i,j,eigval(i,j)
             if ((i.ne.n).or.(j.ne.k)) then
                call io_error('param_read: mismatch in '//trim(seedname)//'.eig')
             end if
          enddo
       end do
       close(unit)
    end if

    dis_win_min=-1.0_dp;dis_win_max=0.0_dp
    if(eig_found) dis_win_min = minval(eigval)       
    call param_get_keyword('dis_win_min',found,r_value=dis_win_min)

    if(eig_found) dis_win_max = maxval(eigval)
    call param_get_keyword('dis_win_max',found,r_value=dis_win_max)

    if ( dis_win_max.lt.dis_win_min ) &
         call io_error('Error: param_read: check disentanglement windows')

    dis_froz_min=-1.0_dp;dis_froz_max=0.0_dp
    ! no default for dis_froz_max
    frozen_states=.false.
    call param_get_keyword('dis_froz_max',found,r_value=dis_froz_max)
    if (found) then 
       frozen_states=.true.
       dis_froz_min = dis_win_min ! default value for the bottom of frozen window
    end if
    call param_get_keyword('dis_froz_min',found2,r_value=dis_froz_min)
    if ( dis_froz_max.lt.dis_froz_min ) &
       call io_error('Error: param_read: check disentanglement frozen windows')
    if(found2 .and. .not. found) &
           call io_error('Error: param_read: found dis_froz_min but not dis_froz_max')

    dis_num_iter      = 200
    call param_get_keyword('dis_num_iter',found,i_value=dis_num_iter)
    if (dis_num_iter<0) call io_error('Error: dis_num_iter must be positive')       

    dis_mix_ratio     = 0.5d0
    call param_get_keyword('dis_mix_ratio',found,r_value=dis_mix_ratio)
    if (dis_mix_ratio<0.d0) call io_error('Error: dis_mix_ratio must be positive')

    dis_conv_tol      = 0.0d0       
    call param_get_keyword('dis_conv_tol',found,r_value=dis_conv_tol)
    if (dis_conv_tol<0.d0) call io_error('Error: dis_conv_tol must be positive')

    dis_conv_window=3   
    call param_get_keyword('dis_conv_window',found,i_value=dis_conv_window)
    if (dis_conv_window<0) call io_error('Error: dis_conv_window must be positive')       

    !%%%%%%%%%%%%%%%%
    !  Other Stuff
    !%%%%%%%%%%%%%%%%

    num_nnmax                 = 12
    call param_get_keyword('num_nnmax',found,i_value=num_nnmax)
    if (num_nnmax<0) call io_error('Error: num_nnmax must be positive')       

    num_shells                   = 0 
    call param_get_keyword('num_shells',found,i_value=num_shells)
    if(num_shells<0 .or. num_shells>6) call io_error('Error: num_shells must be between zero and six')
    if (num_shells==0) then
       allocate( shell_list(max_shells))
    else
       allocate( shell_list(num_shells))
    end if
    call param_get_keyword_vector('shell_list',found,num_shells,i_value=shell_list)
    if(num_shells==0 .and. found) call io_error('Error: shell_list has no effect when num_shells=0')
    if(num_shells/=0 .and. .not. found) call io_error('Error: shell_list must be set when when num_shells>0')
    if(num_shells/=0 .and. any(shell_list<1)) call io_error('Error: shell_list must be positive')

    insign=1

    call param_get_keyword_block('unit_cell_cart',found,3,3,r_value=real_lattice)
    !This is a hack. I must workout what is the sensible way to read and store this jry
    real_lattice=transpose(real_lattice)

    if(.not. found) call io_error('Error: Did not find the cell information in the input file')

    call utility_recip_lattice(real_lattice,recip_lattice,cell_volume)
    call utility_compute_metric(real_lattice,recip_lattice,real_metric,recip_metric)

    allocate ( kpt_cart(3,num_kpts) )
    allocate ( kpt_latt(3,num_kpts) )
    allocate ( wtkpt(num_kpts)  )

    call param_get_keyword_block('kpoints',found,num_kpts,3,r_value=kpt_latt)
    if(.not. found) call io_error('Error: Did not find the kpoint information in the input file')

    do i=1,num_kpts
       wtkpt(i)=1.d0/real(num_kpts,dp)
    end do

    ! Calculate the kpoints in cartesian coordinates
    do nkp=1,num_kpts
       kpt_cart(:,nkp)=matmul(kpt_latt(:,nkp),recip_lattice(:,:))
    end do

    ! Atoms
    num_atoms=0
    call param_get_block_length('atoms_frac',found,i_temp)
    call param_get_block_length('atoms_cart',found2,i_temp2,lunits)
    if (found .and. found2) call io_error('Error: Cannot specify both atoms_frac and atoms_cart')
    if (found .and. i_temp>0) then
       lunits=.false.
       num_atoms=i_temp
    elseif (found2 .and. i_temp2>0) then
       num_atoms=i_temp2
       if (lunits) num_atoms=num_atoms-1
    end if
    if(num_atoms>0) then
       call param_get_atoms(lunits)
    end if

    call param_get_block_length('projections',found,i_temp)
    if (found) call param_get_projections
    if (guiding_centres .and. .not. found) &
       call io_error('param_read: Guiding centres requested, but no projection block found')

    ! check to see that there are no unrecognised keywords

    if ( any(len_trim(in_data(:))>0 )) then
       write(stdout,'(1x,a)') 'The following section of file '//trim(seedname)//'.win contained unrecognised keywords'
       write(stdout,*) 
       do loop=1,num_lines
          if (len_trim(in_data(loop))>0) then
             write(stdout,'(1x,a)') trim(in_data(loop))
          end if
       end do
       write(stdout,*) 
       call io_error('Unrecognised keyword(s) in input file')
    end if



    ! For aesthetic purposes, convert some things to uppercase
    call param_uppercase()


    deallocate(in_data)

    ! Some checks 
    if (restart.ne.' ') disentanglement=.false.

    if ( (restart.ne.' ').and.(restart.ne.'default') & 
         .and.(restart.ne.'wannierise').and.(restart.ne.'plot') ) &
         call io_error('Error in input file: value of restart not recognised')

    if (disentanglement) allocate(ndimwin(num_kpts))
    if (disentanglement) allocate(lwindow(num_bands,num_kpts))

    ! Initialise
    omega_invariant = -999.0_dp
    have_disentangled = .false.

    return

105 call io_error('Error: Problem opening eigenvalue file '//trim(seedname)//'.eig')

  end subroutine param_read



  !===================================================================
  subroutine param_uppercase
    !===================================================================
    !                                                                  !
    ! Convert a few things to uppercase to look nice in the output     !
    !                                                                  !
    !===================================================================  

    implicit none

    integer :: nsp,ic,loop

    ! Atom labels (eg, si --> Si)
    do nsp=1,num_species
       ic=ichar(atoms_label(nsp)(1:1))                           
       if ((ic.ge.ichar('a')).and.(ic.le.ichar('z'))) &
            atoms_label(nsp)(1:1) = char(ic+ichar('Z')-ichar('z'))
    enddo

    ! Bands labels (eg, x --> X)
    do loop=1,bands_num_spec_points
       ic=ichar(bands_label(loop))                           
       if ((ic.ge.ichar('a')).and.(ic.le.ichar('z'))) &
            bands_label(loop) = char(ic+ichar('Z')-ichar('z'))
    enddo

    ! Length unit (ang --> Ang, bohr --> Bohr)
    ic=ichar(length_unit(1:1))
    if ((ic.ge.ichar('a')).and.(ic.le.ichar('z'))) &
         length_unit(1:1) = char(ic+ichar('Z')-ichar('z'))

    return

  end subroutine param_uppercase


  !===================================================================
  subroutine param_write
    !==================================================================!
    !                                                                  !
    ! write parameters to stdout                                       !
    !                                                                  !
    !===================================================================  

    implicit none

    integer :: i,nkp,loop,nat,nsp,ic

    ! System
    write(stdout,*)
    write(stdout,'(36x,a6)') '------' 
    write(stdout,'(36x,a6)') 'SYSTEM' 
    write(stdout,'(36x,a6)') '------' 
    write(stdout,*)
    if (lenconfac.eq.1.0_dp) then
       write(stdout,'(30x,a21)') 'Lattice Vectors (Ang)' 
    else
       write(stdout,'(28x,a22)') 'Lattice Vectors (Bohr)' 
    endif
    write(stdout,101) 'a_1',(real_lattice(1,I)*lenconfac, i=1,3)
    write(stdout,101) 'a_2',(real_lattice(2,I)*lenconfac, i=1,3)
    write(stdout,101) 'a_3',(real_lattice(3,I)*lenconfac, i=1,3)
    write(stdout,*)   
    write(stdout,'(19x,a17,3x,f11.5))',advance='no') &
         'Unit Cell Volume:',cell_volume*lenconfac**3
    if (lenconfac.eq.1.0_dp) then
       write(stdout,'(2x,a7)') '(Ang^3)'
    else
       write(stdout,'(2x,a8)') '(Bohr^3)'
    endif
    write(stdout,*)   
    if (lenconfac.eq.1.0_dp) then
       write(stdout,'(24x,a33)') 'Reciprocal-Space Vectors (Ang^-1)'
    else
       write(stdout,'(22x,a34)') 'Reciprocal-Space Vectors (Bohr^-1)'
    endif
    write(stdout,101) 'b_1',(recip_lattice(1,I)/lenconfac, i=1,3)
    write(stdout,101) 'b_2',(recip_lattice(2,I)/lenconfac, i=1,3)
    write(stdout,101) 'b_3',(recip_lattice(3,I)/lenconfac, i=1,3)
    write(stdout,*)   ' '
    ! Atoms
    if(num_atoms>0) then
       write(stdout,'(1x,a)') '*----------------------------------------------------------------------------*'
       if (lenconfac.eq.1.0_dp) then
          write(stdout,'(1x,a)') '|   Site       Fractional Coordinate          Cartesian Coordinate (Ang)     |'
       else
          write(stdout,'(1x,a)') '|   Site       Fractional Coordinate          Cartesian Coordinate (Bohr)    |'
       endif
       write(stdout,'(1x,a)') '+----------------------------------------------------------------------------+'
       do nsp=1,num_species
          do nat=1,atoms_species_num(nsp)
             write(stdout,'(1x,a1,1x,a2,1x,i3,3F10.5,3x,a1,1x,3F10.5,4x,a1)') '|',atoms_label(nsp),nat,atoms_pos_frac(:,nat,nsp),&
                  '|',atoms_pos_cart(:,nat,nsp)*lenconfac,'|'
          end do
       end do
       write(stdout,'(1x,a)') '*----------------------------------------------------------------------------*'
    else
       write(stdout,'(25x,a)') 'No atom positions specified'
    end if
    write(stdout,*) ' '
    ! Projections
    if(iprint>1 .and. allocated(proj_site) ) then
       write(stdout,'(32x,a)') '-----------'
       write(stdout,'(32x,a)') 'PROJECTIONS'
       write(stdout,'(32x,a)') '-----------'
       write(stdout,*) ' '
       write(stdout,'(1x,a)') '+----------------------------------------------------------------------------+'
       write(stdout,'(1x,a)') '|    Frac. Coord.     l mr r        z-axis            x-axis        diff box |'
       write(stdout,'(1x,a)') '+----------------------------------------------------------------------------+'
       do nsp=1,num_wann
          write(stdout,'(1x,a1,3(1x,f5.2),2x,i2,1x,i2,1x,i1,1x,3(1x,f5.2),1x,3(1x,f5.2),&
               & 2x,f4.1,1x,f3.1,1x,a1)')  '|',proj_site(1,nsp),proj_site(2,nsp),&
               proj_site(3,nsp),proj_l(nsp), proj_m(nsp),proj_radial(nsp),&
               proj_z(1,nsp),proj_z(2,nsp),proj_z(3,nsp),proj_x(1,nsp),&
               proj_x(2,nsp),proj_x(3,nsp),proj_zona(nsp), proj_box(nsp),'|'
       end do
       write(stdout,'(1x,a)') '+----------------------------------------------------------------------------+'
       write(stdout,*) ' '
    end if

    ! K-points
    write(stdout,'(32x,a)') '------------'
    write(stdout,'(32x,a)') 'K-POINT GRID'
    write(stdout,'(32x,a)') '------------'
    write(stdout,*) ' '
    write(stdout,'(13x,a,i3,1x,a1,i3,1x,a1,i3,6x,a,i5)') 'Grid size =',mp_grid(1),'x',mp_grid(2),'x',mp_grid(3),&
         'Total points =',num_kpts
    write(stdout,*) ' '
    if(iprint>1) then
       write(stdout,'(1x,a)') '*----------------------------------------------------------------------------*'
       if (lenconfac.eq.1.0_dp) then
          write(stdout,'(1x,a)') '| k-point      Fractional Coordinate        Cartesian Coordinate (Ang^-1)    |'
       else
          write(stdout,'(1x,a)') '| k-point      Fractional Coordinate        Cartesian Coordinate (Bohr^-1)   |'
       endif
       write(stdout,'(1x,a)') '+----------------------------------------------------------------------------+'
       do nkp=1,num_kpts
          write(stdout,'(1x,a1,i6,1x,3F10.5,3x,a1,1x,3F10.5,4x,a1)') '|',nkp,kpt_latt(:,nkp),'|',kpt_cart(:,nkp)/lenconfac,'|'
       end do
       write(stdout,'(1x,a)') '*----------------------------------------------------------------------------*'
       write(stdout,*) ' '
    end if
    ! Main
    write(stdout,*) ' '
    write(stdout,'(1x,a78)') '*---------------------------------- MAIN ------------------------------------*'
    write(stdout,'(1x,a46,10x,I8,13x,a1)') '|  Output verbosity (1=low, 5=high)          :',iprint,'|'
    write(stdout,'(1x,a46,10x,a8,13x,a1)') '|  Length Unit                               :',length_unit,'|'  
    write(stdout,'(1x,a78)') '*----------------------------------------------------------------------------*'

    ! Wannierise
    write(stdout,*) ' '
    write(stdout,'(1x,a78)') '*------------------------------- WANNIERISE ---------------------------------*'
    write(stdout,'(1x,a46,10x,I8,13x,a1)')   '|  Total number of iterations                :',num_iter,'|'
    write(stdout,'(1x,a46,10x,I8,13x,a1)')   '|  Number of CG steps before reset           :',num_cg_steps,'|'
    write(stdout,'(1x,a46,8x,E10.3,13x,a1)') '|  Convergence tolerence                     :',conv_tol,'|'
    write(stdout,'(1x,a46,10x,I8,13x,a1)')   '|  Convergence window                        :',conv_window,'|'
    write(stdout,'(1x,a46,10x,I8,13x,a1)')   '|  Iterations between writing output         :',num_print_cycles,'|'
    write(stdout,'(1x,a46,10x,I8,13x,a1)')   '|  Iterations between backing up to disk     :',num_dump_cycles,'|'
    write(stdout,'(1x,a78)') '*----------------------------------------------------------------------------*'
    !
    ! Disentanglement
    !
    if (disentanglement .or. iprint>2) then
       write(stdout,'(1x,a78)') '*------------------------------- DISENTANGLE --------------------------------*'
       write(stdout,'(1x,a46,10x,L8,13x,a1)')   '|  Using band disentanglement                :',disentanglement,'|'
       write(stdout,'(1x,a46,10x,I8,13x,a1)')   '|  Total number of iterations                :',dis_num_iter,'|'
       write(stdout,'(1x,a46,10x,F8.3,13x,a1)') '|  Mixing ratio                              :',dis_mix_ratio,'|'
       write(stdout,'(1x,a46,8x,E10.3,13x,a1)') '|  Convergence tolerence                     :',dis_conv_tol,'|'
       write(stdout,'(1x,a46,10x,I8,13x,a1)')   '|  Convergence window                        :',dis_conv_window,'|'
       write(stdout,'(1x,a78)') '*----------------------------------------------------------------------------*'
    end if
    !
    ! Plotting
    !
    if (wannier_plot .or. bands_plot .or. fermi_surface_plot .or. slice_plot &
         .or. dos_plot .or. iprint>2) then
       !
       write(stdout,'(1x,a78)') '*-------------------------------- PLOTTING ----------------------------------*'
       !
       if (wannier_plot .or. iprint>2) then
          write(stdout,'(1x,a46,10x,L8,13x,a1)') '|  Plotting Wannier functions                :',wannier_plot,'|'
          write(stdout,'(1x,a46,10x,I8,13x,a1)') '|   Size of supercell for plotting           :',wannier_plot_supercell,'|'
          write(stdout,'(1x,a46,10x,a8,13x,a1)') '|   Plotting format                          :',wannier_plot_format,'|'
          write(stdout,'(1x,a78)') '*----------------------------------------------------------------------------*'
       end if
       !
       if (fermi_surface_plot .or. iprint>2) then
          write(stdout,'(1x,a46,10x,L8,13x,a1)') '|  Plotting Fermi surface                    :',fermi_surface_plot,'|'
          write(stdout,'(1x,a46,10x,I8,13x,a1)') '|   Number of plotting points (along b_1)    :',fermi_surface_num_points,'|'
          write(stdout,'(1x,a46,10x,a8,13x,a1)') '|   Plotting format                          :',fermi_surface_plot_format,'|'
          write(stdout,'(1x,a78)') '*----------------------------------------------------------------------------*'
       end if
       !
       if (bands_plot .or. iprint>2) then
          write(stdout,'(1x,a46,10x,L8,13x,a1)') '|  Plotting interpolated bandstructure       :',bands_plot,'|'
          write(stdout,'(1x,a46,10x,I8,13x,a1)') '|   Number of K-path sections                :',bands_num_spec_points/2,'|'
          write(stdout,'(1x,a46,10x,I8,13x,a1)') '|   Divisions along first K-path section     :',bands_num_points,'|'
          write(stdout,'(1x,a46,10x,a8,13x,a1)') '|   Output format                            :',bands_plot_format,'|'
          write(stdout,'(1x,a78)') '|   K-space path sections:                                                   |'
          if(bands_num_spec_points==0) then
             write(stdout,'(1x,a78)') '|     None defined                                                           |'
          else
             do loop=1,bands_num_spec_points,2
                write(stdout,'(1x,a10,2x,a1,2x,3F7.3,5x,a3,2x,a1,2x,3F7.3,7x,a1)') '|    From:',bands_label(loop),&
                     (bands_spec_points(i,loop),i=1,3),'To:',bands_label(loop+1),(bands_spec_points(i,loop+1),i=1,3),'|'
             end do
          end if
          write(stdout,'(1x,a78)') '*----------------------------------------------------------------------------*'
       end if
       !
    endif


101 format(20x,a3,2x,3F11.6)

  end subroutine param_write

  subroutine param_write_header
    use io, only : io_date
    implicit none


    character (len=9) :: cdate, ctime

    call io_date(cdate, ctime)

    write(stdout,*)
    write(stdout,*)  '            +---------------------------------------------------+'
    write(stdout,*)  '            |                                                   |'
    write(stdout,*)  '            |                    WANNIER                        |'
    write(stdout,*)  '            |                                                   |'
    write(stdout,*)  '            +---------------------------------------------------+'
    write(stdout,*)  '            |                                                   |'
    write(stdout,*)  '            |        Welcome to the Maximally-Localized         |'
    write(stdout,*)  '            |        Generalized Wannier Functions code         |'
    write(stdout,*)  '            |            http://www.wannier.org                 |'
    write(stdout,*)  '            |                                                   |'
    write(stdout,*)  '            |  Authors:                                         |'
    write(stdout,*)  '            |           Arash Mostofi   (MIT)                   |'
    write(stdout,*)  '            |           Jonathan Yates  (LBNL and UC Berkeley)  |'
    write(stdout,*)  '            |                                                   |'
    write(stdout,*)  '            |  Based on the codes written by                    |'
    write(stdout,*)  '            |   Nicola Marzari, Ivo Souza and David Vanderbilt  |'
    write(stdout,*)  '            |  and the methods detailed in                      |'
    write(stdout,*)  '            |                                                   |'
    write(stdout,*)  '            |  "Maximally Localised Generalised Wannier         |'
    write(stdout,*)  '            |   Functions for Composite Energy Bands"           |'
    write(stdout,*)  '            |  N. Marzari and D. Vanderbilt PRB 56 12847 (1997) |'
    write(stdout,*)  '            |                                                   |'
    write(stdout,*)  '            |  "Maximally Localised Wannier                     |'
    write(stdout,*)  '            |   Functions for Entangled Energy Bands"           |'
    write(stdout,*)  '            |  I. Souza, N. Marzari and D. Vanderbilt           |'
    write(stdout,*)  '            |  PRB 65 035109 (2001)                             |'
    write(stdout,*)  '            |                                                   |'
    write(stdout,*)  '            | Copyright (c) 1997-2006 J. Yates, A. Mostofi,     |'
    write(stdout,*)  '            |          N. Marzari, I. Souza, D. Vanderbilt      |'
    write(stdout,*)  '            |                                                   |'
    write(stdout,*)  '            |       Release: 1.91  18th November 2005           |'
    write(stdout,*)  '            |                                                   |'
    write(stdout,*)  '            | This program is free software; you can            |'
    write(stdout,*)  '            | redistribute it and/or modify it under the terms  |'
    write(stdout,*)  '            | of the GNU General Public License as published by |'
    write(stdout,*)  '            | the Free Software Foundation; either version 2 of |'
    write(stdout,*)  '            | the License, or (at your option) any later version|'
    write(stdout,*)  '            |                                                   |'
    write(stdout,*)  '            | This program is distributed in the hope that it   |'
    write(stdout,*)  '            | will be useful, but WITHOUT ANY WARRANTY; without |'
    write(stdout,*)  '            | even the implied warranty of MERCHANTABILITY or   |'
    write(stdout,*)  '            | FITNESS FOR A PARTICULAR PURPOSE. See the GNU     |'
    write(stdout,*)  '            | General Public License for more details.          |'
    write(stdout,*)  '            |                                                   |'
    write(stdout,*)  '            | You should have received a copy of the GNU General|'
    write(stdout,*)  '            | Public License along with this program; if not,   |'
    write(stdout,*)  '            | write to the Free Software Foundation, Inc.,      |'
    write(stdout,*)  '            | 675 Mass Ave, Cambridge, MA 02139, USA.           |'
    write(stdout,*)  '            |                                                   |'
    write(stdout,*)  '            +---------------------------------------------------+'
    write(stdout,*)  '            |    Execution started on ',cdate,' at ',ctime,'    |'
    write(stdout,*)  '            +---------------------------------------------------+'

  end subroutine param_write_header


  !==================================================================!
  subroutine param_dealloc
    !==================================================================!
    !                                                                  !
    ! release memory from allocated parameters                         !
    !                                                                  !
    !===================================================================  
    use io, only : io_error

    implicit none
    integer :: ierr

    if ( allocated ( ndimwin ) ) then
       deallocate (  ndimwin, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating ndimwin in param_dealloc')
    end if
    if ( allocated ( lwindow ) ) then
       deallocate (  lwindow, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating lwindow in param_dealloc')
    end if
    deallocate ( eigval, stat=ierr  )
    if (ierr/=0) call io_error('Error in deallocating eigval in param_dealloc')
    deallocate ( shell_list, stat=ierr  )
    if (ierr/=0) call io_error('Error in deallocating shell_list in param_dealloc')
    deallocate ( wtkpt , stat=ierr  )
    if (ierr/=0) call io_error('Error in deallocating wtkpt in param_dealloc')
    deallocate ( kpt_latt, stat=ierr  )
    if (ierr/=0) call io_error('Error in deallocating kpt_latt in param_dealloc')
    deallocate ( kpt_cart, stat=ierr  )
    if (ierr/=0) call io_error('Error in deallocating kpt_cart in param_dealloc')
    if ( allocated ( bands_label ) ) then
       deallocate (  bands_label, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating bands_label in param_dealloc')
    end if
    if ( allocated ( bands_spec_points ) ) then
       deallocate (  bands_spec_points, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating bands_spec_points in param_dealloc')
    end if
    if ( allocated ( atoms_label ) ) then
       deallocate (  atoms_label, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating atoms_label in param_dealloc')
    end if
    if ( allocated ( atoms_pos_frac ) ) then
       deallocate (  atoms_pos_frac, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating atom_pos_frac in param_dealloc')
    end if
    if ( allocated ( atoms_pos_cart ) ) then
       deallocate (  atoms_pos_cart, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating atoms_pos_cart in param_dealloc')
    end if
    if ( allocated ( atoms_species_num ) ) then
       deallocate (atoms_species_num, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating atoms_species_num in param_dealloc')
    end if
    if ( allocated( proj_site ) ) then
       deallocate( proj_site, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating proj_site in param_dealloc')
    end if
    if ( allocated( proj_l ) ) then
       deallocate( proj_l, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating proj_l in param_dealloc')
    end if
    if ( allocated( proj_m ) ) then
       deallocate( proj_m, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating proj_m in param_dealloc')
    end if
    if ( allocated( proj_z ) ) then
       deallocate( proj_z, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating proj_z in param_dealloc')
    end if
    if ( allocated( proj_x ) ) then
       deallocate( proj_x, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating proj_x in param_dealloc') 
    end if
    if ( allocated( proj_radial ) ) then
       deallocate( proj_radial, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating proj_radial in param_dealloc')
    end if
    if ( allocated( proj_zona ) ) then
       deallocate( proj_zona, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating proj_zona in param_dealloc')
    end if
    if ( allocated( proj_box ) ) then
       deallocate( proj_box, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating proj_box in param_dealloc')
    end if

    return

  end subroutine param_dealloc


  !================================!
  subroutine param_write_um
    !================================!
    !                                !
    ! Dump the U and M to *_um.dat   !
    !                                !
    !================================!


    use io,        only : io_file_unit,io_error,seedname,io_date
    implicit none

    integer :: i,j,k,l,um_unit
    character (len=9) :: cdate, ctime
    character(len=33) :: header

    call io_date(cdate, ctime)
    header='written on '//cdate//' at '//ctime

    um_unit=io_file_unit()
    open(unit=um_unit,file=trim(seedname)//'_um.dat',form='unformatted')
    write(um_unit) header
    write(um_unit) omega_invariant
    write(um_unit) num_wann,num_kpts,num_nnmax    
    write(um_unit) (((u_matrix(i,j,k),i=1,num_wann),j=1,num_wann),k=1,num_kpts)
    write(um_unit) ((((m_matrix(i,j,k,l),i=1,num_wann),j=1,num_wann),k=1,nntot),l=1,num_kpts)
    close(um_unit)

    return

  end subroutine param_write_um


  !================================!
  subroutine param_read_um
    !================================!
    !                                !
    ! Restore U and M from file      !
    !                                !
    !================================!

    use io,        only : io_file_unit,io_error,seedname
    implicit none

    integer       :: tmp_num_wann,tmp_num_kpts,tmp_num_nnmax    
    integer       :: i,j,k,l,um_unit,ierr
    character(len=33) :: header
    real(kind=dp) :: tmp_omi

    um_unit=io_file_unit()
    open(unit=um_unit,file=trim(seedname)//'_um.dat',status="old",form='unformatted',err=105)
    read(um_unit) header
    write(stdout,'(1x,4(a))') 'Reading U and M from file ',trim(seedname),'_um.dat ', header 
    read(um_unit) tmp_omi
    if ( have_disentangled ) then
       if ( abs(tmp_omi-omega_invariant).gt.1.0e-10_dp )  &
            call io_error('Error in restart: omega_invariant in .chk and um.dat files do not match')
    endif
    read(um_unit) tmp_num_wann,tmp_num_kpts,tmp_num_nnmax    
    if(tmp_num_wann/=num_wann) call io_error('Error in param_read_um: num_wann mismatch')
    if(tmp_num_kpts/=num_kpts) call io_error('Error in param_read_um: num_kpts mismatch')
    if(tmp_num_nnmax/=num_nnmax) call io_error('Error in param_read_um: num_nnmax mismatch')
    if (.not.allocated(u_matrix)) then
       allocate(u_matrix(num_wann,num_wann,num_kpts),stat=ierr)
       if (ierr/=0) call io_error('Error allocating u_matrix in param_read_um')
    endif
    read(um_unit) (((u_matrix(i,j,k),i=1,num_wann),j=1,num_wann),k=1,num_kpts)
    if (.not.allocated(m_matrix)) then
       allocate(m_matrix(num_wann,num_wann,nntot,num_kpts),stat=ierr)
       if (ierr/=0) call io_error('Error allocating m_matrix in param_read_um')
    endif
    read(um_unit) ((((m_matrix(i,j,k,l),i=1,num_wann),j=1,num_wann),k=1,nntot),l=1,num_kpts)
    close(um_unit)

    return

105 call io_error('Error: Problem opening file '//trim(seedname)//'_um.dat in param_read_um')

  end subroutine param_read_um



  !=================================================!
  subroutine param_write_chkpt(chkpt)
    !=================================================!
    ! Write checkpoint file                           !
    !=================================================!

    use io, only : io_file_unit,io_date,seedname

    implicit none

    character(len=*), intent(in) :: chkpt

    integer :: chk_unit,nkp,nn,i,j
    character (len=9) :: cdate,ctime
    character (len=33) :: header
    character (len=20) :: chkpt1

    write(stdout,'(/1x,3a)',advance='no') 'Writing checkpoint file ',trim(seedname),'.chk...'

    call io_date(cdate,ctime)
    header='written on '//cdate//' at '//ctime

    chk_unit=io_file_unit()
    open(unit=chk_unit,file=trim(seedname)//'.chk',form='unformatted')

    write(chk_unit) header                                   ! Date and time
    write(chk_unit) ((real_lattice(i,j),i=1,3),j=1,3)        ! Real lattice
    write(chk_unit) ((recip_lattice(i,j),i=1,3),j=1,3)       ! Reciprocal lattice
    write(chk_unit) num_kpts
    write(chk_unit) ((kpt_latt(i,nkp),i=1,3),nkp=1,num_kpts) ! K-points
    write(chk_unit) nntot                  ! Number of nearest k-point neighbours
    write(chk_unit) num_wann               ! Number of wannier functions
    chkpt1 = adjustl(trim(chkpt))
    write(chk_unit) chkpt1                 ! Position of checkpoint
    write(chk_unit) have_disentangled      ! Whether a disentanglement has been performed
    if (have_disentangled) then
       write(chk_unit) omega_invariant     ! Omega invariant
       ! lwindow, ndimwin and U_matrix_opt 
       write(chk_unit) ((lwindow(i,nkp),i=1,num_bands),nkp=1,num_kpts)
       write(chk_unit) (ndimwin(nkp),nkp=1,num_kpts)
       write(chk_unit) (((u_matrix_opt(i,j,nkp),i=1,num_bands),j=1,num_wann),nkp=1,num_kpts)
    endif

    close(chk_unit)

    write(stdout,'(a/)') ' done'

    return

  end subroutine param_write_chkpt


  !=======================================!
  subroutine param_read_chkpt
    !=======================================!
    ! Read checkpoint file                  !
    !=======================================!

    use io,      only : io_error,io_file_unit,stdout,seedname
    use utility, only : utility_strip

    implicit none

    integer :: chk_unit,nkp,nn,i,j,ntmp,ierr
    character(len=33) :: header
    real(kind=dp) :: tmp_latt(3,3), tmp_kpt_latt(3,num_kpts)

    write(stdout,'(1x,3a)') 'Reading restart information from file ',trim(seedname),'.chk :'

    chk_unit=io_file_unit()
    open(unit=chk_unit,file=trim(seedname)//'.chk',status='old',form='unformatted',err=121)

    ! Read comment line
    read(chk_unit) header
    write(stdout,'(1x,a)',advance='no') trim(header)

    ! Consistency checks
    read(chk_unit) ((tmp_latt(i,j),i=1,3),j=1,3)  ! Real lattice
    do j=1,3
       do i=1,3
          if (abs(tmp_latt(i,j)-real_lattice(i,j)).gt.1.0e-6_dp) &
               call io_error('param_read_chk: Mismatch in real_lattice')
       enddo
    enddo
    read(chk_unit) ((tmp_latt(i,j),i=1,3),j=1,3)  ! Reciprocal lattice
    do j=1,3
       do i=1,3
          if (abs(tmp_latt(i,j)-recip_lattice(i,j)).gt.1.0e-6_dp) &
               call io_error('param_read_chk: Mismatch in recip_lattice')
       enddo
    enddo
    read(chk_unit) ntmp                ! K-points
    if (ntmp.ne.num_kpts) &
         call io_error('param_read_chk: Mismatch in num_kpts')
    read(chk_unit) ((tmp_kpt_latt(i,nkp),i=1,3),nkp=1,num_kpts)
    do nkp=1,num_kpts
       do i=1,3
          if (abs(tmp_kpt_latt(i,nkp)-kpt_latt(i,nkp)).gt.1.0e-6_dp) &
               call io_error('param_read_chk: Mismatch in kpt_latt')
       enddo
    enddo
    read(chk_unit) ntmp                ! nntot
    if (ntmp.ne.nntot) &
         call io_error('param_read_chk: Mismatch in nntot')
    read(chk_unit) ntmp                ! num_wann
    if (ntmp.ne.num_wann) &
         call io_error('param_read_chk: Mismatch in num_wann')

    ! End of consistency checks

    read(chk_unit) checkpoint             ! checkpoint
    checkpoint=adjustl(trim(checkpoint))

    read(chk_unit) have_disentangled      ! whether a disentanglement has been performed

    if (have_disentangled) then

       read(chk_unit) omega_invariant     ! omega invariant

       ! Allocate matrices if required
       if (.not.allocated(u_matrix_opt)) then
          allocate(u_matrix_opt(num_bands,num_wann,num_kpts),stat=ierr)
          if (ierr/=0) call io_error('Error allocating u_matrix_opt in param_read_chkpt')
       endif

       if (.not.allocated(lwindow)) then
          allocate(lwindow(num_bands,num_kpts),stat=ierr)
          if (ierr/=0) call io_error('Error allocating lwindow in param_read_chkpt')
       endif

       if (.not.allocated(ndimwin)) then
          allocate(ndimwin(num_kpts),stat=ierr)
          if (ierr/=0) call io_error('Error allocating ndimwin in param_read_chkpt')
       endif

       ! U matrix opt
       read(chk_unit,err=122) ((lwindow(i,nkp),i=1,num_bands),nkp=1,num_kpts)
       read(chk_unit,err=123) (ndimwin(nkp),nkp=1,num_kpts)
       read(chk_unit,err=124) (((u_matrix_opt(i,j,nkp),i=1,num_bands),j=1,num_wann),nkp=1,num_kpts)

    endif

    close(chk_unit)

    write(stdout,'(a/)') ' ... done'

    return

121 call io_error('Error opening '//trim(seedname)//'.chk in param_read_chkpt')
122 call io_error('Error reading lwindow from '//trim(seedname)//'.chk in param_read_chkpt')
123 call io_error('Error reading ndimwin from '//trim(seedname)//'.chk in param_read_chkpt')
124 call io_error('Error reading u_matrix_opt from '//trim(seedname)//'.chk in param_read_chkpt')

  end subroutine param_read_chkpt


  !=======================================!
  subroutine param_in_file
    !=======================================!
    ! Load the *.win file into a character  !
    ! array in_file, ignoring comments and  !
    ! blank lines and converting everything !
    ! to lowercase characters               !
    !=======================================!

    use io,        only : io_file_unit,io_error,seedname,stdout
    use utility,   only : utility_lowercase

    implicit none

    integer           :: in_unit,tot_num_lines,ierr,line_counter,loop,in1,in2
    character(len=maxlen) :: dummy

    in_unit=io_file_unit( )
    open (in_unit, file=trim(seedname)//'.win',form='formatted',status='old',err=101)

    num_lines=0;tot_num_lines=0
    do
       read(in_unit, '(a)', iostat = ierr, err= 200, end =210 ) dummy
       dummy=adjustl(dummy)
       tot_num_lines=tot_num_lines+1
       if( .not.dummy(1:1)=='!'  .and. .not. dummy(1:1)=='#' ) then
          if(len(trim(dummy)) > 0 ) num_lines=num_lines+1
       endif

    end do

101 call io_error('Error: Problem opening input file '//trim(seedname)//'.win')
200 call io_error('Error: Problem reading input file '//trim(seedname)//'.win')
210 continue
    rewind(in_unit)

    allocate(in_data(num_lines))

    line_counter=0
    do loop=1,tot_num_lines
       read(in_unit, '(a)', iostat = ierr, err= 200 ) dummy
       dummy=utility_lowercase(dummy)
       dummy=adjustl(dummy)
       if( dummy(1:1)=='!' .or.  dummy(1:1)=='#' ) cycle
       if(len(trim(dummy)) == 0 ) cycle
       line_counter=line_counter+1
       in1=index(dummy,'!')
       in2=index(dummy,'#')
       if(in1==0 .and. in2==0)  in_data(line_counter)=dummy
       if(in1==0 .and. in2>0 )  in_data(line_counter)=dummy(:in2-1)
       if(in2==0 .and. in1>0 )  in_data(line_counter)=dummy(:in1-1)
       if(in2>0 .and. in1>0 )   in_data(line_counter)=dummy(:min(in1,in2)-1)
    end do

    close(in_unit)

  end subroutine param_in_file


  !===========================================================================!
  subroutine param_get_keyword(keyword,found,c_value,l_value,i_value,r_value)
    !===========================================================================!
    !                                                                           !
    !             Finds the value of the required keyword.                      !
    !                                                                           !
    !===========================================================================!

    use io,        only : io_error

    implicit none

    character(*),      intent(in)  :: keyword
    logical          , intent(out) :: found
    character(*)     ,optional, intent(inout) :: c_value
    logical          ,optional, intent(inout) :: l_value
    integer          ,optional, intent(inout) :: i_value
    real(kind=dp)    ,optional, intent(inout) :: r_value

    integer           :: kl, in,loop,itmp
    character(len=maxlen) :: dummy

    kl=len_trim(keyword)

    found=.false.

    do loop=1,num_lines
       in=index(in_data(loop),trim(keyword))
       if (in==0 .or. in>1 ) cycle
       itmp=in+len(trim(keyword))
       if (in_data(loop)(itmp:itmp)/='=' &
            .and. in_data(loop)(itmp:itmp)/=':' &
            .and. in_data(loop)(itmp:itmp)/=' ') cycle
       if (found) then
          call io_error('Error: Found keyword '//trim(keyword)//' more than once in input file')
       endif
       found=.true.
       dummy=in_data(loop)(kl+1:)
       in_data(loop)(1:maxlen) = ' '
       dummy=adjustl(dummy)
       if( dummy(1:1)=='=' .or. dummy(1:1)==':') then
          dummy=dummy(2:)
          dummy=adjustl(dummy)
       end if
    end do

    if(found) then
       if( present(c_value) ) c_value=dummy
       if( present(l_value) ) then
          if (index(dummy,'t') > 0) then
             l_value=.true.
          elseif (index(dummy,'f') > 0) then
             l_value=.false.
          else
             call io_error('Error: Problem reading logical keyword '//trim(keyword))
          endif
       endif
       if( present(i_value) ) read(dummy,*,err=220) i_value
       if( present(r_value) ) read(dummy,*,err=220) r_value
    end if

    return

220 call io_error('Error: Problem reading keyword '//trim(keyword))


  end subroutine param_get_keyword


  !=========================================================================================!
  subroutine param_get_keyword_vector(keyword,found,length,c_value,l_value,i_value,r_value)
    !=========================================================================================!
    !                                                                                         !
    !                  Finds the values of the required keyword vector                        !
    !                                                                                         !
    !=========================================================================================!

    use io,        only : io_error

    implicit none

    character(*),      intent(in)  :: keyword
    logical          , intent(out) :: found
    integer,           intent(in)  :: length
    character(*)     ,optional, intent(inout) :: c_value(length)
    logical          ,optional, intent(inout) :: l_value(length)
    integer          ,optional, intent(inout) :: i_value(length)
    real(kind=dp)    ,optional, intent(inout) :: r_value(length)

    integer           :: kl, in,loop,i
    character(len=maxlen) :: dummy

    kl=len_trim(keyword)

    found=.false.



    do loop=1,num_lines
       in=index(in_data(loop),trim(keyword))
       if (in==0 .or. in>1 ) cycle
       if (found) then
          call io_error('Error: Found keyword '//trim(keyword)//' more than once in input file')
       endif
       found=.true.
       dummy=in_data(loop)(kl+1:)
       in_data(loop)(1:maxlen) = ' '
       dummy=adjustl(dummy)
       if( dummy(1:1)=='=' .or. dummy(1:1)==':') then
          dummy=dummy(2:)
          dummy=adjustl(dummy)
       end if
    end do

    if(found) then
       if( present(c_value) ) read(dummy,*,err=230) (c_value(i),i=1,length)
       if( present(l_value) ) then
          ! I don't think we need this. Maybe read into a dummy charater
          ! array and convert each element to logical
       endif
       if( present(i_value) ) read(dummy,*,err=230) (i_value(i),i=1,length)
       if( present(r_value) ) read(dummy,*,err=230) (r_value(i),i=1,length)
    end if



    return

230 call io_error('Error: Problem reading keyword '//trim(keyword)//' in param_get_keyword_vector')


  end subroutine param_get_keyword_vector



  !========================================================!
  subroutine param_get_vector_length(keyword,found,length)
    !======================================================!
    !                                                      !
    !        Returns the length of a keyword vector        !
    !                                                      !
    !======================================================!

    use io,        only : io_error

    implicit none

    character(*),      intent(in)  :: keyword
    logical          , intent(out) :: found
    integer,           intent(out)  :: length

    integer           :: kl, in,loop,i,pos
    character(len=maxlen) :: dummy

    kl=len_trim(keyword)

    found=.false.



    do loop=1,num_lines
       in=index(in_data(loop),trim(keyword))
       if (in==0 .or. in>1 ) cycle
       if (found) then
          call io_error('Error: Found keyword '//trim(keyword)//' more than once in input file')
       endif
       found=.true.
       dummy=in_data(loop)(kl+1:)
       dummy=adjustl(dummy)
       if( dummy(1:1)=='=' .or. dummy(1:1)==':') then
          dummy=dummy(2:)
          dummy=adjustl(dummy)
       end if
    end do

    length=0
    if(found) then
       if (len_trim(dummy)==0) call io_error('Error: keyword '//trim(keyword)//' is blank')
       length=1
       dummy=adjustl(dummy)
       do 
          pos=index(dummy,' ')
          dummy=dummy(pos+1:)
          dummy=adjustl(dummy)
          if(len_trim(dummy)>0) then
             length=length+1
          else
             exit
          endif

       end do

    end if



    return

230 call io_error('Error: Problem reading keyword '//trim(keyword)//' in param_get_keyword_vector')


  end subroutine param_get_vector_length


  !==============================================================================================!
  subroutine param_get_keyword_block(keyword,found,rows,columns,c_value,l_value,i_value,r_value)
    !==============================================================================================!
    !                                                                                              !
    !                           Finds the values of the required data block                        !
    !                                                                                              !
    !==============================================================================================!

    use constants, only : bohr
    use io,        only : io_error

    implicit none

    character(*),      intent(in)  :: keyword
    logical          , intent(out) :: found
    integer,           intent(in)  :: rows
    integer,           intent(in)  :: columns
    character(*)     ,optional, intent(inout) :: c_value(columns,rows)
    logical          ,optional, intent(inout) :: l_value(columns,rows)
    integer          ,optional, intent(inout) :: i_value(columns,rows)
    real(kind=dp)    ,optional, intent(inout) :: r_value(columns,rows)

    integer           :: in,ins,ine,loop,i,line_e,line_s,counter,blen
    logical           :: found_e,found_s,lconvert
    character(len=maxlen) :: dummy,end_st,start_st

    found_s=.false.
    found_e=.false.

    start_st='begin '//trim(keyword)
    end_st='end '//trim(keyword)


    do loop=1,num_lines
       ins=index(in_data(loop),trim(keyword))
       if (ins==0 ) cycle
       in=index(in_data(loop),'begin')
       if (in==0 .or. in>1) cycle
       line_s=loop
       if (found_s) then
          call io_error('Error: Found '//trim(start_st)//' more than once in input file')
       endif
       found_s=.true.
    end do

    if(.not. found_s) then
       found=.false.
       return
    end if


    do loop=1,num_lines
       ine=index(in_data(loop),trim(keyword))
       if (ine==0 ) cycle
       in=index(in_data(loop),'end')
       if (in==0 .or. in>1) cycle
       line_e=loop
       if (found_e) then
          call io_error('Error: Found '//trim(end_st)//' more than once in input file')
       endif
       found_e=.true.
    end do

    if(.not. found_e) then
       call io_error('Error: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file')
    end if

    if(line_e<=line_s) then
       call io_error('Error: '//trim(end_st)//' comes before '//trim(start_st)//' in input file')
    end if

    ! number of lines of data in block
    blen = line_e-line_s-1

    !    if( blen /= rows) then
    !       if ( index(trim(keyword),'unit_cell_cart').ne.0 ) then
    !          if ( blen /= rows+1 ) call io_error('Error: Wrong number of lines in block '//trim(keyword))
    !       else
    !          call io_error('Error: Wrong number of lines in block '//trim(keyword))          
    !       endif
    !    endif

    if ( (blen.ne.rows) .and. (blen.ne.rows+1) ) &
         call io_error('Error: Wrong number of lines in block '//trim(keyword))          

    if ( (blen.eq.rows+1) .and. (index(trim(keyword),'unit_cell_cart').eq.0) ) &
         call io_error('Error: Wrong number of lines in block '//trim(keyword))          


    found=.true.

    lconvert=.false.
    if (blen==rows+1) then
       dummy=in_data(line_s+1)
       if ( index(dummy,'ang').ne.0 ) then
          lconvert=.false.
       elseif ( index(dummy,'bohr').ne.0 ) then
          lconvert=.true.
       else
          call io_error('Error: Units in block '//trim(keyword)//' not recognised')
       endif
       in_data(line_s)(1:maxlen) = ' '
       line_s=line_s+1
    endif

!    r_value=1.0_dp
    counter=0
    do loop=line_s+1,line_e-1
       dummy=in_data(loop)
       counter=counter+1
       if( present(c_value) ) read(dummy,*,err=240) (c_value(i,counter),i=1,columns)
       if( present(l_value) ) then
          ! I don't think we need this. Maybe read into a dummy charater
          ! array and convert each element to logical
       endif
       if( present(i_value) ) read(dummy,*,err=240) (i_value(i,counter),i=1,columns)
       if( present(r_value) ) read(dummy,*,err=240) (r_value(i,counter),i=1,columns)
    end do

    if (lconvert) then
       if (present(r_value)) then
          r_value=r_value*bohr
       endif
    endif

    in_data(line_s:line_e)(1:maxlen) = ' '


    return

240 call io_error('Error: Problem reading block keyword '//trim(keyword))


  end subroutine param_get_keyword_block

  !=====================================================!
  subroutine param_get_block_length(keyword,found,rows,lunits)
    !=====================================================!
    !                                                     !
    !       Finds the length of the data block            !
    !                                                     !
    !=====================================================!

    use io,        only : io_error

    implicit none

    character(*),      intent(in)  :: keyword
    logical,           intent(out) :: found
    integer,           intent(out) :: rows
    logical, optional, intent(out) :: lunits

    integer           :: i,in,ins,ine,loop,line_e,line_s
    logical           :: found_e,found_s
    character(len=maxlen) :: end_st,start_st,dummy
    character(len=2)  :: atsym
    real(kind=dp)     :: atpos(3)

    found_s=.false.
    found_e=.false.

    start_st='begin '//trim(keyword)
    end_st='end '//trim(keyword)


    do loop=1,num_lines
       ins=index(in_data(loop),trim(keyword))
       if (ins==0 ) cycle
       in=index(in_data(loop),'begin')
       if (in==0 .or. in>1) cycle
       line_s=loop
       if (found_s) then
          call io_error('Error: Found '//trim(start_st)//' more than once in input file')
       endif
       found_s=.true.
    end do

    if(.not. found_s) then
       found=.false.
       return
    end if


    do loop=1,num_lines
       ine=index(in_data(loop),trim(keyword))
       if (ine==0 ) cycle
       in=index(in_data(loop),'end')
       if (in==0 .or. in>1) cycle
       line_e=loop
       if (found_e) then
          call io_error('Error: Found '//trim(end_st)//' more than once in input file')
       endif
       found_e=.true.
    end do

    if(.not. found_e) then
       call io_error('Error: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file')
    end if

    if(line_e<=line_s) then
       call io_error('Error: '//trim(end_st)//' comes before '//trim(start_st)//' in input file')
    end if

    rows=line_e-line_s-1

    found=.true.

    if (present(lunits)) then
       dummy=in_data(line_s+1)
       !       write(stdout,*) dummy
       !       write(stdout,*) trim(dummy)
       read(dummy,*,end=555) atsym, (atpos(i),i=1,3)
       lunits=.false.
    endif

    return

555 lunits=.true.

    return

  end subroutine param_get_block_length


  !===================================!
  subroutine param_get_atoms(lunits)
    !===================================!
    !                                   !
    !   Fills the atom data block       !
    !                                   !
    !===================================!

    use constants, only : bohr
    use utility,   only : utility_frac_to_cart,utility_cart_to_frac
    use io,        only : io_error
    implicit none

    logical, intent(in) :: lunits

    real(kind=dp)     :: atoms_pos_frac_tmp(3,num_atoms)
    real(kind=dp)     :: atoms_pos_cart_tmp(3,num_atoms)
    character(len=20) :: keyword
    integer           :: in,ins,ine,loop,i,line_e,line_s,counter
    integer           :: i_temp,loop2,max_sites
    logical           :: found_e,found_s,found,frac
    character(len=maxlen) :: dummy,end_st,start_st
    character(len=2)  :: ctemp(num_atoms)
    character(len=2)  :: atoms_label_tmp(num_atoms)
    logical           :: lconvert

    keyword="atoms_cart"
    frac=.false.
    call param_get_block_length("atoms_frac",found,i_temp)
    if (found) then
       keyword="atoms_frac"
       frac=.true.
    end if


    found_s=.false.
    found_e=.false.

    start_st='begin '//trim(keyword)
    end_st='end '//trim(keyword)


    do loop=1,num_lines
       ins=index(in_data(loop),trim(keyword))
       if (ins==0 ) cycle
       in=index(in_data(loop),'begin')
       if (in==0 .or. in>1) cycle
       line_s=loop
       if (found_s) then
          call io_error('Error: Found '//trim(start_st)//' more than once in input file')
       endif
       found_s=.true.
    end do



    do loop=1,num_lines
       ine=index(in_data(loop),trim(keyword))
       if (ine==0 ) cycle
       in=index(in_data(loop),'end')
       if (in==0 .or. in>1) cycle
       line_e=loop
       if (found_e) then
          call io_error('Error: Found '//trim(end_st)//' more than once in input file')
       endif
       found_e=.true.
    end do

    if(.not. found_e) then
       call io_error('Error: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file')
    end if

    if(line_e<=line_s) then
       call io_error('Error: '//trim(end_st)//' comes before '//trim(start_st)//' in input file')
    end if

    lconvert=.false.
    if (lunits) then
       dummy=in_data(line_s+1)
       if ( index(dummy,'ang').ne.0 ) then
          lconvert=.false.
       elseif ( index(dummy,'bohr').ne.0 ) then
          lconvert=.true.
       else
          call io_error('Error: Units in block atoms_cart not recognised in param_get_atoms')
       endif
       in_data(line_s)(1:maxlen) = ' '
       line_s=line_s+1
    endif

    counter=0
    do loop=line_s+1,line_e-1
       dummy=in_data(loop)
       counter=counter+1
       if(frac) then
          read(dummy,*,err=240) atoms_label_tmp(counter),(atoms_pos_frac_tmp(i,counter),i=1,3)
       else
          read(dummy,*,err=240) atoms_label_tmp(counter),(atoms_pos_cart_tmp(i,counter),i=1,3)
       end if
    end do

    if (lconvert) atoms_pos_cart_tmp = atoms_pos_cart_tmp*bohr

    in_data(line_s:line_e)(1:maxlen) = ' '



    if(frac) then
       do loop=1,num_atoms
          call utility_frac_to_cart (atoms_pos_frac_tmp(:,loop),atoms_pos_cart_tmp(:,loop),real_lattice)
       end do
    else
       do loop=1,num_atoms
          call utility_cart_to_frac (atoms_pos_cart_tmp(:,loop),atoms_pos_frac_tmp(:,loop),recip_lattice)
       end do
    end if


    ! Now we sort the data into the proper structures
    num_species=1
    ctemp(1)=atoms_label_tmp(1)
    do loop=2,num_atoms
       do loop2=1,loop-1
          if( trim(atoms_label_tmp(loop))==trim(atoms_label_tmp(loop2) )) exit
          if (loop2==loop-1) then 
             num_species=num_species+1
             ctemp(num_species)=atoms_label_tmp(loop)
          end if
       end do
    end do

    allocate(atoms_species_num(num_species))
    allocate(atoms_label(num_atoms))
    atoms_species_num(:)=0

    do loop=1,num_species
       atoms_label(loop)=ctemp(loop)
       do loop2=1,num_atoms
          if( trim(atoms_label(loop))==trim(atoms_label_tmp(loop2) )) then
             atoms_species_num(loop)=atoms_species_num(loop)+1
          end if
       end do
    end do

    max_sites=maxval(atoms_species_num)
    allocate(atoms_pos_frac(3,max_sites,num_species))
    allocate(atoms_pos_cart(3,max_sites,num_species))

    do loop=1,num_species
       counter=0
       do loop2=1,num_atoms
          if( trim(atoms_label(loop))==trim(atoms_label_tmp(loop2) )) then
             counter=counter+1
             atoms_pos_frac(:,counter,loop)=atoms_pos_frac_tmp(:,loop2)
             atoms_pos_cart(:,counter,loop)=atoms_pos_cart_tmp(:,loop2)
          end if
       end do
    end do


    return

240 call io_error('Error: Problem reading block keyword '//trim(keyword))

  end subroutine param_get_atoms


  !===================================!
  subroutine param_get_projections
    !===================================!
    !                                   !
    !  Fills the projection data block  !
    !                                   !
    !===================================!

    use utility,   only : utility_frac_to_cart,utility_cart_to_frac,&
         utility_string_to_coord,utility_strip
    use io,        only : io_error
    implicit none


    real(kind=dp)     :: pos_frac(3)
    real(kind=dp)     :: pos_cart(3)
    character(len=20) :: keyword
    integer           :: in,ins,ine,loop,line_e,line_s,counter
    integer           :: sites,species,line,pos1,pos2,pos3,m_tmp,l_tmp,mstate
    integer           :: loop_l,loop_m,loop_sites
    logical           :: found_e,found_s
    character(len=maxlen) :: dummy,end_st,start_st
    character(len=maxlen) :: ctemp,ctemp2,ctemp3,ctemp4,ctemp5,m_string
    !
    integer, parameter :: min_l=-5
    integer, parameter :: max_l=3
    integer, parameter :: min_m=1
    integer, parameter :: max_m=7
    integer            :: ang_states(min_m:max_m,min_l:max_l)
    ! default values for the optional part of the projection definitions
    real(kind=dp), parameter :: proj_z_def(3)=(/0.d0,0.d0,1.d0/)
    real(kind=dp), parameter :: proj_x_def(3)=(/1.d0,0.d0,0.d0/)
    real(kind=dp), parameter :: proj_zona_def=1.d0
    real(kind=dp), parameter :: proj_box_def=1.d0
    integer, parameter       :: proj_radial_def=1.d0
    !
    real(kind=dp) :: proj_z_tmp(3)
    real(kind=dp) :: proj_x_tmp(3)
    real(kind=dp) :: proj_zona_tmp
    real(kind=dp) :: proj_box_tmp
    integer       :: proj_radial_tmp

    keyword="projections"

    found_s=.false.
    found_e=.false.

    start_st='begin '//trim(keyword)
    end_st='end '//trim(keyword)

    allocate( proj_site(3,num_wann) )
    allocate( proj_l(num_wann)  )
    allocate( proj_m(num_wann)  )
    allocate( proj_z(3,num_wann) )
    allocate( proj_x(3,num_wann) )
    allocate( proj_radial(num_wann) )  
    allocate( proj_zona(num_wann) )
    allocate( proj_box(num_wann) )




    do loop=1,num_lines
       ins=index(in_data(loop),trim(keyword))
       if (ins==0 ) cycle
       in=index(in_data(loop),'begin')
       if (in==0 .or. in>1) cycle
       line_s=loop
       if (found_s) then
          call io_error('Error: Found '//trim(start_st)//' more than once in input file')
       endif
       found_s=.true.
    end do



    do loop=1,num_lines
       ine=index(in_data(loop),trim(keyword))
       if (ine==0 ) cycle
       in=index(in_data(loop),'end')
       if (in==0 .or. in>1) cycle
       line_e=loop
       if (found_e) then
          call io_error('param_get_projections: Found '//trim(end_st)//' more than once in input file')
       endif
       found_e=.true.
    end do

    if(.not. found_e) then
       call io_error('param_get_projections: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file')
    end if

    if(line_e<=line_s) then
       call io_error('param_get_projections: '//trim(end_st)//' comes before '//trim(start_st)//' in input file')
    end if

    counter=0
    do line=line_s+1,line_e-1
       ang_states=0
       !Assume the default values
       proj_z_tmp      = proj_z_def  
       proj_x_tmp      = proj_x_def  
       proj_zona_tmp   = proj_zona_def  
       proj_box_tmp    = proj_box_def   
       proj_radial_tmp =  proj_radial_def
       ! Strip input line of all spaces
       dummy=utility_strip(in_data(line))
       dummy=adjustl(dummy)
       pos1=index(dummy,':')
       if(pos1==0) call io_error('param_read_projection: malformed projection &
            &definition: '//trim(dummy))
       sites=0
       ctemp=dummy(:pos1-1)
       ! Read the atomic site
       if(index(ctemp,'c=')>0) then
          sites=-1
          ctemp=ctemp(3:)
          call utility_string_to_coord(ctemp,pos_cart)
          call utility_cart_to_frac (pos_cart(:),pos_frac(:),recip_lattice)          
       elseif(index(ctemp,'f=')>0) then
          sites=-1
          ctemp=ctemp(3:)
          call utility_string_to_coord(ctemp,pos_frac)
       else
          do loop=1,num_species
             if(trim(ctemp)==atoms_label(loop)) then
                species=loop
                sites=atoms_species_num(loop)
                exit
             end if
             if(loop==num_species) call io_error('param_get_projection: Atom site not recognised '//trim(ctemp))
          end do
       end if
       !Now we know the sites for this line. Get the angular momentum states
       dummy=dummy(pos1+1:)
       pos1=index(dummy,':')
       if (pos1>0) then
          ctemp=dummy(:pos1-1)
       else
          ctemp=dummy
       end if
       ctemp2=ctemp
       do
          pos2=index(ctemp2,';')
          if (pos2==0) then
             ctemp3=ctemp2
          else
             ctemp3=ctemp2(:pos2-1)
          endif
          if (index(ctemp3,'l=')==1) then
             mstate=index(ctemp3,',')
             if(mstate>0) then
                read(ctemp3(3:mstate-1),*,err=101) l_tmp
             else
                read(ctemp3(3:),*,err=101) l_tmp
             end if
             if (l_tmp<-5 .or. l_tmp>3) call io_error('param_get_projection: Incorrect l state requested')
             if (mstate==0) then
                if(l_tmp>=0) then
                   do loop_m=1,2*l_tmp+1
                      ang_states(loop_m,l_tmp)=1
                   end do
                elseif (l_tmp==-1) then !sp
                   ang_states(1:2,l_tmp)=1
                elseif (l_tmp==-2) then !sp2
                   ang_states(1:3,l_tmp)=1
                elseif (l_tmp==-3) then !sp3
                   ang_states(1:4,l_tmp)=1
                elseif (l_tmp==-4) then !sp3d
                   ang_states(1:5,l_tmp)=1
                elseif (l_tmp==-5) then !sp3d2
                   ang_states(1:6,l_tmp)=1
                endif
             else
                if (index(ctemp3,'mr=')/=mstate+1) &
                     call io_error('param_get_projection: Problem reading m state')
                ctemp4=ctemp3(mstate+4:)
                do
                   pos3=index(ctemp4,',')
                   if (pos3==0) then
                      ctemp5=ctemp4
                   else
                      ctemp5=ctemp4(:pos3-1)
                   endif
                   read(ctemp5(1:),*,err=102) m_tmp
                   if(l_tmp>=0) then
                      if ((m_tmp>2*l_tmp+1) .or. (m_tmp<=0)) call io_error('param_get_projection: m is > l !')
                   elseif (l_tmp==-1 .and. (m_tmp>2 .or. m_tmp<=0)) then
                      call io_error('param_get_projection: m has incorrect value (1)')
                   elseif (l_tmp==-2 .and. (m_tmp>3 .or. m_tmp<=0)) then
                      call io_error('param_get_projection: m has incorrect value (2)')
                   elseif (l_tmp==-3 .and. (m_tmp>4 .or. m_tmp<=0)) then
                      call io_error('param_get_projection: m has incorrect value (3)')
                   elseif (l_tmp==-4 .and. (m_tmp>5 .or. m_tmp<=0)) then
                      call io_error('param_get_projection: m has incorrect value (4)')
                   elseif (l_tmp==-5 .and. (m_tmp>6 .or. m_tmp<=0)) then
                      call io_error('param_get_projection: m has incorrect value (5)')
                   endif
                   ang_states(m_tmp,l_tmp)=1
                   if (pos3==0) exit
                   ctemp4=ctemp4(pos3+1:)
                enddo
             end if
          else
             do
                pos3=index(ctemp3,',')
                if (pos3==0) then
                   ctemp4=ctemp3
                else
                   ctemp4=ctemp3(:pos3-1)
                endif
                read(ctemp4(1:),*,err=106) m_string
                select case (trim(adjustl(m_string)))
                case ('s')
                   ang_states(1,0)=1
                case ('pz')
                   ang_states(1,1)=1
                case ('px')
                   ang_states(2,1)=1
                case ('py')
                   ang_states(3,1)=1
                case ('dz2')
                   ang_states(1,2)=1
                case ('dxz')
                   ang_states(2,2)=1
                case ('dyz')
                   ang_states(3,2)=1
                case ('dx2-y2')
                   ang_states(4,2)=1
                case ('dxy')
                   ang_states(5,2)=1
                case ('fz2')
                   ang_states(1,3)=1
                case ('fxz2')
                   ang_states(2,3)=1
                case ('fyz2')
                   ang_states(3,3)=1
                case ('fxyz')
                   ang_states(4,3)=1
                case ('fz(x2-y2)')
                   ang_states(5,3)=1
                case ('fx(x2-3y2)')
                   ang_states(6,3)=1
                case ('fy(3x2-y2)')
                   ang_states(7,3)=1
                case ('sp')
                   ang_states(1:2,-1)=1
                case ('sp-1')
                   ang_states(1,-1)=1
                case ('sp-2')
                   ang_states(2,-1)=1
                case ('sp2')
                   ang_states(1:3,-2)=1
                case ('sp2-1')
                   ang_states(1,-2)=1
                case ('sp2-2')
                   ang_states(2,-2)=1
                case ('sp2-3')
                   ang_states(3,-2)=1
                case ('sp3')
                   ang_states(1:4,-3)=1
                case ('sp3-1')
                   ang_states(1,-3)=1
                case ('sp3-2')
                   ang_states(2,-3)=1
                case ('sp3-3')
                   ang_states(3,-3)=1
                case ('sp3-4')
                   ang_states(4,-3)=1
                case ('sp3d')
                   ang_states(1:5,-4)=1
                case ('sp3d-1')
                   ang_states(1,-4)=1
                case ('sp3d-2')
                   ang_states(2,-4)=1
                case ('sp3d-3')
                   ang_states(3,-4)=1
                case ('sp3d-4')
                   ang_states(4,-4)=1
                case ('sp3d-5')
                   ang_states(5,-4)=1
                case ('sp3d2')
                   ang_states(1:6,-5)=1
                case ('sp3d2-1')
                   ang_states(1,-5)=1
                case ('sp3d2-2')
                   ang_states(2,-5)=1
                case ('sp3d2-3')
                   ang_states(3,-5)=1
                case ('sp3d2-4')
                   ang_states(4,-5)=1
                case ('sp3d2-5')
                   ang_states(5,-5)=1
                case ('sp3d2-6')
                   ang_states(6,-5)=1
                case default
                   call io_error('param_get_projection: Problem reading l state '//trim(ctemp3))
                end select
                if (pos3==0) exit
                ctemp3=ctemp3(pos3+1:)
             enddo
          endif
          if(pos2==0) exit
          ctemp2=ctemp2(pos2+1:)
       enddo
       ! check for non-default values
       if(pos1>0) then
          dummy=dummy(pos1+1:)
          ! z axis
          pos1=index(dummy,'z=')
          if(pos1>0) then
             ctemp=(dummy(pos1+2:))
             pos2=index(ctemp,':')
             if(pos2>0) ctemp=ctemp(:pos2-1)
             call utility_string_to_coord(ctemp,proj_z_tmp)
          endif
          ! x axis
          pos1=index(dummy,'x=')
          if(pos1>0) then
             ctemp=(dummy(pos1+2:))
             pos2=index(ctemp,':')
             if(pos2>0) ctemp=ctemp(:pos2-1)
             call utility_string_to_coord(ctemp,proj_x_tmp)
          endif
          ! projection box
          pos1=index(dummy,'b=')
          if(pos1>0) then
             ctemp=(dummy(pos1+2:))
             pos2=index(ctemp,':')
             if(pos2>0) ctemp=ctemp(:pos2-1)
             read(ctemp,*,err=103) proj_box_tmp
          endif
          ! diffusivity of orbital
          pos1=index(dummy,'zona=')
          if(pos1>0) then
             ctemp=(dummy(pos1+2:))
             pos2=index(ctemp,':')
             if(pos2>0) ctemp=ctemp(:pos2-1)
             read(ctemp,*,err=104) proj_zona_tmp
          endif
          ! nodes for the radial part
          pos1=index(dummy,'r=')
          if(pos1>0) then
             ctemp=(dummy(pos1+2:))
             pos2=index(ctemp,':')
             if(pos2>0) ctemp=ctemp(:pos2-1)
             read(ctemp,*,err=105) proj_radial_tmp
          endif
       end if
       if(sites==-1) then
          if(counter+sum(ang_states) > num_wann) call io_error('param_get_projection: &
               &too many projections defined')
       else
          if(counter+sites*sum(ang_states) > num_wann) call io_error('param_get_projection:&
               & too many projections defined')
       end if
       !
       if(sites==-1) then
          do loop_l=min_l,max_l
             do loop_m=min_m,max_m
                if(ang_states(loop_m,loop_l)==1) then
                   counter=counter+1
                   proj_site(:,counter) = pos_frac
                   proj_l(counter)      = loop_l
                   proj_m(counter)      = loop_m
                   proj_z(:,counter)    = proj_z_tmp
                   proj_x(:,counter)    = proj_x_tmp
                   proj_radial(counter) = proj_radial_tmp
                   proj_zona(counter)   = proj_zona_tmp
                   proj_box(counter)    = proj_box_tmp
                end if
             end do
          end do
       else
          do loop_sites=1,sites
             do loop_l=min_l,max_l
                do loop_m=min_m,max_m
                   if(ang_states(loop_m,loop_l)==1) then
                      counter=counter+1
                      proj_site(:,counter) = atoms_pos_frac(:,loop_sites,species)
                      proj_l(counter)      = loop_l
                      proj_m(counter)      = loop_m
                      proj_z(:,counter)    = proj_z_tmp
                      proj_x(:,counter)    = proj_x_tmp
                      proj_radial(counter) = proj_radial_tmp
                      proj_zona(counter)   = proj_zona_tmp
                      proj_box(counter)    = proj_box_tmp
                   end if
                end do
             end do
          end do
       end if

    end do !end loop over projection block
    if (counter.ne.num_wann) call io_error('param_get_projections:&
         & Fewer projections defined than the number of Wannier functions requested')


    in_data(line_s:line_e)(1:maxlen) = ' '

    ! Check
    do loop=1,num_wann
       if ( sum(proj_z(:,loop)*proj_x(:,loop)).gt.1.0e-6_dp ) then
          write(stdout,*) ' Projection:',loop
          call io_error(' Error in projections: z and x axes are not orthogonal')
       endif
    enddo

    return


101 call io_error('param_get_projection: Problem reading l state into integer '//trim(ctemp3))
102 call io_error('param_get_projection: Problem reading m state into integer '//trim(ctemp3))
103 call io_error('param_get_projection: Problem reading box size into real '//trim(ctemp))
104 call io_error('param_get_projection: Problem reading zona into real '//trim(ctemp))
105 call io_error('param_get_projection: Problem reading radial state into integer '//trim(ctemp))
106 call io_error('param_get_projection: Problem reading m state into string '//trim(ctemp3))



  end subroutine param_get_projections

  !===================================!
  subroutine param_get_keyword_kpath
    !===================================!
    !                                   !
    !  Fills the kpath data block       !
    !                                   !
    !===================================!
    use io,        only : io_error

    implicit none

    character(len=20) :: keyword
    integer           :: in,ins,ine,loop,i,line_e,line_s,counter
    logical           :: found_e,found_s
    character(len=maxlen) :: dummy,end_st,start_st

    keyword="kpoint_path"

    found_s=.false.
    found_e=.false.

    start_st='begin '//trim(keyword)
    end_st='end '//trim(keyword)


    do loop=1,num_lines
       ins=index(in_data(loop),trim(keyword))
       if (ins==0 ) cycle
       in=index(in_data(loop),'begin')
       if (in==0 .or. in>1) cycle
       line_s=loop
       if (found_s) then
          call io_error('Error: Found '//trim(start_st)//' more than once in input file')
       endif
       found_s=.true.
    end do



    do loop=1,num_lines
       ine=index(in_data(loop),trim(keyword))
       if (ine==0 ) cycle
       in=index(in_data(loop),'end')
       if (in==0 .or. in>1) cycle
       line_e=loop
       if (found_e) then
          call io_error('Error: Found '//trim(end_st)//' more than once in input file')
       endif
       found_e=.true.
    end do

    if(.not. found_e) then
       call io_error('Error: Found '//trim(start_st)//' but no '//trim(end_st)//' in input file')
    end if

    if(line_e<=line_s) then
       call io_error('Error: '//trim(end_st)//' comes before '//trim(start_st)//' in input file')
    end if

    counter=0
    do loop=line_s+1,line_e-1

       counter=counter+2
       dummy=in_data(loop)
       read(dummy,*,err=240) bands_label(counter-1),(bands_spec_points(i,counter-1),i=1,3)&
            ,bands_label(counter),(bands_spec_points(i,counter),i=1,3)
    end do


    in_data(line_s:line_e)(1:maxlen) = ' '

    return


240 call io_error('param_get_keyword_kpath: Problem reading kpath '//trim(dummy))

  end subroutine param_get_keyword_kpath

end module parameters
