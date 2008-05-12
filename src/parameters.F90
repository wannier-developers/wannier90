!-*- mode: F90; mode: font-lock; column-number-mode: true -*-!
!                                                            !
! Copyright (C) 2007 Jonathan Yates, Arash Mostofi,          !
!  Young-Su Lee, Nicola Marzari, Ivo Souza, David Vanderbilt !
!                                                            !
! This file is distributed under the terms of the GNU        !
! General Public License. See the file `LICENSE' in          !
! the root directory of the present distribution, or         !
! http://www.gnu.org/copyleft/gpl.txt .                      !
!                                                            !
!------------------------------------------------------------!

module w90_parameters

  use w90_constants, only : dp
  use w90_io,        only : stdout,maxlen

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
  logical,           public, save :: gamma_only  
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
  character(len=20), public, save :: wannier_plot_mode
  logical,           public, save :: bands_plot
  integer,           public, save :: bands_num_points
  character(len=20), public, save :: bands_plot_format
  character(len=20), public, save :: bands_plot_mode
  integer, allocatable, public, save :: bands_plot_project(:)
  integer,           public, save :: bands_plot_dim         
  logical,           public, save :: hr_plot
  real(kind=dp),     public, save :: hr_cutoff
  real(kind=dp),     public, save :: dist_cutoff
  character(len=20), public, save :: dist_cutoff_mode
  character(len=20), public, save :: one_dim_axis
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
  logical,           public, save :: transport
  character(len=20), public, save :: transport_mode
  real(kind=dp),     public, save :: tran_win_min
  real(kind=dp),     public, save :: tran_win_max
  real(kind=dp),     public, save :: tran_energy_step
  integer,           public, save :: tran_num_bb
  integer,           public, save :: tran_num_ll
  integer,           public, save :: tran_num_rr
  integer,           public, save :: tran_num_cc
  integer,           public, save :: tran_num_lc
  integer,           public, save :: tran_num_cr
  integer,           public, save :: tran_num_bandc
  logical,           public, save :: tran_write_ht
  logical,           public, save :: tran_read_ht 
  logical,           public, save :: tran_use_same_lead
  real(kind=dp),     public, save :: translation_centre_frac(3)
  integer,           public, save :: num_shells    !no longer an input keyword
  integer, allocatable, public,save :: shell_list(:)
  real(kind=dp), allocatable,    public, save :: kpt_latt(:,:) !kpoints in lattice vecs
  real(kind=dp),     public, save :: real_lattice(3,3)
  logical,           public, save :: postproc_setup
  logical,           public, save :: cp_pp ! Car-Parinello post-proc flag
  logical,           public, save :: calc_only_A
  logical,           public, save :: use_bloch_phases
  character(len=20), public, save :: restart
  logical,           public, save :: write_r2mn
  logical,           public, save :: guiding_centres
  integer,           public, save :: num_guide_cycles
  integer,           public, save :: num_no_guide_iter
  real(kind=dp),     public, save :: fixed_step
  real(kind=dp),     public, save :: trial_step
  logical,           public, save :: write_proj
  integer,           public, save :: timing_level
  logical,           public, save :: spinors   !are our WF spinors?
  logical,           public, save :: translate_home_cell
  logical,           public, save :: write_xyz
  real(kind=dp),     public, save :: conv_noise_amp
  integer,           public, save :: conv_noise_num
  real(kind=dp),     public, save :: wannier_plot_radius
  integer,           public, save :: search_shells   !for kmesh
  real(kind=dp),     public, save :: kmesh_tol

  ! Restarts
  real(kind=dp),     public, save :: omega_invariant
  character(len=20), public, save :: checkpoint
  logical,           public, save :: have_disentangled

  ! Atom sites
  real(kind=dp), allocatable,     public, save :: atoms_pos_frac(:,:,:)
  real(kind=dp), allocatable,     public, save :: atoms_pos_cart(:,:,:)
  integer, allocatable,           public, save :: atoms_species_num(:)  
  character(len=maxlen), allocatable,  public, save :: atoms_label(:)
  character(len=2), allocatable,  public, save :: atoms_symbol(:)
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
  integer,                        public, save :: num_proj

  !parameters dervied from input
  integer,           public, save :: num_kpts
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
  integer,           public, save :: num_bands_project
  integer,           public, save :: num_exclude_bands
  logical,           public, save :: lfixstep

  ! kmesh parameters (set in kmesh)

  integer,       public, save              :: nnh           ! the number of b-directions (bka)
  integer,       public, save              :: nntot         ! total number of neighbours for each k-point
  integer,       public, save, allocatable :: nnlist(:,:)   ! list of neighbours for each k-point
  integer,       public, save, allocatable :: neigh(:,:)    
  integer,       public, save, allocatable :: nncell(:,:,:) ! gives BZ of each neighbour of each k-point
  real(kind=dp), public, save              :: wbtot
  real(kind=dp), public, save, allocatable :: wb(:)         ! weights associated with neighbours of each k-point
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

!!$![ysl-b]
!!$  ! ph_g = phase factor of Bloch functions at Gamma
!!$  !  assuming that Bloch functions at Gamma are real except this phase factor
!!$  complex(kind=dp), allocatable, save, public :: ph_g(:)
!!$![ysl-e]

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
  integer, parameter, public :: num_nnmax=12

  ! Are we running as a libarary
  logical, save, public :: library

  ! Wannier centres and spreads
  real(kind=dp), public, save, allocatable :: wannier_centres(:,:)
  real(kind=dp), public, save, allocatable :: wannier_spreads(:)
  real(kind=dp), public, save :: omega_total
  real(kind=dp), public, save :: omega_tilde
  ! [ omega_invariant is declared above ]

  ! For Hamiltonian matrix in WF representation
  logical,          public, save              :: automatic_translation
  integer,          public, save              :: one_dim_dir

  ! Private data
  integer                            :: num_lines
  character(len=maxlen), allocatable :: in_data(:)
  logical                            :: ltmp

  public :: param_read
  public :: param_write
  public :: param_dealloc
  public :: param_write_header
  public :: param_write_chkpt
  public :: param_read_chkpt
  public :: param_lib_set_atoms
  public :: param_memory_estimate

contains



  !==================================================================!
  subroutine param_read ( )
  !==================================================================!
  !                                                                  !
  ! Read parameters and calculate derived values                     !
  !                                                                  !
  !===================================================================  
    use w90_constants, only : bohr, eps6
    use w90_utility,   only : utility_recip_lattice,utility_compute_metric
    use w90_io,        only : io_error,io_file_unit,seedname,post_proc_flag
    implicit none

    !local variables
    real(kind=dp)  :: real_lattice_tmp(3,3)
    integer :: nkp,i,j,n,k,itmp,i_temp,i_temp2,eig_unit,loop,ierr,iv_temp(3)
    logical :: found,found2,eig_found,lunits,chk_found
    character(len=6) :: spin_str
    real(kind=dp) :: cosa(3),rv_temp(3)

    call param_in_file

    !%%%%%%%%%%%%%%%%
    ! Transport 
    !%%%%%%%%%%%%%%%%

    transport         = .false.
    call param_get_keyword('transport',found,l_value=transport)    

    tran_read_ht           = .false. 
    call param_get_keyword('tran_read_ht',found,l_value=tran_read_ht)

    if (transport .and. tran_read_ht) restart = ' '

    !%%%%%%%%%%%%%%%%
    !System variables
    !%%%%%%%%%%%%%%%%

    timing_level    =  1             ! Verbosity of timing output info
    call param_get_keyword('timing_level',found,i_value=timing_level)

    iprint          =  1             ! Verbosity
    call param_get_keyword('iprint',found,i_value=iprint)

    if (transport .and. tran_read_ht) goto 301 

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
    if(.not. found) call io_error('Error: You must specify num_wann')
    if(num_wann<=0) then
       call io_error('Error: num_wann must be greater than zero')
    endif

    num_exclude_bands=0
    call param_get_range_vector('exclude_bands',found,num_exclude_bands,lcount=.true.)
    if(found) then
       if(num_exclude_bands<1) call io_error('Error: problem reading exclude_bands')
       allocate(exclude_bands(num_exclude_bands),stat=ierr)
       if (ierr/=0) call io_error('Error allocating exclude_bands in param_read')
       call param_get_range_vector('exclude_bands',found,num_exclude_bands,.false.,exclude_bands)
       if (any(exclude_bands<1)  ) &
            call io_error('Error: exclude_bands must contain positive numbers')
    end if

!    num_bands       =   -1   
    call param_get_keyword('num_bands',found,i_value=i_temp)
    if(found.and.library) write(stdout,'(/a)') ' Ignoring <num_bands> in input file'
    if (.not. library) then
       if(found) num_bands=i_temp
       if(.not.found) num_bands=num_wann
       if(found .and. num_bands<num_wann) then
          call io_error('Error: num_bands must be greater than or equal to num_wann')
       endif
    end if

    num_dump_cycles =   100          ! frequency to write backups at
    call param_get_keyword('num_dump_cycles',found,i_value=num_dump_cycles)
    if (num_dump_cycles<0) call io_error('Error: num_dump_cycles must be positive')       

    num_print_cycles =   1          ! frequency to write at
    call param_get_keyword('num_print_cycles',found,i_value=num_print_cycles)
    if (num_print_cycles<0) call io_error('Error: num_print_cycles must be positive')       

    devel_flag      =   ' '          !       
    call param_get_keyword('devel_flag',found,c_value=devel_flag)

!    mp_grid=-99
    call param_get_keyword_vector('mp_grid',found,3,i_value=iv_temp)
    if(found.and.library) write(stdout,'(a)') ' Ignoring <mp_grid> in input file'
    if(.not.library) then
       if(found) mp_grid=iv_temp
       if (.not. found) then
          call io_error('Error: You must specify dimensions of the Monkhorst-Pack grid by setting mp_grid')
       elseif (any(mp_grid<1)) then
          call io_error('Error: mp_grid must be greater than zero')
       end if
       num_kpts= mp_grid(1)*mp_grid(2)*mp_grid(3)
    end if

![ysl-b]
    ltmp=.false.
    call param_get_keyword('gamma_only',found,l_value=ltmp)
    if (.not.library) then
       gamma_only=ltmp
       if ( gamma_only .and. (num_kpts.ne.1) ) &
            call io_error('Error: gamma_only is true, but num_kpts > 1')
    else
       if (found) write(stdout,'(a)') ' Ignoring <gamma_only> in input file'
    endif
![ysl-e]

    automatic_mp_grid = .false.
    call param_get_keyword('automatic_mp_grid',found,l_value=automatic_mp_grid)

    postproc_setup = .false.            ! set to true to write .nnkp file and exit
    call param_get_keyword('postproc_setup',found,l_value=postproc_setup)
    ! We allow this keyword to be overriden by a command line arg -pp
    if(post_proc_flag) postproc_setup=.true.      


    cp_pp = .false.                  ! set to true if doing CP post-processing
    call param_get_keyword('cp_pp',found,l_value=cp_pp)

    calc_only_A = .false.
    call param_get_keyword('calc_only_A',found,l_value=calc_only_A)

    restart = ' '
    call param_get_keyword('restart',found,c_value=restart)
    if (found) then
       if ( (restart.ne.'default').and.(restart.ne.'wannierise') &
          .and.(restart.ne.'plot').and.(restart.ne.'transport') ) then 
          call io_error('Error in input file: value of restart not recognised')
       else
          inquire(file=trim(seedname)//'.chk',exist=chk_found)
          if (.not. chk_found) &
               call io_error('Error: restart requested but '//trim(seedname)//'.chk file not found')
       endif
    endif
    !post processing takes priority (must warn user of this)
    if (postproc_setup) restart = ' '

    write_r2mn = .false.
    call param_get_keyword('write_r2mn',found,l_value=write_r2mn)

    write_proj = .false.
    call param_get_keyword('write_proj',found,l_value=write_proj)


    ltmp=.false.  ! by default our WF are not spinors
    call param_get_keyword('spinors',found,l_value=ltmp)
    if (.not.library) then
       spinors=ltmp
    else
       if (found) write(stdout,'(a)') ' Ignoring <spinors> in input file'
    endif
    if(spinors .and. (2*(num_wann/2))/=num_wann) &
       call io_error('Error: For spinor WF num_wann ust be even')

    translate_home_cell = .false.
    call param_get_keyword('translate_home_cell',found,l_value=translate_home_cell)

    write_xyz = .false.
    call param_get_keyword('write_xyz',found,l_value=write_xyz)

    !%%%%%%%%%%%
    ! Wannierise
    !%%%%%%%%%%%

    num_iter          = 100    
    call param_get_keyword('num_iter',found,i_value=num_iter)
    if (num_iter<0) call io_error('Error: num_iter must be positive')       

    num_cg_steps      =   5
    call param_get_keyword('num_cg_steps',found,i_value=num_cg_steps)
    if (num_cg_steps<0) call io_error('Error: num_cg_steps must be positive')       

    conv_tol=1.0e-10_dp
    call param_get_keyword('conv_tol',found,r_value=conv_tol)
    if (conv_tol<0.0_dp) call io_error('Error: conv_tol must be positive')

    conv_noise_amp=-1.0_dp
    call param_get_keyword('conv_noise_amp',found,r_value=conv_noise_amp)

    conv_window=-1
    if (conv_noise_amp>0.0_dp) conv_window=5
    call param_get_keyword('conv_window',found,i_value=conv_window)

    conv_noise_num=3
    call param_get_keyword('conv_noise_num',found,i_value=conv_noise_num)
    if (conv_noise_num<0) call io_error('Error: conv_noise_num must be positive')

    guiding_centres=.false.
    call param_get_keyword('guiding_centres',found,l_value=guiding_centres)

    num_guide_cycles=1
    call param_get_keyword('num_guide_cycles',found,i_value=num_guide_cycles)
    if (num_guide_cycles<0) call io_error('Error: num_guide_cycles must be >= 0')

    num_no_guide_iter=0
    call param_get_keyword('num_no_guide_iter',found,i_value=num_no_guide_iter)
    if (num_no_guide_iter<0) call io_error('Error: num_no_guide_iter must be >= 0')

    fixed_step=-999.0_dp ; lfixstep=.false.
    call param_get_keyword('fixed_step',found,r_value=fixed_step)
    if ( found.and.(fixed_step<0.0_dp) ) call io_error('Error: fixed_step must be > 0')
    if ( fixed_step > 0.0_dp ) lfixstep=.true.

    trial_step=2.0_dp
    call param_get_keyword('trial_step',found,r_value=trial_step)
    if ( found.and.lfixstep ) call io_error('Error: cannot specify both fixed_step and trial_step')

    !%%%%%%%%%
    ! Plotting
    !%%%%%%%%%

    wannier_plot              = .false.
    call param_get_keyword('wannier_plot',found,l_value=wannier_plot)

    wannier_plot_supercell    = 2
    call param_get_keyword('wannier_plot_supercell',found,i_value=wannier_plot_supercell)

    wannier_plot_format       = 'xcrysden'
    call param_get_keyword('wannier_plot_format',found,c_value=wannier_plot_format)

    wannier_plot_mode       = 'crystal'
    call param_get_keyword('wannier_plot_mode',found,c_value=wannier_plot_mode)
    
    call param_get_range_vector('wannier_plot_list',found,num_wannier_plot,lcount=.true.)
    if(found) then
       if(num_wannier_plot<1) call io_error('Error: problem reading wannier_plot_list')
       allocate(wannier_plot_list(num_wannier_plot),stat=ierr)
       if (ierr/=0) call io_error('Error allocating wannier_plot_list in param_read')
       call param_get_range_vector('wannier_plot_list',found,num_wannier_plot,.false.,wannier_plot_list)
       if (any(wannier_plot_list<1) .or. any(wannier_plot_list>num_wann) ) &
            call io_error('Error: wannier_plot_list asks for a non-valid wannier function to be plotted')
    else
       ! we plot all wannier functions
       num_wannier_plot=num_wann
       allocate(wannier_plot_list(num_wannier_plot),stat=ierr)
       if (ierr/=0) call io_error('Error allocating wannier_plot_list in param_read')
       do loop=1,num_wann
          wannier_plot_list(loop)=loop
       end do
    end if

    wannier_plot_radius = 3.5_dp
    call param_get_keyword('wannier_plot_radius',found,r_value=wannier_plot_radius)

    ! checks
    if (wannier_plot) then
       if ( (index(wannier_plot_format,'xcrys').eq.0) .and. (index(wannier_plot_format,'cub').eq.0) ) &
            call io_error('Error: wannier_plot_format not recognised')
       if ( (index(wannier_plot_mode,'crys').eq.0) .and. (index(wannier_plot_mode,'mol').eq.0) ) &
            call io_error('Error: wannier_plot_mode not recognised')
       if ( wannier_plot_radius < 0.0_dp ) call io_error('Error: wannier_plot_radius must be positive')
       if ( wannier_plot_supercell < 0 ) call io_error('Error: wannier_plot_supercell must be positive')       
    endif

    bands_plot                = .false.
    call param_get_keyword('bands_plot',found,l_value=bands_plot)

    bands_num_points          = 100
    call param_get_keyword('bands_num_points',found,i_value=bands_num_points)

    bands_plot_format         = 'gnuplot'
    call param_get_keyword('bands_plot_format',found,c_value=bands_plot_format)

    bands_plot_mode             = 's-k'
    call param_get_keyword('bands_plot_mode',found,c_value=bands_plot_mode)

    bands_plot_dim             = 3
    call param_get_keyword('bands_plot_dim',found,i_value=bands_plot_dim)

    num_bands_project=0
    call param_get_range_vector('bands_plot_project',found,num_bands_project,lcount=.true.)
    if(found) then
       if(num_bands_project<1) call io_error('Error: problem reading bands_plot_project')
       allocate(bands_plot_project(num_bands_project),stat=ierr)
       if (ierr/=0) call io_error('Error allocating bands_plot_project in param_read')
       call param_get_range_vector('bands_plot_project',found,num_bands_project,.false.,bands_plot_project)
       if (any(bands_plot_project<1) .or. any(bands_plot_project>num_wann) ) &
            call io_error('Error: bands_plot_project asks for a non-valid wannier function to be projected')
    endif 

    bands_num_spec_points=0
    call param_get_block_length('kpoint_path',found,i_temp)
    if (found) then
       bands_num_spec_points=i_temp*2
       allocate(bands_label(bands_num_spec_points),stat=ierr)
       if (ierr/=0) call io_error('Error allocating bands_label in param_read')
       allocate(bands_spec_points(3,bands_num_spec_points),stat=ierr)
       if (ierr/=0) call io_error('Error allocating bands_spec_points in param_read')
       call param_get_keyword_kpath
    end if
    if(.not.found .and. bands_plot) &
         call io_error('A bandstructure plot has been requested but there is no kpoint_path block') 

    ! checks
    if (bands_plot) then
       if ( (index(bands_plot_format,'gnu').eq.0) .and. (index(bands_plot_format,'xmgr').eq.0) ) &
            call io_error('Error: bands_plot_format not recognised')
       if ( (index(bands_plot_mode,'s-k').eq.0) .and. (index(bands_plot_mode,'cut').eq.0) ) &
            call io_error('Error: bands_plot_mode not recognised')
       if ( bands_num_points < 0 ) call io_error('Error: bands_num_points must be positive')       
    endif

    fermi_surface_plot        =  .false.
    call param_get_keyword('fermi_surface_plot',found,l_value=fermi_surface_plot)

    fermi_surface_num_points  = 50
    call param_get_keyword('fermi_surface_num_points',found,i_value=fermi_surface_num_points)
 
    fermi_surface_plot_format = 'xcrysden'
    call param_get_keyword('fermi_surface_plot_format',found,c_value=fermi_surface_plot_format)

    fermi_energy=0.0_dp
    call param_get_keyword('fermi_energy',found,r_value=fermi_energy)

    ! checks
    if (fermi_surface_plot) then
       if ( (index(fermi_surface_plot_format,'xcrys').eq.0) ) &
            call io_error('Error: fermi_surface_plot_format not recognised')    
       if ( fermi_surface_num_points < 0 ) call io_error('Error: fermi_surface_num_points must be positive')
    endif

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

    hr_plot                    = .false.
    call param_get_keyword('hr_plot',found,l_value=hr_plot)
                                                                                           
    hr_cutoff                 = 0.0_dp
    call param_get_keyword('hr_cutoff',found,r_value=hr_cutoff)

    dist_cutoff_mode             = 'three_dim'
    call param_get_keyword('dist_cutoff_mode',found,c_value=dist_cutoff_mode)
    if (      (index(dist_cutoff_mode,'three_dim').eq.0)  &
        .and. (index(dist_cutoff_mode,'two_dim')  .eq.0)  &
        .and. (index(dist_cutoff_mode,'one_dim')  .eq.0)) &
       call io_error('Error: dist_cutoff_mode not recognised')

    dist_cutoff                 = 1000.0_dp
    call param_get_keyword('dist_cutoff',found,r_value=dist_cutoff)

    one_dim_axis                = 'none'
    call param_get_keyword('one_dim_axis',found,c_value=one_dim_axis)
    one_dim_dir=0
    if (index(one_dim_axis,'x')>0 ) one_dim_dir = 1
    if (index(one_dim_axis,'y')>0 ) one_dim_dir = 2
    if (index(one_dim_axis,'z')>0 ) one_dim_dir = 3
    if ( transport.and. .not.tran_read_ht .and.(one_dim_dir.eq.0) ) call io_error('Error: one_dim_axis not recognised')
    if ( bands_plot.and.(index(bands_plot_mode,'cut').ne.0).and.(one_dim_dir.eq.0) ) &
         call io_error('Error: one_dim_axis not recognised')

301  continue

    !%%%%%%%%%%%%%%%%
    ! Transport 
    !%%%%%%%%%%%%%%%%

    transport_mode    = 'bulk'
    call param_get_keyword('transport_mode',found,c_value=transport_mode)

    if ( .not.tran_read_ht  .and. (index(transport_mode,'lcr').ne.0) ) &
       call io_error('Error: transport_mode.eq.lcr not compatible with tran_read_ht.eq.false')

    tran_win_min      = -3.0_dp
    call param_get_keyword('tran_win_min',found,r_value=tran_win_min)

    tran_win_max      =  3.0_dp
    call param_get_keyword('tran_win_max',found,r_value=tran_win_max)

    tran_energy_step  =  0.01_dp
    call param_get_keyword('tran_energy_step',found,r_value=tran_energy_step)

    tran_num_bb            = 0
    call param_get_keyword('tran_num_bb',found,i_value=tran_num_bb)

    tran_num_ll            = 0
    call param_get_keyword('tran_num_ll',found,i_value=tran_num_ll)

    tran_num_rr            = 0
    call param_get_keyword('tran_num_rr',found,i_value=tran_num_rr)

    tran_num_cc            = 0
    call param_get_keyword('tran_num_cc',found,i_value=tran_num_cc)

    tran_num_lc            = 0
    call param_get_keyword('tran_num_lc',found,i_value=tran_num_lc)

    tran_num_cr            = 0
    call param_get_keyword('tran_num_cr',found,i_value=tran_num_cr)

    tran_num_bandc         = 0
    call param_get_keyword('tran_num_bandc',found,i_value=tran_num_bandc)

    tran_write_ht          = .false. 
    call param_get_keyword('tran_write_ht',found,l_value=tran_write_ht)

    tran_use_same_lead     = .true.
    call param_get_keyword('tran_use_same_lead',found,l_value=tran_use_same_lead)

    ! checks
    if (transport) then
       if ( (index(transport_mode,'bulk').eq.0) .and. (index(transport_mode,'lcr').eq.0) ) &
            call io_error('Error: transport_mode not recognised')
       if ( tran_num_bb < 0 )    call io_error('Error: tran_num_bb < 0')
       if ( tran_num_ll < 0 )    call io_error('Error: tran_num_ll < 0')
       if ( tran_num_rr < 0 )    call io_error('Error: tran_num_rr < 0')
       if ( tran_num_cc < 0 )    call io_error('Error: tran_num_cc < 0')
       if ( tran_num_lc < 0 )    call io_error('Error: tran_num_lc < 0')
       if ( tran_num_cr < 0 )    call io_error('Error: tran_num_cr < 0')
       if ( tran_num_bandc < 0 ) call io_error('Error: tran_num_bandc < 0')
    endif

    if (transport .and. tran_read_ht) goto 302 

    !%%%%%%%%%%%%%%%%
    ! Disentanglement
    !%%%%%%%%%%%%%%%%

    disentanglement=.false.
    if(num_bands>num_wann) disentanglement=.true.

    ! Read the eigenvalues from wannier.eig
    eig_found=.false.
    if(.not. library) then
       allocate(eigval(num_bands,num_kpts),stat=ierr)
       if (ierr/=0) call io_error('Error allocating eigval in param_read')
       
       if(.not.postproc_setup)  then
          inquire(file=trim(seedname)//'.eig',exist=eig_found)
          if(.not. eig_found) then
             if ( disentanglement) then
                call io_error('No '//trim(seedname)//'.eig file found. Needed for disentanglement')
             else if ((bands_plot .or. dos_plot .or. fermi_surface_plot) ) then
                call io_error('No '//trim(seedname)//'.eig file found. Needed for interpolation')
             end if
          else
             eig_unit=io_file_unit()
             open(unit=eig_unit,file=trim(seedname)//'.eig',form='formatted',status='old',err=105)
             do k=1,num_kpts
                do n=1,num_bands
                   read(eig_unit,*,err=106,end=106) i,j,eigval(n,k)
                   if ((i.ne.n).or.(j.ne.k)) then
                      write(stdout,'(a)') 'Found a mismatch in '//trim(seedname)//'.eig' 
                      write(stdout,'(a,i0,a,i0)') 'Wanted band  : ',n,' found band  : ',i
                      write(stdout,'(a,i0,a,i0)') 'Wanted kpoint: ',k,' found kpoint: ',j
                      write(stdout,'(a)') ' '
                      write(stdout,'(a)') 'A common cause of this error is using the wrong'
                      write(stdout,'(a)') 'number of bands. Check your input files.'
                      write(stdout,'(a)') 'If your pseudopotentials have shallow core states remember'
                      write(stdout,'(a)') 'to account for these electrons.'   
                      write(stdout,'(a)') ' '
                      call io_error('param_read: mismatch in '//trim(seedname)//'.eig')
                   end if
                enddo
             end do
             close(eig_unit)
          end if
       end if
    end if

    if(library .and. allocated(eigval) ) eig_found=.true.

    dis_win_min=-1.0_dp;dis_win_max=0.0_dp
    if(eig_found) dis_win_min = minval(eigval)       
    call param_get_keyword('dis_win_min',found,r_value=dis_win_min)

    if(eig_found) dis_win_max = maxval(eigval)
    call param_get_keyword('dis_win_max',found,r_value=dis_win_max)
    if ( eig_found .and. (dis_win_max.lt.dis_win_min) ) &
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
    if (eig_found) then
       if ( dis_froz_max.lt.dis_froz_min ) &
            call io_error('Error: param_read: check disentanglement frozen windows')
       if(found2 .and. .not. found) &
            call io_error('Error: param_read: found dis_froz_min but not dis_froz_max')
    endif

    dis_num_iter      = 200
    call param_get_keyword('dis_num_iter',found,i_value=dis_num_iter)
    if (dis_num_iter<0) call io_error('Error: dis_num_iter must be positive')       

    dis_mix_ratio     = 0.5_dp
    call param_get_keyword('dis_mix_ratio',found,r_value=dis_mix_ratio)
    if (dis_mix_ratio<=0.0_dp .or. dis_mix_ratio>1.0_dp) &
         call io_error('Error: dis_mix_ratio must be greater than 0.0 but not greater than 1.0')


    dis_conv_tol      = 1.0e-10_dp      
    call param_get_keyword('dis_conv_tol',found,r_value=dis_conv_tol)
    if (dis_conv_tol<0.0_dp) call io_error('Error: dis_conv_tol must be positive')

    dis_conv_window=3   
    call param_get_keyword('dis_conv_window',found,i_value=dis_conv_window)
    if (dis_conv_window<0) call io_error('Error: dis_conv_window must be positive')       

    !%%%%%%%%%%%%%%%%
    !  Other Stuff
    !%%%%%%%%%%%%%%%%

    automatic_translation=.true.
    translation_centre_frac=0.0_dp
    call param_get_keyword_vector('translation_centre_frac',found,3,r_value=rv_temp)
    if (found) then
       translation_centre_frac=rv_temp
       automatic_translation=.false.
    endif

    use_bloch_phases = .false.
    call param_get_keyword('use_bloch_phases',found,l_value=use_bloch_phases)
    if(disentanglement .and. use_bloch_phases) &
         call io_error('Error: Cannot use bloch phases for disentanglement')

    search_shells                 = 12
    call param_get_keyword('search_shells',found,i_value=search_shells)
    if (search_shells<0) call io_error('Error: search_shells must be positive')       

    kmesh_tol=0.000001_dp
    call param_get_keyword('kmesh_tol',found,r_value=kmesh_tol)
    if (kmesh_tol<0.0_dp) call io_error('Error: kmesh_tol must be positive')

    num_shells  = 0
    call param_get_range_vector('shell_list',found,num_shells,lcount=.true.)
    if(found) then
       if(num_shells<0 .or. num_shells>max_shells) &
            call io_error('Error: number of shell in shell_list must be between zero and six')
       allocate(shell_list(num_shells),stat=ierr)
       if (ierr/=0) call io_error('Error allocating shell_list in param_read')
       call param_get_range_vector('shell_list',found,num_shells,.false.,shell_list)
       if (any(shell_list<1)  ) &
            call io_error('Error: shell_list must contain positive numbers')
    else
       allocate( shell_list(max_shells),stat=ierr)
       if (ierr/=0) call io_error('Error allocating shell_list in param_read')
    end if

    call param_get_keyword('num_shells',found,i_value=itmp)
    if(found .and. (itmp/=num_shells)) &
        call io_error('Error: Found obsolete keyword num_shells. Its value does not agree with shell_list')


    call param_get_keyword_block('unit_cell_cart',found,3,3,r_value=real_lattice_tmp)
    if(found.and.library) write(stdout,'(a)') ' Ignoring <unit_cell_cart> in input file'
    if (.not. library) then
       real_lattice=transpose(real_lattice_tmp)
       if(.not. found) call io_error('Error: Did not find the cell information in the input file')
    end if

    if(.not. library) &
         call utility_recip_lattice(real_lattice,recip_lattice,cell_volume)
    call utility_compute_metric(real_lattice,recip_lattice,real_metric,recip_metric)

    allocate ( kpt_cart(3,num_kpts) ,stat=ierr)
    if (ierr/=0) call io_error('Error allocating kpt_cart in param_read')
    if(.not. library) then
       allocate ( kpt_latt(3,num_kpts) ,stat=ierr)
       if (ierr/=0) call io_error('Error allocating kpt_latt in param_read')
    end if

    call param_get_keyword_block('kpoints',found,num_kpts,3,r_value=kpt_cart)
    if(found.and.library) write(stdout,'(a)') ' Ignoring <kpoints> in input file'
    if (.not. library) then
       kpt_latt=kpt_cart
       if(.not. found) call io_error('Error: Did not find the kpoint information in the input file')
    end if

    ! Calculate the kpoints in cartesian coordinates
    do nkp=1,num_kpts
       kpt_cart(:,nkp)=matmul(kpt_latt(:,nkp),recip_lattice(:,:))
    end do

    ! Atoms
    if (.not.library) num_atoms=0
    call param_get_block_length('atoms_frac',found,i_temp)
    if (found.and.library) write(stdout,'(a)') ' Ignoring <atoms_frac> in input file'
    call param_get_block_length('atoms_cart',found2,i_temp2,lunits)
    if (found2.and.library) write(stdout,'(a)') ' Ignoring <atoms_cart> in input file'
    if (.not.library) then
       if (found.and.found2) call io_error('Error: Cannot specify both atoms_frac and atoms_cart')
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
    endif

    ! Projections
    call param_get_block_length('projections',found,i_temp)
    if (found) call param_get_projections
    if (guiding_centres .and. .not. found .and. .not.(gamma_only.and.use_bloch_phases)) & 
       call io_error('param_read: Guiding centres requested, but no projection block found')

    ! check to see that there are no unrecognised keywords

302  continue

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

    if (transport .and. tran_read_ht) goto 303 

    ! For aesthetic purposes, convert some things to uppercase
    call param_uppercase()

303  continue

    deallocate(in_data,stat=ierr)
    if (ierr/=0) call io_error('Error deallocating in_data in param_read')

    if (transport .and. tran_read_ht) return 

    ! =============================== !
    ! Some checks and initialisations !
    ! =============================== !

!    if (restart.ne.' ') disentanglement=.false.

    if (disentanglement) then 
       allocate(ndimwin(num_kpts),stat=ierr)
       if (ierr/=0) call io_error('Error allocating ndimwin in param_read')
       allocate(lwindow(num_bands,num_kpts),stat=ierr)
       if (ierr/=0) call io_error('Error allocating lwindow in param_read')
    endif
    
    if ( wannier_plot .and. (index(wannier_plot_format,'cub').ne.0) ) then
       cosa(1)=dot_product(real_lattice(1,:),real_lattice(2,:))
       cosa(2)=dot_product(real_lattice(1,:),real_lattice(3,:))
       cosa(3)=dot_product(real_lattice(2,:),real_lattice(3,:))
       cosa = abs(cosa)
       if (any(cosa.gt.eps6)) &
            call io_error('Error: plotting in cube format requires orthogonal lattice vectors')
    endif

    ! Initialise
    omega_total = -999.0_dp
    omega_tilde = -999.0_dp
    omega_invariant = -999.0_dp
    have_disentangled = .false.

    allocate(wannier_centres(3,num_wann),stat=ierr)
    if (ierr/=0) call io_error('Error allocating wannier_centres in param_read')
    wannier_centres=0.0_dp
    allocate(wannier_spreads(num_wann),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating wannier_spreads in param_read')
    wannier_spreads=0.0_dp

    return

105 call io_error('Error: Problem opening eigenvalue file '//trim(seedname)//'.eig')
106 call io_error('Error: Problem reading eigenvalue file '//trim(seedname)//'.eig')

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

    do nsp=1,num_species
       ic=ichar(atoms_symbol(nsp)(1:1))
       if ((ic.ge.ichar('a')).and.(ic.le.ichar('z'))) &
            atoms_symbol(nsp)(1:1) = char(ic+ichar('Z')-ichar('z'))
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

    integer :: i,nkp,loop,nat,nsp

    if (transport .and. tran_read_ht) goto 401

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
    write(stdout,'(19x,a17,3x,f11.5)',advance='no') &
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
             write(stdout,'(1x,a1,1x,a2,1x,i3,3F10.5,3x,a1,1x,3F10.5,4x,a1)') '|',atoms_symbol(nsp),nat,atoms_pos_frac(:,nat,nsp),&
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
       write(stdout,'(1x,a)') '|     Frac. Coord.   l mr  r        z-axis               x-axis          Z/a |'
       write(stdout,'(1x,a)') '+----------------------------------------------------------------------------+'
       do nsp=1,num_proj
          write(stdout,'(1x,a1,3(1x,f5.2),1x,i2,1x,i2,1x,i2,3(1x,f6.3),3(1x,f6.3),&
               & 2x,f4.1,1x,a1)')  '|',proj_site(1,nsp),proj_site(2,nsp),&
               proj_site(3,nsp),proj_l(nsp), proj_m(nsp),proj_radial(nsp),&
               proj_z(1,nsp),proj_z(2,nsp),proj_z(3,nsp),proj_x(1,nsp),&
               proj_x(2,nsp),proj_x(3,nsp),proj_zona(nsp),'|'
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
    write(stdout,'(1x,a46,10x,I8,13x,a1)') '|  Number of Wannier Functions               :',num_wann,'|'
    write(stdout,'(1x,a46,10x,I8,13x,a1)') '|  Number of input Bloch states              :',num_bands,'|'
    write(stdout,'(1x,a46,10x,I8,13x,a1)') '|  Output verbosity (1=low, 5=high)          :',iprint,'|'
    write(stdout,'(1x,a46,10x,I8,13x,a1)') '|  Timing Level (1=low, 5=high)              :',timing_level,'|'
    write(stdout,'(1x,a46,10x,a8,13x,a1)') '|  Length Unit                               :',trim(length_unit),'|'  
    write(stdout,'(1x,a46,10x,L8,13x,a1)') '|  Post-processing setup (write *.nnkp)      :',postproc_setup,'|'
    write(stdout,'(1x,a46,10x,L8,13x,a1)') '|  Using Gamma-only branch of algorithms     :',gamma_only,'|'
    if(cp_pp .or. iprint>2) &
                  write(stdout,'(1x,a46,10x,L8,13x,a1)') '|  CP code post-processing                   :',cp_pp,'|'
    if(wannier_plot .or. iprint>2) then
       if(wvfn_formatted) then
             write(stdout,'(1x,a46,9x,a9,13x,a1)') '|  Wavefunction (UNK) file-type              :','formatted','|'
        else
         write(stdout,'(1x,a46,7x,a11,13x,a1)') '|  Wavefunction (UNK) file-type              :','unformatted','|'
        endif
       if(spin==1) then
         write(stdout,'(1x,a46,16x,a2,13x,a1)') '|  Wavefunction spin channel                 :','up','|' 
       else
         write(stdout,'(1x,a46,14x,a4,13x,a1)') '|  Wavefunction spin channel                 :','down','|'  
       endif
     endif  

    write(stdout,'(1x,a78)') '*----------------------------------------------------------------------------*'

    ! Wannierise
    write(stdout,'(1x,a78)') '*------------------------------- WANNIERISE ---------------------------------*'
    write(stdout,'(1x,a46,10x,I8,13x,a1)')   '|  Total number of iterations                :',num_iter,'|'
    write(stdout,'(1x,a46,10x,I8,13x,a1)')   '|  Number of CG steps before reset           :',num_cg_steps,'|'
    if (lfixstep) then
       write(stdout,'(1x,a46,10x,f8.3,13x,a1)')   '|  Fixed step length for minimisation        :',fixed_step,'|'
    else
       write(stdout,'(1x,a46,10x,f8.3,13x,a1)')   '|  Trial step length for line search         :',trial_step,'|'       
    endif
    write(stdout,'(1x,a46,8x,E10.3,13x,a1)') '|  Convergence tolerence                     :',conv_tol,'|'
    write(stdout,'(1x,a46,10x,I8,13x,a1)')   '|  Convergence window                        :',conv_window,'|'
    write(stdout,'(1x,a46,10x,I8,13x,a1)')   '|  Iterations between writing output         :',num_print_cycles,'|'
    write(stdout,'(1x,a46,10x,I8,13x,a1)')   '|  Iterations between backing up to disk     :',num_dump_cycles,'|'
    write(stdout,'(1x,a46,10x,L8,13x,a1)')   '|  Write r^2_nm to file                      :',write_r2mn,'|'
    write(stdout,'(1x,a46,10x,L8,13x,a1)')   '|  Write xyz WF centres to file              :',write_xyz,'|'
    write(stdout,'(1x,a46,10x,L8,13x,a1)')   '|  Use guiding centre to control phases      :',guiding_centres,'|'
    if(guiding_centres .or. iprint>2) then
    write(stdout,'(1x,a46,10x,I8,13x,a1)')   '|  Iterations before starting guiding centres:',num_no_guide_iter,'|'
    write(stdout,'(1x,a46,10x,I8,13x,a1)')   '|  Iterations between using guiding centres  :',num_guide_cycles,'|'
    end if
    write(stdout,'(1x,a78)') '*----------------------------------------------------------------------------*'
    !
    ! Disentanglement
    !
    if (disentanglement .or. iprint>2) then
       write(stdout,'(1x,a78)') '*------------------------------- DISENTANGLE --------------------------------*'
       write(stdout,'(1x,a46,10x,L8,13x,a1)')   '|  Using band disentanglement                :',disentanglement,'|'
       write(stdout,'(1x,a46,10x,I8,13x,a1)')   '|  Total number of iterations                :',dis_num_iter,'|'
       write(stdout,'(1x,a46,10x,F8.3,13x,a1)') '|  Mixing ratio                              :',dis_mix_ratio,'|'
       write(stdout,'(1x,a46,8x,ES10.3,13x,a1)') '|  Convergence tolerence                     :',dis_conv_tol,'|'
       write(stdout,'(1x,a46,10x,I8,13x,a1)')   '|  Convergence window                        :',dis_conv_window,'|'
       write(stdout,'(1x,a78)') '*----------------------------------------------------------------------------*'
    end if
    !
    ! Plotting
    !
    if (wannier_plot .or. bands_plot .or. fermi_surface_plot .or. slice_plot &
         .or. dos_plot .or. hr_plot .or. iprint>2) then
       !
       write(stdout,'(1x,a78)') '*-------------------------------- PLOTTING ----------------------------------*'
       !
       if (wannier_plot .or. iprint>2) then
          write(stdout,'(1x,a46,10x,L8,13x,a1)') '|  Plotting Wannier functions                :',wannier_plot,'|'
          write(stdout,'(1x,a46,10x,I8,13x,a1)') '|   Size of supercell for plotting           :',wannier_plot_supercell,'|'
          write(stdout,'(1x,a46,10x,a8,13x,a1)') '|   Plotting mode (molecule or crystal)      :',trim(wannier_plot_mode),'|'
          write(stdout,'(1x,a46,10x,a8,13x,a1)') '|   Plotting format                          :',trim(wannier_plot_format),'|'
          write(stdout,'(1x,a78)') '*----------------------------------------------------------------------------*'
       end if
       !
       if (fermi_surface_plot .or. iprint>2) then
          write(stdout,'(1x,a46,10x,L8,13x,a1)') '|  Plotting Fermi surface                    :',fermi_surface_plot,'|'
          write(stdout,'(1x,a46,10x,I8,13x,a1)') '|   Number of plotting points (along b_1)    :',fermi_surface_num_points,'|'
          write(stdout,'(1x,a46,10x,a8,13x,a1)') '|   Plotting format                          :' &
                                                                                 ,trim(fermi_surface_plot_format),'|'
          write(stdout,'(1x,a78)') '*----------------------------------------------------------------------------*'
       end if
       !
       if (bands_plot .or. iprint>2) then
          write(stdout,'(1x,a46,10x,L8,13x,a1)') '|  Plotting interpolated bandstructure       :',bands_plot,'|'
          write(stdout,'(1x,a46,10x,I8,13x,a1)') '|   Number of K-path sections                :',bands_num_spec_points/2,'|'
          write(stdout,'(1x,a46,10x,I8,13x,a1)') '|   Divisions along first K-path section     :',bands_num_points,'|'
          write(stdout,'(1x,a46,10x,a8,13x,a1)') '|   Output format                            :',trim(bands_plot_format),'|'
          write(stdout,'(1x,a46,10x,a8,13x,a1)') '|   Output mode                              :',trim(bands_plot_mode),'|'
          if(index(bands_plot_mode,'cut').ne.0) then
             write(stdout,'(1x,a46,10x,I8,13x,a1)')   '|   Dimension of the system                  :',bands_plot_dim,'|'
             if (bands_plot_dim .eq. 1) &
             write(stdout,'(1x,a46,10x,a8,13x,a1)')   '|   System extended in                       :',trim(one_dim_axis),'|'
             if (bands_plot_dim .eq. 2) &
             write(stdout,'(1x,a46,10x,a8,13x,a1)')   '|   System confined in                       :',trim(one_dim_axis),'|'
             write(stdout,'(1x,a46,10x,F8.3,13x,a1)') '|   Hamiltonian cut-off value                :',hr_cutoff,'|'
             write(stdout,'(1x,a46,10x,F8.3,13x,a1)') '|   Hamiltonian cut-off distance             :',dist_cutoff,'|'
             write(stdout,'(1x,a46,10x,a8,13x,a1)')   '|   Hamiltonian cut-off distance mode        :',trim(dist_cutoff_mode),'|'
          endif
          write(stdout,'(1x,a78)') '*----------------------------------------------------------------------------*'
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
       if (hr_plot .or. iprint>2) then
          write(stdout,'(1x,a46,10x,L8,13x,a1)')   '|  Plotting Hamiltonian in WF basis          :',hr_plot,'|'
          write(stdout,'(1x,a78)') '*----------------------------------------------------------------------------*'
       endif
       !
    endif

401  continue
    !
    ! Transport
    !
    if (transport.or.iprint>2 ) then
       !
       write(stdout,'(1x,a78)') '*------------------------------- TRANSPORT ----------------------------------*'
       !
       write(stdout,'(1x,a46,10x,a8,13x,a1)') '|  Transport mode                            :',trim(transport_mode),'|'
       !
       if (tran_read_ht) then
       !
       write(stdout,'(1x,a46,10x,a8,13x,a1)') '|   Hamiltonian from external files          :','T','|'
       !
       else
       !
       write(stdout,'(1x,a46,10x,a8,13x,a1)') '|   Hamiltonian from external files          :','F','|'
       write(stdout,'(1x,a46,10x,a8,13x,a1)') '|   System extended in                       :',trim(one_dim_axis),'|'
       !
       end if
       !
       write(stdout,'(1x,a78)') '*----------------------------------------------------------------------------*'
       !
    endif


101 format(20x,a3,2x,3F11.6)

  end subroutine param_write

  subroutine param_write_header
    use w90_io, only : io_date
    implicit none


    character (len=9) :: cdate, ctime

    call io_date(cdate, ctime)

    write(stdout,*)
    write(stdout,*)  '            +---------------------------------------------------+'
    write(stdout,*)  '            |                                                   |'
    write(stdout,*)  '            |                   WANNIER90                       |'
    write(stdout,*)  '            |                                                   |'
    write(stdout,*)  '            +---------------------------------------------------+'
    write(stdout,*)  '            |                                                   |'
    write(stdout,*)  '            |        Welcome to the Maximally-Localized         |'
    write(stdout,*)  '            |        Generalized Wannier Functions code         |'
    write(stdout,*)  '            |            http://www.wannier.org                 |'
    write(stdout,*)  '            |                                                   |'
    write(stdout,*)  '            |  Authors:                                         |'
    write(stdout,*)  '            |    Arash A. Mostofi   (Imperial College London)   |'
    write(stdout,*)  '            |    Jonathan R. Yates  (University of Cambridge)   |'
    write(stdout,*)  '            |    Young-Su Lee       (KIST, S. Korea)            |'
    write(stdout,*)  '            |                                                   |'
    write(stdout,*)  '            |                                                   |'
    write(stdout,*)  '            |  Please cite                                      |'
    write(stdout,*)  '            |                                                   |'
    write(stdout,*)  '            |  [ref] "Wannier90: A Tool for Obtaining Maximally |'
    write(stdout,*)  '            |         Localised Wannier Functions"              |'
    write(stdout,*)  '            |        A. A. Mostofi, J. R. Yates, Y.-S. Lee,     |'
    write(stdout,*)  '            |        I. Souza, D. Vanderbilt and N. Marzari     |'
    write(stdout,*)  '            |        Comput. Phys. Commun. 178, 685 (2008)      |'
    write(stdout,*)  '            |                                                   |'
    write(stdout,*)  '            |  in any publications arising from the use of      |'
    write(stdout,*)  '            |  this code.                                       |'
    write(stdout,*)  '            |                                                   |'
    write(stdout,*)  '            |                                                   |'
    write(stdout,*)  '            |  Wannier90 is based on routines written by        |'
    write(stdout,*)  '            |  Nicola Marzari, Ivo Souza and David Vanderbilt.  |'
    write(stdout,*)  '            |  For the method please cite                       |'
    write(stdout,*)  '            |                                                   |'
    write(stdout,*)  '            |  [ref] "Maximally Localized Generalised Wannier   |'
    write(stdout,*)  '            |         Functions for Composite Energy Bands"     |'
    write(stdout,*)  '            |         N. Marzari and D. Vanderbilt              |'
    write(stdout,*)  '            |         Phys. Rev. B 56 12847 (1997)              |'
    write(stdout,*)  '            |                                                   |'
    write(stdout,*)  '            |  [ref] "Maximally Localized Wannier Functions     |'
    write(stdout,*)  '            |         for Entangled Energy Bands"               |'
    write(stdout,*)  '            |         I. Souza, N. Marzari and D. Vanderbilt    |'
    write(stdout,*)  '            |         Phys. Rev. B 65 035109 (2001)             |'
    write(stdout,*)  '            |                                                   |'
    write(stdout,*)  '            |                                                   |'
    write(stdout,*)  '            | Copyright (c) 1997-2007 J. Yates, A. Mostofi,     |'
    write(stdout,*)  '            |   Y.-S Lee, N. Marzari, I. Souza, D. Vanderbilt   |'
    write(stdout,*)  '            |                                                   |'
    write(stdout,*)  '            |       Release: 1.1.1        10th May 2008         |'
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
    use w90_io, only : io_error

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
    if ( allocated (eigval) ) then
       deallocate ( eigval, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating eigval in param_dealloc')
    endif
    if ( allocated (shell_list) ) then
       deallocate ( shell_list, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating shell_list in param_dealloc')
    endif
    if ( allocated(kpt_latt) ) then
       deallocate ( kpt_latt, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating kpt_latt in param_dealloc')
    endif
    if ( allocated(kpt_cart) ) then
       deallocate ( kpt_cart, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating kpt_cart in param_dealloc')
    endif
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
    if ( allocated ( atoms_symbol ) ) then
       deallocate (  atoms_symbol, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating atoms_symbol in param_dealloc')
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
    if ( allocated( wannier_plot_list )  ) then
       deallocate( wannier_plot_list, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating wannier_plot_list in param_dealloc')
    end if
    if( allocated( exclude_bands ) ) then
       deallocate( exclude_bands, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating exclude_bands in param_dealloc')
    end if
    if( allocated( wannier_centres ) ) then
       deallocate( wannier_centres, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating wannier_centres in param_dealloc')
    end if
    if( allocated( wannier_spreads ) ) then
       deallocate( wannier_spreads, stat=ierr  )
       if (ierr/=0) call io_error('Error in deallocating wannier_spreads in param_dealloc')
    endif
    return

  end subroutine param_dealloc


!!$  !================================!
!!$  subroutine param_write_um
!!$    !================================!
!!$    !                                !
!!$    ! Dump the U and M to *_um.dat   !
!!$    !                                !
!!$    !================================!
!!$
!!$
!!$    use w90_io,        only : io_file_unit,io_error,seedname,io_date
!!$    implicit none
!!$
!!$    integer :: i,j,k,l,um_unit
!!$    character (len=9) :: cdate, ctime
!!$    character(len=33) :: header
!!$
!!$    call io_date(cdate, ctime)
!!$    header='written on '//cdate//' at '//ctime
!!$
!!$    um_unit=io_file_unit()
!!$    open(unit=um_unit,file=trim(seedname)//'_um.dat',form='unformatted')
!!$    write(um_unit) header
!!$    write(um_unit) omega_invariant
!!$    write(um_unit) num_wann,num_kpts,num_nnmax    
!!$    write(um_unit) (((u_matrix(i,j,k),i=1,num_wann),j=1,num_wann),k=1,num_kpts)
!!$    write(um_unit) ((((m_matrix(i,j,k,l),i=1,num_wann),j=1,num_wann),k=1,nntot),l=1,num_kpts)
!!$    close(um_unit)
!!$
!!$    return
!!$
!!$  end subroutine param_write_um


!!$  !================================!
!!$  subroutine param_read_um
!!$    !================================!
!!$    !                                !
!!$    ! Restore U and M from file      !
!!$    !                                !
!!$    !================================!
!!$
!!$    use w90_io,        only : io_file_unit,io_error,seedname
!!$    implicit none
!!$
!!$    integer       :: tmp_num_wann,tmp_num_kpts,tmp_num_nnmax    
!!$    integer       :: i,j,k,l,um_unit,ierr
!!$    character(len=33) :: header
!!$    real(kind=dp) :: tmp_omi
!!$
!!$    um_unit=io_file_unit()
!!$    open(unit=um_unit,file=trim(seedname)//'_um.dat',status="old",form='unformatted',err=105)
!!$    read(um_unit) header
!!$    write(stdout,'(1x,4(a))') 'Reading U and M from file ',trim(seedname),'_um.dat ', header 
!!$    read(um_unit) tmp_omi
!!$    if ( have_disentangled ) then
!!$       if ( abs(tmp_omi-omega_invariant).gt.1.0e-10_dp )  &
!!$            call io_error('Error in restart: omega_invariant in .chk and um.dat files do not match')
!!$    endif
!!$    read(um_unit) tmp_num_wann,tmp_num_kpts,tmp_num_nnmax    
!!$    if(tmp_num_wann/=num_wann) call io_error('Error in param_read_um: num_wann mismatch')
!!$    if(tmp_num_kpts/=num_kpts) call io_error('Error in param_read_um: num_kpts mismatch')
!!$    if(tmp_num_nnmax/=num_nnmax) call io_error('Error in param_read_um: num_nnmax mismatch')
!!$    if (.not.allocated(u_matrix)) then
!!$       allocate(u_matrix(num_wann,num_wann,num_kpts),stat=ierr)
!!$       if (ierr/=0) call io_error('Error allocating u_matrix in param_read_um')
!!$    endif
!!$    read(um_unit) (((u_matrix(i,j,k),i=1,num_wann),j=1,num_wann),k=1,num_kpts)
!!$    if (.not.allocated(m_matrix)) then
!!$       allocate(m_matrix(num_wann,num_wann,nntot,num_kpts),stat=ierr)
!!$       if (ierr/=0) call io_error('Error allocating m_matrix in param_read_um')
!!$    endif
!!$    read(um_unit) ((((m_matrix(i,j,k,l),i=1,num_wann),j=1,num_wann),k=1,nntot),l=1,num_kpts)
!!$    close(um_unit)
!!$
!!$    return
!!$
!!$105 call io_error('Error: Problem opening file '//trim(seedname)//'_um.dat in param_read_um')
!!$
!!$  end subroutine param_read_um



  !=================================================!
  subroutine param_write_chkpt(chkpt)
    !=================================================!
    ! Write checkpoint file                           !
    !=================================================!

    use w90_io, only : io_file_unit,io_date,seedname

    implicit none

    character(len=*), intent(in) :: chkpt

    integer :: chk_unit,nkp,i,j,k,l
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
    write(chk_unit) (((u_matrix(i,j,k),i=1,num_wann),j=1,num_wann),k=1,num_kpts)               ! U_matrix
    write(chk_unit) ((((m_matrix(i,j,k,l),i=1,num_wann),j=1,num_wann),k=1,nntot),l=1,num_kpts) ! M_matrix
    write(chk_unit) ((wannier_centres(i,j),i=1,3),j=1,num_wann)

    close(chk_unit)

    write(stdout,'(a/)') ' done'

    return

  end subroutine param_write_chkpt


  !=======================================!
  subroutine param_read_chkpt
    !=======================================!
    ! Read checkpoint file                  !
    !=======================================!

    use w90_constants, only : eps6
    use w90_io,        only : io_error,io_file_unit,stdout,seedname

    implicit none

    integer :: chk_unit,nkp,i,j,k,l,ntmp,ierr
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
          if (abs(tmp_latt(i,j)-real_lattice(i,j)).gt.eps6) &
               call io_error('param_read_chk: Mismatch in real_lattice')
       enddo
    enddo
    read(chk_unit) ((tmp_latt(i,j),i=1,3),j=1,3)  ! Reciprocal lattice
    do j=1,3
       do i=1,3
          if (abs(tmp_latt(i,j)-recip_lattice(i,j)).gt.eps6) &
               call io_error('param_read_chk: Mismatch in recip_lattice')
       enddo
    enddo
    read(chk_unit) ntmp                ! K-points
    if (ntmp.ne.num_kpts) &
         call io_error('param_read_chk: Mismatch in num_kpts')
    read(chk_unit) ((tmp_kpt_latt(i,nkp),i=1,3),nkp=1,num_kpts)
    do nkp=1,num_kpts
       do i=1,3
          if (abs(tmp_kpt_latt(i,nkp)-kpt_latt(i,nkp)).gt.eps6) &
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

       ! lwindow
       if (.not.allocated(lwindow)) then
          allocate(lwindow(num_bands,num_kpts),stat=ierr)
          if (ierr/=0) call io_error('Error allocating lwindow in param_read_chkpt')
       endif
       read(chk_unit,err=122) ((lwindow(i,nkp),i=1,num_bands),nkp=1,num_kpts)

       ! ndimwin
       if (.not.allocated(ndimwin)) then
          allocate(ndimwin(num_kpts),stat=ierr)
          if (ierr/=0) call io_error('Error allocating ndimwin in param_read_chkpt')
       endif
       read(chk_unit,err=123) (ndimwin(nkp),nkp=1,num_kpts)

       ! U_matrix_opt
       if (.not.allocated(u_matrix_opt)) then
          allocate(u_matrix_opt(num_bands,num_wann,num_kpts),stat=ierr)
          if (ierr/=0) call io_error('Error allocating u_matrix_opt in param_read_chkpt')
       endif
       read(chk_unit,err=124) (((u_matrix_opt(i,j,nkp),i=1,num_bands),j=1,num_wann),nkp=1,num_kpts)

    endif

    ! U_matrix
    if (.not.allocated(u_matrix)) then
       allocate(u_matrix(num_wann,num_wann,num_kpts),stat=ierr)
       if (ierr/=0) call io_error('Error allocating u_matrix in param_read_chkpt')
    endif
    read(chk_unit,err=125) (((u_matrix(i,j,k),i=1,num_wann),j=1,num_wann),k=1,num_kpts)

    ! M_matrix
    if (.not.allocated(m_matrix)) then
       allocate(m_matrix(num_wann,num_wann,nntot,num_kpts),stat=ierr)
       if (ierr/=0) call io_error('Error allocating m_matrix in param_read_chkpt')
    endif
    read(chk_unit,err=126) ((((m_matrix(i,j,k,l),i=1,num_wann),j=1,num_wann),k=1,nntot),l=1,num_kpts)

    ! wannier_centres
    read(chk_unit,err=127) ((wannier_centres(i,j),i=1,3),j=1,num_wann)

    close(chk_unit)

    write(stdout,'(a/)') ' ... done'

    return

121 call io_error('Error opening '//trim(seedname)//'.chk in param_read_chkpt')
122 call io_error('Error reading lwindow from '//trim(seedname)//'.chk in param_read_chkpt')
123 call io_error('Error reading ndimwin from '//trim(seedname)//'.chk in param_read_chkpt')
124 call io_error('Error reading u_matrix_opt from '//trim(seedname)//'.chk in param_read_chkpt')
125 call io_error('Error reading u_matrix from '//trim(seedname)//'.chk in param_read_chkpt')
126 call io_error('Error reading m_matrix from '//trim(seedname)//'.chk in param_read_chkpt')
127 call io_error('Error reading wannier_centres from '//trim(seedname)//'.chk in param_read_chkpt')

  end subroutine param_read_chkpt


  !=======================================!
  subroutine param_in_file
    !=======================================!
    ! Load the *.win file into a character  !
    ! array in_file, ignoring comments and  !
    ! blank lines and converting everything !
    ! to lowercase characters               !
    !=======================================!

    use w90_io,        only : io_file_unit,io_error,seedname
    use w90_utility,   only : utility_lowercase

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

    allocate(in_data(num_lines),stat=ierr)
    if (ierr/=0) call io_error('Error allocating in_data in param_in_file')

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

    use w90_io,        only : io_error

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
       if( present(i_value) ) read(dummy,*,err=220,end=220) i_value
       if( present(r_value) ) read(dummy,*,err=220,end=220) r_value
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

    use w90_io,        only : io_error

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
       if( present(c_value) ) read(dummy,*,err=230,end=230) (c_value(i),i=1,length)
       if( present(l_value) ) then
          ! I don't think we need this. Maybe read into a dummy charater
          ! array and convert each element to logical
       endif
       if( present(i_value) ) read(dummy,*,err=230,end=230) (i_value(i),i=1,length)
       if( present(r_value) ) read(dummy,*,err=230,end=230) (r_value(i),i=1,length)
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

    use w90_io,        only : io_error

    implicit none

    character(*),      intent(in)  :: keyword
    logical          , intent(out) :: found
    integer,           intent(out)  :: length

    integer           :: kl, in,loop,pos
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


  end subroutine param_get_vector_length


  !==============================================================================================!
  subroutine param_get_keyword_block(keyword,found,rows,columns,c_value,l_value,i_value,r_value)
    !==============================================================================================!
    !                                                                                              !
    !                           Finds the values of the required data block                        !
    !                                                                                              !
    !==============================================================================================!

    use w90_constants, only : bohr
    use w90_io,        only : io_error

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
       if( present(c_value) ) read(dummy,*,err=240,end=240) (c_value(i,counter),i=1,columns)
       if( present(l_value) ) then
          ! I don't think we need this. Maybe read into a dummy charater
          ! array and convert each element to logical
       endif
       if( present(i_value) ) read(dummy,*,err=240,end=240) (i_value(i,counter),i=1,columns)
       if( present(r_value) ) read(dummy,*,err=240,end=240) (r_value(i,counter),i=1,columns)
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

    use w90_io,        only : io_error

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

    ! Ignore atoms_cart and atoms_frac blocks if running in library mode
    if (library) then
       if (trim(keyword).eq.'atoms_cart' .or. trim(keyword).eq.'atoms_frac') then
          in_data(line_s:line_e)(1:maxlen) = ' '
       endif
    endif

    if (present(lunits)) then
       dummy=in_data(line_s+1)
       !       write(stdout,*) dummy
       !       write(stdout,*) trim(dummy)
       read(dummy,*,end=555) atsym, (atpos(i),i=1,3)
       lunits=.false.
    endif

    if(rows<=0) then !cope with empty blocks
       found=.false.
       in_data(line_s:line_e)(1:maxlen) = ' '
    end if


    return

555 lunits=.true.

    if(rows<=1) then !cope with empty blocks
       found=.false.
       in_data(line_s:line_e)(1:maxlen) = ' '
    end if


    return

  end subroutine param_get_block_length


  !===================================!
  subroutine param_get_atoms(lunits)
    !===================================!
    !                                   !
    !   Fills the atom data block       !
    !                                   !
    !===================================!

    use w90_constants, only : bohr
    use w90_utility,   only : utility_frac_to_cart,utility_cart_to_frac
    use w90_io,        only : io_error
    implicit none

    logical, intent(in) :: lunits

    real(kind=dp)     :: atoms_pos_frac_tmp(3,num_atoms)
    real(kind=dp)     :: atoms_pos_cart_tmp(3,num_atoms)
    character(len=20) :: keyword
    integer           :: in,ins,ine,loop,i,line_e,line_s,counter
    integer           :: i_temp,loop2,max_sites,ierr,ic
    logical           :: found_e,found_s,found,frac
    character(len=maxlen) :: dummy,end_st,start_st
    character(len=maxlen) :: ctemp(num_atoms)
    character(len=maxlen) :: atoms_label_tmp(num_atoms)
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
          read(dummy,*,err=240,end=240) atoms_label_tmp(counter),(atoms_pos_frac_tmp(i,counter),i=1,3)
       else
          read(dummy,*,err=240,end=240) atoms_label_tmp(counter),(atoms_pos_cart_tmp(i,counter),i=1,3)
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

    allocate(atoms_species_num(num_species),stat=ierr)
       if (ierr/=0) call io_error('Error allocating atoms_species_num in param_get_atoms')
    allocate(atoms_label(num_species),stat=ierr)
       if (ierr/=0) call io_error('Error allocating atoms_label in param_get_atoms')
    allocate(atoms_symbol(num_species),stat=ierr)
       if (ierr/=0) call io_error('Error allocating atoms_symbol in param_get_atoms')
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
    allocate(atoms_pos_frac(3,max_sites,num_species),stat=ierr)
       if (ierr/=0) call io_error('Error allocating atoms_pos_frac in param_get_atoms')
    allocate(atoms_pos_cart(3,max_sites,num_species),stat=ierr)
       if (ierr/=0) call io_error('Error allocating atoms_pos_cart in param_get_atoms')

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

    ! Strip any numeric characters from atoms_label to get atoms_symbol
    do loop=1,num_species    
       atoms_symbol(loop)(1:2)=atoms_label(loop)(1:2)
       ic=ichar(atoms_symbol(loop)(2:2))
       if ((ic.lt.ichar('a')).or.(ic.gt.ichar('z'))) &
         atoms_symbol(loop)(2:2)=' '
    end do

    return

240 call io_error('Error: Problem reading block keyword '//trim(keyword))

  end subroutine param_get_atoms

    !=====================================================!
     subroutine param_lib_set_atoms(atoms_label_tmp,atoms_pos_cart_tmp)
    !=====================================================!
    !                                                     !
    !   Fills the atom data block during a library call   !
    !                                                     !
    !=====================================================!

    use w90_utility,   only : utility_cart_to_frac, utility_lowercase
    use w90_io,        only : io_error

    implicit none

    character(len=*), intent(in) :: atoms_label_tmp(num_atoms)
    real(kind=dp), intent(in)      :: atoms_pos_cart_tmp(3,num_atoms)

    real(kind=dp)     :: atoms_pos_frac_tmp(3,num_atoms)
    integer           :: loop2,max_sites,ierr,ic,loop,counter
    character(len=maxlen) :: ctemp(num_atoms)
    character(len=maxlen) :: tmp_string


    do loop=1,num_atoms
       call utility_cart_to_frac(atoms_pos_cart_tmp(:,loop),&
            atoms_pos_frac_tmp(:,loop),recip_lattice)
    enddo
       
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

    allocate(atoms_species_num(num_species),stat=ierr)
       if (ierr/=0) call io_error('Error allocating atoms_species_num in param_lib_set_atoms')
    allocate(atoms_label(num_species),stat=ierr)
       if (ierr/=0) call io_error('Error allocating atoms_label in param_lib_set_atoms')
    allocate(atoms_symbol(num_species),stat=ierr)
       if (ierr/=0) call io_error('Error allocating atoms_symbol in param_lib_set_atoms')
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
    allocate(atoms_pos_frac(3,max_sites,num_species),stat=ierr)
    if (ierr/=0) call io_error('Error allocating atoms_pos_frac in param_lib_set_atoms')
    allocate(atoms_pos_cart(3,max_sites,num_species),stat=ierr)
    if (ierr/=0) call io_error('Error allocating atoms_pos_cart in param_lib_set_atoms')
    
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

    ! Strip any numeric characters from atoms_label to get atoms_symbol
    do loop=1,num_species    
       atoms_symbol(loop)(1:2)=atoms_label(loop)(1:2)
       ic=ichar(atoms_symbol(loop)(2:2))
       if ((ic.lt.ichar('a')).or.(ic.gt.ichar('z'))) &
         atoms_symbol(loop)(2:2)=' '
       tmp_string = trim(adjustl(utility_lowercase(atoms_symbol(loop))))
       atoms_symbol(loop)(1:2)=tmp_string(1:2)
       tmp_string = trim(adjustl(utility_lowercase(atoms_label(loop))))
       atoms_label(loop)(1:2)=tmp_string(1:2)
    end do

    return


  end subroutine param_lib_set_atoms


    !====================================================================!
    subroutine param_get_range_vector(keyword,found,length,lcount,i_value)
    !====================================================================!
    !   Read a range vector eg. 1,2,3,4-10  or 1 3 400:100               !
    !   if(lcount) we return the number of states in length              !
    !====================================================================!
    use w90_io,        only : io_error

    implicit none

    character(*),      intent(in)    :: keyword
    logical          , intent(out)   :: found
    integer,           intent(inout) :: length
    logical,           intent(in)    :: lcount
    integer, optional, intent(out)   :: i_value(length)

    integer   :: kl, in,loop,num1,num2,i_punc
    integer   :: counter,i_digit,loop_r,range_size
    character(len=maxlen) :: dummy
    character(len=10), parameter :: c_digit="0123456789"
    character(len=2) , parameter :: c_range="-:"
    character(len=3) , parameter :: c_sep=" ,;"
    character(len=5) , parameter :: c_punc=" ,;-:"
    character(len=5)  :: c_num1,c_num2

    
    if(lcount .and. present(i_value) ) call io_error('param_get_range_vector: incorrect call')

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
       if(.not. lcount) in_data(loop)(1:maxlen) = ' '
       if( dummy(1:1)=='=' .or. dummy(1:1)==':') then
          dummy=dummy(2:)
          dummy=adjustl(dummy)
       end if
    end do

    if(.not. found) return

    counter=0
    if (len_trim(dummy)==0) call io_error('Error: keyword '//trim(keyword)//' is blank')
    dummy=adjustl(dummy)
    do 
       i_punc=scan(dummy,c_punc)
       if(i_punc==0) call io_error('Error parsing keyword '//trim(keyword)) 
       c_num1=dummy(1:i_punc-1)
       read(c_num1,*,err=101,end=101) num1
       dummy=adjustl(dummy(i_punc:))
       !look for range
       if(scan(dummy,c_range)==1) then
          i_digit=scan(dummy,c_digit)
          dummy=adjustl(dummy(i_digit:))
          i_punc=scan(dummy,c_punc)
          c_num2=dummy(1:i_punc-1)
          read(c_num2,*,err=101,end=101) num2
          dummy=adjustl(dummy(i_punc:))
          range_size=abs(num2-num1)+1
          do loop_r=1,range_size
             counter=counter+1
             if(.not. lcount) i_value(counter)=min(num1,num2)+loop_r-1
          end do
       else
          counter=counter+1 
          if(.not. lcount) i_value(counter)=num1
       end if

       if(scan(dummy,c_sep)==1) dummy=adjustl(dummy(2:))
       if(scan(dummy,c_range)==1) call io_error('Error parsing keyword '//trim(keyword)//' incorrect range') 
       if(index(dummy,' ')==1) exit
    end do

    if(lcount) length=counter
    if(.not.lcount) then
       do loop=1,counter-1
          do loop_r=loop+1,counter 
             if(i_value(loop)==i_value(loop_r)) &
                call io_error('Error parsing keyword '//trim(keyword)//' duplicate values')
          end do
        end do
    end if

    return

101 call io_error('Error parsing keyword '//trim(keyword))


   end  subroutine param_get_range_vector

  !===================================!
   subroutine param_get_projections
     !===================================!
     !                                   !
     !  Fills the projection data block  !
     !                                   !
     !===================================!

     use w90_constants, only : bohr,eps6,eps2
     use w90_utility,   only : utility_cart_to_frac,&
          utility_string_to_coord,utility_strip
     use w90_io,        only : io_error

     implicit none


     real(kind=dp)     :: pos_frac(3)
     real(kind=dp)     :: pos_cart(3)
     character(len=20) :: keyword
     integer           :: in,ins,ine,loop,line_e,line_s,counter
     integer           :: sites,species,line,pos1,pos2,pos3,m_tmp,l_tmp,mstate
     integer           :: loop_l,loop_m,loop_sites,ierr
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
     real(kind=dp), parameter :: proj_z_def(3)=(/0.0_dp,0.0_dp,1.0_dp/)
     real(kind=dp), parameter :: proj_x_def(3)=(/1.0_dp,0.0_dp,0.0_dp/)
     real(kind=dp), parameter :: proj_zona_def=1.0_dp
     integer, parameter       :: proj_radial_def=1
     !
     real(kind=dp) :: proj_z_tmp(3)
     real(kind=dp) :: proj_x_tmp(3)
     real(kind=dp) :: proj_zona_tmp
     integer       :: proj_radial_tmp
     logical       :: lconvert,lrandom
     logical       :: lpartrandom
     !
     real(kind=dp) :: xnorm,znorm,cosphi,sinphi,xnorm_new,cosphi_new

     keyword="projections"

     found_s=.false.
     found_e=.false.

     start_st='begin '//trim(keyword)
     end_st='end '//trim(keyword)

     num_proj=num_wann
     if(spinors) num_proj=num_wann/2


     allocate( proj_site(3,num_proj),stat=ierr)
     if (ierr/=0) call io_error('Error allocating proj_site in param_get_projections') 
     allocate( proj_l(num_proj) ,stat=ierr)
     if (ierr/=0) call io_error('Error allocating proj_l in param_get_projections') 
     allocate( proj_m(num_proj)  ,stat=ierr)
     if (ierr/=0) call io_error('Error allocating proj_m in param_get_projections')
     allocate( proj_z(3,num_proj) ,stat=ierr)
     if (ierr/=0) call io_error('Error allocating proj_z in param_get_projections')
     allocate( proj_x(3,num_proj) ,stat=ierr)
     if (ierr/=0) call io_error('Error allocating proj_x in param_get_projections')
     allocate( proj_radial(num_proj)   ,stat=ierr)
     if (ierr/=0) call io_error('Error allocating proj_radial in param_get_projections')
     allocate( proj_zona(num_proj) ,stat=ierr)
     if (ierr/=0) call io_error('Error allocating proj_zona in param_get_projections')




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

     dummy=in_data(line_s+1)
     lconvert=.false.
     lrandom=.false.
     lpartrandom=.false.
     if ( index(dummy,'ang').ne.0 ) then
        in_data(line_s)(1:maxlen) = ' '
        line_s = line_s + 1
     elseif ( index(dummy,'bohr').ne.0 ) then
        in_data(line_s)(1:maxlen) = ' '
        line_s = line_s + 1
        lconvert=.true.
     elseif ( index(dummy,'random').ne.0 ) then
        in_data(line_s)(1:maxlen) = ' '
        line_s = line_s + 1
        if (index(in_data(line_s+1),end_st).ne.0) then
           lrandom=.true.     ! all projections random
        else
           lpartrandom=.true. ! only some projections random
           if ( index(in_data(line_s+1),'ang').ne.0 ) then
              in_data(line_s)(1:maxlen) = ' '
              line_s = line_s + 1
           elseif ( index(in_data(line_s+1),'bohr').ne.0 ) then
              in_data(line_s)(1:maxlen) = ' '
              line_s = line_s + 1
              lconvert=.true.
           endif
        endif
     endif

     if(.not. lrandom) then
        counter=0
        do line=line_s+1,line_e-1
           ang_states=0
           !Assume the default values
           proj_z_tmp      = proj_z_def  
           proj_x_tmp      = proj_x_def  
           proj_zona_tmp   = proj_zona_def  
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
              if (lconvert) pos_cart = pos_cart * bohr
              call utility_cart_to_frac (pos_cart(:),pos_frac(:),recip_lattice)          
           elseif(index(ctemp,'f=')>0) then
              sites=-1
              ctemp=ctemp(3:)
              call utility_string_to_coord(ctemp,pos_frac)
           else
              if(num_species==0) &
                   call io_error('param_get_projection: Atom centred projection requested but no atoms defined')
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
                    read(ctemp3(3:mstate-1),*,err=101,end=101) l_tmp
                 else
                    read(ctemp3(3:),*,err=101,end=101) l_tmp
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
                       read(ctemp5(1:),*,err=102,end=102) m_tmp
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
                    read(ctemp4(1:),*,err=106,end=106) m_string
                    select case (trim(adjustl(m_string)))
                    case ('s')
                       ang_states(1,0)=1
                    case ('p')
                       ang_states(1:3,1)=1
                    case ('pz')
                       ang_states(1,1)=1
                    case ('px')
                       ang_states(2,1)=1
                    case ('py')
                       ang_states(3,1)=1
                    case ('d')
                       ang_states(1:5,2)=1
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
                    case ('f')
                       ang_states(1:7,3)=1
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
              ! diffusivity of orbital
              pos1=index(dummy,'zona=')
              if(pos1>0) then
                 ctemp=(dummy(pos1+5:))
                 pos2=index(ctemp,':')
                 if(pos2>0) ctemp=ctemp(:pos2-1)
                 read(ctemp,*,err=104,end=104) proj_zona_tmp
              endif
              ! nodes for the radial part
              pos1=index(dummy,'r=')
              if(pos1>0) then
                 ctemp=(dummy(pos1+2:))
                 pos2=index(ctemp,':')
                 if(pos2>0) ctemp=ctemp(:pos2-1)
                 read(ctemp,*,err=105,end=105) proj_radial_tmp
              endif
           end if
           if(sites==-1) then
              if(counter+sum(ang_states) > num_proj) call io_error('param_get_projection: &
                   &too many projections defined')
           else
              if(counter+sites*sum(ang_states) > num_proj) call io_error('param_get_projection:&
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
                       end if
                    end do
                 end do
              end do
           end if

        end do !end loop over projection block

        ! check there are enough projections and add random projections if required
        if (.not. lpartrandom) then
              if (counter.ne.num_proj) call io_error('param_get_projections:&
                   & Fewer projections defined than the number of Wannier functions requested')
        else
           call random_seed()
           do loop=counter+1,num_proj
              call random_number(proj_site(:,loop))
              proj_l(loop)      = 0
              proj_m(loop)      = 1
              proj_z(:,loop)    = proj_z_def  
              proj_x(:,loop)    = proj_x_def  
              proj_zona(loop)   = proj_zona_def  
              proj_radial(loop) = proj_radial_def             
           enddo
        endif

     elseif(lrandom) then

        call random_seed() ! comment out this line for reproducible random positions!
        do loop=1,num_proj
           call random_number(proj_site(:,loop))
           proj_l(loop)      = 0
           proj_m(loop)      = 1
           proj_z(:,loop)    = proj_z_def  
           proj_x(:,loop)    = proj_x_def  
           proj_zona(loop)   = proj_zona_def  
           proj_radial(loop) = proj_radial_def
        end do

     end if

     in_data(line_s:line_e)(1:maxlen) = ' '

!!$     ! Check
!!$     do loop=1,num_wann
!!$        if ( abs(sum(proj_z(:,loop)*proj_x(:,loop))).gt.1.0e-6_dp ) then
!!$           write(stdout,*) ' Projection:',loop
!!$           call io_error(' Error in projections: z and x axes are not orthogonal')
!!$        endif
!!$     enddo

     ! Normalise z-axis and x-axis and check/fix orthogonality
     do loop=1,num_proj

        znorm=sqrt(sum(proj_z(:,loop)*proj_z(:,loop)))
        xnorm=sqrt(sum(proj_x(:,loop)*proj_x(:,loop)))
        proj_z(:,loop)=proj_z(:,loop)/znorm             ! normalise z
        proj_x(:,loop)=proj_x(:,loop)/xnorm             ! normalise x
        cosphi=sum(proj_z(:,loop)*proj_x(:,loop))       

        ! Check whether z-axis and z-axis are orthogonal
        if ( abs(cosphi).gt.eps6 ) then

           ! Special case of circularly symmetric projections (pz, dz2, fz3)
           ! just choose an x-axis that is perpendicular to the given z-axis
           if ( (proj_l(loop).ge.0) .and. (proj_m(loop).eq.1) ) then
              proj_x_tmp(:) = proj_x(:,loop)            ! copy of original x-axis
              call random_seed()
              call random_number(proj_z_tmp(:))         ! random vector
              ! calculate new x-axis as the cross (vector) product of random vector with z-axis
              proj_x(1,loop)=proj_z_tmp(2)*proj_z(3,loop) - proj_z_tmp(3)*proj_z(2,loop)
              proj_x(2,loop)=proj_z_tmp(3)*proj_z(1,loop) - proj_z_tmp(1)*proj_z(3,loop)
              proj_x(3,loop)=proj_z_tmp(1)*proj_z(2,loop) - proj_z_tmp(2)*proj_z(1,loop)
              xnorm_new=sqrt(sum(proj_x(:,loop)*proj_x(:,loop)))
              proj_x(:,loop)=proj_x(:,loop)/xnorm_new   ! normalise
              goto 555
           endif

           ! If projection axes non-orthogonal enough, then
           ! user may have made a mistake and should check
           if ( abs(cosphi).gt.eps2 ) then 
              write(stdout,*) ' Projection:',loop
              call io_error(' Error in projections: z and x axes are not orthogonal')
           endif

           ! If projection axes are "reasonably orthogonal", project x-axis
           ! onto plane perpendicular to z-axis to make them more so
           sinphi=sqrt(1-cosphi*cosphi)
           proj_x_tmp(:) = proj_x(:,loop)               ! copy of original x-axis
           ! calculate new x-axis:
           ! x = z \cross (x_tmp \cross z) / sinphi = ( x_tmp - z(z.x_tmp) ) / sinphi
           proj_x(:,loop) = ( proj_x_tmp(:) - cosphi*proj_z(:,loop) ) / sinphi

           ! Final check
555        cosphi_new=sum(proj_z(:,loop)*proj_x(:,loop))
           if ( abs(cosphi_new).gt.eps6 ) then
              write(stdout,*) ' Projection:',loop
              call io_error(' Error: z and x axes are still not orthogonal after projection')
           endif

        endif

     enddo

     return


101  call io_error('param_get_projection: Problem reading l state into integer '//trim(ctemp3))
102  call io_error('param_get_projection: Problem reading m state into integer '//trim(ctemp3))
104  call io_error('param_get_projection: Problem reading zona into real '//trim(ctemp))
105  call io_error('param_get_projection: Problem reading radial state into integer '//trim(ctemp))
106  call io_error('param_get_projection: Problem reading m state into string '//trim(ctemp3))



   end subroutine param_get_projections

  !===================================!
  subroutine param_get_keyword_kpath
    !===================================!
    !                                   !
    !  Fills the kpath data block       !
    !                                   !
    !===================================!
    use w90_io,        only : io_error

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
       read(dummy,*,err=240,end=240) bands_label(counter-1),(bands_spec_points(i,counter-1),i=1,3)&
            ,bands_label(counter),(bands_spec_points(i,counter),i=1,3)
    end do


    in_data(line_s:line_e)(1:maxlen) = ' '

    return


240 call io_error('param_get_keyword_kpath: Problem reading kpath '//trim(dummy))

  end subroutine param_get_keyword_kpath


    !===========================================!
    subroutine param_memory_estimate
    !===========================================!
    !                                           !
    ! Estimate how much memory we will allocate !
    !                                           !
    !===========================================!

    implicit none

    real(kind=dp), parameter :: size_log=1.0_dp
    real(kind=dp), parameter :: size_int=4.0_dp
    real(kind=dp), parameter :: size_real=8.0_dp
    real(kind=dp), parameter :: size_cmplx=16.0_dp
    real(kind=dp) :: mem_wan,mem_wan1,mem_param,mem_dis,mem_dis2,mem_dis1

    mem_param=0    
    mem_dis=0    
    mem_dis1=0    
    mem_dis2=0    
    mem_wan=0    
    mem_wan1=0    


    ! First the data stored in the parameters module
    mem_param= mem_param+num_wann*num_wann*num_kpts*size_cmplx                   !u_matrix
    if(.not.disentanglement) &
         mem_param= mem_param+ num_wann*num_wann*nntot*num_kpts*size_cmplx       !m_matrix
       
    if (disentanglement) then
       mem_param= mem_param+  num_bands*num_wann*num_kpts*size_cmplx             ! u_matrix_opt
    endif

    if(allocated(atoms_species_num))then 
       mem_param= mem_param+(num_species)*size_int                               !atoms_species_num
       mem_param= mem_param+(num_species)*size_real                              !atoms_label
       mem_param= mem_param+(num_species)*size_real                              !atoms_symbol
       mem_param= mem_param+(3*maxval(atoms_species_num)*num_species)*size_real  !atoms_pos_frac
       mem_param= mem_param+(3*maxval(atoms_species_num)*num_species)*size_real  !atoms_pos_cart
    endif

    if(allocated(proj_site))then 
       mem_param= mem_param+ (3*num_proj)*size_real              !proj_site
       mem_param= mem_param+ (num_proj) *size_int                !proj_l
       mem_param= mem_param+ (num_proj)*size_int                 !proj_m
       mem_param= mem_param+ (3*num_proj) *size_real             !proj_z
       mem_param= mem_param+ (3*num_proj) *size_real             !proj_x
       mem_param= mem_param+ (num_proj)*size_real                !proj_radial
       mem_param= mem_param+ (num_proj)*size_real                !proj_zona
    endif

    mem_param=mem_param+num_kpts*nntot*size_int                  !nnlist
    mem_param=mem_param+num_kpts*nntot/2*size_int                !neigh
    mem_param=mem_param+3*num_kpts*nntot*size_int                !nncell
    mem_param=mem_param+nntot*size_real                          !wb
    mem_param=mem_param+3*nntot/2*size_real                      !bka
    mem_param=mem_param+3*nntot*num_kpts*size_real               !bk

    mem_param=mem_param+num_bands*num_kpts*size_real             !eigval
    mem_param=mem_param+3*num_kpts*size_real                     !kpt_cart
    mem_param=mem_param+3*num_kpts*size_real                     !kpt_latt
    if (disentanglement) then
       mem_param=mem_param+num_kpts*size_int                     !ndimwin
       mem_param=mem_param+num_bands*num_kpts*size_log           !lwindow
    endif
    mem_param=mem_param+3*num_wann*size_real                     !wannier_centres
    mem_param=mem_param+num_wann*size_real                       !wannier_spreads

    if (disentanglement) then
       ! Module vars
       mem_dis=mem_dis+num_bands*num_kpts*size_real              !eigval_opt
       mem_dis=mem_dis+num_kpts*size_int                         !nfirstwin
       mem_dis=mem_dis+num_kpts*size_int                         !ndimfroz
       mem_dis=mem_dis+num_bands*num_kpts*size_int               !indxfroz
       mem_dis=mem_dis+num_bands*num_kpts*size_int               !indxnfroz
       mem_dis=mem_dis+num_bands*num_kpts*size_log               !lfrozen

       !the memory high-water wiil occur in dis_extract or when we allocate m_matrix

       mem_dis1=mem_dis1+num_wann*num_bands*size_cmplx              !cwb
       mem_dis1=mem_dis1+num_wann*num_wann*size_cmplx               !cww
       mem_dis1=mem_dis1+num_bands*num_wann*size_cmplx              !cbw
       mem_dis1=mem_dis1+5*num_bands*size_int                       !iwork
       mem_dis1=mem_dis1+num_bands*size_int                         !ifail
       mem_dis1=mem_dis1+num_bands*size_real                        !w
       mem_dis1=mem_dis1+7*num_bands*size_real                      !rwork
       mem_dis1=mem_dis1+(num_bands*(num_bands+1))/2*size_cmplx     !cap
       mem_dis1=mem_dis1+2*num_bands*size_cmplx                     !cwork
       mem_dis1=mem_dis1+num_bands*num_bands*size_cmplx             !cz
       mem_dis1=mem_dis1+num_kpts*size_real                         !wkomegai1
       mem_dis1=mem_dis1+num_bands*num_bands*num_kpts*size_cmplx    !ceamp
       mem_dis1=mem_dis1+num_bands*num_bands*num_kpts*size_cmplx    !cham
       mem_dis2=mem_dis2+num_wann*num_wann*nntot*num_kpts*size_cmplx!m_matrix

       if(index(devel_flag,'memory')>0) then
          mem_dis=mem_dis+mem_dis1
       else
          mem_dis=mem_dis+max(mem_dis1,mem_dis2)
       endif

       mem_dis= mem_dis+num_bands*num_bands*nntot*num_kpts*size_cmplx      ! m_matrix_orig
       mem_dis= mem_dis+num_bands*num_wann*num_kpts*size_cmplx             ! a_matrix

    endif

    !Wannierise

    mem_wan1=mem_wan1+(num_wann*num_wann*nntot*num_kpts)*size_cmplx     !  'm0' 
    if(index(devel_flag,'memory')==0) then
       mem_wan=mem_wan+mem_wan1
    endif
    mem_wan=mem_wan+(num_wann* num_wann* num_kpts)*size_cmplx           !  'u0' 
    mem_wan=mem_wan+(num_wann* nntot* num_kpts)*size_real               !  'rnkb' 
    mem_wan=mem_wan+(num_wann* nntot* num_kpts)*size_real               !  'ln_tmp' 
    mem_wan=mem_wan+(num_wann* nntot* num_kpts)*size_cmplx              !  'csheet' 
    mem_wan=mem_wan+(num_wann* num_wann* num_kpts) *size_cmplx          !  'cdodq' 
    mem_wan=mem_wan+(num_wann* nntot* num_kpts)*size_real               !  'sheet' 
    mem_wan=mem_wan+(3* num_wann)*size_real                             !  'rave' 
    mem_wan=mem_wan+(num_wann)  *size_real                              !  'r2ave' 
    mem_wan=mem_wan+(num_wann) *size_real                               !  'rave2' 
    mem_wan=mem_wan+(3* num_wann) *size_real                            !  'rguide' 
    mem_wan=mem_wan+ (num_wann)*size_cmplx                              !  'cwschur1'
    mem_wan=mem_wan+ (10 * num_wann) *size_cmplx                        !  'cwschur2'
    mem_wan=mem_wan+ (num_wann)*size_cmplx                              !  'cwschur3'
    mem_wan=mem_wan+ (num_wann) *size_cmplx                             !  'cwschur4'
    mem_wan=mem_wan+ (num_wann* num_wann* num_kpts) *size_cmplx         !  'cdq'
    mem_wan=mem_wan+ (num_wann* num_wann)  *size_cmplx                  !  'cz'
    mem_wan=mem_wan+ (num_wann* num_wann) *size_cmplx                   !  'cmtmp'
    mem_wan=mem_wan+ (num_wann* num_wann* num_kpts) *size_cmplx         !  'cdqkeep'
    mem_wan=mem_wan+(num_wann*num_wann)*size_cmplx                      !  'tmp_cdq'
    mem_wan=mem_wan+ (num_wann)*size_real                               !  'evals'
    mem_wan=mem_wan+ (4*num_wann)*size_cmplx                            !  'cwork'
    mem_wan=mem_wan+ (3*num_wann-2)*size_real                           !  'rwork'
    !d_omega
    mem_wan=mem_wan+(num_wann* num_wann) *size_cmplx   !  'cr' 
    mem_wan=mem_wan+(num_wann* num_wann)*size_cmplx    !  'crt' 
    if(disentanglement) &
         mem_wan= mem_wan+ num_wann*num_wann*nntot*num_kpts*size_cmplx       !m_matrix

     write(stdout,'(1x,a)') '*============================================================================*'
     write(stdout,'(1x,a)')  '|                              MEMORY ESTIMATE                               |'
     write(stdout,'(1x,a)')  '|         Maximum RAM allocated during each phase of the calculation         |'
     write(stdout,'(1x,a)')  '*============================================================================*'
     if(disentanglement) &
          write(stdout,'(1x,"|",24x,a15,f16.2,a,18x,"|")') 'Disentanglement:',(mem_param+mem_dis)/(1024**2),' Mb'
     write(stdout,'(1x,"|",24x,a15,f16.2,a,18x,"|")') 'Wannierise:',(mem_param+mem_wan)/(1024**2),' Mb'
     if(index(devel_flag,'memory')==0.and. iprint>1) then
        write(stdout,'(1x,a)')  '|                                                                            |'
        write(stdout,'(1x,a)')  '| N.B. by setting the page file option memory usage will be reduced to:      |'
        write(stdout,'(1x,"|",24x,a15,f16.2,a,18x,"|")') 'Disentanglement:',(mem_param+mem_dis- &
             max(mem_dis1,mem_dis2)+mem_dis1)/(1024**2),' Mb'
        write(stdout,'(1x,"|",24x,a15,f16.2,a,18x,"|")') 'Wannierise:',(mem_param+mem_wan-mem_wan1)/(1024**2),' Mb'
     endif

     write(stdout,'(1x,a)')  '*----------------------------------------------------------------------------*'
     write(stdout,*) ' '


!    if(disentanglement) then
!       write(*,'(a12,f12.4,a)') 'Disentangle',(mem_param+mem_dis)/(1024**2),' Mb'
!    end if
!    write(*,'(a12,f12.4,a)') 'Wannierise ',(mem_wan+mem_param)/(1024**2),' Mb'
!    write(*,'(a12,f12.4,a)') 'Module',(mem_param)/(1024**2),' Mb'

    return
  end subroutine param_memory_estimate


end module w90_parameters
