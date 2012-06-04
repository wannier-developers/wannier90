!-*- mode: F90; mode: font-lock -*-!

program postw90

  use w90_constants, only : dp
  use w90_parameters
  use w90_io

  use w90_kmesh
  use w90_comms, only : on_root,num_nodes, comms_setup, comms_end
  use w90_wanint_common

  ! These modules deal with the interpolation of specific physical properties
  !
  use w90_dos_wanint
  use w90_berry_wanint
  use w90_spin_wanint
  use w90_kpath_plot
  use w90_slice_plot

  use w90_boltzwann
  use w90_geninterp

  implicit none

  integer       :: nkp
  logical       :: have_gamma
  real(kind=dp) :: time0,time1,time2

  ! Put some descriptive comments here
  !
  call comms_setup

  library = .false.
  ispostw90 = .true.
  
  if(on_root) then
     time0=io_time()
     call io_get_seedname
     stdout=io_file_unit()
     open(unit=stdout,file=trim(seedname)//'.wpout')
     call param_write_header
     if(num_nodes==1) then
#ifdef MPI        
        write(stdout,'(/,1x,a)') 'Running in serial (with parallel executable)'
#else
        write(stdout,'(/,1x,a)') 'Running in serial (with serial executable)'
#endif
     else
        write(stdout,'(/,1x,a,i3,a/)')&
             'Running in parallel on ',num_nodes,' CPUs'
     endif
  end if

  ! Read onto the root node all the input parameters from seendame.wanint, 
  ! as well as the energy eigenvalues on the ab-initio q-mesh from seedname.eig 
  !
  if(on_root) then 
       call param_read
       time1=io_time()
       write(stdout,'(1x,a25,f11.3,a)')&
            'Time to read parameters  ',time1-time0,' (sec)'

       ! Check if the q-mesh includes the gamma point
       !
       have_gamma=.false.
       do nkp=1,num_kpts
          if (all(kpt_latt(:,nkp)<0.000001_dp)) have_gamma=.true.       
       end do
       if(.not. have_gamma) write(stdout,'(1x,a)')&
           'Ab-initio does not include Gamma. Interpolation may be incorrect!!!'
       !
       ! Need nntot, wb, and bk to evaluate WF matrix elements of the position
       ! operator in reciprocal space. Also need nnlist to compute the
       ! additional matrix elements entering the orbital magnetization
       !
       call kmesh_get
       time2=io_time()
       write(stdout,'(1x,a25,f11.3,a)')&
            'Time to get kmesh        ',time2-time1,' (sec)'

       ! GP, May 10, 2012: for the moment I leave this commented 
       ! since we need first to tune that routine so that it doesn't
       ! print the memory information related to wannier90.x.
       ! Note that the code for the memory estimation for the
       ! Boltzwann routine is already there.
       !       call param_memory_estimate
  end if

  ! We now distribute a subset of the parameters to the other nodes
  !
  call wanint_param_dist

  ! Read files seedname.chk (overlap matrices, unitary matrices for both
  ! disentanglement and maximal localization, etc.)
  !
  if(on_root) then 
     call param_read_chkpt()
  end if

  ! Distribute the information in the um and chk files to the other nodes
  ! Ivo: For interpolation purposes we do not need u_matrix_opt and u_matrix
  !      separately, only their product v_matrix, and this is what is 
  !      distributed now
  !
  call wanint_data_dist

  ! Read list of k-points in irreducible BZ and their weights
  !
  ! Should this be done on root node only?
  !
  if(wanint_kpoint_file) call wanint_get_kpoint_file

  ! Setup a number of common variables for all interpolation tasks
  !
  call wanint_setup

  if(on_root) then
     time1=io_time()
     write(stdout,'(/1x,a25,f11.3,a)')&
          'Time to read and process .chk    ',time1-time2,' (sec)'
  endif

  ! Now perform one or more of the following tasks

  ! ---------------------------------------------------------------
  ! Density of states calculated using a uniform interpolation mesh
  ! ---------------------------------------------------------------
  !
  if(do_dos .and. index(dos_task,'dos_plot')>0) call dos

  if(do_dos .and. index(dos_task,'find_fermi_energy')>0) call find_fermi_level

  ! --------------------------------------------------------------------
  ! Bands, Berry curvature, or orbital magnetization plot along a k-path
  ! --------------------------------------------------------------------
  !
  if(kpath_plot) call k_path

  ! ---------------------------------------------------------------------------
  ! Bands, Berry curvature, or orbital magnetization plot on a slice in k-space
  ! ---------------------------------------------------------------------------
  !
  if(slice_plot) call k_slice

  ! --------------------
  ! Spin magnetic moment
  ! --------------------
  !
  if(evaluate_spin_moment) call spin_moment

  ! -------------------------------------------------------------------
  ! dc Anomalous Hall conductivity and eventually (if 'mcd' string also 
  ! present in addition to 'ahe', e.g., 'ahe+mcd') dichroic optical
  ! conductivity, both calculated on the same (adaptively-refined) mesh
  ! -------------------------------------------------------------------
  !
  ! ---------------------------------------------------------------
  ! Absorptive dichroic optical conductivity & JDOS on uniform mesh
  ! ---------------------------------------------------------------
  !
  ! -----------------------------------------------------------------
  ! Absorptive ordinary optical conductivity & JDOS on a uniform mesh
  ! -----------------------------------------------------------------
  !
  ! -----------------------------------------------------------------
  ! Orbital magnetization
  ! -----------------------------------------------------------------
  !
  if(optics_plot) call berry
  ! -----------------------------------------------------------------
  ! Boltzmann transport coefficients (BoltzWann module)
  ! -----------------------------------------------------------------
  !
  if(on_root) then
     time1=io_time()
  endif

  if(geninterp) call geninterp_main

  if(boltzwann) call boltzwann_main

  if(on_root.and.boltzwann) then
     time2=io_time()
     write(stdout,'(/1x,a,f11.3,a)')&
          'Time for BoltzWann (Boltzmann transport) ',time2-time1,' (sec)'
  endif

  ! I put a barrier here before calling the final time printing routines,
  ! just to be sure that all processes have arrived here.
  call comms_barrier
  if(on_root) then
     write(stdout,'(/,1x,a25,f11.3,a)')&
          'Total Execution Time     ',io_time(),' (sec)'
     if (timing_level>0) call io_print_timings()
     write(stdout,*) 
     write(stdout,'(/,1x,a)') 'All done: postw90 exiting'
     close(stdout)
  end if

  call comms_end

end program postw90
  

