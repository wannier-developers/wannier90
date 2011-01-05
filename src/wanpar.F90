!-*- mode: F90; mode: font-lock -*-!

program wannier_parint

  use w90_constants, only : dp
  use w90_parameters
  use w90_io

  use w90_kmesh
  use w90_comms
  use w90_wanint_common



  implicit none
 
  integer       :: ierr,nkp
  logical       :: have_gamma

  ! Put some descriptive comments here
  !
  call comms_setup

  if(on_root) then
     call io_get_seedname
     stdout=io_file_unit()
     open(unit=stdout,file=trim(seedname)//'.wpout')
     call param_write_header
  end if

  ! Read onto the root node all the input parameters from seendame.win, 
  ! as well as the energy eigenvalues on the ab-initio q-mesh from seedname.eig 
  !
  if(on_root) then 
       call param_read
       !
       ! Need nntot, wb, and bk to evaluate WF matrix elements of the position
       ! operator in reciprocal space. Also need nnlist to compute the
       ! additional matrix elements entering the orbital magnetization
       !
       call kmesh_get
  end if

  ! We now distribute a subset of the parameters to the other nodes
  ! The subset is defined in the USE statement
  !
  call wanint_param_dist

  ! Read files seedname.chk and seedname_um.dat (overlap matrices, unitary
  ! matrices for both disentanglement and maximal localization, etc.)
  !
  !
  if(on_root) then 
     call param_read_chkpt
  end if

  ! Distribute the information in the chk file to the other nodes
  !
  call wanint_data_dist

  ! Read list of k-points in irreducible BZ and their weights
  !
  ! Should this be done on root node only?
  !
!  if(wanint_kpoint_file) call wanint_get_kpoint_file


!  All interesting stuff happens here
!  call tran_main

  if(on_root) then
     write(stdout,'(/,1x,a)') 'All done: wanint exiting'
     close(stdout)
  end if

  call comms_end

100 format(1x,'Execution started on ',a9,' at ',a9)
101 format(1x,'Reading parameters')

end program wannier_parint
  

