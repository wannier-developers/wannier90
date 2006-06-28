!-*- mode: F90; mode: font-lock; column-number-mode: true -*-!
!                                                            !
! Copyright (C) 2004,2006 Jonathan Yates, Arash Mostofi,     !
!            Nicola Marzari, Ivo Souza, David Vanderbilt     !
!                                                            !
! This file is distributed under the terms of the GNU        !
! General Public License. See the file `LICENSE' in          !
! the root directory of the present distribution, or         !
! http://www.gnu.org/copyleft/gpl.txt .                      !
!                                                            !
!------------------------------------------------------------!

subroutine wannier_setup(mp_grid_loc,num_kpts_loc,real_lattice_loc,kpt_latt_loc, &
     num_bands_loc,num_atoms_loc,atom_symbols_loc,atoms_cart_loc,num_nnmax_loc, &
     nntot_loc,nnlist_loc,nncell_loc, &
     proj_site_loc,proj_l_loc,proj_m_loc,proj_radial_loc,proj_z_loc,proj_x_loc,proj_zona_loc, &
     proj_box_loc,exclude_bands_loc)

  use constants
  use parameters
  use io
  use kmesh
  use disentangle
  use overlap
  use wannierise
  use plot
 
  implicit none

  integer, dimension(3), intent(in) :: mp_grid_loc
  integer, intent(in) :: num_kpts_loc
  real(kind=dp), dimension(3,3), intent(in) :: real_lattice_loc
  real(kind=dp), dimension(3,num_kpts), intent(in) :: kpt_latt_loc
  integer, intent(in) :: num_bands_loc
  integer, intent(in) :: num_atoms_loc
  character(len=20), dimension(num_atoms_loc), intent(in) :: atom_symbols_loc
  real(kind=dp), dimension(3,num_atoms_loc), intent(in) :: atoms_cart_loc
  integer, intent(in) :: num_nnmax_loc
  integer, intent(out) :: nntot_loc
  integer, dimension(num_kpts_loc,num_nnmax_loc), intent(out) :: nnlist_loc
  integer,dimension(3,num_kpts_loc,num_nnmax_loc), intent(out) :: nncell_loc
  real(kind=dp), dimension(3,num_bands_loc), intent(out) :: proj_site_loc
  integer, dimension(num_bands_loc), intent(out) :: proj_l_loc
  integer, dimension(num_bands_loc), intent(out) :: proj_m_loc
  integer, dimension(num_bands_loc), intent(out) :: proj_radial_loc
  real(kind=dp), dimension(3,num_bands_loc), intent(out) :: proj_z_loc
  real(kind=dp), dimension(3,num_bands_loc), intent(out) :: proj_x_loc
  real(kind=dp), dimension(num_bands_loc), intent(out) :: proj_zona_loc
  real(kind=dp), dimension(num_bands_loc), intent(out) :: proj_box_loc
  integer, dimension(num_bands_loc), intent(out) :: exclude_bands_loc

  real(kind=dp) time0,time1,time2
  character(len=9) :: stat,pos,cdate,ctime
  integer :: ierr
  logical :: wout_found

  time0=io_time()

  library=.true.
  seedname="wannier"

  stdout=io_file_unit()
  open(unit=stdout,file=trim(seedname)//'.wout')
  call param_write_header

  ! copy local data into module variables
  num_bands=num_bands_loc
  mp_grid=mp_grid_loc
  num_kpts=num_kpts_loc
  real_lattice=real_lattice_loc
  allocate ( kpt_latt(3,num_kpts) ,stat=ierr)
  if (ierr/=0) call io_error('Error allocating kpt_latt in wannier_setup')
  kpt_latt=kpt_latt_loc
  num_atoms=num_atoms_loc
  call param_lib_set_atoms(atom_symbols_loc,atoms_cart_loc)
  
  call param_read
  call param_write

  time1=io_time()
  write(stdout,'(1x,a25,f11.3,a)') 'Time to read parameters  ',time1-time0,' (sec)'

  call io_stopwatch('kmesh_get',1)
  call kmesh_get
  call io_stopwatch('kmesh_get',2)


  ! Now we zero all of the local output data, then copy in the data
  ! from the parameters module

  nntot_loc       = 0
  nnlist_loc      = 0   
  nncell_loc      = 0     
  proj_site_loc      = 0.0_dp   
  proj_l_loc           = 0
  proj_m_loc          = 0
  proj_z_loc      = 0.0_dp
  proj_x_loc      = 0.0_dp
  proj_radial_loc           = 0
  proj_zona_loc        = 0.0_dp
  proj_box_loc    =0.0_dp
  exclude_bands_loc = 0
  

  nntot_loc       = nntot        
  nnlist_loc(:,1:nntot)   =  nnlist_loc(:,1:nntot)       
  nncell_loc(:,:,1:nntot) =  nncell(:,:,1:nntot)       
  if(allocated(proj_site)) then
     proj_site_loc(:,1:nntot)      = proj_site(:,1:nntot)    
     proj_l_loc(1:num_wann)        = proj_l(1:num_wann)          
     proj_m_loc(1:num_wann)        = proj_m(1:num_wann)           
     proj_z_loc(:,1:num_wann)      = proj_z(:,1:num_wann)     
     proj_x_loc(:,1:num_wann)      = proj_x(:,1:num_wann)       
     proj_radial_loc(1:num_wann)   = proj_radial(1:num_wann)            
     proj_zona_loc(1:num_wann)     = proj_zona(1:num_wann)         
     proj_box_loc(1:num_wann)      = proj_box(1:num_wann)   
  end if
  if(allocated(exclude_bands)) then
     exclude_bands_loc(1:num_exclude_bands) = exclude_bands(1:num_exclude_bands)
  end if

  call kmesh_dealloc
  call param_dealloc
  write(stdout,'(1x,a25,f11.3,a)') 'Time to write kmesh      ',io_time(),' (sec)'
  close(stdout)

  
end subroutine wannier_setup


subroutine wannier_run(mp_grid_loc,num_kpts_loc,real_lattice_loc,kpt_latt_loc, &
     num_bands_loc,num_atoms_loc,atom_symbols_loc, &
     atoms_cart_loc,M_matrix_orig_loc, &
     A_matrix_loc,eigenvalues_loc,U_matrix_loc,U_matrix_opt_loc,lwindow_loc, &
     wann_centres_loc,wann_spreads_loc,spread_loc)


  use constants
  use parameters
  use io
  use kmesh
  use disentangle
  use overlap
  use wannierise
  use plot

  implicit none

  integer, dimension(3), intent(in) :: mp_grid_loc
  integer, intent(in) :: num_kpts_loc
  real(kind=dp), dimension(3,3), intent(in) :: real_lattice_loc
  real(kind=dp), dimension(3,num_kpts), intent(in) :: kpt_latt_loc
  integer, intent(in) :: num_bands_loc
  integer, intent(in) :: num_atoms_loc
  character(len=20), dimension(num_atoms_loc), intent(in) :: atom_symbols_loc
  real(kind=dp), dimension(3,num_atoms_loc), intent(in) :: atoms_cart_loc
  complex(kind=dp), dimension(num_bands,num_bands,num_kpts,nntot), intent(in) :: M_matrix_orig_loc
  complex(kind=dp), dimension(num_bands,num_wann,num_kpts), intent(in) :: A_matrix_loc
  real(kind=dp), dimension(num_bands,num_kpts), intent(in) :: eigenvalues_loc
  complex(kind=dp), dimension(num_wann,num_wann,num_kpts), intent(out) :: U_matrix_loc
  complex(kind=dp), dimension(num_bands,num_wann,num_kpts), optional, intent(out) :: U_matrix_opt_loc
  logical, dimension(num_bands,num_kpts), optional, intent(out) :: lwindow_loc
  real(kind=dp), dimension(3,num_wann), optional, intent(out) :: wann_centres_loc
  real(kind=dp), dimension(num_wann), optional, intent(out) :: wann_spreads_loc
  real(kind=dp), dimension(3), optional, intent(out) :: spread_loc

  real(kind=dp) time0,time1,time2
  character(len=9) :: stat,pos,cdate,ctime
  integer :: ierr
  logical :: wout_found

  time0=io_time()

  seedname="wannier"

  write(*,*) seedname  

  ! copy local data into module variables
  num_bands=num_bands_loc
  mp_grid=mp_grid_loc
  num_kpts=num_kpts_loc
  real_lattice=real_lattice_loc
  allocate ( kpt_latt(3,num_kpts) ,stat=ierr)
  if (ierr/=0) call io_error('Error allocating kpt_latt in wannier_setup')
  kpt_latt=kpt_latt_loc
  num_atoms=num_atoms_loc
  call param_lib_set_atoms(atom_symbols_loc,atoms_cart_loc)
  


  call param_read
  call param_write


end subroutine wannier_run
