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

module w90_hamiltonian

  use w90_constants, only : dp

  implicit none

  private

  public :: hamiltonian_get_hr
  public :: hamiltonian_write_hr

contains

 !============================================!
  subroutine hamiltonian_get_hr()
  !============================================!
  !                                            !
  !  Calculate the Hamiltonian in the WF basis !
  !                                            !
  !============================================!

    use w90_constants,  only : cmplx_0,cmplx_i,twopi
    use w90_io,         only : io_error,io_stopwatch,io_file_unit, &
                               stdout, seedname, io_date
    use w90_parameters, only : num_bands,num_kpts,num_wann,u_matrix,eigval,kpt_latt, &
         u_matrix_opt,lwindow,ndimwin,have_disentangled,ham_r,ham_k,irvec,nrpts, &
         have_ham_k,have_ham_r,have_translated,use_translation,hr_plot,timing_level

    implicit none
  
    integer              :: file_unit
    integer, allocatable :: shift_vec(:,:)
    complex(kind=dp)     :: fac
    real(kind=dp)        :: rdotk
    real(kind=dp)        :: eigval_opt(num_bands,num_kpts)
    real(kind=dp)        :: eigval2(num_wann,num_kpts)
    real(kind=dp)        :: irvec_tmp(3)
    integer              :: loop_kpt,i,j,m,loop_rpt,ierr,counter
    character (len=33)   :: header
    character (len=9)    :: cdate,ctime

    if (timing_level>1) call io_stopwatch('hamiltonian: get_hr',1)

    if(have_ham_r) then 
      if (have_translated .eqv. use_translation) then
         goto 200
      else 
         goto 100
      endif
    end if

    if(have_ham_k) go to 100

    if (.not. allocated(ham_k)) then
       allocate(ham_k(num_wann,num_wann,num_kpts),stat=ierr)
       if (ierr/=0) call io_error('Error in allocating ham_k in hamiltonian_get_hr')
    end if

    ham_k=cmplx_0
    eigval_opt=0.0_dp
    eigval2=0.0_dp

    if(have_disentangled) then

       ! slim down eigval to contain states within the outer window

       do loop_kpt=1,num_kpts
          counter=0
          do j=1,num_bands
             if(lwindow(j,loop_kpt)) then
                counter=counter+1
                eigval_opt(counter,loop_kpt)=eigval(j,loop_kpt)
             end if
          end do
       end do
       
       ! rotate eigval into the optimal subspace
       ! in general eigval would be a matrix at each kpoints
       ! but we choose u_matrix_opt such that the Hamiltonian is
       ! diagonal at each kpoint. (I guess we should check it here)
       
       do loop_kpt=1,num_kpts
          do j=1,num_wann
             do m=1,ndimwin(loop_kpt)
                eigval2(j,loop_kpt)=eigval2(j,loop_kpt)+eigval_opt(m,loop_kpt)* &
                     real(conjg(u_matrix_opt(m,j,loop_kpt))*u_matrix_opt(m,j,loop_kpt),dp)
             enddo
          enddo
       enddo

    else
       eigval2(1:num_wann,:)=eigval(1:num_wann,:)
    end if


    ! At this point eigval2 contains num_wann values which belong to the wannier subspace.


    ! Rotate Hamiltonian into the basis of smooth bloch states 
    !          H(k)=U^{dagger}(k).H_0(k).U(k)
    ! Note: we enforce hermiticity here

    do loop_kpt=1,num_kpts
       do j=1,num_wann
          do i=1,j
             do m=1,num_wann
                ham_k(i,j,loop_kpt)=ham_k(i,j,loop_kpt)+eigval2(m,loop_kpt)* &
                     conjg(u_matrix(m,i,loop_kpt))*u_matrix(m,j,loop_kpt)
             enddo
             if(i.lt.j) ham_k(j,i,loop_kpt)=conjg(ham_k(i,j,loop_kpt))
          enddo
       enddo
    enddo

    have_ham_k = .true.

100 continue

    ! Fourier transform rotated hamiltonian into WF basis
    ! H_ij(k) --> H_ij(R) = (1/N_kpts) sum_k e^{-ikR} H_ij(k)
    if (.not.allocated(ham_r)) then
      allocate(ham_r(num_wann,num_wann,nrpts),stat=ierr)
      if (ierr/=0) call io_error('Error in allocating ham_r in hamiltonian_get_hr')
    end if
    
    ham_r=cmplx_0

    if (.not. use_translation) then

       do loop_rpt=1,nrpts
          do loop_kpt=1,num_kpts
             rdotk=twopi*dot_product(kpt_latt(:,loop_kpt),real(irvec(:,loop_rpt),dp))
             fac=exp(-cmplx_i*rdotk)/real(num_kpts,dp)
             ham_r(:,:,loop_rpt)=ham_r(:,:,loop_rpt)+fac*ham_k(:,:,loop_kpt)
          enddo
       enddo
    
       have_translated = .false.

    else

       allocate(shift_vec(3,num_wann),stat=ierr)
       if (ierr/=0) call io_error('Error in allocating shift_vec in hamiltonian_get_hr')
       call internal_translate_wannier_centres()

       do loop_rpt=1,nrpts
          do loop_kpt=1,num_kpts
             do i=1,num_wann
                do j=1,num_wann
                   ! ham_r(j,i,loop_rpt)
                   ! interaction btw j at 0 and i at irvec(:,loop_rpt)
                   irvec_tmp(:)=irvec(:,loop_rpt)+shift_vec(:,i)-shift_vec(:,j)   
                   rdotk=twopi*dot_product(kpt_latt(:,loop_kpt),real(irvec_tmp(:),dp))
                   fac=exp(-cmplx_i*rdotk)/real(num_kpts,dp)
                   ham_r(j,i,loop_rpt)=ham_r(j,i,loop_rpt)+fac*ham_k(j,i,loop_kpt)
                end do
             end do
          enddo
       enddo

       have_translated = .true.

    end if

    have_ham_r = .true.

200 continue

    if (allocated(shift_vec)) then
       deallocate(shift_vec,stat=ierr)
       if (ierr/=0) call io_error('Error in deallocating shift_vec in hamiltonian_get_hr')
    end if

    if (timing_level>1) call io_stopwatch('hamiltonian: get_hr',2)

    return

    contains

    !====================================================!
      subroutine internal_translate_wannier_centres()
    !====================================================!

      use w90_parameters, only : num_wann,real_lattice,recip_lattice,wannier_centres, &
                                 num_atoms,atoms_pos_cart,translation_centre_frac, &
                                 automatic_translation,wannier_centres_translated, &
                                 num_species,atoms_species_num,lenconfac
      use w90_io,         only : stdout,io_error
      use w90_utility,    only : utility_cart_to_frac,utility_frac_to_cart

      implicit none

      ! <<<local variables>>>
      integer :: iw,loop,ierr,nat,nsp,ind
      real(kind=dp), allocatable :: r_home(:,:),r_frac(:,:)
      real(kind=dp) :: c_pos_cart(3), c_pos_frac(3)
      real(kind=dp) :: r_frac_min(3)

      if (.not.allocated(wannier_centres_translated)) then
         allocate(wannier_centres_translated(3,num_wann),stat=ierr)
         if (ierr/=0) call io_error('Error in allocating wannier_centres_translated &
              &in internal_translate_wannier_centres')
      end if
      
      allocate(r_home(3,num_wann),stat=ierr)
      if (ierr/=0) call io_error('Error in allocating r_home in internal_translate_wannier_centres')
      allocate(r_frac(3,num_wann),stat=ierr)
      if (ierr/=0) call io_error('Error in allocating r_frac in internal_translate_wannier_centres')
      r_home=0.0_dp;r_frac=0.0_dp

      if (automatic_translation) then
         ! Calculate centre of atomic positions
         c_pos_cart=0.0_dp;c_pos_frac=0.0_dp
         do nsp=1,num_species
            do nat=1,atoms_species_num(nsp)
               c_pos_cart(:) = c_pos_cart(:) + atoms_pos_cart(:,nat,nsp)
            enddo
         enddo
!!@         do loop=1,num_atoms
!!@            c_pos_cart(:) = c_pos_cart(:) + atoms_pos_cart(:,loop)
!!@         end do
         c_pos_cart = c_pos_cart / num_atoms
         ! Cartesian --> fractional
         call utility_cart_to_frac(c_pos_cart,c_pos_frac,recip_lattice)      
         ! Wannier function centres will be in [c_pos_frac-0.5,c_pos_frac+0.5]
         r_frac_min(:)=c_pos_frac(:)-0.5_dp
      else
         r_frac_min(:)=translation_centre_frac(:)-0.5_dp
      end if
       
      ! Cartesian --> fractional
      do iw=1,num_wann
         call utility_cart_to_frac(wannier_centres(:,iw),r_frac(:,iw),recip_lattice)
         ! Rationalise r_frac - r_frac_min to interval [0,1]
         !  by applying shift of -floor(r_frac - r_frac_min)
         shift_vec(:,iw)=-floor(r_frac(:,iw)-r_frac_min(:))
         r_frac(:,iw)=r_frac(:,iw)+real(shift_vec(:,iw),kind=dp)
         ! Fractional --> Cartesian
         call utility_frac_to_cart(r_frac(:,iw),r_home(:,iw),real_lattice)
      end do

      ! NEVER overwrite wannier_centres
      !wannier_centres = r_home

      write(stdout,'(1x,a)') 'Translated centres (translated to the center of atomic position)'
      do iw=1,num_wann
         write(stdout,888) iw,(r_home(ind,iw)*lenconfac,ind=1,3)
      end do
      write(stdout,'(1x,a78)') repeat('-',78)
      write(stdout,*)

      wannier_centres_translated = r_home

      deallocate(r_frac,stat=ierr)
      if (ierr/=0) call io_error('Error in deallocating r_frac in internal_translate_wannier_centres')
      deallocate(r_home,stat=ierr)
      if (ierr/=0) call io_error('Error in deallocating r_home in internal_translate_wannier_centres')

      return

888   format(2x,'WF centre and spread',i5,2x,'(',f10.6,',',f10.6,',',f10.6,' )',f15.8)

    end subroutine internal_translate_wannier_centres

  end subroutine hamiltonian_get_hr


  !============================================!
  subroutine hamiltonian_write_hr()
  !============================================!
  !  Write the Hamiltonian in the WF basis     !
  !============================================!

    use w90_io,        only : io_error,io_stopwatch,io_file_unit, &
                              stdout,seedname,io_date
    use w90_parameters, only : num_wann,timing_level,ham_r,nrpts,irvec,hr_written

    integer            :: i,j,loop_rpt,file_unit
    character (len=33) :: header
    character (len=9)  :: cdate,ctime

    if (hr_written) return

    if (timing_level>1) call io_stopwatch('hamiltonian: write_hr',1)

    ! write the  whole matrix with all the indices 
 
    file_unit=io_file_unit()
    open(file_unit,file=trim(seedname)//'_hr.dat',form='formatted',status='unknown',err=101)

    call io_date(cdate,ctime)
    header='written on '//cdate//' at '//ctime

    write(file_unit,*) header ! Date and time

    do loop_rpt=1,nrpts
       do i=1,num_wann
          do j=1,num_wann
             write( file_unit,'(5I5,2F12.6)') irvec(:,loop_rpt), j, i, ham_r(j,i,loop_rpt)
          end do
       end do
    end do

    close(file_unit)

    hr_written=.true.

    if (timing_level>1) call io_stopwatch('hamiltonian: write_hr',2)

    return

101 call io_error('Error: hamiltonian_write_hr: problem opening file '//trim(seedname)//'_hr.dat')

  end subroutine hamiltonian_write_hr
  

end module w90_hamiltonian
