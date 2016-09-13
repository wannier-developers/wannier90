!-*- mode: F90; mode: font-lock; column-number-mode: true -*-!
!                                                            !
! This module calculates and prints r_ab for a and b Wannier !
! functions.                                                 !
!                                                            !
! The module is heavily based on w90_hamiltonian substituting!
! the hamiltonian calculation by                             !
! Eq. C16 of Marzari and Vanderbilt PRB 56, 12847 (1997)     !
! This is called from plot.F90 to get the file at the end of !
! the run if wannier90_hr.dat is requested.                  !
!                                                            !
! To install this module it is necessary to add the line:    !
!     call position_setup()                                  !
! to plot.F90 right after                                    !
!     if(hr_plot) call hamiltonian_write_hr()                !
! And force the compilation by adding the lines:             !
! position.o: position.F90 constants.o io.o utility.o parameters.o
!        $(F90) $(FCOPTS) -c position.F90
! in Makefile after the definition of hamiltonian.o          !
!                                                            !
!------------------------------------------------------------!

module w90_position

  use w90_constants, only : dp

  implicit none

  private

  public :: position_setup

contains

  !============================================!
  subroutine position_setup()
    !============================================!
    use w90_parameters, only : m_matrix, wb, bk, num_wann, num_kpts, kpt_latt,&
                               nntot
    use w90_hamiltonian, only: irvec, nrpts
    use w90_constants, only  : twopi, cmplx_i
    use w90_io,         only : io_error, io_file_unit, seedname,io_date

    implicit none

    complex(kind=dp)     :: fac
    real(kind=dp)        :: rdotk
    real(kind=dp)        :: delta
    integer              :: loop_rpt, m, n, nkp, ind, nn, file_unit
    complex(kind=dp)     :: position(3)
    character (len=33) :: header
    character (len=9)  :: cdate,ctime

    file_unit=io_file_unit()
    open(file_unit,file=trim(seedname)//'_pos.dat',form='formatted',status='unknown',err=101)
    call io_date(cdate,ctime)

    header='written on '//cdate//' at '//ctime
    write(file_unit,*) header ! Date and time
    write(file_unit,*) num_wann

    do loop_rpt=1,nrpts
       do m=1,num_wann
          do n=1,num_wann
             delta=0._dp
             if (m.eq.n) delta=1._dp
             position(:)=0._dp
             do nkp=1,num_kpts
                rdotk=twopi*dot_product(kpt_latt(:,nkp),real(irvec(:,loop_rpt),dp))
                fac=exp(-cmplx_i*rdotk)/real(num_kpts,dp)
                do ind = 1, 3  
                   do nn = 1, nntot 
                      ! Eq. C16 of Marzari and Vanderbilt PRB 56, 12847 (1997)     !
                      position(ind) = position(ind) + &
                                      wb(nn) * bk(ind,nn,nkp) * (m_matrix(n,m,nn,nkp) - delta) * fac
                   end do   
                end do
             end do
             write( file_unit ,'(5I5,6F12.6)') irvec(:,loop_rpt),n,m,cmplx_i*position(:)  
          end do
       end do
    end do
    
    close(file_unit)

    return

101 call io_error('Error: position_write_pos: problem opening file '//trim(seedname)//'_pos.dat')

  end subroutine position_setup
  
end module w90_position
