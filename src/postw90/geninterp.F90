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
!============================================================!
!                                                            !
! Generic Interpolation Routine written by Giovanni Pizzi    !
! THEOS, EPFL, Station 12, 1015 Lausanne (Switzerland)       !
! June, 2012                                                 !
!                                                            !
!============================================================!

module w90_geninterp

  use w90_constants
  use w90_parameters, only : geninterp_alsofirstder, num_wann, recip_lattice, real_lattice
  use w90_io, only         : io_error,stdout,io_stopwatch,io_file_unit,seedname  
  use w90_get_oper, only      : get_HH_R, HH_R
  use w90_comms
  use w90_utility, only       : utility_diagonalize
  use w90_wanint_common, only : fourier_R_to_k
  use w90_wan_ham, only       : get_deleig_a
  implicit none

  private 
  public :: geninterp_main

contains 

  subroutine geninterp_main()
    integer :: kpt_unit, outdat_unit, num_kpts, ierr, i, j, enidx
    character(len=500) :: commentline
    character(len=50) :: cdum
    integer, dimension(:), allocatable :: kpointidx
    real(kind=dp), dimension(:,:), allocatable :: kpoints
    complex(kind=dp), dimension(:,:), allocatable :: HH
    complex(kind=dp), dimension(:,:), allocatable :: UU
    complex(kind=dp), dimension(:,:,:), allocatable :: delHH
    real(kind=dp), dimension(3) :: kpt, frac
    real(kind=dp) :: del_eig(num_wann,3)
    real(kind=dp) :: eig(num_wann), levelspacing_k(num_wann)
    logical :: absoluteCoords

    kpt_unit=io_file_unit()
    open(unit=kpt_unit,file=trim(seedname)//'_geninterp.kpt',form='formatted',status='old',err=105)

    ! First line: comment (e.g. creation date, author, ...)
    read(kpt_unit,'(A500)',err=106,end=106) commentline
    read(kpt_unit,*,err=106,end=106) cdum

    if (index(cdum,'crystal')>0) then
       absoluteCoords = .false.
    elseif (index(cdum,'rel')>0) then
       absoluteCoords = .false.       
    elseif (index(cdum,'frac')>0) then
       absoluteCoords = .true.
    elseif (index(cdum,'abs')>0) then
       absoluteCoords = .true.
    else
       call io_error('Error on second line of file '//trim(seedname)//'_geninterp.kpt: ' // &
            'unable to recognize keyword')
    end if

    ! Third line: number of following kpoints
    read(kpt_unit,*,err=106,end=106) num_kpts

    allocate(kpointidx(num_kpts),stat=ierr)
    if (ierr/=0) call io_error('Error allocating kpointidx in geinterp_main.')

    allocate(kpoints(3,num_kpts),stat=ierr)
    if (ierr/=0) call io_error('Error allocating kpoints in geinterp_main.')

    ! Lines with integer identifier and three coordinates
    ! (in crystallographic coordinates relative to the reciprocal lattice vectors)
    do i=1,num_kpts
       read(kpt_unit,*,err=106,end=106) kpointidx(i), kpt
       ! Internally, I need the relative (fractional) coordinates in units of the reciprocal-lattice vectors
       if (absoluteCoords.eqv..false.) then
          kpoints(:,i) = kpt
       else
          kpoints(:,i) = 0._dp
          ! I use the real_lattice (transposed) and a factor of 2pi instead of inverting again recip_lattice
          do j=1,3
             kpoints(j,i)=real_lattice(j,1)*kpt(1) + real_lattice(j,2)*kpt(2) + real_lattice(j,3)*kpt(3) 
          end do
          kpoints(:,i) = kpoints(:,i) / (2._dp * pi)
       end if
    end do

    close(kpt_unit)

    allocate(HH(num_wann,num_wann),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating HH in calcTDF')
    allocate(UU(num_wann,num_wann),stat=ierr)
    if (ierr/=0) call io_error('Error in allocating UU in calcTDF')    
    if (geninterp_alsofirstder) then
       allocate(delHH(num_wann,num_wann,3),stat=ierr)
       if (ierr/=0) call io_error('Error in allocating delHH in calcTDF')
    end if
    
    ! I call once the routine to calculate the Hamiltonian in real-space <0n|H|Rm>
    call get_HH_R

    ! outuput band structure
    outdat_unit=io_file_unit()
    open(unit=outdat_unit,file=trim(seedname)//'_geninterp.dat',form='formatted',err=107)

    ! I rewrite the comment line on the output
    write(outdat_unit, '(A)') "# Input file comment: " // trim(commentline)

    if (geninterp_alsofirstder) then
       write(outdat_unit, '(A)') "#  Kpt_idx   K_x (2pi/ang)     K_y (2pi/ang)     K_z (2pi/ang)     Energy (eV)" // &
            "      EnergyDer_x       EnergyDer_y       EnergyDer_z"
    else
       write(outdat_unit, '(A)') "#  Kpt_idx   K_x (2pi/ang)     K_y (2pi/ang)     K_z (2pi/ang)     Energy (eV)"
    end if

    do i=1, num_kpts
       kpt = kpoints(:,i)
       frac = 0._dp
       do j=1,3
          frac(j)=recip_lattice(1,j)*kpt(1) + recip_lattice(2,j)*kpt(2) + recip_lattice(3,j)*kpt(3) 
       end do
       ! Here I get the band energies and the velocities (if required)
       call fourier_R_to_k(kpt,HH_R,HH,0) 
       call utility_diagonalize(HH,num_wann,eig,UU) 
       if (geninterp_alsofirstder) then
          call fourier_R_to_k(kpt,HH_R,delHH(:,:,1),1) 
          call fourier_R_to_k(kpt,HH_R,delHH(:,:,2),2) 
          call fourier_R_to_k(kpt,HH_R,delHH(:,:,3),3) 
          call get_deleig_a(del_eig(:,1),eig,delHH(:,:,1),UU)
          call get_deleig_a(del_eig(:,2),eig,delHH(:,:,2),UU)
          call get_deleig_a(del_eig(:,3),eig,delHH(:,:,3),UU)
          do enidx=1,num_wann
             write(outdat_unit, '(I10,7G18.10)') kpointidx(i), frac, eig(enidx), del_eig(enidx,:)
          end do
       else
          do enidx=1,num_wann
             write(outdat_unit, '(I10,4G18.10)') kpointidx(i), frac, eig(enidx)
          end do
       end if
    end do

    close(outdat_unit)

    ! Final deallocations
    if (allocated(kpointidx)) deallocate(kpointidx)
    if (allocated(kpoints)) deallocate(kpoints)
    if (allocated(HH)) deallocate(HH)
    if (allocated(UU)) deallocate(UU)
    if (allocated(delHH)) deallocate(delHH)

    return

105 call io_error('Error: Problem opening k-point file '//trim(seedname)//'_geninterp.kpt')
106 call io_error('Error: Problem reading k-point file '//trim(seedname)//'_geninterp.kpt')
107 call io_error('Error: Problem opening output file '//trim(seedname)//'_geninterp.dat')
  end subroutine geninterp_main

end module w90_geninterp
