!-*- mode: F90; mode: font-lock -*-!

module w90_wan_ham

  use w90_constants, only : dp

  implicit none

  contains


  subroutine get_D_h_a(delHH_a,UU,eig,ef,D_h_a)
  !===============================================!
  !                                               !
  ! Compute D^H_a=UU^dag.del_a UU (a=alpha,beta), !
  ! using Eq.(24) of WYSV06                       !
  !                                               !
  !===============================================!

    use w90_constants, only      : dp,cmplx_0
    use w90_parameters, only     : num_wann
    use w90_utility, only        : utility_rotate
    use w90_postw90_common, only : get_occ

    ! Arguments
    !
    complex(kind=dp), dimension(:,:), intent(in)  :: delHH_a
    complex(kind=dp), dimension(:,:), intent(in)  :: UU
    real(kind=dp),    dimension(:),   intent(in)  :: eig
    real(kind=dp),                    intent(in)  :: ef
    complex(kind=dp), dimension(:,:), intent(out) :: D_h_a

    complex(kind=dp), allocatable :: delHH_a_bar(:,:)
    real(kind=dp)                 :: occ(num_wann)
    integer                       :: n,m

    call get_occ(eig,occ,ef)

    allocate(delHH_a_bar(num_wann,num_wann))
    delHH_a_bar=utility_rotate(delHH_a,UU,num_wann)
    do m=1,num_wann
       do n=1,num_wann
          if(occ(n)>0.999_dp.and.occ(m)<0.001_dp) then
             D_h_a(n,m)=delHH_a_bar(n,m)/(eig(m)-eig(n))
          else
             D_h_a(n,m)=cmplx_0
          end if
       end do
    end do
    D_h_a=D_h_a-conjg(transpose(D_h_a))

  end subroutine get_D_h_a

  subroutine get_D_h(delHH,UU,eig,D_h)
  !=========================================!
  !                                         !
  ! Compute D^H_a=UU^dag.del_a UU (a=x,y,z) !
  ! using Eq.(24) of WYSV06                 !
  !                                         !
  !=========================================!

    ! TO DO: Implement version where energy denominators only connect
    !        occupied and empty states. In this case probably do not need
    !        to worry about avoiding small energy denominators

    use w90_constants, only     : dp,cmplx_0
    use w90_parameters, only    : num_wann
    use w90_utility, only       : utility_rotate

    ! Arguments
    !
    complex(kind=dp), dimension(:,:,:), intent(in)  :: delHH
    complex(kind=dp), dimension(:,:), intent(in)    :: UU
    real(kind=dp),    dimension(:),   intent(in)    :: eig
    complex(kind=dp), dimension(:,:,:), intent(out) :: D_h

    complex(kind=dp), allocatable :: delHH_bar_i(:,:)
    integer                       :: n,m,i

    allocate(delHH_bar_i(num_wann,num_wann))
    D_h=cmplx_0
    do i=1,3
       delHH_bar_i(:,:)=utility_rotate(delHH(:,:,i),UU,num_wann)
       do m=1,num_wann
          do n=1,num_wann
             if(n==m .or. abs(eig(m)-eig(n))<1.0e-7_dp) cycle
             D_h(n,m,i)=delHH_bar_i(n,m)/(eig(m)-eig(n))
          end do
       end do
    enddo

  end subroutine get_D_h


  subroutine get_JJp_list(delHH,UU,eig,JJp_list)
  !====================================!
  !                                    !
  ! Compute JJ^+_a (a=Cartesian index) !
  ! for a list of Fermi energies       !
  !                                    !
  !====================================!

    use w90_constants, only  : dp,cmplx_0,cmplx_i
    use w90_parameters, only : num_wann,nfermi,fermi_energy_list
    use w90_utility, only    : utility_rotate

    complex(kind=dp), dimension(:,:), intent(in)    :: delHH
    complex(kind=dp), dimension(:,:), intent(in)    :: UU
    real(kind=dp),    dimension(:),   intent(in)    :: eig
    complex(kind=dp), dimension(:,:,:), intent(out) :: JJp_list

    complex(kind=dp), allocatable :: delHH_bar(:,:)
    integer                       :: n,m,if

    allocate(delHH_bar(num_wann,num_wann))
    delHH_bar=utility_rotate(delHH,UU,num_wann)
    do if=1,nfermi
       do m=1,num_wann
          do n=1,num_wann
             if(eig(n)>fermi_energy_list(if) .and.&
                  eig(m)<fermi_energy_list(if)) then
                JJp_list(n,m,if)=cmplx_i*delHH_bar(n,m)/(eig(m)-eig(n))
             else
                JJp_list(n,m,if)=cmplx_0
             end if
          enddo
       end do
    end do
    do if=1,nfermi
       JJp_list(:,:,if)=&
            utility_rotate(JJp_list(:,:,if),conjg(transpose(UU)),num_wann)
    enddo

  end subroutine get_JJp_list


  subroutine get_JJm_list(delHH,UU,eig,JJm_list)
  !====================================!
  !                                    !
  ! Compute JJ^-_a (a=Cartesian index) !
  ! for a list of Fermi energies       !
  !                                    !
  !====================================!

    use w90_constants, only  : dp,cmplx_0,cmplx_i
    use w90_parameters, only : num_wann,nfermi,fermi_energy_list
    use w90_utility, only    : utility_rotate

    complex(kind=dp), dimension(:,:), intent(in)   :: delHH
    complex(kind=dp), dimension(:,:), intent(in)   :: UU
    real(kind=dp),    dimension(:),   intent(in)   :: eig
    complex(kind=dp), dimension(:,:,:), intent(out) :: JJm_list

    complex(kind=dp), allocatable :: delHH_bar(:,:)
    integer                       :: n,m,if

    allocate(delHH_bar(num_wann,num_wann))
    delHH_bar=utility_rotate(delHH,UU,num_wann)
    do if=1,nfermi
       do m=1,num_wann
          do n=1,num_wann
             if(eig(m)>fermi_energy_list(if) .and.&
                  eig(n)<fermi_energy_list(if)) then
                JJm_list(n,m,if)=cmplx_i*delHH_bar(n,m)/(eig(m)-eig(n))
             else
                JJm_list(n,m,if)=cmplx_0
             end if
          enddo
       end do
    end do
    do if=1,nfermi
       JJm_list(:,:,if)=&
            utility_rotate(JJm_list(:,:,if),conjg(transpose(UU)),num_wann)
    enddo

  end subroutine get_JJm_list


  subroutine get_occ_mat_list(eig,UU,f_list,g_list)
  !================================!
  !                                !
  ! Occupation matrix f, and g=1-f !
  ! for a list of Fermi energies   !
  !                                !
  !================================!
    
    use w90_constants, only      : dp,cmplx_0,cmplx_1
    use w90_parameters, only     : num_wann,nfermi,fermi_energy_list
    use w90_postw90_common, only : get_occ

    ! Arguments
    !
    real(kind=dp),    dimension(:),     intent(in)  :: eig
    complex(kind=dp), dimension(:,:),   intent(in)  :: UU
    complex(kind=dp), dimension(:,:,:), intent(out) :: f_list
    complex(kind=dp), dimension(:,:,:), intent(out) :: g_list

    real(kind=dp) :: occ_list(num_wann,nfermi)
    integer       :: n,m,i,if

    do if=1,nfermi
       call get_occ(eig,occ_list(:,if),fermi_energy_list(if))
    enddo
    f_list=cmplx_0
    do if=1,nfermi
       do n=1,num_wann
          do m=1,num_wann
             do i=1,num_wann
                f_list(n,m,if)=f_list(n,m,if)&
                     +UU(n,i)*occ_list(i,if)*conjg(UU(m,i))
             enddo
             g_list(n,m,if)=-f_list(n,m,if)
             if(m==n) g_list(n,n,if)=g_list(n,n,if)+cmplx_1
          enddo
       enddo
    enddo

  end subroutine get_occ_mat_list


  subroutine get_deleig_a(deleig_a,eig,delHH_a,UU)
  !==========================!
  !                          !
  ! Band derivatives dE/dk_a !
  !                          !
  !==========================!

    use w90_constants, only  : dp,cmplx_0,cmplx_i
    use w90_utility, only    : utility_diagonalize,utility_rotate,&
                               utility_rotate_diag
    use w90_parameters, only : num_wann,use_degen_pert,degen_thr

    ! Arguments
    !
    real(kind=dp),                    intent(out) :: deleig_a(num_wann)
    real(kind=dp),                    intent(in)  :: eig(num_wann)
    complex(kind=dp), dimension(:,:), intent(in)  :: delHH_a
    complex(kind=dp), dimension(:,:), intent(in)  :: UU

    ! Misc/Dummy
    !
    integer                       :: i,degen_min,degen_max,dim
    real(kind=dp)                 :: diff
    complex(kind=dp), allocatable :: delHH_bar_a(:,:),U_deg(:,:)

    allocate(delHH_bar_a(num_wann,num_wann))
    allocate(U_deg(num_wann,num_wann))
    
    if(use_degen_pert) then
       
       delHH_bar_a=utility_rotate(delHH_a,UU,num_wann)
       
       ! Assuming that the energy eigenvalues are stored in eig(:) in
       ! increasing order (diff >= 0)
       
       i=0
       do 
          i=i+1
          if(i>num_wann) exit
          if(i+1 <= num_wann) then
             diff=eig(i+1)-eig(i)
          else
             !
             ! i-th is the highest band, and it is non-degenerate
             !
             diff =degen_thr+1.0_dp
          end if
          if(diff<degen_thr) then
             !
             ! Bands i and i+1 are degenerate 
             !
             degen_min=i
             degen_max=degen_min+1
             !
             ! See if any higher bands are in the same degenerate group
             !
             do
                if(degen_max+1>num_wann) exit
                diff=eig(degen_max+1)-eig(degen_max)
                if(diff<degen_thr) then
                   degen_max=degen_max+1
                else
                   exit
                end if
             end do
             !
             ! Bands from degen_min to degen_max are degenerate. Diagonalize 
             ! the submatrix in Eq.(31) YWVS07 over this degenerate subspace.
             ! The eigenvalues are the band gradients
             !
             !
             dim=degen_max-degen_min+1
             call utility_diagonalize(delHH_bar_a(degen_min:degen_max,&
                  degen_min:degen_max),dim,&
                  deleig_a(degen_min:degen_max),U_deg(1:dim,1:dim))
             !
             ! Scanned bands up to degen_max
             !
             i=degen_max
          else
             !
             ! Use non-degenerate form [Eq.(27) YWVS07] for current (i-th) band
             !
             deleig_a(i)=aimag(cmplx_i*delHH_bar_a(i,i))
          end if
       end do
       
    else
       
       ! Use non-degenerate form for all bands
       !
       deleig_a(:)=aimag(cmplx_i*utility_rotate_diag(delHH_a(:,:),UU,num_wann))

    end if
    
  end subroutine get_deleig_a

  !> Given a k point, this function returns eigenvalues E and
  !> derivatives of the eigenvalues dE/dk_a, using get_deleig_a
  !> \param kpt the three coordinates of the k point vector
  !>            (in relative coordinates)
  !> \param eig the calculated eigenvalues at kpt
  !> \param deleig the calculated derivatives of the eigenvalues
  !>        at kpt [first component: band; second component: 1,2,3
  !>        for the derivatives along the three k directions]
  !> \param HH the Hamiltonian matrix at kpt
  !> \param delHH the delHH matrix (derivative of H) at kpt
  !> \param UU the rotation matrix that gives the eigenvectors of HH
  subroutine get_eig_deleig(kpt,eig,del_eig,HH,delHH,UU)
    use w90_parameters, only: num_wann
    use w90_get_oper, only: HH_R, get_HH_R
    use w90_postw90_common, only : fourier_R_to_k
    use w90_utility, only : utility_diagonalize

    real(kind=dp), dimension(3), intent(in)         :: kpt
    real(kind=dp), intent(out)                      :: eig(num_wann)
    real(kind=dp), intent(out)                      :: del_eig(num_wann,3)
    complex(kind=dp), dimension(:,:), intent(out)   :: HH
    complex(kind=dp), dimension(:,:,:), intent(out) :: delHH
    complex(kind=dp), dimension(:,:), intent(out)   :: UU
   
    ! I call it to be sure that it has been called already once, 
    ! and that HH_R contains the actual matrix. 
    ! Further calls should return very fast.
    call get_HH_R

    call fourier_R_to_k(kpt,HH_R,HH,0) 
    call utility_diagonalize(HH,num_wann,eig,UU) 
    call fourier_R_to_k(kpt,HH_R,delHH(:,:,1),1) 
    call fourier_R_to_k(kpt,HH_R,delHH(:,:,2),2) 
    call fourier_R_to_k(kpt,HH_R,delHH(:,:,3),3) 
    call get_deleig_a(del_eig(:,1),eig,delHH(:,:,1),UU)
    call get_deleig_a(del_eig(:,2),eig,delHH(:,:,2),UU)
    call get_deleig_a(del_eig(:,3),eig,delHH(:,:,3),UU)

  end subroutine get_eig_deleig

  
  subroutine get_eig_UU_HH_JJlist(kpt,eig,UU,HH,JJp_list,JJm_list)
  !========================================================!
  !                                                        ! 
  ! Wrapper routine used to reduce number of Fourier calls !
  !                                                        ! 
  !========================================================!

    use w90_parameters, only     : num_wann
    use w90_get_oper, only       : HH_R,get_HH_R
    use w90_postw90_common, only : fourier_R_to_k_new
    use w90_utility, only        : utility_diagonalize

    real(kind=dp), dimension(3), intent(in)           :: kpt
    real(kind=dp), intent(out)                        :: eig(num_wann)
    complex(kind=dp), dimension(:,:), intent(out)     :: UU
    complex(kind=dp), dimension(:,:), intent(out)     :: HH
    complex(kind=dp), dimension(:,:,:,:), intent(out) :: JJp_list
    complex(kind=dp), dimension(:,:,:,:), intent(out) :: JJm_list

    integer                       :: i
    complex(kind=dp), allocatable :: delHH(:,:,:)

    call get_HH_R

    allocate(delHH(num_wann,num_wann,3))
    call fourier_R_to_k_new(kpt,HH_R,OO=HH,&
                                     OO_dx=delHH(:,:,1),&
                                     OO_dy=delHH(:,:,2),&
                                     OO_dz=delHH(:,:,3))
    call utility_diagonalize(HH,num_wann,eig,UU) 
    do i=1,3
       call get_JJp_list(delHH(:,:,i),UU,eig,JJp_list(:,:,:,i))
       call get_JJm_list(delHH(:,:,i),UU,eig,JJm_list(:,:,:,i))
    enddo

  end subroutine get_eig_UU_HH_JJlist

end module w90_wan_ham
