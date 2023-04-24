!
! Copyright (C) 2003-2013 Quantum ESPRESSO and Wannier90 groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! pw2wannier was written by Stefano de Gironcoli
! with later additions by
! Jonathan Yates - spinors
! Arash Mostofi - gamma point and transport things
! Timo Thonhauser, Graham Lopez, Ivo Souza
!         uHu, uIu terms for orbital magnetisation
! please send bugs and comments to
! Jonathan Yates and Arash Mostofi
!
! Known limitations:
!  spinors and orbital magnetisation term are not
!  yet inplimented for ultrasofts or PAW
!
!
!
module wannier
   USE kinds, only : DP
   !integer, allocatable :: nnb(:)       ! #b  (ik)
   integer              :: nnb          ! #b
   integer, allocatable :: kpb(:,:)     ! k+b (ik,ib)
   integer, allocatable :: g_kpb(:,:,:) ! G_k+b (ipol,ik,ib)
   integer, allocatable :: ig_(:,:)     ! G_k+b (ipol,ik,ib)
   integer, allocatable :: lw(:,:), mw(:,:) ! l and m of wannier (16,n_wannier)
   integer, allocatable :: num_sph(:)   ! num. func. in lin. comb., (n_wannier)
   logical, allocatable :: excluded_band(:)
   ! begin change Lopez, Thonhauser, Souza
   integer  :: iun_nnkp,iun_mmn,iun_amn,iun_band,iun_spn,iun_plot,iun_parity,&
        nnbx,nexband,iun_uhu,&
        iun_uIu !ivo
   ! end change Lopez, Thonhauser, Souza
   integer  :: n_wannier !number of WF
   integer  :: n_proj    !number of projection
   complex(DP), allocatable :: gf(:,:)  ! guding_function(npwx,n_wannier)
   integer               :: ispinw, ikstart, ikstop, iknum
   character(LEN=15)     :: wan_mode    ! running mode
   logical               :: logwann, wvfn_formatted, write_unk, &
   ! begin change Lopez, Thonhauser, Souza
                            write_amn,write_mmn,reduce_unk,write_spn,&
                            write_unkg,write_uhu,&
                            write_uIu, spn_formatted, uHu_formatted, uIu_formatted !ivo
   ! end change Lopez, Thonhauser, Souza
   ! input data from nnkp file
   real(DP), allocatable :: center_w(:,:)     ! center_w(3,n_wannier)
   integer,  allocatable :: spin_eig(:)
   real(DP), allocatable :: spin_qaxis(:,:)
   integer, allocatable  :: l_w(:), mr_w(:) ! l and mr of wannier (n_wannier) as from table 3.1,3.2 of spec.
   integer, allocatable  :: r_w(:)      ! index of radial function (n_wannier) as from table 3.3 of spec.
   real(DP), allocatable :: xaxis(:,:),zaxis(:,:) ! xaxis and zaxis(3,n_wannier)
   real(DP), allocatable :: alpha_w(:)  ! alpha_w(n_wannier) ( called zona in wannier spec)
   !
   real(DP), allocatable :: csph(:,:)    ! expansion coefficients of gf on QE ylm function (16,n_wannier)
   CHARACTER(len=256) :: seedname  = 'wannier'  ! prepended to file names in wannier90
   ! For implementation of wannier_lib
   integer               :: mp_grid(3)            ! dimensions of MP k-point grid
   real(DP)              :: rlatt(3,3),glatt(3,3) ! real and recip lattices (Cartesian co-ords, units of Angstrom)
   real(DP), allocatable :: kpt_latt(:,:)  ! k-points in crystal co-ords. kpt_latt(3,iknum)
   real(DP), allocatable :: atcart(:,:)    ! atom centres in Cartesian co-ords and Angstrom units. atcart(3,nat)
   integer               :: num_bands      ! number of bands left after exclusions
   character(len=3), allocatable :: atsym(:) ! atomic symbols. atsym(nat)
   integer               :: num_nnmax=12
   complex(DP), allocatable :: m_mat(:,:,:,:), a_mat(:,:,:)
   complex(DP), allocatable :: u_mat(:,:,:), u_mat_opt(:,:,:)
   logical, allocatable     :: lwindow(:,:)
   real(DP), allocatable    :: wann_centers(:,:),wann_spreads(:)
   real(DP)                 :: spreads(3)
   real(DP), allocatable    :: eigval(:,:)
   logical                  :: old_spinor_proj  ! for compatability for nnkp files prior to W90v2.0
end module wannier
!


!------------------------------------------------------------------------
PROGRAM pw2wannier90
  ! This is the interface to the Wannier90 code: see http://www.wannier.org
  !------------------------------------------------------------------------
  !
  USE io_global,  ONLY : stdout, ionode, ionode_id
  USE mp_global,  ONLY : mp_startup, mpime, kunit
  USE mp,         ONLY : mp_bcast
  USE cell_base,  ONLY : at, bg
  USE lsda_mod,   ONLY : nspin, isk
  USE klist,      ONLY : nkstot
  USE io_files,   ONLY : prefix, tmp_dir
  USE noncollin_module, ONLY : noncolin
  USE control_flags,    ONLY : gamma_only
  USE environment,ONLY : environment_start
  USE wvfct,      ONLY : ecutwfc
  USE wannier
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  INTEGER :: ios
  CHARACTER(len=4) :: spin_component
  CHARACTER(len=256) :: outdir

  ! these are in wannier module.....-> integer :: ispinw, ikstart, ikstop, iknum
  NAMELIST / inputpp / outdir, prefix, spin_component, wan_mode, &
       seedname, write_unk, write_amn, write_mmn, write_spn, &
   ! begin change Lopez, Thonhauser, Souza
       wvfn_formatted, reduce_unk, write_unkg, write_uhu,&
       write_uIu, spn_formatted, uHu_formatted, uIu_formatted !ivo
   ! end change Lopez, Thonhauser, Souza
  !
  ! initialise environment
  !
#ifdef __MPI
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'PW2WANNIER' )
  !
  ! Read input on i/o node and broadcast to the rest
  !
  ios = 0
  IF(ionode) THEN
     !
     ! Check to see if we are reading from a file
     !
     CALL input_from_file()
     !
     !   set default values for variables in namelist
     !
     CALL get_env( 'ESPRESSO_TMPDIR', outdir )
     IF ( trim( outdir ) == ' ' ) outdir = './'
     prefix = ' '
     seedname = 'wannier'
     spin_component = 'none'
     wan_mode = 'standalone'
     wvfn_formatted = .false.
     spn_formatted=.false.
     uHu_formatted=.false.
     uIu_formatted=.false.
     write_unk = .false.
     write_amn = .true.
     write_mmn = .true.
     write_spn = .false.
     ! begin change Lopez, Thonhauser, Souza
     write_uhu = .false.
     write_uIu = .false. !ivo
     ! end change Lopez, Thonhauser, Souza
     reduce_unk= .false.
     write_unkg= .false.
     !
     !     reading the namelist inputpp
     !
     READ (5, inputpp, iostat=ios)
     !
     !     Check of namelist variables
     !
     tmp_dir = trimcheck(outdir)
     ! back to all nodes
  ENDIF
  !
  CALL mp_bcast(ios,ionode_id)
  IF (ios /= 0) CALL errore( 'pw2wannier90', 'reading inputpp namelist', abs(ios))
  !
  ! broadcast input variable to all nodes
  !
  CALL mp_bcast(outdir,ionode_id)
  CALL mp_bcast(tmp_dir,ionode_id)
  CALL mp_bcast(prefix,ionode_id)
  CALL mp_bcast(seedname,ionode_id)
  CALL mp_bcast(spin_component,ionode_id)
  CALL mp_bcast(wan_mode,ionode_id)
  CALL mp_bcast(wvfn_formatted,ionode_id)
  CALL mp_bcast(write_unk,ionode_id)
  CALL mp_bcast(write_amn,ionode_id)
  CALL mp_bcast(write_mmn,ionode_id)
  ! begin change Lopez, Thonhauser, Souza
  CALL mp_bcast(write_uhu,ionode_id)
  CALL mp_bcast(write_uIu,ionode_id) !ivo
  ! end change Lopez, Thonhauser, Souza
  CALL mp_bcast(write_spn,ionode_id)
  CALL mp_bcast(reduce_unk,ionode_id)
  CALL mp_bcast(write_unkg,ionode_id)
  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  logwann = .true.
  WRITE(stdout,*)
  WRITE(stdout,*) ' Reading nscf_save data'
  CALL read_file
  WRITE(stdout,*)
  !
  IF (noncolin.and.gamma_only) CALL errore('pw2wannier90',&
       'Non-collinear and gamma_only not implemented',1)

  ! Here we should trap restarts from a different number of nodes.
  ! or attempts at kpoint distribution
  !
  SELECT CASE ( trim( spin_component ) )
  CASE ( 'up' )
     WRITE(stdout,*) ' Spin CASE ( up )'
     ispinw  = 1
     ikstart = 1
     ikstop  = nkstot/2
     iknum   = nkstot/2
  CASE ( 'down' )
     WRITE(stdout,*) ' Spin CASE ( down )'
     ispinw = 2
     ikstart = nkstot/2 + 1
     ikstop  = nkstot
     iknum   = nkstot/2
  CASE DEFAULT
     IF(noncolin) THEN
        WRITE(stdout,*) ' Spin CASE ( non-collinear )'
     ELSE
        WRITE(stdout,*) ' Spin CASE ( default = unpolarized )'
     ENDIF
     ispinw = 0
     ikstart = 1
     ikstop  = nkstot
     iknum   = nkstot
  END SELECT
  !
  WRITE(stdout,*)
  WRITE(stdout,*) ' Wannier mode is: ',wan_mode
  WRITE(stdout,*)
  !
  IF(wan_mode=='standalone') THEN
     !
     WRITE(stdout,*) ' -----------------'
     WRITE(stdout,*) ' *** Reading nnkp '
     WRITE(stdout,*) ' -----------------'
     WRITE(stdout,*)
     CALL read_nnkp
     WRITE(stdout,*) ' Opening pp-files '
     CALL openfil_pp
     CALL ylm_expansion
     WRITE(stdout,*)
     WRITE(stdout,*)
     IF(write_amn) THEN
        WRITE(stdout,*) ' ---------------'
        WRITE(stdout,*) ' *** Compute  A '
        WRITE(stdout,*) ' ---------------'
        WRITE(stdout,*)
        CALL compute_amn
        WRITE(stdout,*)
     ELSE
        WRITE(stdout,*) ' -----------------------------'
        WRITE(stdout,*) ' *** A matrix is not computed '
        WRITE(stdout,*) ' -----------------------------'
        WRITE(stdout,*)
     ENDIF
    IF(write_mmn) THEN
        WRITE(stdout,*) ' ---------------'
        WRITE(stdout,*) ' *** Compute  M '
        WRITE(stdout,*) ' ---------------'
        WRITE(stdout,*)
        CALL compute_mmn
        WRITE(stdout,*)
     ELSE
        WRITE(stdout,*) ' -----------------------------'
        WRITE(stdout,*) ' *** M matrix is not computed '
        WRITE(stdout,*) ' -----------------------------'
        WRITE(stdout,*)
     ENDIF
     if(noncolin) then
        IF(write_spn) THEN
           WRITE(stdout,*) ' ------------------'
           WRITE(stdout,*) ' *** Compute  Spin '
           WRITE(stdout,*) ' ------------------'
           WRITE(stdout,*)
           CALL compute_spin
           WRITE(stdout,*)
        ELSE
           WRITE(stdout,*) ' --------------------------------'
           WRITE(stdout,*) ' *** Spin matrix is not computed '
           WRITE(stdout,*) ' --------------------------------'
           WRITE(stdout,*)
        ENDIF
     elseif(write_spn) then
        write(stdout,*) ' -----------------------------------'
        write(stdout,*) ' *** Non-collinear calculation is   '
        write(stdout,*) '     required for spin              '
        write(stdout,*) '     term  to be computed           '
        write(stdout,*) ' -----------------------------------'
     endif
     IF(write_uHu.or.write_uIu) THEN
        WRITE(stdout,*) ' ----------------'
        WRITE(stdout,*) ' *** Compute Orb '
        WRITE(stdout,*) ' ----------------'
        WRITE(stdout,*)
        CALL compute_orb
        WRITE(stdout,*)
     ELSE
        WRITE(stdout,*) ' -----------------------------------'
        WRITE(stdout,*) ' *** Orbital terms are not computed '
        WRITE(stdout,*) ' -----------------------------------'
        WRITE(stdout,*)
     ENDIF
     WRITE(stdout,*) ' ----------------'
     WRITE(stdout,*) ' *** Write bands '
     WRITE(stdout,*) ' ----------------'
     WRITE(stdout,*)
     CALL write_band
     WRITE(stdout,*)
     IF(write_unk) THEN
        WRITE(stdout,*) ' --------------------'
        WRITE(stdout,*) ' *** Write plot info '
        WRITE(stdout,*) ' --------------------'
        WRITE(stdout,*)
        CALL write_plot
        WRITE(stdout,*)
     ELSE
        WRITE(stdout,*) ' -----------------------------'
        WRITE(stdout,*) ' *** Plot info is not printed '
        WRITE(stdout,*) ' -----------------------------'
        WRITE(stdout,*)
     ENDIF
     IF(write_unkg) THEN
        WRITE(stdout,*) ' --------------------'
        WRITE(stdout,*) ' *** Write parity info '
        WRITE(stdout,*) ' --------------------'
        WRITE(stdout,*)
        CALL write_parity
        WRITE(stdout,*)
     ELSE
        WRITE(stdout,*) ' -----------------------------'
        WRITE(stdout,*) ' *** Parity info is not printed '
        WRITE(stdout,*) ' -----------------------------'
        WRITE(stdout,*)
     ENDIF
     WRITE(stdout,*) ' ------------'
     WRITE(stdout,*) ' *** Stop pp '
     WRITE(stdout,*) ' ------------'
     WRITE(stdout,*)
     CALL stop_pp
     !
  ENDIF
  !
  IF(wan_mode=='library') THEN
     !
!     seedname='wannier'
     WRITE(stdout,*) ' Setting up...'
     CALL setup_nnkp
     WRITE(stdout,*)
     WRITE(stdout,*) ' Opening pp-files '
     CALL openfil_pp
     WRITE(stdout,*)
     WRITE(stdout,*) ' Ylm expansion'
     CALL ylm_expansion
     WRITE(stdout,*)
     CALL compute_amn
     CALL compute_mmn
     CALL write_band
     CALL run_wannier
     IF(write_unk) CALL write_plot
     CALL lib_dealloc
     CALL stop_pp
     !
  ENDIF
  !
  IF(wan_mode=='wannier2sic') THEN
     !
     CALL read_nnkp
     CALL wan2sic
     !
  ENDIF
  !
  STOP
END PROGRAM pw2wannier90
!
!-----------------------------------------------------------------------
SUBROUTINE lib_dealloc
  !-----------------------------------------------------------------------
  !
  USE wannier

  IMPLICIT NONE

  DEALLOCATE(m_mat,u_mat,u_mat_opt,a_mat,eigval)

  RETURN
END SUBROUTINE lib_dealloc
!
!-----------------------------------------------------------------------
SUBROUTINE setup_nnkp
  !-----------------------------------------------------------------------
  !
  USE io_global, ONLY : stdout, ionode, ionode_id
  USE kinds,     ONLY : DP
  USE constants, ONLY : eps6, tpi, bohr => BOHR_RADIUS_ANGS
  USE cell_base, ONLY : at, bg, alat
  USE gvect,     ONLY : g, gg
  USE ions_base, ONLY : nat, tau, ityp, atm
  USE klist,     ONLY : xk
  USE mp,        ONLY : mp_bcast, mp_sum
  USE mp_global, ONLY : intra_pool_comm
  USE wvfct,     ONLY : nbnd,npwx
  USE control_flags,    ONLY : gamma_only
  USE noncollin_module, ONLY : noncolin
  USE wannier

  IMPLICIT NONE
  real(DP) :: g_(3), gg_
  INTEGER  :: ik, ib, ig, iw, ia, indexb, TYPE
  INTEGER, ALLOCATABLE :: ig_check(:,:)
  real(DP) :: xnorm, znorm, coseno
  INTEGER  :: exclude_bands(nbnd)

  ! aam: translations between PW2Wannier90 and Wannier90
  ! pw2wannier90   <==>   Wannier90
  !    nbnd                num_bands_tot
  !    n_wannier           num_wann
  !    num_bands           num_bands
  !    nat                 num_atoms
  !    iknum               num_kpts
  !    rlatt               transpose(real_lattice)
  !    glatt               transpose(recip_lattice)
  !    kpt_latt            kpt_latt
  !    nnb                 nntot
  !    kpb                 nnlist
  !    g_kpb               nncell
  !    mp_grid             mp_grid
  !    center_w            proj_site
  !    l_w,mr_w,r_w        proj_l,proj_m,proj_radial
  !    xaxis,zaxis         proj_x,proj_z
  !    alpha_w             proj_zona
  !    exclude_bands       exclude_bands
  !    atcart              atoms_cart
  !    atsym               atom_symbols

  ALLOCATE( kpt_latt(3,iknum) )
  ALLOCATE( atcart(3,nat), atsym(nat) )
  ALLOCATE( kpb(iknum,num_nnmax), g_kpb(3,iknum,num_nnmax) )
  ALLOCATE( center_w(3,nbnd), alpha_w(nbnd), l_w(nbnd), &
       mr_w(nbnd), r_w(nbnd), zaxis(3,nbnd), xaxis(3,nbnd) )
  ALLOCATE( excluded_band(nbnd) )

  ! real lattice (Cartesians, Angstrom)
  rlatt(:,:) = transpose(at(:,:))*alat*bohr
  ! reciprocal lattice (Cartesians, Angstrom)
  glatt(:,:) = transpose(bg(:,:))*tpi/(alat*bohr)
  ! convert Cartesian k-points to crystallographic co-ordinates
  kpt_latt(:,1:iknum)=xk(:,1:iknum)
  CALL cryst_to_cart(iknum,kpt_latt,at,-1)
  ! atom co-ordinates in Cartesian co-ords and Angstrom units
  atcart(:,:) = tau(:,:)*bohr*alat
  ! atom symbols
  DO ia=1,nat
     TYPE=ityp(ia)
     atsym(ia)=atm(TYPE)
  ENDDO

  ! MP grid dimensions
  CALL find_mp_grid()

  WRITE(stdout,'("  - Number of atoms is (",i3,")")') nat

#ifdef __WANLIB
  IF (ionode) THEN
     CALL wannier_setup(seedname,mp_grid,iknum,rlatt, &               ! input
          glatt,kpt_latt,nbnd,nat,atsym,atcart,gamma_only,noncolin, & ! input
          nnb,kpb,g_kpb,num_bands,n_wannier,center_w, &               ! output
          l_w,mr_w,r_w,zaxis,xaxis,alpha_w,exclude_bands)             ! output
  ENDIF
#endif

  CALL mp_bcast(nnb,ionode_id)
  CALL mp_bcast(kpb,ionode_id)
  CALL mp_bcast(g_kpb,ionode_id)
  CALL mp_bcast(num_bands,ionode_id)
  CALL mp_bcast(n_wannier,ionode_id)
  CALL mp_bcast(center_w,ionode_id)
  CALL mp_bcast(l_w,ionode_id)
  CALL mp_bcast(mr_w,ionode_id)
  CALL mp_bcast(r_w,ionode_id)
  CALL mp_bcast(zaxis,ionode_id)
  CALL mp_bcast(xaxis,ionode_id)
  CALL mp_bcast(alpha_w,ionode_id)
  CALL mp_bcast(exclude_bands,ionode_id)

  IF(noncolin) THEN
     n_proj=n_wannier/2
  ELSE
     n_proj=n_wannier
  ENDIF

  ALLOCATE( gf(npwx,n_proj), csph(16,n_proj) )

  WRITE(stdout,'("  - Number of wannier functions is (",i3,")")') n_wannier

  excluded_band(1:nbnd)=.false.
  nexband=0
  band_loop: DO ib=1,nbnd
     indexb=exclude_bands(ib)
     IF (indexb>nbnd .or. indexb<0) THEN
        CALL errore('setup_nnkp',' wrong excluded band index ', 1)
     ELSEIF (indexb==0) THEN
        exit band_loop
     ELSE
        nexband=nexband+1
        excluded_band(indexb)=.true.
     ENDIF
  ENDDO band_loop

  IF ( (nbnd-nexband)/=num_bands ) &
       CALL errore('setup_nnkp',' something wrong with num_bands',1)

  DO iw=1,n_proj
     xnorm = sqrt(xaxis(1,iw)*xaxis(1,iw) + xaxis(2,iw)*xaxis(2,iw) + &
          xaxis(3,iw)*xaxis(3,iw))
     IF (xnorm < eps6) CALL errore ('setup_nnkp',' |xaxis| < eps ',1)
     znorm = sqrt(zaxis(1,iw)*zaxis(1,iw) + zaxis(2,iw)*zaxis(2,iw) + &
          zaxis(3,iw)*zaxis(3,iw))
     IF (znorm < eps6) CALL errore ('setup_nnkp',' |zaxis| < eps ',1)
     coseno = (xaxis(1,iw)*zaxis(1,iw) + xaxis(2,iw)*zaxis(2,iw) + &
          xaxis(3,iw)*zaxis(3,iw))/xnorm/znorm
     IF (abs(coseno) > eps6) &
          CALL errore('setup_nnkp',' xaxis and zaxis are not orthogonal !',1)
     IF (alpha_w(iw) < eps6) &
          CALL errore('setup_nnkp',' zona value must be positive', 1)
     ! convert wannier center in cartesian coordinates (in unit of alat)
     CALL cryst_to_cart( 1, center_w(:,iw), at, 1 )
  ENDDO
  WRITE(stdout,*) ' - All guiding functions are given '

  nnbx=0
  nnb=max(nnbx,nnb)

  ALLOCATE( ig_(iknum,nnb), ig_check(iknum,nnb) )

  DO ik=1, iknum
     DO ib = 1, nnb
        g_(:) = REAL( g_kpb(:,ik,ib) )
        CALL cryst_to_cart (1, g_, bg, 1)
        gg_ = g_(1)*g_(1) + g_(2)*g_(2) + g_(3)*g_(3)
        ig_(ik,ib) = 0
        ig = 1
        DO WHILE  (gg(ig) <= gg_ + eps6)
           IF ( (abs(g(1,ig)-g_(1)) < eps6) .and.  &
                (abs(g(2,ig)-g_(2)) < eps6) .and.  &
                (abs(g(3,ig)-g_(3)) < eps6)  ) ig_(ik,ib) = ig
           ig= ig +1
        ENDDO
     ENDDO
  ENDDO

  ig_check(:,:) = ig_(:,:)
  CALL mp_sum( ig_check, intra_pool_comm )
  DO ik=1, iknum
     DO ib = 1, nnb
        IF (ig_check(ik,ib) ==0) &
          CALL errore('setup_nnkp', &
                      ' g_kpb vector is not in the list of Gs', 100*ik+ib )
     ENDDO
  ENDDO
  DEALLOCATE (ig_check)

  WRITE(stdout,*) ' - All neighbours are found '
  WRITE(stdout,*)

  RETURN
END SUBROUTINE setup_nnkp
 !
 !-----------------------------------------------------------------------
SUBROUTINE run_wannier
  !-----------------------------------------------------------------------
  !
  USE io_global, ONLY : ionode, ionode_id
  USE ions_base, ONLY : nat
  USE mp,        ONLY : mp_bcast
  USE control_flags, ONLY : gamma_only
  USE wannier

  IMPLICIT NONE

  ALLOCATE(u_mat(n_wannier,n_wannier,iknum))
  ALLOCATE(u_mat_opt(num_bands,n_wannier,iknum))
  ALLOCATE(lwindow(num_bands,iknum))
  ALLOCATE(wann_centers(3,n_wannier))
  ALLOCATE(wann_spreads(n_wannier))

#ifdef __WANLIB
  IF (ionode) THEN
     CALL wannier_run(seedname,mp_grid,iknum,rlatt, &                ! input
          glatt,kpt_latt,num_bands,n_wannier,nnb,nat, &              ! input
          atsym,atcart,gamma_only,m_mat,a_mat,eigval, &              ! input
          u_mat,u_mat_opt,lwindow,wann_centers,wann_spreads,spreads) ! output
  ENDIF
#endif

  CALL mp_bcast(u_mat,ionode_id)
  CALL mp_bcast(u_mat_opt,ionode_id)
  CALL mp_bcast(lwindow,ionode_id)
  CALL mp_bcast(wann_centers,ionode_id)
  CALL mp_bcast(wann_spreads,ionode_id)
  CALL mp_bcast(spreads,ionode_id)

  RETURN
END SUBROUTINE run_wannier
!-----------------------------------------------------------------------
!
SUBROUTINE find_mp_grid()
  !-----------------------------------------------------------------------
  !
  USE io_global, ONLY : stdout
  USE kinds,     ONLY: DP
  USE wannier

  IMPLICIT NONE

  ! <<<local variables>>>
  INTEGER  :: ik,ntemp,ii
  real(DP) :: min_k,temp(3,iknum),mpg1

  min_k=minval(kpt_latt(1,:))
  ii=0
  DO ik=1,iknum
     IF (kpt_latt(1,ik)==min_k) THEN
        ii=ii+1
        temp(:,ii)=kpt_latt(:,ik)
     ENDIF
  ENDDO
  ntemp=ii

  min_k=minval(temp(2,1:ntemp))
  ii=0
  DO ik=1,ntemp
     IF (temp(2,ik)==min_k) THEN
        ii=ii+1
     ENDIF
  ENDDO
  mp_grid(3)=ii

  min_k=minval(temp(3,1:ntemp))
  ii=0
  DO ik=1,ntemp
     IF (temp(3,ik)==min_k) THEN
        ii=ii+1
     ENDIF
  ENDDO
  mp_grid(2)=ii

  IF ( (mp_grid(2)==0) .or. (mp_grid(3)==0) ) &
       CALL errore('find_mp_grid',' one or more mp_grid dimensions is zero', 1)

  mpg1=iknum/(mp_grid(2)*mp_grid(3))

  mp_grid(1) = nint(mpg1)

  WRITE(stdout,*)
  WRITE(stdout,'(3(a,i3))') '  MP grid is ',mp_grid(1),' x',mp_grid(2),' x',mp_grid(3)

  IF (real(mp_grid(1),kind=DP)/=mpg1) &
       CALL errore('find_mp_grid',' determining mp_grid failed', 1)

  RETURN
END SUBROUTINE find_mp_grid
!-----------------------------------------------------------------------
!
SUBROUTINE read_nnkp
  !-----------------------------------------------------------------------
  !
  USE io_global, ONLY : stdout, ionode, ionode_id
  USE kinds,     ONLY: DP
  USE constants, ONLY : eps6, tpi, bohr => BOHR_RADIUS_ANGS
  USE cell_base, ONLY : at, bg, alat
  USE gvect,     ONLY : g, gg
  USE klist,     ONLY : nkstot, xk
  USE mp,        ONLY : mp_bcast, mp_sum
  USE mp_global, ONLY : intra_pool_comm
  USE wvfct,     ONLY : npwx, nbnd
  USE noncollin_module, ONLY : noncolin
  USE wannier

  IMPLICIT NONE
  !
  INTEGER, EXTERNAL :: find_free_unit
  !
  real(DP) :: g_(3), gg_
  INTEGER :: ik, ib, ig, ipol, iw, idum, indexb
  INTEGER numk, i, j
  INTEGER, ALLOCATABLE :: ig_check(:,:)
  real(DP) :: xx(3), xnorm, znorm, coseno
  LOGICAL :: have_nnkp,found

  IF (ionode) THEN  ! Read nnkp file on ionode only

     INQUIRE(file=trim(seedname)//".nnkp",exist=have_nnkp)
     IF(.not. have_nnkp) THEN
        WRITE(stdout,*) ' Could not find the file '//trim(seedname)//'.nnkp'
        STOP
     ENDIF

     iun_nnkp = find_free_unit()
     OPEN (unit=iun_nnkp, file=trim(seedname)//".nnkp",form='formatted')

  ENDIF

  nnbx=0

  !   check the information from *.nnkp with the nscf_save data
  WRITE(stdout,*) ' Checking info from wannier.nnkp file'
  WRITE(stdout,*)

  IF (ionode) THEN   ! read from ionode only

     CALL scan_file_to('real_lattice',found)
     if(.not.found) then
        WRITE(stdout,*) ' Could not find real_lattice block in '//trim(seedname)//'.nnkp'
        STOP
     endif
     DO j=1,3
        READ(iun_nnkp,*) (rlatt(i,j),i=1,3)
        DO i = 1,3
           rlatt(i,j) = rlatt(i,j)/(alat*bohr)
        ENDDO
     ENDDO
     DO j=1,3
        DO i=1,3
           IF(abs(rlatt(i,j)-at(i,j))>eps6) THEN
              WRITE(stdout,*)  ' Something wrong! '
              WRITE(stdout,*)  ' rlatt(i,j) =',rlatt(i,j),  ' at(i,j)=',at(i,j)
              STOP
           ENDIF
        ENDDO
     ENDDO
     WRITE(stdout,*) ' - Real lattice is ok'

     CALL scan_file_to('recip_lattice',found)
     if(.not.found) then
        WRITE(stdout,*) ' Could not find recip_lattice block in '//trim(seedname)//'.nnkp'
        STOP
     endif
     DO j=1,3
        READ(iun_nnkp,*) (glatt(i,j),i=1,3)
        DO i = 1,3
           glatt(i,j) = (alat*bohr)*glatt(i,j)/tpi
        ENDDO
     ENDDO
     DO j=1,3
        DO i=1,3
           IF(abs(glatt(i,j)-bg(i,j))>eps6) THEN
              WRITE(stdout,*)  ' Something wrong! '
              WRITE(stdout,*)  ' glatt(i,j)=',glatt(i,j), ' bg(i,j)=',bg(i,j)
              STOP
           ENDIF
        ENDDO
     ENDDO
     WRITE(stdout,*) ' - Reciprocal lattice is ok'

     CALL scan_file_to('kpoints',found)
     if(.not.found) then
        WRITE(stdout,*) ' Could not find kpoints block in '//trim(seedname)//'.nnkp'
        STOP
     endif
     READ(iun_nnkp,*) numk
     IF(numk/=iknum) THEN
        WRITE(stdout,*)  ' Something wrong! '
        WRITE(stdout,*)  ' numk=',numk, ' iknum=',iknum
        STOP
     ENDIF
     DO i=1,numk
        READ(iun_nnkp,*) xx(1), xx(2), xx(3)
        CALL cryst_to_cart( 1, xx, bg, 1 )
        IF(abs(xx(1)-xk(1,i))>eps6.or. &
             abs(xx(2)-xk(2,i))>eps6.or. &
             abs(xx(3)-xk(3,i))>eps6) THEN
           WRITE(stdout,*)  ' Something wrong! '
           WRITE(stdout,*) ' k-point ',i,' is wrong'
           WRITE(stdout,*) xx(1), xx(2), xx(3)
           WRITE(stdout,*) xk(1,i), xk(2,i), xk(3,i)
           STOP
        ENDIF
     ENDDO
     WRITE(stdout,*) ' - K-points are ok'

  ENDIF ! ionode

  ! Broadcast
  CALL mp_bcast(rlatt,ionode_id)
  CALL mp_bcast(glatt,ionode_id)

  IF (ionode) THEN   ! read from ionode only
     if(noncolin) then
        old_spinor_proj=.false.
        CALL scan_file_to('spinor_projections',found)
        if(.not.found) then
           !try old style projections
           CALL scan_file_to('projections',found)
           if(found) then
              old_spinor_proj=.true.
           else
              WRITE(stdout,*) ' Could not find projections block in '//trim(seedname)//'.nnkp'
              STOP
           endif
        end if
     else
        CALL scan_file_to('projections',found)
        if(.not.found) then
           WRITE(stdout,*) ' Could not find projections block in '//trim(seedname)//'.nnkp'
           STOP
        endif
     endif
     READ(iun_nnkp,*) n_proj
  ENDIF

  ! Broadcast
  CALL mp_bcast(n_proj,ionode_id)
  CALL mp_bcast(old_spinor_proj,ionode_id)


  IF(old_spinor_proj) THEN
     n_wannier=n_proj*2
  ELSE
     n_wannier=n_proj
  ENDIF



  ALLOCATE( center_w(3,n_proj), alpha_w(n_proj), gf(npwx,n_proj), &
       l_w(n_proj), mr_w(n_proj), r_w(n_proj), &
       zaxis(3,n_proj), xaxis(3,n_proj), csph(16,n_proj) )
  if(noncolin.and..not.old_spinor_proj) then
     ALLOCATE( spin_eig(n_proj),spin_qaxis(3,n_proj) )
  endif

  WRITE(stdout,'("  - Number of wannier functions is ok (",i3,")")') n_wannier

  IF (ionode) THEN   ! read from ionode only
     DO iw=1,n_proj
        READ(iun_nnkp,*) (center_w(i,iw), i=1,3), l_w(iw), mr_w(iw), r_w(iw)
        READ(iun_nnkp,*) (zaxis(i,iw),i=1,3),(xaxis(i,iw),i=1,3),alpha_w(iw)
        xnorm = sqrt(xaxis(1,iw)*xaxis(1,iw) + xaxis(2,iw)*xaxis(2,iw) + &
             xaxis(3,iw)*xaxis(3,iw))
        IF (xnorm < eps6) CALL errore ('read_nnkp',' |xaxis| < eps ',1)
        znorm = sqrt(zaxis(1,iw)*zaxis(1,iw) + zaxis(2,iw)*zaxis(2,iw) + &
             zaxis(3,iw)*zaxis(3,iw))
        IF (znorm < eps6) CALL errore ('read_nnkp',' |zaxis| < eps ',1)
        coseno = (xaxis(1,iw)*zaxis(1,iw) + xaxis(2,iw)*zaxis(2,iw) + &
             xaxis(3,iw)*zaxis(3,iw))/xnorm/znorm
        IF (abs(coseno) > eps6) &
             CALL errore('read_nnkp',' xaxis and zaxis are not orthogonal !',1)
        IF (alpha_w(iw) < eps6) &
             CALL errore('read_nnkp',' zona value must be positive', 1)
        ! convert wannier center in cartesian coordinates (in unit of alat)
        CALL cryst_to_cart( 1, center_w(:,iw), at, 1 )
        if(noncolin.and..not.old_spinor_proj) then
           READ(iun_nnkp,*) spin_eig(iw),(spin_qaxis(i,iw),i=1,3)
           xnorm = sqrt(spin_qaxis(1,iw)*spin_qaxis(1,iw) + spin_qaxis(2,iw)*spin_qaxis(2,iw) + &
             spin_qaxis(3,iw)*spin_qaxis(3,iw))
           IF (xnorm < eps6) CALL errore ('read_nnkp',' |xaxis| < eps ',1)
           spin_qaxis(:,iw)=spin_qaxis(:,iw)/xnorm
        endif
     ENDDO
  ENDIF

  WRITE(stdout,*) ' - All guiding functions are given '

  ! Broadcast
  CALL mp_bcast(center_w,ionode_id)
  CALL mp_bcast(l_w,ionode_id)
  CALL mp_bcast(mr_w,ionode_id)
  CALL mp_bcast(r_w,ionode_id)
  CALL mp_bcast(zaxis,ionode_id)
  CALL mp_bcast(xaxis,ionode_id)
  CALL mp_bcast(alpha_w,ionode_id)
  CALL mp_bcast(spin_eig,ionode_id)
  CALL mp_bcast(spin_qaxis,ionode_id)
  !
  WRITE(stdout,*)
  WRITE(stdout,*) 'Projections:'
  DO iw=1,n_proj
     WRITE(stdout,'(3f12.6,3i3,f12.6)') &
          center_w(1:3,iw),l_w(iw),mr_w(iw),r_w(iw),alpha_w(iw)
  ENDDO

  IF (ionode) THEN   ! read from ionode only
     CALL scan_file_to('nnkpts',found)
     if(.not.found) then
        WRITE(stdout,*) ' Could not find nnkpts block in '//trim(seedname)//'.nnkp'
        STOP
     endif
     READ (iun_nnkp,*) nnb
  ENDIF

  ! Broadcast
  CALL mp_bcast(nnb,ionode_id)
  !
  nnbx = max (nnbx, nnb )
  !
  ALLOCATE ( kpb(iknum,nnbx), g_kpb(3,iknum,nnbx),&
             ig_(iknum,nnbx), ig_check(iknum,nnbx) )

  !  read data about neighbours
  WRITE(stdout,*)
  WRITE(stdout,*) ' Reading data about k-point neighbours '
  WRITE(stdout,*)

  IF (ionode) THEN
     DO ik=1, iknum
        DO ib = 1, nnb
           READ(iun_nnkp,*) idum, kpb(ik,ib), (g_kpb(ipol,ik,ib), ipol =1,3)
        ENDDO
     ENDDO
  ENDIF

  ! Broadcast
  CALL mp_bcast(kpb,ionode_id)
  CALL mp_bcast(g_kpb,ionode_id)

  DO ik=1, iknum
     DO ib = 1, nnb
        g_(:) = REAL( g_kpb(:,ik,ib) )
        CALL cryst_to_cart (1, g_, bg, 1)
        gg_ = g_(1)*g_(1) + g_(2)*g_(2) + g_(3)*g_(3)
        ig_(ik,ib) = 0
        ig = 1
        DO WHILE  (gg(ig) <= gg_ + eps6)
           IF ( (abs(g(1,ig)-g_(1)) < eps6) .and.  &
                (abs(g(2,ig)-g_(2)) < eps6) .and.  &
                (abs(g(3,ig)-g_(3)) < eps6)  ) ig_(ik,ib) = ig
           ig= ig +1
        ENDDO
     ENDDO
  ENDDO
  ig_check(:,:) = ig_(:,:)
  CALL mp_sum( ig_check, intra_pool_comm )
  DO ik=1, iknum
     DO ib = 1, nnb
        IF (ig_check(ik,ib) ==0) &
          CALL errore('read_nnkp', &
                      ' g_kpb vector is not in the list of Gs', 100*ik+ib )
     ENDDO
  ENDDO
  DEALLOCATE (ig_check)

  WRITE(stdout,*) ' All neighbours are found '
  WRITE(stdout,*)

  ALLOCATE( excluded_band(nbnd) )

  IF (ionode) THEN     ! read from ionode only
     CALL scan_file_to('exclude_bands',found)
     if(.not.found) then
        WRITE(stdout,*) ' Could not find exclude_bands block in '//trim(seedname)//'.nnkp'
        STOP
     endif
     READ (iun_nnkp,*) nexband
     excluded_band(1:nbnd)=.false.
     DO i=1,nexband
        READ(iun_nnkp,*) indexb
        IF (indexb<1 .or. indexb>nbnd) &
             CALL errore('read_nnkp',' wrong excluded band index ', 1)
        excluded_band(indexb)=.true.
     ENDDO
  ENDIF
  num_bands=nbnd-nexband

  ! Broadcast
  CALL mp_bcast(nexband,ionode_id)
  CALL mp_bcast(excluded_band,ionode_id)
  CALL mp_bcast(num_bands,ionode_id)

  IF (ionode) CLOSE (iun_nnkp)   ! ionode only

  RETURN
END SUBROUTINE read_nnkp
!
!-----------------------------------------------------------------------
SUBROUTINE scan_file_to (keyword,found)
   !-----------------------------------------------------------------------
   !
   USE wannier, ONLY :iun_nnkp
   USE io_global,  ONLY : stdout
   IMPLICIT NONE
   CHARACTER(len=*), intent(in) :: keyword
   logical, intent(out) :: found
   CHARACTER(len=80) :: line1, line2
!
! by uncommenting the following line the file scan restarts every time
! from the beginning thus making the reading independent on the order
! of data-blocks
!   rewind (iun_nnkp)
!
10 CONTINUE
   READ(iun_nnkp,*,end=20) line1, line2
   IF(line1/='begin')  GOTO 10
   IF(line2/=keyword) GOTO 10
   found=.true.
   RETURN
20 found=.false.
   rewind (iun_nnkp)

END SUBROUTINE scan_file_to
!
!-----------------------------------------------------------------------
SUBROUTINE compute_mmn
   !-----------------------------------------------------------------------
   !
   USE io_global,  ONLY : stdout, ionode
   USE kinds,           ONLY: DP
   USE wvfct,           ONLY : nbnd, npw, npwx, igk, g2kin
   USE control_flags,   ONLY : gamma_only
   USE wavefunctions_module, ONLY : evc, psic, psic_nc
   USE fft_base,        ONLY : dffts, dfftp
   USE fft_interfaces,  ONLY : fwfft, invfft
   USE gvecs,         ONLY : nls, nlsm
   USE klist,           ONLY : nkstot, xk
   USE io_files,        ONLY : nwordwfc, iunwfc
   USE gvect,           ONLY : g, ngm, gstart
   USE cell_base,       ONLY : tpiba2, omega, alat, tpiba, at, bg
   USE ions_base,       ONLY : nat, ntyp => nsp, ityp, tau
   USE constants,       ONLY : tpi
   USE uspp,            ONLY : nkb, vkb
   USE uspp_param,      ONLY : upf, nh, lmaxq
   USE becmod,          ONLY : bec_type, becp, calbec, &
                               allocate_bec_type, deallocate_bec_type
   USE mp_global,       ONLY : intra_pool_comm
   USE mp,              ONLY : mp_sum
   USE noncollin_module,ONLY : noncolin, npol
   USE wvfct,           ONLY : ecutwfc
   USE wannier

   IMPLICIT NONE
   !
   INTEGER, EXTERNAL :: find_free_unit
   !
   complex(DP), parameter :: cmplx_i=(0.0_DP,1.0_DP)
   !
   INTEGER :: mmn_tot, ik, ikp, ipol, ib, npwq, i, m, n
   INTEGER :: ikb, jkb, ih, jh, na, nt, ijkb0, ind, nbt
   INTEGER :: ikevc, ikpevcq, s, counter
   COMPLEX(DP), ALLOCATABLE :: phase(:), aux(:), aux2(:), evcq(:,:), &
                               becp2(:,:), Mkb(:,:), aux_nc(:,:)
   real(DP), ALLOCATABLE    :: rbecp2(:,:)
   COMPLEX(DP), ALLOCATABLE :: qb(:,:,:,:), qgm(:)
   real(DP), ALLOCATABLE    :: qg(:), ylm(:,:), dxk(:,:)
   INTEGER, ALLOCATABLE     :: igkq(:)
   COMPLEX(DP)              :: mmn, zdotc, phase1
   real(DP)                 :: arg, g_(3)
   CHARACTER (len=9)        :: cdate,ctime
   CHARACTER (len=60)       :: header
   LOGICAL                  :: any_uspp
   INTEGER                  :: nn,inn,loop,loop2
   LOGICAL                  :: nn_found
   INTEGER                  :: istart,iend

   any_uspp = any(upf(1:ntyp)%tvanp)

   IF(any_uspp .and. noncolin) CALL errore('pw2wannier90',&
       'NCLS calculation not implimented with USP',1)

   ALLOCATE( phase(dffts%nnr), igkq(npwx) )
   ALLOCATE( evcq(npol*npwx,nbnd) )

   IF(noncolin) THEN
      ALLOCATE( aux_nc(npwx,npol) )
   ELSE
      ALLOCATE( aux(npwx) )
   ENDIF

   IF (gamma_only) ALLOCATE(aux2(npwx))

   IF (wan_mode=='library') ALLOCATE(m_mat(num_bands,num_bands,nnb,iknum))

   IF (wan_mode=='standalone') THEN
      iun_mmn = find_free_unit()
      IF (ionode) OPEN (unit=iun_mmn, file=trim(seedname)//".mmn",form='formatted')
      CALL date_and_tim( cdate, ctime )
      header='Created on '//cdate//' at '//ctime
      IF (ionode) THEN
         WRITE (iun_mmn,*) header
         WRITE (iun_mmn,*) nbnd-nexband, iknum, nnb
      ENDIF
   ENDIF

   !
   !   USPP
   !
   !
   IF(any_uspp) THEN
      CALL init_us_1
      CALL allocate_bec_type ( nkb, nbnd, becp )
      IF (gamma_only) THEN
         ALLOCATE ( rbecp2(nkb,nbnd))
      ELSE
         ALLOCATE ( becp2(nkb,nbnd) )
      ENDIF
   ENDIF
   !
   !     qb is  FT of Q(r)
   !
   nbt = nnb * iknum
   !
   ALLOCATE( qg(nbt) )
   ALLOCATE (dxk(3,nbt))
   !
   ind = 0
   DO ik=1,iknum
      DO ib=1,nnb
         ind = ind + 1
         ikp = kpb(ik,ib)
         !
         g_(:) = REAL( g_kpb(:,ik,ib) )
         CALL cryst_to_cart (1, g_, bg, 1)
         dxk(:,ind) = xk(:,ikp) +g_(:) - xk(:,ik)
         qg(ind) = dxk(1,ind)*dxk(1,ind)+dxk(2,ind)*dxk(2,ind)+dxk(3,ind)*dxk(3,ind)
      ENDDO
!      write (stdout,'(i3,12f8.4)')  ik, qg((ik-1)*nnb+1:ik*nnb)
   ENDDO
   !
   !  USPP
   !
   IF(any_uspp) THEN

      ALLOCATE( ylm(nbt,lmaxq*lmaxq), qgm(nbt) )
      ALLOCATE( qb (nkb, nkb, ntyp, nbt) )
      !
      CALL ylmr2 (lmaxq*lmaxq, nbt, dxk, qg, ylm)
      qg(:) = sqrt(qg(:)) * tpiba
      !
      DO nt = 1, ntyp
         IF (upf(nt)%tvanp ) THEN
            DO ih = 1, nh (nt)
               DO jh = 1, nh (nt)
                  CALL qvan2 (nbt, ih, jh, nt, qg, qgm, ylm)
                  qb (ih, jh, nt, 1:nbt) = omega * qgm(1:nbt)
               ENDDO
            ENDDO
         ENDIF
      ENDDO
      !
      DEALLOCATE (qg, qgm, ylm )
      !
   ENDIF
   WRITE (stdout,*) "MMN"
   !
   ALLOCATE( Mkb(nbnd,nbnd) )
   !
   WRITE(stdout,'(a,i8)') ' iknum = ',iknum

   ind = 0
   DO ik=1,iknum
      WRITE (stdout,'(i8)') ik
      ikevc = ik + ikstart - 1
      CALL davcio (evc, nwordwfc, iunwfc, ikevc, -1 )
      CALL gk_sort (xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
      !
      !  USPP
      !
      IF(any_uspp) THEN
         CALL init_us_2 (npw, igk, xk(1,ik), vkb)
         ! below we compute the product of beta functions with |psi>
         CALL calbec (npw, vkb, evc, becp)
      ENDIF
      !
      !
      !do ib=1,nnb(ik)
      DO ib=1,nnb
         ind = ind + 1
         ikp = kpb(ik,ib)
! read wfc at k+b
         ikpevcq = ikp + ikstart - 1
!         if(noncolin) then
!            call davcio (evcq_nc, nwordwfc, iunwfc, ikpevcq, -1 )
!         else
            CALL davcio (evcq, nwordwfc, iunwfc, ikpevcq, -1 )
!         end if
         CALL gk_sort (xk(1,ikp), ngm, g, ecutwfc / tpiba2, npwq, igkq, g2kin)
! compute the phase
         phase(:) = (0.d0,0.d0)
         IF ( ig_(ik,ib)>0) phase( nls(ig_(ik,ib)) ) = (1.d0,0.d0)
         CALL invfft ('Wave', phase, dffts)
         !
         !  USPP
         !
         IF(any_uspp) THEN
            CALL init_us_2 (npwq, igkq, xk(1,ikp), vkb)
            ! below we compute the product of beta functions with |psi>
            IF (gamma_only) THEN
               CALL calbec ( npwq, vkb, evcq, rbecp2 )
            ELSE
               CALL calbec ( npwq, vkb, evcq, becp2 )
            ENDIF
         ENDIF
         !
         !
         Mkb(:,:) = (0.0d0,0.0d0)
         !
         IF (any_uspp) THEN
            ijkb0 = 0
            DO nt = 1, ntyp
               IF ( upf(nt)%tvanp ) THEN
                  DO na = 1, nat
                     !
                     arg = dot_product( dxk(:,ind), tau(:,na) ) * tpi
                     phase1 = cmplx( cos(arg), -sin(arg) ,kind=DP)
                     !
                     IF ( ityp(na) == nt ) THEN
                        DO jh = 1, nh(nt)
                           jkb = ijkb0 + jh
                           DO ih = 1, nh(nt)
                              ikb = ijkb0 + ih
                              !
                              DO m = 1,nbnd
                                 IF (excluded_band(m)) CYCLE
                                 IF (gamma_only) THEN
                                    DO n=1,m ! Mkb(m,n) is symmetric in m and n for gamma_only case
                                       IF (excluded_band(n)) CYCLE
                                       Mkb(m,n) = Mkb(m,n) + &
                                            phase1 * qb(ih,jh,nt,ind) * &
                                            becp%r(ikb,m) * rbecp2(jkb,n)
                                    ENDDO
                                 ELSE
                                    DO n=1,nbnd
                                       IF (excluded_band(n)) CYCLE
                                       Mkb(m,n) = Mkb(m,n) + &
                                            phase1 * qb(ih,jh,nt,ind) * &
                                            conjg( becp%k(ikb,m) ) * becp2(jkb,n)
                                    ENDDO
                                 ENDIF
                              ENDDO ! m
                           ENDDO !ih
                        ENDDO !jh
                        ijkb0 = ijkb0 + nh(nt)
                     ENDIF  !ityp
                  ENDDO  !nat
               ELSE  !tvanp
                  DO na = 1, nat
                     IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
                  ENDDO
               ENDIF !tvanp
            ENDDO !ntyp
         ENDIF ! any_uspp
         !
         !
! loops on bands
         !
         IF (wan_mode=='standalone') THEN
            IF (ionode) WRITE (iun_mmn,'(7i5)') ik, ikp, (g_kpb(ipol,ik,ib), ipol=1,3)
         ENDIF
         !
         DO m=1,nbnd
            IF (excluded_band(m)) CYCLE
            !
            IF(noncolin) THEN
               psic_nc(:,:) = (0.d0, 0.d0)
               DO ipol=1,2!npol
                  istart=(ipol-1)*npwx+1
                  iend=istart+npw-1
                  psic_nc(nls (igk (1:npw) ),ipol ) = evc(istart:iend, m)
                  CALL invfft ('Wave', psic_nc(:,ipol), dffts)
                  psic_nc(1:dffts%nnr,ipol) = psic_nc(1:dffts%nnr,ipol) * &
                                                 phase(1:dffts%nnr)
                  CALL fwfft ('Wave', psic_nc(:,ipol), dffts)
                  aux_nc(1:npwq,ipol) = psic_nc(nls (igkq(1:npwq) ),ipol )
               ENDDO
            ELSE
               psic(:) = (0.d0, 0.d0)
               psic(nls (igk (1:npw) ) ) = evc (1:npw, m)
               IF(gamma_only) psic(nlsm(igk (1:npw) ) ) = conjg(evc (1:npw, m))
               CALL invfft ('Wave', psic, dffts)
               psic(1:dffts%nnr) = psic(1:dffts%nnr) * phase(1:dffts%nnr)
               CALL fwfft ('Wave', psic, dffts)
               aux(1:npwq)  = psic(nls (igkq(1:npwq) ) )
            ENDIF
            IF(gamma_only) THEN
               IF (gstart==2) psic(nlsm(1)) = (0.d0,0.d0)
               aux2(1:npwq) = conjg(psic(nlsm(igkq(1:npwq) ) ) )
            ENDIF
            !
            !  Mkb(m,n) = Mkb(m,n) + \sum_{ijI} qb_{ij}^I * e^-i(b*tau_I)
            !             <psi_m,k1| beta_i,k1 > < beta_j,k2 | psi_n,k2 >
            !
            IF (gamma_only) THEN
               DO n=1,m ! Mkb(m,n) is symmetric in m and n for gamma_only case
                  IF (excluded_band(n)) CYCLE
                  mmn = zdotc (npwq, aux,1,evcq(1,n),1) &
                       + conjg(zdotc(npwq,aux2,1,evcq(1,n),1))
                  CALL mp_sum(mmn, intra_pool_comm)
                  Mkb(m,n) = mmn + Mkb(m,n)
                  IF (m/=n) Mkb(n,m) = Mkb(m,n) ! fill other half of matrix by symmetry
               ENDDO
            ELSEIF(noncolin) THEN
               DO n=1,nbnd
                  IF (excluded_band(n)) CYCLE
                  mmn=(0.d0, 0.d0)
!                  do ipol=1,2
!                     mmn = mmn+zdotc (npwq, aux_nc(1,ipol),1,evcq_nc(1,ipol,n),1)
                  mmn = mmn + zdotc (npwq, aux_nc(1,1),1,evcq(1,n),1) &
                       + zdotc (npwq, aux_nc(1,2),1,evcq(npwx+1,n),1)
!                  end do
                  CALL mp_sum(mmn, intra_pool_comm)
                  Mkb(m,n) = mmn + Mkb(m,n)
               ENDDO
            ELSE
               DO n=1,nbnd
                  IF (excluded_band(n)) CYCLE
                  mmn = zdotc (npwq, aux,1,evcq(1,n),1)
                  CALL mp_sum(mmn, intra_pool_comm)
                  Mkb(m,n) = mmn + Mkb(m,n)
               ENDDO
            ENDIF
         ENDDO   ! m


         DO n=1,nbnd
            IF (excluded_band(n)) CYCLE
            DO m=1,nbnd
               IF (excluded_band(m)) CYCLE
               IF (wan_mode=='standalone') THEN
                  IF (ionode) WRITE (iun_mmn,'(2f18.12)') Mkb(m,n)
               ELSEIF (wan_mode=='library') THEN
                  m_mat(m,n,ib,ik)=Mkb(m,n)  !! this is wrong!
               ELSE
                  CALL errore('compute_mmn',' value of wan_mode not recognised',1)
               ENDIF
            ENDDO
         ENDDO

      ENDDO !ib
   ENDDO  !ik


   IF (ionode .and. wan_mode=='standalone') CLOSE (iun_mmn)

   IF (gamma_only) DEALLOCATE(aux2)
   DEALLOCATE (Mkb, dxk, phase, igkq)
   IF(noncolin) THEN
      DEALLOCATE(aux_nc)
   ELSE
      DEALLOCATE(aux)
   ENDIF
   DEALLOCATE(evcq)

   IF(any_uspp) THEN
      DEALLOCATE (  qb)
      CALL deallocate_bec_type (becp)
      IF (gamma_only) THEN
          DEALLOCATE (rbecp2)
       ELSE
          DEALLOCATE (becp2)
       ENDIF
    ENDIF
!
   WRITE(stdout,*)
   WRITE(stdout,*) ' MMN calculated'

   RETURN
END SUBROUTINE compute_mmn

!-----------------------------------------------------------------------
SUBROUTINE compute_spin
   !-----------------------------------------------------------------------
   !
   USE io_global,  ONLY : stdout, ionode
   USE kinds,           ONLY: DP
   USE wvfct,           ONLY : nbnd, npw, npwx, igk, g2kin
   USE control_flags,   ONLY : gamma_only
   USE wavefunctions_module, ONLY : evc, psic, psic_nc
   USE fft_base,        ONLY : dffts, dfftp
   USE fft_interfaces,  ONLY : fwfft, invfft
   USE gvecs,         ONLY : nls, nlsm
   USE klist,           ONLY : nkstot, xk
   USE io_files,        ONLY : nwordwfc, iunwfc
   USE gvect,           ONLY : g, ngm, gstart
   USE cell_base,       ONLY : tpiba2, omega, alat, tpiba, at, bg
   USE ions_base,       ONLY : nat, ntyp => nsp, ityp, tau
   USE constants,       ONLY : tpi
   USE uspp,            ONLY : nkb, vkb
   USE uspp_param,      ONLY : upf, nh, lmaxq
   USE becmod,          ONLY : bec_type, becp, calbec, &
                               allocate_bec_type, deallocate_bec_type
   USE mp_global,       ONLY : intra_pool_comm
   USE mp,              ONLY : mp_sum
   USE noncollin_module,ONLY : noncolin, npol
   USE wvfct,           ONLY : ecutwfc
   USE wannier
   ! begin change Lopez, Thonhauser, Souza
   USE mp,              ONLY : mp_barrier
   USE scf,             ONLY : vrs, vltot, v, kedtau
   USE gvecs,           ONLY : doublegrid
   USE lsda_mod,        ONLY : nspin
   USE constants,       ONLY : rytoev

   IMPLICIT NONE
   !
   INTEGER, EXTERNAL :: find_free_unit
   !
   complex(DP), parameter :: cmplx_i=(0.0_DP,1.0_DP)
   !
   INTEGER :: mmn_tot, ik, ikp, ipol, ib, npwq, i, m, n
   INTEGER :: ikb, jkb, ih, jh, na, nt, ijkb0, ind, nbt
   INTEGER :: ikevc, ikpevcq, s, counter
   COMPLEX(DP)              :: mmn, zdotc, phase1
   real(DP)                 :: arg, g_(3)
   CHARACTER (len=9)        :: cdate,ctime
   CHARACTER (len=60)       :: header
   LOGICAL                  :: any_uspp
   INTEGER                  :: nn,inn,loop,loop2
   LOGICAL                  :: nn_found
   INTEGER                  :: istart,iend
   COMPLEX(DP)              :: sigma_x,sigma_y,sigma_z,cdum1,cdum2
   complex(DP), allocatable :: spn(:,:)

   any_uspp = any(upf(1:ntyp)%tvanp)

   IF(any_uspp .and. noncolin) CALL errore('pw2wannier90',&
       'NCLS calculation not implimented with USP',1)

   if (write_spn) allocate(spn(3,(num_bands*(num_bands+1))/2))

!ivo
! not sure this is really needed
   if((write_spn.or.write_uhu.or.write_uIu).and.wan_mode=='library')&
        call errore('pw2wannier90',&
        'write_spn, write_uhu, and write_uIu not meant to work library mode',1)
!endivo

   IF(write_spn.and.noncolin) THEN
      IF (ionode) then
         iun_spn = find_free_unit()
         CALL date_and_tim( cdate, ctime )
         header='Created on '//cdate//' at '//ctime
         if(spn_formatted) then
            OPEN (unit=iun_spn, file=trim(seedname)//".spn",form='formatted')
            WRITE (iun_spn,*) header !ivo
            WRITE (iun_spn,*) nbnd-nexband,iknum
         else
            OPEN (unit=iun_spn, file=trim(seedname)//".spn",form='unformatted')
            WRITE (iun_spn) header !ivo
            WRITE (iun_spn) nbnd-nexband,iknum
         endif
      ENDIF
   ENDIF
   !
   WRITE(stdout,'(a,i8)') ' iknum = ',iknum

   ind = 0
   DO ik=1,iknum
      WRITE (stdout,'(i8)') ik
      ikevc = ik + ikstart - 1
      CALL davcio (evc, nwordwfc, iunwfc, ikevc, -1 )
      CALL gk_sort (xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
      IF(write_spn.and.noncolin) THEN
         counter=0
         DO m=1,nbnd
            if(excluded_band(m)) cycle !ivo
            DO n=1,m
               if(excluded_band(n)) cycle !ivo
               cdum1=ZDOTC(npw,evc(1,n),1,evc(npwx+1,m),1)
               call mp_sum(cdum1,intra_pool_comm)
               cdum2=ZDOTC(npw,evc(npwx+1,n),1,evc(1,m),1)
               call mp_sum(cdum2,intra_pool_comm)
               sigma_x=cdum1+cdum2
               sigma_y=cmplx_i*(cdum2-cdum1)
               sigma_z=ZDOTC(npw,evc(1,n),1,evc(1,m),1)&
                    -ZDOTC(npw,evc(npwx+1,n),1,evc(npwx+1,m),1)
               call mp_sum(sigma_z,intra_pool_comm)
               counter=counter+1
               spn(1,counter)=sigma_x
               spn(2,counter)=sigma_y
               spn(3,counter)=sigma_z
            ENDDO
         ENDDO
         if(ionode) then ! root node for i/o
            if(spn_formatted) then ! slow formatted way
               counter=0
               do m=1,num_bands
                  do n=1,m
                     counter=counter+1
                     do s=1,3
                         write(iun_spn,'(2es26.16)') spn(s,counter)
                      enddo
                   enddo
                enddo
             else ! fast unformatted way
                write(iun_spn) ((spn(s,m),s=1,3),m=1,((num_bands*(num_bands+1))/2))
             endif
          endif ! end of root activity


      ENDIF

   end DO

   IF (ionode .and. write_spn .and. noncolin) CLOSE (iun_spn)

   if(write_spn.and.noncolin) deallocate(spn)

   WRITE(stdout,*)
   WRITE(stdout,*) ' SPIN calculated'

   RETURN
END SUBROUTINE compute_spin

!-----------------------------------------------------------------------
SUBROUTINE compute_orb
   !-----------------------------------------------------------------------
   !
   USE io_global,  ONLY : stdout, ionode
   USE kinds,           ONLY: DP
   USE wvfct,           ONLY : nbnd, npw, npwx, igk, g2kin
   USE control_flags,   ONLY : gamma_only
   USE wavefunctions_module, ONLY : evc, psic, psic_nc
   USE fft_base,        ONLY : dffts, dfftp
   USE fft_interfaces,  ONLY : fwfft, invfft
   USE gvecs,         ONLY : nls, nlsm
   USE klist,           ONLY : nkstot, xk
   USE io_files,        ONLY : nwordwfc, iunwfc
   USE gvect,           ONLY : g, ngm, gstart
   USE cell_base,       ONLY : tpiba2, omega, alat, tpiba, at, bg
   USE ions_base,       ONLY : nat, ntyp => nsp, ityp, tau
   USE constants,       ONLY : tpi
   USE uspp,            ONLY : nkb, vkb
   USE uspp_param,      ONLY : upf, nh, lmaxq
   USE becmod,          ONLY : bec_type, becp, calbec, &
                               allocate_bec_type, deallocate_bec_type
   USE mp_global,       ONLY : intra_pool_comm
   USE mp,              ONLY : mp_sum
   USE noncollin_module,ONLY : noncolin, npol
   USE wvfct,           ONLY : ecutwfc
   USE wannier
   ! begin change Lopez, Thonhauser, Souza
   USE mp,              ONLY : mp_barrier
   USE scf,             ONLY : vrs, vltot, v, kedtau
   USE gvecs,           ONLY : doublegrid
   USE lsda_mod,        ONLY : nspin
   USE constants,       ONLY : rytoev

   IMPLICIT NONE
   !
   INTEGER, EXTERNAL :: find_free_unit
   !
   complex(DP), parameter :: cmplx_i=(0.0_DP,1.0_DP)
   !
   INTEGER :: mmn_tot, ik, ikp, ipol, ib, npwq, i, m, n
   INTEGER :: ikb, jkb, ih, jh, na, nt, ijkb0, ind, nbt
   INTEGER :: ikevc, ikpevcq, s, counter
   COMPLEX(DP), ALLOCATABLE :: phase(:), aux(:), aux2(:), evcq(:,:), &
                               becp2(:,:), Mkb(:,:), aux_nc(:,:)
   real(DP), ALLOCATABLE    :: rbecp2(:,:)
   COMPLEX(DP), ALLOCATABLE :: qb(:,:,:,:), qgm(:)
   real(DP), ALLOCATABLE    :: qg(:), ylm(:,:)
   INTEGER, ALLOCATABLE     :: igkq(:)
   COMPLEX(DP)              :: mmn, zdotc, phase1
   real(DP)                 :: arg, g_(3)
   CHARACTER (len=9)        :: cdate,ctime
   CHARACTER (len=60)       :: header
   LOGICAL                  :: any_uspp
   INTEGER                  :: nn,inn,loop,loop2
   LOGICAL                  :: nn_found
   INTEGER                  :: istart,iend
   ! begin change Lopez, Thonhauser, Souza
   COMPLEX(DP)              :: sigma_x,sigma_y,sigma_z,cdum1,cdum2
   integer                  :: npw_b1, npw_b2, i_b1, i_b2, ikp_b1, ikp_b2
   integer, allocatable     :: igk_b1(:), igk_b2(:)
   complex(DP), allocatable :: evc_b1(:,:),evc_b2(:,:),evc_aux(:,:),H_evc(:,:)
   complex(DP), allocatable :: uHu(:,:),uIu(:,:),spn(:,:)
   ! end change Lopez, Thonhauser, Souza

   any_uspp = any(upf(1:ntyp)%tvanp)

   IF(any_uspp .and. noncolin) CALL errore('pw2wannier90',&
       'NCLS calculation not implimented with USP',1)

   ALLOCATE( phase(dffts%nnr) )
   ALLOCATE( evcq(npol*npwx,nbnd) )

   IF(noncolin) THEN
      ALLOCATE( aux_nc(npwx,npol) )
   ELSE
      ALLOCATE( aux(npwx) )
   ENDIF

   IF (gamma_only) ALLOCATE(aux2(npwx))

   IF (wan_mode=='library') ALLOCATE(m_mat(num_bands,num_bands,nnb,iknum))

   if (write_uHu) allocate(uhu(num_bands,num_bands))
   if (write_uIu) allocate(uIu(num_bands,num_bands))


!ivo
! not sure this is really needed
   if((write_uhu.or.write_uIu).and.wan_mode=='library')&
        call errore('pw2wannier90',&
        'write_uhu, and write_uIu not meant to work library mode',1)
!endivo


   !
   !
   ! begin change Lopez, Thonhauser, Souza
   !
   !====================================================================
   !
   ! The following code was inserted by Timo Thonhauser, Ivo Souza, and
   ! Graham Lopez in order to calculate the matrix elements
   ! <u_n(q+b1)|H(q)|u_m(q+b2)> necessary for the Wannier interpolation
   ! of the orbital magnetization
   !
   !====================================================================
   !
   !
   !
   if(write_uHu.or.write_uIu) then !ivo
     !
     if(gamma_only) call errore('pw2wannier90',&
      'write_uHu and write_uIu not yet implemented for gamma_only case',1) !ivo
     if(any_uspp) call errore('pw2wannier90',&
      'write_uHu and write_uIu not yet implemented with USP',1) !ivo
     !
     !
     allocate(igk_b1(npwx),igk_b2(npwx),evc_b1(npol*npwx,nbnd),&
          evc_b2(npol*npwx,nbnd),&
          evc_aux(npol*npwx,nbnd))
     !
     if(write_uHu) then
        allocate(H_evc(npol*npwx,nbnd))
        write(stdout,*)
        write(stdout,*) ' -----------------'
        write(stdout,*) ' *** Compute  uHu '
        write(stdout,*) ' -----------------'
        write(stdout,*)
        iun_uhu = find_free_unit()
        if (ionode) then
           CALL date_and_tim( cdate, ctime )
           header='Created on '//cdate//' at '//ctime
           if(uHu_formatted) then
              open  (unit=iun_uhu, file=TRIM(seedname)//".uHu",form='FORMATTED')
              write (iun_uhu,*) header
              write (iun_uhu,*) nbnd, iknum, nnb
           else
              open  (unit=iun_uhu, file=TRIM(seedname)//".uHu",form='UNFORMATTED')
              write (iun_uhu) header
              write (iun_uhu) nbnd, iknum, nnb
           endif
        endif
     endif
     if(write_uIu) then
        write(stdout,*)
        write(stdout,*) ' -----------------'
        write(stdout,*) ' *** Compute  uIu '
        write(stdout,*) ' -----------------'
        write(stdout,*)
        iun_uIu = find_free_unit()
        if (ionode) then
           CALL date_and_tim( cdate, ctime )
           header='Created on '//cdate//' at '//ctime
           if(uIu_formatted) then
              open  (unit=iun_uIu, file=TRIM(seedname)//".uIu",form='FORMATTED')
              write (iun_uIu,*) header
              write (iun_uIu,*) nbnd, iknum, nnb
           else
              open  (unit=iun_uIu, file=TRIM(seedname)//".uIu",form='UNFORMATTED')
              write (iun_uIu) header
              write (iun_uIu) nbnd, iknum, nnb
           endif
        endif
     endif

     CALL set_vrs(vrs,vltot,v%of_r,kedtau,v%kin_r,dfftp%nnr,nspin,doublegrid)
     call allocate_bec_type ( nkb, nbnd, becp )

     write(stdout,'(a,i8)') ' iknum = ',iknum
     do ik = 1, iknum ! loop over k points
        !
        write (stdout,'(i8)') ik
        !
        ! sort the wfc at k and set up stuff for h_psi
        CALL gk_sort(xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
        CALL init_us_2(npw,igk,xk(1,ik),vkb)
        !
        ! compute  " H | u_n,k+b2 > "
        !
        do i_b2 = 1, nnb ! nnb = # of nearest neighbors
           !
           ! read wfc at k+b2
           ikp_b2 = kpb(ik,i_b2) ! for kpoint 'ik', index of neighbor 'i_b2'
           !
!           call davcio  (evc_b2, nwordwfc, iunwfc, ikp_b2, -1 ) !ivo
           call davcio  (evc_b2, nwordwfc, iunwfc, ikp_b2+ikstart-1, -1 ) !ivo
!           call gk_sort (xk(1,ikp_b2), ngm, g, ecutwfc/tpiba2, npw_b1, igk_b1, g2kin) !ivo
           call gk_sort (xk(1,ikp_b2), ngm, g, ecutwfc/tpiba2, npw_b2, igk_b2, g2kin) !ivo
           !
           ! compute the phase
           phase(:) = ( 0.0D0, 0.0D0 )
           if (ig_(ik,i_b2)>0) phase( nls(ig_(ik,i_b2)) ) = ( 1.0D0, 0.0D0 )
           call invfft('Wave', phase, dffts)
           !
           ! loop on bands
           evc_aux = ( 0.0D0, 0.0D0 )
           do n = 1, nbnd !ivo replaced dummy m --> n everywhere on this do loop, for consistency w/ band indices in comments
              if (excluded_band(n)) cycle
              if(noncolin) then
                 psic_nc = ( 0.0D0, 0.0D0 ) !ivo
                 do ipol = 1, 2
!                    psic_nc = ( 0.0D0, 0.0D0 ) !ivo
                    istart=(ipol-1)*npwx+1
                    iend=istart+npw_b2-1 !ivo npw_b1 --> npw_b2
                    psic_nc(nls (igk_b2(1:npw_b2) ),ipol ) = evc_b2(istart:iend, n) !ivo igk_b1, npw_b1 --> igk_b2, npw_b2
                    ! multiply by phase in real space - '1' unless neighbor is in a bordering BZ
                    call invfft ('Wave', psic_nc(:,ipol), dffts)
                    psic_nc(1:dffts%nnr,ipol) = psic_nc(1:dffts%nnr,ipol) * conjg(phase(1:dffts%nnr))
                    call fwfft ('Wave', psic_nc(:,ipol), dffts)
                    ! save the result
                    iend=istart+npw-1
                    evc_aux(istart:iend,n) = psic_nc(nls (igk(1:npw) ),ipol )
                 end do
              else ! this is modeled after the pre-existing code at 1162
                 psic = ( 0.0D0, 0.0D0 )
                 ! Graham, changed npw --> npw_b2 on RHS. Do you agree?!
                 psic(nls (igk_b2(1:npw_b2) ) ) = evc_b2(1:npw_b2, n)
                 call invfft ('Wave', psic, dffts)
                 psic(1:dffts%nnr) = psic(1:dffts%nnr) * conjg(phase(1:dffts%nnr))
                 call fwfft ('Wave', psic, dffts)
                 evc_aux(1:npw,n) = psic(nls (igk(1:npw) ) )
              end if
           end do !n

           if(write_uHu) then !ivo
              !
              ! gk_sort overwrites the kinetic energy - recalculate at ik
              g2kin(1:npw) = ( ( xk(1,ik) + g(1,igk(1:npw)) )**2 + &
                   ( xk(2,ik) + g(2,igk(1:npw)) )**2 + &
                   ( xk(3,ik) + g(3,igk(1:npw)) )**2 ) * tpiba2
              !
              CALL h_psi(npwx, npw, nbnd, evc_aux, H_evc)
           endif
           !
           ! compute  " < u_m,k+b1 | "
           !
           do i_b1 = 1, nnb
              !
              ! read wfc at k+b1 !ivo replaced k+b2 --> k+b1
              ikp_b1 = kpb(ik,i_b1)
!              call davcio  (evc_b1, nwordwfc, iunwfc, ikp_b1, -1 ) !ivo
              call davcio  (evc_b1, nwordwfc, iunwfc, ikp_b1+ikstart-1, -1 ) !ivo

!              call gk_sort (xk(1,ikp_b1), ngm, g, ecutwfc/tpiba2, npw_b2, igk_b2, g2kin) !ivo
              call gk_sort (xk(1,ikp_b1), ngm, g, ecutwfc/tpiba2, npw_b1, igk_b1, g2kin) !ivo
              !
              ! compute the phase
              phase(:) = ( 0.0D0, 0.0D0 )
              if (ig_(ik,i_b1)>0) phase( nls(ig_(ik,i_b1)) ) = ( 1.0D0, 0.0D0 )
              !call cft3s (phase, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, +2)
              call invfft('Wave', phase, dffts)
              !
              ! loop on bands
              do m = 1, nbnd
                 if (excluded_band(m)) cycle
                 if(noncolin) then
                    aux_nc  = ( 0.0D0, 0.0D0 )
                    psic_nc = ( 0.0D0, 0.0D0 ) !ivo
                    do ipol = 1, 2
!                      psic_nc = ( 0.0D0, 0.0D0 ) !ivo
                       istart=(ipol-1)*npwx+1
                       iend=istart+npw_b1-1  !ivo npw_b2 --> npw_b1
                       psic_nc(nls (igk_b1(1:npw_b1) ),ipol ) = evc_b1(istart:iend, m) !ivo igk_b2,npw_b2 --> igk_b1,npw_b1
                       ! multiply by phase in real space - '1' unless neighbor is in a different BZ
                       call invfft ('Wave', psic_nc(:,ipol), dffts)
                       !psic_nc(1:nrxxs,ipol) = psic_nc(1:nrxxs,ipol) * conjg(phase(1:nrxxs))
                       psic_nc(1:dffts%nnr,ipol) = psic_nc(1:dffts%nnr,ipol) * conjg(phase(1:dffts%nnr))
                       call fwfft ('Wave', psic_nc(:,ipol), dffts)
                       ! save the result
                       aux_nc(1:npw,ipol) = psic_nc(nls (igk(1:npw) ),ipol )
                    end do
                 else ! this is modeled after the pre-existing code at 1162
                    aux  = ( 0.0D0 )
                    psic = ( 0.0D0, 0.0D0 )
                    ! Graham, changed npw --> npw_b1 on RHS. Do you agree?!
                    psic(nls (igk_b1(1:npw_b1) ) ) = evc_b1(1:npw_b1, m) !ivo igk_b2 --> igk_b1
                    call invfft ('Wave', psic, dffts)
                    !psic(1:nrxxs) = psic(1:nrxxs) * conjg(phase(1:nrxxs))
                    psic(1:dffts%nnr) = psic(1:dffts%nnr) * conjg(phase(1:dffts%nnr))
                    call fwfft ('Wave', psic, dffts)
                    aux(1:npw) = psic(nls (igk(1:npw) ) )
                 end if

                !
                !
                if(write_uHu) then !ivo
                   do n = 1, nbnd  ! loop over bands of already computed ket
                      if (excluded_band(n)) cycle
                      if(noncolin) then
                         mmn = ZDOTC (npw, aux_nc(1,1),1,H_evc(1,n),1) + &
                              ZDOTC (npw, aux_nc(1,2),1,H_evc(1+npwx,n),1)
                      else
                         mmn = ZDOTC (npw, aux,1,H_evc(1,n),1)
                      end if
                      mmn = mmn * rytoev ! because wannier90 works in eV
                      call mp_sum(mmn, intra_pool_comm)
!                      if (ionode) write (iun_uhu) mmn
                      uHu(n,m)=mmn
                      !
                   end do !n
                endif
                if(write_uIu) then !ivo
                   do n = 1, nbnd  ! loop over bands of already computed ket
                      if (excluded_band(n)) cycle
                      if(noncolin) then
                         mmn = ZDOTC (npw, aux_nc(1,1),1,evc_aux(1,n),1) + &
                              ZDOTC (npw, aux_nc(1,2),1,evc_aux(1+npwx,n),1)
                      else
                         mmn = ZDOTC (npw, aux,1,evc_aux(1,n),1)
                      end if
                      call mp_sum(mmn, intra_pool_comm)
!                      if (ionode) write (iun_uIu) mmn
                      uIu(n,m)=mmn
                      !
                   end do !n
                endif
                !
             end do ! m = 1, nbnd
             if (ionode) then  ! write the files out to disk
                if(write_uhu) then
                   if(uHu_formatted) then ! slow bulky way for transferable files
                      do n=1,num_bands
                         do m=1,num_bands
                            write(iun_uHu,'(2ES20.10)') uHu(m,n)
                         enddo
                      enddo
                   else  ! the fast way
                      write(iun_uHu) ((uHu(n,m),n=1,num_bands),m=1,num_bands)
                   endif
                endif
                if(write_uiu) then
                   if(uIu_formatted) then ! slow bulky way for transferable files
                      do n=1,num_bands
                         do m=1,num_bands
                            write(iun_uIu,'(2ES20.10)') uIu(m,n)
                         enddo
                      enddo
                   else ! the fast way
                      write(iun_uIu) ((uIu(n,m),n=1,num_bands),m=1,num_bands)
                   endif
                endif
             endif ! end of io
          end do ! i_b1
       end do ! i_b2
    end do ! ik
    !
    deallocate(igk_b1,igk_b2,evc_b1,evc_b2,evc_aux)
    if(write_uHu) then
       deallocate(H_evc)
       deallocate(uHu)
    end if
    if(write_uIu) deallocate(uIu)
    if (ionode.and.write_uHu) close (iun_uhu) !ivo
    if (ionode.and.write_uIu) close (iun_uIu) !ivo
    !
 else
    if(.not.write_uHu) then
       write(stdout,*)
       write(stdout,*) ' -------------------------------'
       write(stdout,*) ' *** uHu matrix is not computed '
       write(stdout,*) ' -------------------------------'
       write(stdout,*)
    endif
    if(.not.write_uIu) then
       write(stdout,*)
       write(stdout,*) ' -------------------------------'
       write(stdout,*) ' *** uIu matrix is not computed '
       write(stdout,*) ' -------------------------------'
       write(stdout,*)
    endif
 end if
   !
   !
   !
   !
   !
   !
   !====================================================================
   !
   ! END_m_orbit
   !
   !====================================================================
   !
   ! end change Lopez, Thonhauser, Souza
   !
   !
   !

   IF (gamma_only) DEALLOCATE(aux2)
   DEALLOCATE (phase)
   IF(noncolin) THEN
      DEALLOCATE(aux_nc)
   ELSE
      DEALLOCATE(aux)
   ENDIF
   DEALLOCATE(evcq)
   if(write_spn.and.noncolin) deallocate(spn)

   IF(any_uspp) THEN
      DEALLOCATE (  qb)
      CALL deallocate_bec_type (becp)
      IF (gamma_only) THEN
          DEALLOCATE (rbecp2)
       ELSE
          DEALLOCATE (becp2)
       ENDIF
    ENDIF
!
   WRITE(stdout,*)
   WRITE(stdout,*) ' uHu calculated'

   RETURN
END SUBROUTINE compute_orb
!
!-----------------------------------------------------------------------
SUBROUTINE compute_amn
   !-----------------------------------------------------------------------
   !
   USE io_global,  ONLY : stdout, ionode
   USE kinds,           ONLY : DP
   USE klist,           ONLY : nkstot, xk
   USE wvfct,           ONLY : nbnd, npw, npwx, igk, g2kin
   USE control_flags,   ONLY : gamma_only
   USE wavefunctions_module, ONLY : evc
   USE io_files,        ONLY : nwordwfc, iunwfc
   USE gvect,           ONLY : g, ngm, gstart
   USE cell_base,       ONLY : tpiba2
   USE uspp,            ONLY : nkb, vkb
   USE becmod,          ONLY : bec_type, becp, calbec, &
                               allocate_bec_type, deallocate_bec_type
   USE wannier
   USE ions_base,       ONLY : nat, ntyp => nsp, ityp, tau
   USE uspp_param,      ONLY : upf
   USE mp_global,       ONLY : intra_pool_comm
   USE mp,              ONLY : mp_sum
   USE noncollin_module,ONLY : noncolin, npol
   USE wvfct,           ONLY : ecutwfc
   USE constants,       ONLY : eps6

   IMPLICIT NONE
   !
   INTEGER, EXTERNAL :: find_free_unit
   !
   COMPLEX(DP) :: amn, zdotc,amn_tmp,fac(2)
   real(DP):: ddot
   COMPLEX(DP), ALLOCATABLE :: sgf(:,:)
   INTEGER :: ik, ibnd, ibnd1, iw,i, ikevc, nt, ipol
   CHARACTER (len=9)  :: cdate,ctime
   CHARACTER (len=60) :: header
   LOGICAL            :: any_uspp, opnd, exst,spin_z_pos, spin_z_neg
   INTEGER            :: istart

   !nocolin: we have half as many projections g(r) defined as wannier
   !         functions. We project onto (1,0) (ie up spin) and then onto
   !         (0,1) to obtain num_wann projections. jry


   !call read_gf_definition.....>   this is done at the beging

   any_uspp =any (upf(1:ntyp)%tvanp)

   IF(any_uspp .and. noncolin) CALL errore('pw2wannier90',&
       'NCLS calculation not implimented with USP',1)

   IF (wan_mode=='library') ALLOCATE(a_mat(num_bands,n_wannier,iknum))

   IF (wan_mode=='standalone') THEN
      iun_amn = find_free_unit()
      IF (ionode) OPEN (unit=iun_amn, file=trim(seedname)//".amn",form='formatted')
   ENDIF

   WRITE (stdout,*) "AMN"

   IF (wan_mode=='standalone') THEN
      CALL date_and_tim( cdate, ctime )
      header='Created on '//cdate//' at '//ctime
      IF (ionode) THEN
         WRITE (iun_amn,*) header
         WRITE (iun_amn,*) nbnd-nexband,  iknum, n_wannier
      ENDIF
   ENDIF
   !
   ALLOCATE( sgf(npwx,n_proj))
   !
   IF (any_uspp) THEN
      CALL allocate_bec_type ( nkb, n_wannier, becp)
      CALL init_us_1
   ENDIF
   !
   WRITE(stdout,'(a,i8)') ' iknum = ',iknum
   DO ik=1,iknum
      WRITE (stdout,'(i8)') ik
      ikevc = ik + ikstart - 1
!      if(noncolin) then
!         call davcio (evc_nc, nwordwfc, iunwfc, ikevc, -1 )
!      else
         CALL davcio (evc, nwordwfc, iunwfc, ikevc, -1 )
!      end if
      CALL gk_sort (xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
      CALL generate_guiding_functions(ik)   ! they are called gf(npw,n_proj)
      !
      !  USPP
      !
      IF(any_uspp) THEN
         CALL init_us_2 (npw, igk, xk (1, ik), vkb)
         ! below we compute the product of beta functions with trial func.
         IF (gamma_only) THEN
            CALL calbec ( npw, vkb, gf, becp, n_proj )
         ELSE
            CALL calbec ( npw, vkb, gf, becp, n_proj )
         ENDIF
         ! and we use it for the product S|trial_func>
         CALL s_psi (npwx, npw, n_proj, gf, sgf)
      ELSE
         sgf(:,:) = gf(:,:)
      ENDIF
      !
      IF(noncolin) THEN
         if(old_spinor_proj) then
            ! we do the projection as g(r)*a(r) and g(r)*b(r)
            DO ipol=1,npol
               istart = (ipol-1)*npwx + 1
               DO iw = 1,n_proj
                  ibnd1 = 0
                  DO ibnd = 1,nbnd
                     IF (excluded_band(ibnd)) CYCLE
                     amn=(0.0_dp,0.0_dp)
                     !                  amn = zdotc(npw,evc_nc(1,ipol,ibnd),1,sgf(1,iw),1)
                     amn = zdotc(npw,evc(istart,ibnd),1,sgf(1,iw),1)
                     CALL mp_sum(amn, intra_pool_comm)
                     ibnd1=ibnd1+1
                     IF (wan_mode=='standalone') THEN
                        IF (ionode) WRITE(iun_amn,'(3i5,2f18.12)') ibnd1, iw+n_proj*(ipol-1), ik, amn
                     ELSEIF (wan_mode=='library') THEN
                        a_mat(ibnd1,iw+n_proj*(ipol-1),ik) = amn
                     ELSE
                        CALL errore('compute_amn',' value of wan_mode not recognised',1)
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            DO iw = 1,n_proj
               spin_z_pos=.false.;spin_z_neg=.false.
               ! detect if spin quantisation axis is along z
               if((abs(spin_qaxis(1,iw)-0.0d0)<eps6).and.(abs(spin_qaxis(2,iw)-0.0d0)<eps6) &
                    .and.(abs(spin_qaxis(3,iw)-1.0d0)<eps6)) then
                  spin_z_pos=.true.
               elseif(abs(spin_qaxis(1,iw)-0.0d0)<eps6.and.abs(spin_qaxis(2,iw)-0.0d0)<eps6 &
                    .and.abs(spin_qaxis(3,iw)+1.0d0)<eps6) then
                  spin_z_neg=.true.
               endif
               if(spin_z_pos .or. spin_z_neg) then
                  ibnd1 = 0
                  DO ibnd = 1,nbnd
                     IF (excluded_band(ibnd)) CYCLE
                     if(spin_z_pos) then
                        ipol=(3-spin_eig(iw))/2
                     else
                        ipol=(3+spin_eig(iw))/2
                     endif
                     istart = (ipol-1)*npwx + 1
                     amn=(0.0_dp,0.0_dp)
                     amn = zdotc(npw,evc(istart,ibnd),1,sgf(1,iw),1)
                     CALL mp_sum(amn, intra_pool_comm)
                     ibnd1=ibnd1+1
                     IF (wan_mode=='standalone') THEN
                        IF (ionode) WRITE(iun_amn,'(3i5,2f18.12)') ibnd1, iw, ik, amn
                     ELSEIF (wan_mode=='library') THEN
                        a_mat(ibnd1,iw+n_proj*(ipol-1),ik) = amn
                     ELSE
                        CALL errore('compute_amn',' value of wan_mode not recognised',1)
                     ENDIF
                  ENDDO
               else
                  ! general routine
                  ! for quantisation axis (a,b,c)
                  ! 'up'    eigenvector is 1/sqrt(1+c) [c+1,a+ib]
                  ! 'down'  eigenvector is 1/sqrt(1-c) [c-1,a+ib]
                  if(spin_eig(iw)==1) then
                     fac(1)=(1.0_dp/sqrt(1+spin_qaxis(3,iw)))*(spin_qaxis(3,iw)+1)*cmplx(1.0d0,0.0d0,dp)
                     fac(2)=(1.0_dp/sqrt(1+spin_qaxis(3,iw)))*cmplx(spin_qaxis(1,iw),spin_qaxis(2,iw),dp)
                  else
                     fac(1)=(1.0_dp/sqrt(1+spin_qaxis(3,iw)))*(spin_qaxis(3,iw))*cmplx(1.0d0,0.0d0,dp)
                     fac(2)=(1.0_dp/sqrt(1-spin_qaxis(3,iw)))*cmplx(spin_qaxis(1,iw),spin_qaxis(2,iw),dp)
                  endif
                  ibnd1 = 0
                  DO ibnd = 1,nbnd
                     IF (excluded_band(ibnd)) CYCLE
                     amn=(0.0_dp,0.0_dp)
                     DO ipol=1,npol
                        istart = (ipol-1)*npwx + 1
                        amn_tmp=(0.0_dp,0.0_dp)
                        amn_tmp = zdotc(npw,evc(istart,ibnd),1,sgf(1,iw),1)
                        CALL mp_sum(amn_tmp, intra_pool_comm)
                        amn=amn+fac(ipol)*amn_tmp
                     enddo
                     ibnd1=ibnd1+1
                     IF (wan_mode=='standalone') THEN
                        IF (ionode) WRITE(iun_amn,'(3i5,2f18.12)') ibnd1, iw, ik, amn
                     ELSEIF (wan_mode=='library') THEN
                           a_mat(ibnd1,iw+n_proj*(ipol-1),ik) = amn
                        ELSE
                           CALL errore('compute_amn',' value of wan_mode not recognised',1)
                        ENDIF
                     ENDDO
                  endif
               end do
            endif
         ELSE ! scalar wavefunctions
            DO iw = 1,n_proj
               ibnd1 = 0
               DO ibnd = 1,nbnd
                  IF (excluded_band(ibnd)) CYCLE
                  IF (gamma_only) THEN
                     amn = 2.0_dp*ddot(2*npw,evc(1,ibnd),1,sgf(1,iw),1)
                     IF (gstart==2) amn = amn - real(conjg(evc(1,ibnd))*sgf(1,iw))
                  ELSE
                     amn = zdotc(npw,evc(1,ibnd),1,sgf(1,iw),1)
                  ENDIF
                  CALL mp_sum(amn, intra_pool_comm)
                  ibnd1=ibnd1+1
                  IF (wan_mode=='standalone') THEN
                     IF (ionode) WRITE(iun_amn,'(3i5,2f18.12)') ibnd1, iw, ik, amn
                  ELSEIF (wan_mode=='library') THEN
                     a_mat(ibnd1,iw,ik) = amn
                  ELSE
                     CALL errore('compute_amn',' value of wan_mode not recognised',1)
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
      ENDDO  ! k-points
      DEALLOCATE (sgf,csph)
   IF(any_uspp) THEN
     CALL deallocate_bec_type (becp)
   ENDIF
   !
   IF (ionode .and. wan_mode=='standalone') CLOSE (iun_amn)

   WRITE(stdout,*)
   WRITE(stdout,*) ' AMN calculated'

   RETURN
END SUBROUTINE compute_amn
!
!
SUBROUTINE generate_guiding_functions(ik)
   !
   USE io_global,  ONLY : stdout
   USE constants, ONLY : pi, tpi, fpi, eps8
   USE wvfct, ONLY : npw, g2kin, igk
   USE control_flags, ONLY : gamma_only
   USE gvect, ONLY : g, gstart
   USE cell_base,  ONLY : tpiba2, omega, tpiba
   USE wannier
   USE klist,      ONLY : xk
   USE cell_base, ONLY : bg
   USE mp, ONLY : mp_sum
   USE mp_global, ONLY : intra_pool_comm

   IMPLICIT NONE

   INTEGER, PARAMETER :: lmax=3, lmax2=(lmax+1)**2
   INTEGER :: iw, ig, ik, bgtau(3), isph, l, mesh_r
   INTEGER :: lmax_iw, lm, ipol, n1, n2, n3, nr1, nr2, nr3, iig
   real(DP) :: arg, anorm, fac, alpha_w2, yy, alfa, ddot
   COMPLEX(DP) :: zdotc, kphase, lphase, gff, lph
   real(DP), ALLOCATABLE :: gk(:,:), qg(:), ylm(:,:), radial(:,:)
   COMPLEX(DP), ALLOCATABLE :: sk(:)
   !
   ALLOCATE( gk(3,npw), qg(npw), ylm(npw,lmax2), sk(npw), radial(npw,0:lmax) )
   !
   DO ig = 1, npw
      gk (1,ig) = xk(1, ik) + g(1, igk(ig) )
      gk (2,ig) = xk(2, ik) + g(2, igk(ig) )
      gk (3,ig) = xk(3, ik) + g(3, igk(ig) )
      qg(ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
   ENDDO

   CALL ylmr2 (lmax2, npw, gk, qg, ylm)
   ! define qg as the norm of (k+g) in a.u.
   qg(:) = sqrt(qg(:)) * tpiba

   DO iw = 1, n_proj
      !
      gf(:,iw) = (0.d0,0.d0)

      CALL radialpart(npw, qg, alpha_w(iw), r_w(iw), lmax, radial)

      DO lm = 1, lmax2
         IF ( abs(csph(lm,iw)) < eps8 ) CYCLE
         l = int (sqrt( lm-1.d0))
         lphase = (0.d0,-1.d0)**l
         !
         DO ig=1,npw
            gf(ig,iw) = gf(ig,iw) + csph(lm,iw) * ylm(ig,lm) * radial(ig,l) * lphase
         ENDDO !ig
      ENDDO ! lm
      DO ig=1,npw
         iig = igk(ig)
         arg = ( gk(1,ig)*center_w(1,iw) + gk(2,ig)*center_w(2,iw) + &
                                           gk(3,ig)*center_w(3,iw) ) * tpi
         ! center_w are cartesian coordinates in units of alat
         sk(ig) = cmplx(cos(arg), -sin(arg) ,kind=DP)
         gf(ig,iw) = gf(ig,iw) * sk(ig)
      ENDDO
      IF (gamma_only) THEN
          anorm = 2.0_dp*ddot(2*npw,gf(1,iw),1,gf(1,iw),1)
          IF (gstart==2) anorm = anorm - abs(gf(1,iw))**2
      ELSE
          anorm = REAL(zdotc(npw,gf(1,iw),1,gf(1,iw),1))
      ENDIF
      CALL mp_sum(anorm, intra_pool_comm)
!      write (stdout,*) ik, iw, anorm
      gf(:,iw) = gf(:,iw) / dsqrt(anorm)
   ENDDO
   !
   DEALLOCATE ( gk, qg, ylm, sk, radial)
   RETURN
END SUBROUTINE generate_guiding_functions

SUBROUTINE write_band
   USE io_global,  ONLY : stdout, ionode
   USE wvfct, ONLY : nbnd, et
   USE klist, ONLY : nkstot
   USE constants, ONLY: rytoev
   USE wannier

   IMPLICIT NONE
   !
   INTEGER, EXTERNAL :: find_free_unit
   !
   INTEGER ik, ibnd, ibnd1, ikevc

   IF (wan_mode=='standalone') THEN
      iun_band = find_free_unit()
      IF (ionode) OPEN (unit=iun_band, file=trim(seedname)//".eig",form='formatted')
   ENDIF

   IF (wan_mode=='library') ALLOCATE(eigval(num_bands,iknum))

   DO ik=ikstart,ikstop
      ikevc = ik - ikstart + 1
      ibnd1=0
      DO ibnd=1,nbnd
         IF (excluded_band(ibnd)) CYCLE
         ibnd1=ibnd1 + 1
         IF (wan_mode=='standalone') THEN
            IF (ionode) WRITE (iun_band,'(2i5,f18.12)') ibnd1, ikevc, et(ibnd,ik)*rytoev
         ELSEIF (wan_mode=='library') THEN
            eigval(ibnd1,ikevc) = et(ibnd,ik)*rytoev
         ELSE
            CALL errore('write_band',' value of wan_mode not recognised',1)
         ENDIF
      ENDDO
   ENDDO
   RETURN
END SUBROUTINE write_band

SUBROUTINE write_plot
   USE io_global,  ONLY : stdout, ionode
   USE wvfct, ONLY : nbnd, npw, igk, g2kin, ecutwfc
   USE control_flags, ONLY : gamma_only
   USE wavefunctions_module, ONLY : evc, psic
   USE io_files, ONLY : nwordwfc, iunwfc
   USE wannier
   USE gvecs,         ONLY : nls, nlsm
   USE klist,           ONLY : nkstot, xk
   USE gvect,           ONLY : g, ngm
   USE cell_base,       ONLY : tpiba2
   USE fft_base,        ONLY : cgather_smooth, dffts
   USE fft_interfaces,  ONLY : invfft
   USE noncollin_module,ONLY : noncolin

   IMPLICIT NONE
   !
   INTEGER, EXTERNAL :: find_free_unit
   !
   INTEGER ik, ibnd, ibnd1, ikevc, i1, j, spin
   CHARACTER*20 wfnname

   ! aam: 1/5/06: for writing smaller unk files
   INTEGER :: n1by2,n2by2,n3by2,i,k,idx,pos
   COMPLEX(DP),ALLOCATABLE :: psic_small(:)
   !-------------------------------------------!

#ifdef __MPI
   INTEGER nxxs
   COMPLEX(DP),ALLOCATABLE :: psic_all(:)
   nxxs = dffts%nr1x * dffts%nr2x * dffts%nr3x
   ALLOCATE(psic_all(nxxs) )
#endif

   IF(noncolin) CALL errore('pw2wannier90',&
       'write_unk not implemented with ncls',1)

   IF (reduce_unk) THEN
      WRITE(stdout,'(3(a,i5))') 'nr1s =',dffts%nr1,'nr2s=',dffts%nr2,'nr3s=',dffts%nr3
      n1by2=(dffts%nr1+1)/2
      n2by2=(dffts%nr2+1)/2
      n3by2=(dffts%nr3+1)/2
      WRITE(stdout,'(3(a,i5))') 'n1by2=',n1by2,'n2by2=',n2by2,'n3by2=',n3by2
      ALLOCATE(psic_small(n1by2*n2by2*n3by2))
   ENDIF

   DO ik=ikstart,ikstop

      ikevc = ik - ikstart + 1

      iun_plot = find_free_unit()
      !write(wfnname,200) p,spin
      spin=ispinw
      IF(ispinw==0) spin=1
      WRITE(wfnname,200) ikevc, spin
200   FORMAT ('UNK',i5.5,'.',i1)

   IF (ionode) THEN
      IF(wvfn_formatted) THEN
         OPEN (unit=iun_plot, file=wfnname,form='formatted')
         IF (reduce_unk) THEN
            WRITE(iun_plot,*)  n1by2,n2by2,n3by2, ikevc, nbnd-nexband
         ELSE
            WRITE(iun_plot,*)  dffts%nr1,dffts%nr2,dffts%nr3,ikevc,nbnd-nexband
         ENDIF
      ELSE
         OPEN (unit=iun_plot, file=wfnname,form='unformatted')
         IF (reduce_unk) THEN
            WRITE(iun_plot)  n1by2,n2by2,n3by2, ikevc, nbnd-nexband
         ELSE
            WRITE(iun_plot)  dffts%nr1,dffts%nr2,dffts%nr3,ikevc,nbnd-nexband
         ENDIF
      ENDIF
   ENDIF

      CALL davcio (evc, nwordwfc, iunwfc, ik, -1 )
      CALL gk_sort (xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)

      ibnd1 = 0
      DO ibnd=1,nbnd
         IF (excluded_band(ibnd)) CYCLE
         ibnd1=ibnd1 + 1
         psic(:) = (0.d0, 0.d0)
         psic(nls (igk (1:npw) ) ) = evc (1:npw, ibnd)
         IF (gamma_only)  psic(nlsm(igk (1:npw) ) ) = conjg(evc (1:npw, ibnd))
         CALL invfft ('Wave', psic, dffts)
         IF (reduce_unk) pos=0
#ifdef __MPI
         CALL cgather_smooth(psic,psic_all)
         IF (reduce_unk) THEN
            DO k=1,dffts%nr3,2
               DO j=1,dffts%nr2,2
                  DO i=1,dffts%nr1,2
                     idx = (k-1)*dffts%nr2*dffts%nr1 + (j-1)*dffts%nr1 + i
                     pos=pos+1
                     psic_small(pos) = psic_all(idx)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      IF (ionode) THEN
         IF(wvfn_formatted) THEN
            IF (reduce_unk) THEN
               WRITE (iun_plot,'(2ES20.10)') (psic_small(j),j=1,n1by2*n2by2*n3by2)
            ELSE
               WRITE (iun_plot,*) (psic_all(j),j=1,dffts%nr1*dffts%nr2*dffts%nr3)
            ENDIF
         ELSE
            IF (reduce_unk) THEN
               WRITE (iun_plot) (psic_small(j),j=1,n1by2*n2by2*n3by2)
            ELSE
               WRITE (iun_plot) (psic_all(j),j=1,dffts%nr1*dffts%nr2*dffts%nr3)
            ENDIF
         ENDIF
      ENDIF
#else
         IF (reduce_unk) THEN
            DO k=1,dffts%nr3,2
               DO j=1,dffts%nr2,2
                  DO i=1,dffts%nr1,2
                     idx = (k-1)*dffts%nr2*dffts%nr1 + (j-1)*dffts%nr1 + i
                     pos=pos+1
                     psic_small(pos) = psic(idx)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
         IF(wvfn_formatted) THEN
            IF (reduce_unk) THEN
               WRITE (iun_plot,'(2ES20.10)') (psic_small(j),j=1,n1by2*n2by2*n3by2)
            ELSE
               WRITE (iun_plot,*) (psic(j),j=1,dffts%nr1*dffts%nr2*dffts%nr3)
            ENDIF
         ELSE
            IF (reduce_unk) THEN
               WRITE (iun_plot) (psic_small(j),j=1,n1by2*n2by2*n3by2)
            ELSE
               WRITE (iun_plot) (psic(j),j=1,dffts%nr1*dffts%nr2*dffts%nr3)
            ENDIF
         ENDIF
#endif
      ENDDO !ibnd

      IF(ionode) CLOSE (unit=iun_plot)

   ENDDO  !ik

   IF (reduce_unk) DEALLOCATE(psic_small)

#ifdef __MPI
   DEALLOCATE( psic_all )
#endif
   RETURN
END SUBROUTINE write_plot

SUBROUTINE write_parity

   USE mp_global,            ONLY : mpime,nproc,intra_pool_comm
   USE mp,                   ONLY : mp_sum
   USE io_global,            ONLY : stdout, ionode
   USE wvfct,                ONLY : nbnd, npw, igk, g2kin, ecutwfc
   USE control_flags,        ONLY : gamma_only
   USE wavefunctions_module, ONLY : evc
   USE io_files,             ONLY : nwordwfc, iunwfc
   USE wannier
   USE klist,                ONLY : nkstot, xk
   USE gvect,                ONLY : g, ngm
   USE cell_base,            ONLY : tpiba2, at
   USE constants,            ONLY : eps6

   IMPLICIT NONE
   !
   INTEGER, EXTERNAL :: find_free_unit
   !
   INTEGER                      :: ibnd,igv,kgamma,ik,i,ig_idx(32)
   INTEGER,DIMENSION(nproc)     :: num_G,displ

   real(kind=dp)                :: g_abc_1D(32),g_abc_gathered(3,32)
   real(kind=dp),ALLOCATABLE    :: g_abc(:,:),g_abc_pre_gather(:,:,:)
   COMPLEX(kind=dp),ALLOCATABLE :: evc_sub(:,:,:),evc_sub_gathered(:,:)
   COMPLEX(kind=dp)             :: evc_sub_1D(32)
   !
   ! getting the ik index corresponding to the Gamma point
   ! ... and the spin channel (fix due to N Poilvert, Feb 2011)
   !
   IF (.not. gamma_only) THEN
      DO ik=ikstart,ikstop
         IF ( (xk(1,ik)/= 0.d0) .or. (xk(2,ik)/= 0.d0) .or. (xk(3,ik)/= 0.d0) ) THEN
            IF (ik == ikstop) CALL errore('write_parity',&
                 ' parity calculation may only be performed at the gamma point.',1)
            CYCLE
         ELSE
            ! NP: spin unpolarized or "up" component of spin
            IF (ispinw == 0 .or. ispinw == 1) THEN
               kgamma=ik
            ELSE ! NP: "down" component
               kgamma=ik+1
            ENDIF
            exit
         ENDIF
      ENDDO
   ELSE
      ! NP: spin unpolarized or "up" component of spin
      IF (ispinw == 0 .or. ispinw == 1) THEN
         kgamma=1
      ELSE ! NP: "down" component
         kgamma=2
      ENDIF
   ENDIF
   !
   ! building the evc array corresponding to the Gamma point
   !
   CALL davcio (evc, nwordwfc, iunwfc, kgamma, -1 )
   CALL gk_sort (xk(1,kgamma), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
   !
   ! opening the <seedname>.unkg file
   !
   iun_parity = find_free_unit()
   IF (ionode)  THEN
        OPEN (unit=iun_parity, file=trim(seedname)//".unkg",form='formatted')
        WRITE(stdout,*)"Finding the 32 unkg's per band required for parity signature."
   ENDIF
   !
   ! g_abc(:,ipw) are the coordinates of the ipw-th G vector in b1, b2, b3 basis,
   ! we compute them from g(:,ipw) by multiplying : transpose(at) with g(:,ipw)
   !
   ALLOCATE(g_abc(3,npw))
   DO igv=1,npw
       g_abc(:,igk(igv))=matmul(transpose(at),g(:,igk(igv)))
   ENDDO
   !
   ! Count and identify the G vectors we will be extracting for each
   ! cpu.
   !
   ig_idx=0
   num_G = 0
   DO igv=1,npw
       ! 0-th Order
       IF ( (abs(g_abc(1,igv) - 0.d0) <= eps6) .and. &
            (abs(g_abc(2,igv) - 0.d0) <= eps6) .and. &
            (abs(g_abc(3,igv)) - 0.d0 <= eps6) ) THEN ! 1
           num_G(mpime+1) = num_G(mpime+1) + 1
           ig_idx(num_G(mpime+1))=igv
           CYCLE
       ENDIF
       ! 1st Order
       IF ( (abs(g_abc(1,igv) - 1.d0) <= eps6) .and. &
            (abs(g_abc(2,igv) - 0.d0) <= eps6) .and. &
            (abs(g_abc(3,igv)) - 0.d0 <= eps6) ) THEN ! x
           num_G(mpime+1) = num_G(mpime+1) + 1
           ig_idx(num_G(mpime+1))=igv
           CYCLE
       ENDIF
       IF ( (abs(g_abc(1,igv) - 0.d0) <= eps6) .and. &
            (abs(g_abc(2,igv) - 1.d0) <= eps6) .and. &
            (abs(g_abc(3,igv)) - 0.d0 <= eps6) ) THEN ! y
           num_G(mpime+1) = num_G(mpime+1) + 1
           ig_idx(num_G(mpime+1))=igv
           CYCLE
       ENDIF
       IF ( (abs(g_abc(1,igv) - 0.d0) <= eps6) .and. &
            (abs(g_abc(2,igv) - 0.d0) <= eps6) .and. &
            (abs(g_abc(3,igv)) - 1.d0 <= eps6) ) THEN ! z
           num_G(mpime+1) = num_G(mpime+1) + 1
           ig_idx(num_G(mpime+1))=igv
           CYCLE
       ENDIF
       ! 2nd Order
       IF ( (abs(g_abc(1,igv) - 2.d0) <= eps6) .and. &
            (abs(g_abc(2,igv) - 0.d0) <= eps6) .and. &
            (abs(g_abc(3,igv)) - 0.d0 <= eps6) ) THEN ! x^2
           num_G(mpime+1) = num_G(mpime+1) + 1
           ig_idx(num_G(mpime+1))=igv
           CYCLE
       ENDIF
       IF ( (abs(g_abc(1,igv) - 1.d0) <= eps6) .and. &
            (abs(g_abc(2,igv) - 1.d0) <= eps6) .and. &
            (abs(g_abc(3,igv)) - 0.d0 <= eps6) ) THEN ! xy
           num_G(mpime+1) = num_G(mpime+1) + 1
           ig_idx(num_G(mpime+1))=igv
           CYCLE
       ENDIF
       IF ( (abs(g_abc(1,igv) - 1.d0) <= eps6) .and. &
            (abs(g_abc(2,igv) + 1.d0) <= eps6) .and. &
            (abs(g_abc(3,igv)) - 0.d0 <= eps6) ) THEN ! xy
           num_G(mpime+1) = num_G(mpime+1) + 1
           ig_idx(num_G(mpime+1))=igv
           CYCLE
       ENDIF
       IF ( (abs(g_abc(1,igv) - 1.d0) <= eps6) .and. &
            (abs(g_abc(2,igv) - 0.d0) <= eps6) .and. &
            (abs(g_abc(3,igv)) - 1.d0 <= eps6) ) THEN ! xz
           num_G(mpime+1) = num_G(mpime+1) + 1
           ig_idx(num_G(mpime+1))=igv
           CYCLE
       ENDIF
       IF ( (abs(g_abc(1,igv) - 1.d0) <= eps6) .and. &
            (abs(g_abc(2,igv) - 0.d0) <= eps6) .and. &
            (abs(g_abc(3,igv)) + 1.d0 <= eps6) ) THEN ! xz
           num_G(mpime+1) = num_G(mpime+1) + 1
           ig_idx(num_G(mpime+1))=igv
           CYCLE
       ENDIF
       IF ( (abs(g_abc(1,igv) - 0.d0) <= eps6) .and. &
            (abs(g_abc(2,igv) - 2.d0) <= eps6) .and. &
            (abs(g_abc(3,igv)) - 0.d0 <= eps6) ) THEN ! y^2
           num_G(mpime+1) = num_G(mpime+1) + 1
           ig_idx(num_G(mpime+1))=igv
           CYCLE
       ENDIF
       IF ( (abs(g_abc(1,igv) - 0.d0) <= eps6) .and. &
            (abs(g_abc(2,igv) - 1.d0) <= eps6) .and. &
            (abs(g_abc(3,igv)) - 1.d0 <= eps6) ) THEN ! yz
           num_G(mpime+1) = num_G(mpime+1) + 1
           ig_idx(num_G(mpime+1))=igv
           CYCLE
       ENDIF
       IF ( (abs(g_abc(1,igv) - 0.d0) <= eps6) .and. &
            (abs(g_abc(2,igv) - 1.d0) <= eps6) .and. &
            (abs(g_abc(3,igv)) + 1.d0 <= eps6) ) THEN ! yz
           num_G(mpime+1) = num_G(mpime+1) + 1
           ig_idx(num_G(mpime+1))=igv
           CYCLE
       ENDIF
       IF ( (abs(g_abc(1,igv) - 0.d0) <= eps6) .and. &
            (abs(g_abc(2,igv) - 0.d0) <= eps6) .and. &
            (abs(g_abc(3,igv)) - 2.d0 <= eps6) ) THEN ! z^2
           num_G(mpime+1) = num_G(mpime+1) + 1
           ig_idx(num_G(mpime+1))=igv
           CYCLE
       ENDIF
       ! 3rd Order
       IF ( (abs(g_abc(1,igv) - 3.d0) <= eps6) .and. &
            (abs(g_abc(2,igv) - 0.d0) <= eps6) .and. &
            (abs(g_abc(3,igv)) - 0.d0 <= eps6) ) THEN ! x^3
           num_G(mpime+1) = num_G(mpime+1) + 1
           ig_idx(num_G(mpime+1))=igv
           CYCLE
       ENDIF
       IF ( (abs(g_abc(1,igv) - 2.d0) <= eps6) .and. &
            (abs(g_abc(2,igv) - 1.d0) <= eps6) .and. &
            (abs(g_abc(3,igv)) - 0.d0 <= eps6) ) THEN ! x^2y
           num_G(mpime+1) = num_G(mpime+1) + 1
           ig_idx(num_G(mpime+1))=igv
           CYCLE
       ENDIF
       IF ( (abs(g_abc(1,igv) - 2.d0) <= eps6) .and. &
            (abs(g_abc(2,igv) + 1.d0) <= eps6) .and. &
            (abs(g_abc(3,igv)) - 0.d0 <= eps6) ) THEN ! x^2y
           num_G(mpime+1) = num_G(mpime+1) + 1
           ig_idx(num_G(mpime+1))=igv
           CYCLE
       ENDIF
       IF ( (abs(g_abc(1,igv) - 2.d0) <= eps6) .and. &
            (abs(g_abc(2,igv) - 0.d0) <= eps6) .and. &
            (abs(g_abc(3,igv)) - 1.d0 <= eps6) ) THEN ! x^2z
           num_G(mpime+1) = num_G(mpime+1) + 1
           ig_idx(num_G(mpime+1))=igv
           CYCLE
       ENDIF
       IF ( (abs(g_abc(1,igv) - 2.d0) <= eps6) .and. &
            (abs(g_abc(2,igv) - 0.d0) <= eps6) .and. &
            (abs(g_abc(3,igv)) + 1.d0 <= eps6) ) THEN ! x^2z
           num_G(mpime+1) = num_G(mpime+1) + 1
           ig_idx(num_G(mpime+1))=igv
           CYCLE
       ENDIF
       IF ( (abs(g_abc(1,igv) - 1.d0) <= eps6) .and. &
            (abs(g_abc(2,igv) - 2.d0) <= eps6) .and. &
            (abs(g_abc(3,igv)) - 0.d0 <= eps6) ) THEN ! xy^2
           num_G(mpime+1) = num_G(mpime+1) + 1
           ig_idx(num_G(mpime+1))=igv
           CYCLE
       ENDIF
       IF ( (abs(g_abc(1,igv) - 1.d0) <= eps6) .and. &
            (abs(g_abc(2,igv) + 2.d0) <= eps6) .and. &
            (abs(g_abc(3,igv)) - 0.d0 <= eps6) ) THEN ! xy^2
           num_G(mpime+1) = num_G(mpime+1) + 1
           ig_idx(num_G(mpime+1))=igv
           CYCLE
       ENDIF
       IF ( (abs(g_abc(1,igv) - 1.d0) <= eps6) .and. &
            (abs(g_abc(2,igv) - 1.d0) <= eps6) .and. &
            (abs(g_abc(3,igv)) - 1.d0 <= eps6) ) THEN ! xyz
           num_G(mpime+1) = num_G(mpime+1) + 1
           ig_idx(num_G(mpime+1))=igv
           CYCLE
       ENDIF
       IF ( (abs(g_abc(1,igv) - 1.d0) <= eps6) .and. &
            (abs(g_abc(2,igv) - 1.d0) <= eps6) .and. &
            (abs(g_abc(3,igv)) + 1.d0 <= eps6) ) THEN ! xyz
           num_G(mpime+1) = num_G(mpime+1) + 1
           ig_idx(num_G(mpime+1))=igv
           CYCLE
       ENDIF
       IF ( (abs(g_abc(1,igv) - 1.d0) <= eps6) .and. &
            (abs(g_abc(2,igv) + 1.d0) <= eps6) .and. &
            (abs(g_abc(3,igv)) - 1.d0 <= eps6) ) THEN ! xyz
           num_G(mpime+1) = num_G(mpime+1) + 1
           ig_idx(num_G(mpime+1))=igv
           CYCLE
       ENDIF
       IF ( (abs(g_abc(1,igv) - 1.d0) <= eps6) .and. &
            (abs(g_abc(2,igv) + 1.d0) <= eps6) .and. &
            (abs(g_abc(3,igv)) + 1.d0 <= eps6) ) THEN ! xyz
           num_G(mpime+1) = num_G(mpime+1) + 1
           ig_idx(num_G(mpime+1))=igv
           CYCLE
       ENDIF
       IF ( (abs(g_abc(1,igv) - 1.d0) <= eps6) .and. &
            (abs(g_abc(2,igv) - 0.d0) <= eps6) .and. &
            (abs(g_abc(3,igv)) - 2.d0 <= eps6) ) THEN ! xz^2
           num_G(mpime+1) = num_G(mpime+1) + 1
           ig_idx(num_G(mpime+1))=igv
           CYCLE
       ENDIF
       IF ( (abs(g_abc(1,igv) - 1.d0) <= eps6) .and. &
            (abs(g_abc(2,igv) - 0.d0) <= eps6) .and. &
            (abs(g_abc(3,igv)) + 2.d0 <= eps6) ) THEN ! xz^2
           num_G(mpime+1) = num_G(mpime+1) + 1
           ig_idx(num_G(mpime+1))=igv
           CYCLE
       ENDIF
       IF ( (abs(g_abc(1,igv) - 0.d0) <= eps6) .and. &
            (abs(g_abc(2,igv) - 3.d0) <= eps6) .and. &
            (abs(g_abc(3,igv)) - 0.d0 <= eps6) ) THEN ! y^3
           num_G(mpime+1) = num_G(mpime+1) + 1
           ig_idx(num_G(mpime+1))=igv
           CYCLE
       ENDIF
       IF ( (abs(g_abc(1,igv) - 0.d0) <= eps6) .and. &
            (abs(g_abc(2,igv) - 2.d0) <= eps6) .and. &
            (abs(g_abc(3,igv)) - 1.d0 <= eps6) ) THEN ! y^2z
           num_G(mpime+1) = num_G(mpime+1) + 1
           ig_idx(num_G(mpime+1))=igv
           CYCLE
       ENDIF
       IF ( (abs(g_abc(1,igv) - 0.d0) <= eps6) .and. &
            (abs(g_abc(2,igv) - 2.d0) <= eps6) .and. &
            (abs(g_abc(3,igv)) + 1.d0 <= eps6) ) THEN ! y^2z
           num_G(mpime+1) = num_G(mpime+1) + 1
           ig_idx(num_G(mpime+1))=igv
           CYCLE
       ENDIF
       IF ( (abs(g_abc(1,igv) - 0.d0) <= eps6) .and. &
            (abs(g_abc(2,igv) - 1.d0) <= eps6) .and. &
            (abs(g_abc(3,igv)) - 2.d0 <= eps6) ) THEN ! yz^2
           num_G(mpime+1) = num_G(mpime+1) + 1
           ig_idx(num_G(mpime+1))=igv
           CYCLE
       ENDIF
       IF ( (abs(g_abc(1,igv) - 0.d0) <= eps6) .and.&
            (abs(g_abc(2,igv) - 1.d0) <= eps6) .and. &
            (abs(g_abc(3,igv)) + 2.d0 <= eps6) ) THEN ! yz^2
           num_G(mpime+1) = num_G(mpime+1) + 1
           ig_idx(num_G(mpime+1))=igv
           CYCLE
       ENDIF
       IF ( (abs(g_abc(1,igv) - 0.d0) <= eps6) .and. &
            (abs(g_abc(2,igv) - 0.d0) <= eps6) .and. &
            (abs(g_abc(3,igv)) - 3.d0 <= eps6) ) THEN ! z^3
           num_G(mpime+1) = num_G(mpime+1) + 1
           ig_idx(num_G(mpime+1))=igv
           CYCLE
       ENDIF
   ENDDO
   !
   ! Sum laterally across cpus num_G, so it contains
   ! the number of g_vectors on each node, and known to all cpus
   !
   CALL mp_sum(num_G, intra_pool_comm)

   IF (ionode) WRITE(iun_parity,*) sum(num_G)
   IF (sum(num_G) /= 32) CALL errore('write_parity', 'incorrect number of g-vectors extracted',1)
   IF (ionode) THEN
      WRITE(stdout,*)'     ...done'
      WRITE(stdout,*)'G-vector splitting:'
      DO i=1,nproc
         WRITE(stdout,*)' cpu: ',i-1,' number g-vectors: ',num_G(i)
      ENDDO
      WRITE(stdout,*)' Collecting g-vectors and writing to file'
   ENDIF

   !
   ! Define needed intermediate arrays
   !
   ALLOCATE(evc_sub(32,nbnd,nproc))
   ALLOCATE(evc_sub_gathered(32,nbnd))
   ALLOCATE(g_abc_pre_gather(3,32,nproc))
   !
   ! Initialise
   !
   evc_sub=(0.d0,0.d0)
   evc_sub_1D=(0.d0,0.d0)
   evc_sub_gathered=(0.d0,0.d0)
   g_abc_pre_gather=0
   g_abc_1D=0
   g_abc_gathered=0
   !
   ! Compute displacements needed for filling evc_sub
   !
   displ(1)=1
   IF (nproc > 1) THEN
       DO i=2,nproc
           displ(i)=displ(i-1)+num_G(i-1)
       ENDDO
   ENDIF
   !
   ! Fill evc_sub with required fourier component from each cpu dependent evc
   !
   DO i=1,num_G(mpime+1)
       evc_sub(i+displ(mpime+1)-1,:,mpime+1)=evc(ig_idx(i),:)
   ENDDO
   !
   ! g_abc_pre_gather(:,ipw,icpu) are the coordinates of the ipw-th G vector in b1, b2, b3 basis
   ! on icpu and stored sequencially, ready for a lateral mp_sum
   !
   DO igv=1,num_G(mpime+1)
       g_abc_pre_gather(:,igv+displ(mpime+1)-1,mpime+1)=matmul(transpose(at),g(:,ig_idx(igk(igv))))
   ENDDO
   !
   ! Gather evc_sub and  g_abc_pre_gather into common arrays to each cpu
   !
   DO ibnd=1,nbnd
      evc_sub_1D=evc_sub(:,ibnd,mpime+1)
      CALL mp_sum(evc_sub_1D, intra_pool_comm)
      evc_sub_gathered(:,ibnd)=evc_sub_1D
   ENDDO
   !
   DO i=1,3
       g_abc_1D=g_abc_pre_gather(i,:,mpime+1)
       CALL mp_sum(g_abc_1D, intra_pool_comm)
       g_abc_gathered(i,:)=g_abc_1D
   ENDDO
   !
   ! Write to file
   !
   DO ibnd=1,nbnd
      DO igv=1,32
         IF (ionode) WRITE(iun_parity,'(5i5,2f12.7)') ibnd, igv, nint(g_abc_gathered(1,igv)),&
                                                                 nint(g_abc_gathered(2,igv)),&
                                                                 nint(g_abc_gathered(3,igv)),&
                                                                 real(evc_sub_gathered(igv,ibnd)),&
                                                                aimag(evc_sub_gathered(igv,ibnd))
      ENDDO
   ENDDO
   WRITE(stdout,*)'     ...done'
   !
   IF (ionode) CLOSE(unit=iun_parity)
   !
   DEALLOCATE(evc_sub)
   DEALLOCATE(evc_sub_gathered)
   DEALLOCATE(g_abc_pre_gather)

END SUBROUTINE write_parity


SUBROUTINE wan2sic

  USE io_global,  ONLY : stdout
  USE kinds, ONLY : DP
  USE io_files, ONLY : iunwfc, iunatsicwfc, nwordwfc, nwordwann
  USE cell_base, ONLY : omega, tpiba2
  USE gvect, ONLY : g, ngm
  USE gvecs, ONLY: nls
  USE wavefunctions_module, ONLY : evc, psic
  USE wvfct, ONLY : nbnd, npwx, npw, igk, g2kin, ecutwfc
  USE klist, ONLY : nkstot, xk, wk
  USE wannier

  INTEGER :: i, j, nn, ik, ibnd, iw, ikevc
  COMPLEX(DP), ALLOCATABLE :: orbital(:,:), u_matrix(:,:,:)

  OPEN (20, file = trim(seedname)//".dat" , form = 'formatted', status = 'unknown')
  WRITE(stdout,*) ' wannier plot '

  ALLOCATE ( u_matrix( n_wannier, n_wannier, nkstot) )
  ALLOCATE ( orbital( npwx, n_wannier) )

  !
  DO i = 1, n_wannier
     DO j = 1, n_wannier
        DO ik = 1, nkstot
           READ (20, * ) u_matrix(i,j,ik)
           !do nn = 1, nnb(ik)
           DO nn = 1, nnb
              READ (20, * ) ! m_matrix (i,j,nkp,nn)
           ENDDO
        ENDDO  !nkp
     ENDDO !j
  ENDDO !i
  !
  DO ik=1,iknum
     ikevc = ik + ikstart - 1
     CALL davcio (evc, nwordwfc, iunwfc, ikevc, -1)
     CALL gk_sort (xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
     WRITE(stdout,*) 'npw ',npw
     DO iw=1,n_wannier
        DO j=1,npw
           orbital(j,iw) = (0.0d0,0.0d0)
           DO ibnd=1,n_wannier
              orbital(j,iw) = orbital(j,iw) + u_matrix(iw,ibnd,ik)*evc(j,ibnd)
              WRITE(stdout,*) j, iw, ibnd, ik, orbital(j,iw), &
                              u_matrix(iw,ibnd,ik), evc(j,ibnd)
           ENDDO !ibnd
        ENDDO  !j
     ENDDO !wannier
     CALL davcio (orbital, nwordwann, iunatsicwfc, ikevc, +1)
  ENDDO ! k-points

  DEALLOCATE ( u_matrix)
  WRITE(stdout,*) ' dealloc u '
  DEALLOCATE (  orbital)
  WRITE(stdout,*) ' dealloc orbital '
  !
END SUBROUTINE wan2sic

SUBROUTINE ylm_expansion
   USE io_global,  ONLY : stdout
   USE kinds, ONLY :  DP
   USE random_numbers,       ONLY : randy
   USE wannier
   IMPLICIT NONE
   ! local variables
   INTEGER, PARAMETER :: lmax2=16
   INTEGER ::  lm, i, ir, iw, m
   real(DP) :: capel
   real(DP), ALLOCATABLE :: r(:,:), rr(:), rp(:,:), ylm_w(:), ylm(:,:), mly(:,:)
   real(DP) :: u(3,3)

   ALLOCATE (r(3,lmax2), rp(3,lmax2), rr(lmax2), ylm_w(lmax2))
   ALLOCATE (ylm(lmax2,lmax2), mly(lmax2,lmax2) )

   ! generate a set of nr=lmax2 random vectors
   DO ir=1,lmax2
      DO i=1,3
         r(i,ir) = randy() -0.5d0
      ENDDO
   ENDDO
   rr(:) = r(1,:)*r(1,:) + r(2,:)*r(2,:) + r(3,:)*r(3,:)
   !- compute ylm(ir,lm)
   CALL ylmr2(lmax2, lmax2, r, rr, ylm)
   !- store the inverse of ylm(ir,lm) in mly(lm,ir)
   CALL invmat(lmax2, ylm, mly, capel)
   !- check that r points are independent
   CALL check_inverse(lmax2, ylm, mly)

   DO iw=1, n_proj

      !- define the u matrix that rotate the reference frame
      CALL set_u_matrix (xaxis(:,iw),zaxis(:,iw),u)
      !- find rotated r-vectors
      rp(:,:) = matmul ( u(:,:) , r(:,:) )
      !- set ylm funtion according to wannier90 (l,mr) indexing in the rotaterd points
      CALL ylm_wannier(ylm_w,l_w(iw),mr_w(iw),rp,lmax2)

      csph(:,iw) = matmul (mly(:,:), ylm_w(:))

!      write (stdout,*)
!      write (stdout,'(2i4,2(2x,3f6.3))') l_w(iw), mr_w(iw), xaxis(:,iw), zaxis(:,iw)
!      write (stdout,'(16i6)')   (lm, lm=1,lmax2)
!      write (stdout,'(16f6.3)') (csph(lm,iw), lm=1,lmax2)

   ENDDO
   DEALLOCATE (r, rp, rr, ylm_w, ylm, mly )

   RETURN
END SUBROUTINE ylm_expansion

SUBROUTINE check_inverse(lmax2, ylm, mly)
   USE kinds, ONLY :  DP
   USE constants, ONLY :  eps8
   IMPLICIT NONE
   ! I/O variables
   INTEGER :: lmax2
   real(DP) :: ylm(lmax2,lmax2), mly(lmax2,lmax2)
   ! local variables
   real(DP), ALLOCATABLE :: uno(:,:)
   real(DP) :: capel
   INTEGER :: lm
   !
   ALLOCATE (uno(lmax2,lmax2) )
   uno = matmul(mly, ylm)
   capel = 0.d0
   DO lm = 1, lmax2
      uno(lm,lm) = uno(lm,lm) - 1.d0
   ENDDO
   capel = capel + sum ( abs(uno(1:lmax2,1:lmax2) ) )
!   write (stdout,*) "capel = ", capel
   IF (capel > eps8) CALL errore('ylm_expansion', &
                    ' inversion failed: r(*,1:nr) are not all independent !!',1)
   DEALLOCATE (uno)
   RETURN
END SUBROUTINE check_inverse

SUBROUTINE set_u_matrix(x,z,u)
   USE kinds, ONLY :  DP
   USE constants, ONLY : eps6
   IMPLICIT NONE
   ! I/O variables
   real(DP) :: x(3),z(3),u(3,3)
   ! local variables
   real(DP) :: xx, zz, y(3), coseno

   xx = sqrt(x(1)*x(1) + x(2)*x(2) + x(3)*x(3))
   IF (xx < eps6) CALL errore ('set_u_matrix',' |xaxis| < eps ',1)
!   x(:) = x(:)/xx
   zz = sqrt(z(1)*z(1) + z(2)*z(2) + z(3)*z(3))
   IF (zz < eps6) CALL errore ('set_u_matrix',' |zaxis| < eps ',1)
!   z(:) = z(:)/zz

   coseno = (x(1)*z(1) + x(2)*z(2) + x(3)*z(3))/xx/zz
   IF (abs(coseno) > eps6) CALL errore('set_u_matrix',' xaxis and zaxis are not orthogonal !',1)

   y(1) = (z(2)*x(3) - x(2)*z(3))/xx/zz
   y(2) = (z(3)*x(1) - x(3)*z(1))/xx/zz
   y(3) = (z(1)*x(2) - x(1)*z(2))/xx/zz

   u(1,:) = x(:)/xx
   u(2,:) = y(:)
   u(3,:) = z(:)/zz

!   write (stdout,'(3f10.7)') u(:,:)

   RETURN

END SUBROUTINE set_u_matrix

SUBROUTINE ylm_wannier(ylm,l,mr,r,nr)
!
! this routine returns in ylm(r) the values at the nr points r(1:3,1:nr)
! of the spherical harmonic identified  by indices (l,mr)
! in table 3.1 of the wannierf90 specification.
!
! No reference to the particular ylm ordering internal to Quantum ESPRESSO
! is assumed.
!
! If ordering in wannier90 code is changed or extended this should be the
! only place to be modified accordingly
!
   USE kinds, ONLY :  DP
   USE constants, ONLY : pi, fpi, eps8
   IMPLICIT NONE
! I/O variables
!
   INTEGER :: l, mr, nr
   real(DP) :: ylm(nr), r(3,nr)
!
! local variables
!
   real(DP), EXTERNAL :: s, p_z,px,py, dz2, dxz, dyz, dx2my2, dxy
   real(DP), EXTERNAL :: fz3, fxz2, fyz2, fzx2my2, fxyz, fxx2m3y2, fy3x2my2
   real(DP) :: rr, cost, phi
   INTEGER :: ir
   real(DP) :: bs2, bs3, bs6, bs12
   bs2 = 1.d0/sqrt(2.d0)
   bs3=1.d0/sqrt(3.d0)
   bs6 = 1.d0/sqrt(6.d0)
   bs12 = 1.d0/sqrt(12.d0)
!
   IF (l > 3 .or. l < -5 ) CALL errore('ylm_wannier',' l out of range ', 1)
   IF (l>=0) THEN
      IF (mr < 1 .or. mr > 2*l+1) CALL errore('ylm_wannier','mr out of range' ,1)
   ELSE
      IF (mr < 1 .or. mr > abs(l)+1 ) CALL errore('ylm_wannier','mr out of range',1)
   ENDIF

   DO ir=1, nr
      rr = sqrt( r(1,ir)*r(1,ir) +  r(2,ir)*r(2,ir) + r(3,ir)*r(3,ir) )
      IF (rr < eps8) CALL errore('ylm_wannier',' rr too small ',1)

      cost =  r(3,ir) / rr
      !
      !  beware the arc tan, it is defined modulo pi
      !
      IF (r(1,ir) > eps8) THEN
         phi = atan( r(2,ir)/r(1,ir) )
      ELSEIF (r(1,ir) < -eps8 ) THEN
         phi = atan( r(2,ir)/r(1,ir) ) + pi
      ELSE
         phi = sign( pi/2.d0,r(2,ir) )
      ENDIF


      IF (l==0) THEN   ! s orbital
                    ylm(ir) = s(cost,phi)
      ENDIF
      IF (l==1) THEN   ! p orbitals
         IF (mr==1) ylm(ir) = p_z(cost,phi)
         IF (mr==2) ylm(ir) = px(cost,phi)
         IF (mr==3) ylm(ir) = py(cost,phi)
      ENDIF
      IF (l==2) THEN   ! d orbitals
         IF (mr==1) ylm(ir) = dz2(cost,phi)
         IF (mr==2) ylm(ir) = dxz(cost,phi)
         IF (mr==3) ylm(ir) = dyz(cost,phi)
         IF (mr==4) ylm(ir) = dx2my2(cost,phi)
         IF (mr==5) ylm(ir) = dxy(cost,phi)
      ENDIF
      IF (l==3) THEN   ! f orbitals
         IF (mr==1) ylm(ir) = fz3(cost,phi)
         IF (mr==2) ylm(ir) = fxz2(cost,phi)
         IF (mr==3) ylm(ir) = fyz2(cost,phi)
         IF (mr==4) ylm(ir) = fzx2my2(cost,phi)
         IF (mr==5) ylm(ir) = fxyz(cost,phi)
         IF (mr==6) ylm(ir) = fxx2m3y2(cost,phi)
         IF (mr==7) ylm(ir) = fy3x2my2(cost,phi)
      ENDIF
      IF (l==-1) THEN  !  sp hybrids
         IF (mr==1) ylm(ir) = bs2 * ( s(cost,phi) + px(cost,phi) )
         IF (mr==2) ylm(ir) = bs2 * ( s(cost,phi) - px(cost,phi) )
      ENDIF
      IF (l==-2) THEN  !  sp2 hybrids
         IF (mr==1) ylm(ir) = bs3*s(cost,phi)-bs6*px(cost,phi)+bs2*py(cost,phi)
         IF (mr==2) ylm(ir) = bs3*s(cost,phi)-bs6*px(cost,phi)-bs2*py(cost,phi)
         IF (mr==3) ylm(ir) = bs3*s(cost,phi) +2.d0*bs6*px(cost,phi)
      ENDIF
      IF (l==-3) THEN  !  sp3 hybrids
         IF (mr==1) ylm(ir) = 0.5d0*(s(cost,phi)+px(cost,phi)+py(cost,phi)+p_z(cost,phi))
         IF (mr==2) ylm(ir) = 0.5d0*(s(cost,phi)+px(cost,phi)-py(cost,phi)-p_z(cost,phi))
         IF (mr==3) ylm(ir) = 0.5d0*(s(cost,phi)-px(cost,phi)+py(cost,phi)-p_z(cost,phi))
         IF (mr==4) ylm(ir) = 0.5d0*(s(cost,phi)-px(cost,phi)-py(cost,phi)+p_z(cost,phi))
      ENDIF
      IF (l==-4) THEN  !  sp3d hybrids
         IF (mr==1) ylm(ir) = bs3*s(cost,phi)-bs6*px(cost,phi)+bs2*py(cost,phi)
         IF (mr==2) ylm(ir) = bs3*s(cost,phi)-bs6*px(cost,phi)-bs2*py(cost,phi)
         IF (mr==3) ylm(ir) = bs3*s(cost,phi) +2.d0*bs6*px(cost,phi)
         IF (mr==4) ylm(ir) = bs2*p_z(cost,phi)+bs2*dz2(cost,phi)
         IF (mr==5) ylm(ir) =-bs2*p_z(cost,phi)+bs2*dz2(cost,phi)
      ENDIF
      IF (l==-5) THEN  ! sp3d2 hybrids
         IF (mr==1) ylm(ir) = bs6*s(cost,phi)-bs2*px(cost,phi)-bs12*dz2(cost,phi)+.5d0*dx2my2(cost,phi)
         IF (mr==2) ylm(ir) = bs6*s(cost,phi)+bs2*px(cost,phi)-bs12*dz2(cost,phi)+.5d0*dx2my2(cost,phi)
         IF (mr==3) ylm(ir) = bs6*s(cost,phi)-bs2*py(cost,phi)-bs12*dz2(cost,phi)-.5d0*dx2my2(cost,phi)
         IF (mr==4) ylm(ir) = bs6*s(cost,phi)+bs2*py(cost,phi)-bs12*dz2(cost,phi)-.5d0*dx2my2(cost,phi)
         IF (mr==5) ylm(ir) = bs6*s(cost,phi)-bs2*p_z(cost,phi)+bs3*dz2(cost,phi)
         IF (mr==6) ylm(ir) = bs6*s(cost,phi)+bs2*p_z(cost,phi)+bs3*dz2(cost,phi)
      ENDIF

   ENDDO

   RETURN

END SUBROUTINE ylm_wannier

!======== l = 0 =====================================================================
FUNCTION s(cost,phi)
   USE kinds, ONLY :  DP
   USE constants, ONLY : fpi
   IMPLICIT NONE
   real(DP) :: s, cost,phi
   s = 1.d0/ sqrt(fpi)
   RETURN
END FUNCTION s
!======== l = 1 =====================================================================
FUNCTION p_z(cost,phi)
   USE kinds, ONLY :  DP
   USE constants, ONLY : fpi
   IMPLICIT NONE
   real(DP) ::p_z, cost,phi
   p_z =  sqrt(3.d0/fpi) * cost
   RETURN
END FUNCTION p_z
FUNCTION px(cost,phi)
   USE kinds, ONLY :  DP
   USE constants, ONLY : fpi
   IMPLICIT NONE
   real(DP) ::px, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   px =  sqrt(3.d0/fpi) * sint * cos(phi)
   RETURN
END FUNCTION px
FUNCTION py(cost,phi)
   USE kinds, ONLY :  DP
   USE constants, ONLY : fpi
   IMPLICIT NONE
   real(DP) ::py, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   py =  sqrt(3.d0/fpi) * sint * sin(phi)
   RETURN
END FUNCTION py
!======== l = 2 =====================================================================
FUNCTION dz2(cost,phi)
   USE kinds, ONLY :  DP
   USE constants, ONLY : fpi
   IMPLICIT NONE
   real(DP) ::dz2, cost, phi
   dz2 =  sqrt(1.25d0/fpi) * (3.d0* cost*cost-1.d0)
   RETURN
END FUNCTION dz2
FUNCTION dxz(cost,phi)
   USE kinds, ONLY :  DP
   USE constants, ONLY : fpi
   IMPLICIT NONE
   real(DP) ::dxz, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   dxz =  sqrt(15.d0/fpi) * sint*cost * cos(phi)
   RETURN
END FUNCTION dxz
FUNCTION dyz(cost,phi)
   USE kinds, ONLY :  DP
   USE constants, ONLY : fpi
   IMPLICIT NONE
   real(DP) ::dyz, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   dyz =  sqrt(15.d0/fpi) * sint*cost * sin(phi)
   RETURN
END FUNCTION dyz
FUNCTION dx2my2(cost,phi)
   USE kinds, ONLY :  DP
   USE constants, ONLY : fpi
   IMPLICIT NONE
   real(DP) ::dx2my2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   dx2my2 =  sqrt(3.75d0/fpi) * sint*sint * cos(2.d0*phi)
   RETURN
END FUNCTION dx2my2
FUNCTION dxy(cost,phi)
   USE kinds, ONLY :  DP
   USE constants, ONLY : fpi
   IMPLICIT NONE
   real(DP) ::dxy, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   dxy =  sqrt(3.75d0/fpi) * sint*sint * sin(2.d0*phi)
   RETURN
END FUNCTION dxy
!======== l = 3 =====================================================================
FUNCTION fz3(cost,phi)
   USE kinds, ONLY :  DP
   USE constants, ONLY : pi
   IMPLICIT NONE
   real(DP) ::fz3, cost, phi
   fz3 =  0.25d0*sqrt(7.d0/pi) * ( 5.d0 * cost * cost - 3.d0 ) * cost
   RETURN
END FUNCTION fz3
FUNCTION fxz2(cost,phi)
   USE kinds, ONLY :  DP
   USE constants, ONLY : pi
   IMPLICIT NONE
   real(DP) ::fxz2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fxz2 =  0.25d0*sqrt(10.5d0/pi) * ( 5.d0 * cost * cost - 1.d0 ) * sint * cos(phi)
   RETURN
END FUNCTION fxz2
FUNCTION fyz2(cost,phi)
   USE kinds, ONLY :  DP
   USE constants, ONLY : pi
   IMPLICIT NONE
   real(DP) ::fyz2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fyz2 =  0.25d0*sqrt(10.5d0/pi) * ( 5.d0 * cost * cost - 1.d0 ) * sint * sin(phi)
   RETURN
END FUNCTION fyz2
FUNCTION fzx2my2(cost,phi)
   USE kinds, ONLY :  DP
   USE constants, ONLY : pi
   IMPLICIT NONE
   real(DP) ::fzx2my2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fzx2my2 =  0.25d0*sqrt(105d0/pi) * sint * sint * cost * cos(2.d0*phi)
   RETURN
END FUNCTION fzx2my2
FUNCTION fxyz(cost,phi)
   USE kinds, ONLY :  DP
   USE constants, ONLY : pi
   IMPLICIT NONE
   real(DP) ::fxyz, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fxyz =  0.25d0*sqrt(105d0/pi) * sint * sint * cost * sin(2.d0*phi)
   RETURN
END FUNCTION fxyz
FUNCTION fxx2m3y2(cost,phi)
   USE kinds, ONLY :  DP
   USE constants, ONLY : pi
   IMPLICIT NONE
   real(DP) ::fxx2m3y2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fxx2m3y2 =  0.25d0*sqrt(17.5d0/pi) * sint * sint * sint * cos(3.d0*phi)
   RETURN
END FUNCTION fxx2m3y2
FUNCTION fy3x2my2(cost,phi)
   USE kinds, ONLY :  DP
   USE constants, ONLY : pi
   IMPLICIT NONE
   real(DP) ::fy3x2my2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fy3x2my2 =  0.25d0*sqrt(17.5d0/pi) * sint * sint * sint * sin(3.d0*phi)
   RETURN
END FUNCTION fy3x2my2
!
!
!-----------------------------------------------------------------------
SUBROUTINE radialpart(ng, q, alfa, rvalue, lmax, radial)
  !-----------------------------------------------------------------------
  !
  ! This routine computes a table with the radial Fourier transform
  ! of the radial functions.
  !
  USE kinds,      ONLY : dp
  USE constants,  ONLY : fpi
  USE cell_base,  ONLY : omega
  !
  IMPLICIT NONE
  ! I/O
  INTEGER :: ng, rvalue, lmax
  real(DP) :: q(ng), alfa, radial(ng,0:lmax)
  ! local variables
  real(DP), PARAMETER :: xmin=-6.d0, dx=0.025d0, rmax=10.d0

  real(DP) :: rad_int, pref, x
  INTEGER :: l, lp1, ir, ig, mesh_r
  real(DP), ALLOCATABLE :: bes(:), func_r(:), r(:), rij(:), aux(:)

  mesh_r = nint ( ( log ( rmax ) - xmin ) / dx + 1 )
  ALLOCATE ( bes(mesh_r), func_r(mesh_r), r(mesh_r), rij(mesh_r) )
  ALLOCATE ( aux(mesh_r))
  !
  !    compute the radial mesh
  !
  DO ir = 1, mesh_r
     x = xmin  + dble (ir - 1) * dx
     r (ir) = exp (x) / alfa
     rij (ir) = dx  * r (ir)
  ENDDO
  !
  IF (rvalue==1) func_r(:) = 2.d0 * alfa**(3.d0/2.d0) * exp(-alfa*r(:))
  IF (rvalue==2) func_r(:) = 1.d0/sqrt(8.d0) * alfa**(3.d0/2.d0) * &
                     (2.0d0 - alfa*r(:)) * exp(-alfa*r(:)*0.5d0)
  IF (rvalue==3) func_r(:) = sqrt(4.d0/27.d0) * alfa**(3.0d0/2.0d0) * &
                     (1.d0 - 2.0d0/3.0d0*alfa*r(:) + 2.d0*(alfa*r(:))**2/27.d0) * &
                                           exp(-alfa*r(:)/3.0d0)
  pref = fpi/sqrt(omega)
  !
  DO l = 0, lmax
     DO ig=1,ng
       CALL sph_bes (mesh_r, r(1), q(ig), l, bes)
       aux(:) = bes(:) * func_r(:) * r(:) * r(:)
       ! second r factor added upo suggestion by YY Liang
       CALL simpson (mesh_r, aux, rij, rad_int)
       radial(ig,l) = rad_int * pref
     ENDDO
  ENDDO

  DEALLOCATE (bes, func_r, r, rij, aux )
  RETURN
END SUBROUTINE radialpart
