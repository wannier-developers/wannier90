!-*- mode: F90 -*-!
!------------------------------------------------------------!
! This file is distributed as part of the Wannier90 code and !
! under the terms of the GNU General Public License. See the !
! file `LICENSE' in the root directory of the Wannier90      !
! distribution, or http://www.gnu.org/copyleft/gpl.txt       !
!                                                            !
! The webpage of the Wannier90 code is www.wannier.org       !
!                                                            !
! The Wannier90 code is hosted on GitHub:                    !
!                                                            !
! https://github.com/wannier-developers/wannier90            !
!------------------------------------------------------------!
!                                                            !
!  w90_kslice: properties evaluated on isosurface in BZ      !
!                                                            !
!------------------------------------------------------------!

module w90_kslice

  !! Plots the intersections of constant-energy isosurfaces with a BZ
  !! slice, and/or makes a heatmap plot on the slice:
  !!
  !!  - Minus the Berry curvature, summed over occupied bands
  !!
  !!  - The k-integrand of the orbital magnetization formula
  !!
  !!  - The k-integrand of the spin Hall conductivity formula
  !!
  !! The slice is defined in reduced coordinates by three input variables:
  !!
  !!    pw90_kslice_corner(1:3) is the lower left corner
  !!    pw90_kslice_b1(1:3) and pw90_kslice_b2(1:3) are the vectors subtending the slice

  use w90_error, only: w90_error_type, set_error_alloc, set_error_dealloc, set_error_fatal, &
    set_error_input, set_error_fatal, set_error_file

  implicit none

  private

  public :: k_slice

contains

  !================================================!
  !                   PUBLIC PROCEDURES
  !================================================!

  subroutine k_slice(pw90_berry, dis_manifold, fermi_energy_list, kmesh_info, kpt_latt, &
                     pw90_kslice, pw90_oper_read, pw90_band_deriv_degen, pw90_spin, ws_region, &
                     pw90_spin_hall, print_output, wannier_data, ws_distance, wigner_seitz, AA_R, &
                     BB_R, CC_R, HH_R, SH_R, SHR_R, SR_R, SS_R, SAA_R, SBB_R, v_matrix, u_matrix, &
                     bohr, eigval, real_lattice, scissors_shift, mp_grid, fermi_n, num_bands, &
                     num_kpts, num_wann, num_valence_bands, effective_model, have_disentangled, &
                     seedname, stdout, timer, error, comm)
    !================================================!
    !
    !! Main routine
    !
    !================================================!

    use w90_postw90_types, only: pw90_kslice_mod_type, pw90_berry_mod_type, pw90_spin_mod_type, &
      pw90_band_deriv_degen_type, pw90_oper_read_type, pw90_spin_hall_type, wigner_seitz_type
    use w90_berry, only: berry_get_imf_klist, berry_get_imfgh_klist, berry_get_shc_klist
    use w90_comms, only: comms_bcast, w90_comm_type, mpirank, mpisize, comms_gatherv, comms_array_split
    use w90_constants, only: dp, twopi, eps8
    use w90_get_oper, only: get_HH_R, get_AA_R_effective, get_AA_R, get_BB_R, get_CC_R, get_SS_R, get_SHC_R
    use w90_io, only: io_time
    use w90_types, only: dis_manifold_type, kmesh_info_type, print_output_type, &
      wannier_data_type, ws_region_type, ws_distance_type, timer_list_type
    use w90_postw90_common, only: pw90common_fourier_R_to_k
    use w90_spin, only: spin_get_nk
    use w90_utility, only: utility_diagonalize, utility_recip_lattice, utility_recip_lattice_base
    use w90_wan_ham, only: wham_get_eig_deleig

    implicit none

    ! arguments
    type(pw90_berry_mod_type), intent(in) :: pw90_berry
    type(dis_manifold_type), intent(in) :: dis_manifold
    real(kind=dp), allocatable, intent(in) :: fermi_energy_list(:)
    type(kmesh_info_type), intent(in) :: kmesh_info
    real(kind=dp), intent(in) :: kpt_latt(:, :)
    type(pw90_kslice_mod_type), intent(in) :: pw90_kslice
    type(pw90_band_deriv_degen_type), intent(in) :: pw90_band_deriv_degen
    type(pw90_oper_read_type), intent(in) :: pw90_oper_read
    type(pw90_spin_mod_type), intent(in) :: pw90_spin
    type(print_output_type), intent(in) :: print_output
    type(ws_region_type), intent(in) :: ws_region
    type(pw90_spin_hall_type), intent(in) :: pw90_spin_hall
    type(w90_comm_type), intent(in) :: comm
    type(wannier_data_type), intent(in) :: wannier_data
    type(wigner_seitz_type), intent(inout) :: wigner_seitz
    type(ws_distance_type), intent(inout) :: ws_distance
    type(timer_list_type), intent(inout) :: timer
    type(w90_error_type), allocatable, intent(out) :: error

    complex(kind=dp), allocatable, intent(inout) :: AA_R(:, :, :, :)
    complex(kind=dp), allocatable, intent(inout) :: BB_R(:, :, :, :)
    complex(kind=dp), allocatable, intent(inout) :: CC_R(:, :, :, :, :)
    complex(kind=dp), allocatable, intent(inout) :: HH_R(:, :, :)
    complex(kind=dp), allocatable, intent(inout) :: SH_R(:, :, :, :)
    complex(kind=dp), allocatable, intent(inout) :: SHR_R(:, :, :, :, :)
    complex(kind=dp), allocatable, intent(inout) :: SR_R(:, :, :, :, :)
    complex(kind=dp), allocatable, intent(inout) :: SS_R(:, :, :, :)
    complex(kind=dp), allocatable, intent(inout) :: SAA_R(:, :, :, :, :) ! <0n|sigma_x,y,z.(r-R)_alpha|Rm>
    complex(kind=dp), allocatable, intent(inout) :: SBB_R(:, :, :, :, :) ! <0n|sigma_x,y,z.H.(r-R)_alpha|Rm>
    complex(kind=dp), intent(in) :: v_matrix(:, :, :), u_matrix(:, :, :)

    real(kind=dp), intent(in) :: bohr
    real(kind=dp), intent(in) :: eigval(:, :)
    real(kind=dp), intent(in) :: real_lattice(3, 3)
    real(kind=dp), intent(in) :: scissors_shift

    integer, intent(in) :: mp_grid(3)
    integer, intent(in) :: num_bands, num_kpts, num_wann, num_valence_bands, fermi_n
    integer, intent(in) :: stdout

    character(len=50), intent(in) :: seedname
    logical, intent(in) :: have_disentangled
    logical, intent(in) :: effective_model

    ! local variables
    real(kind=dp)     :: recip_lattice(3, 3), volume
    integer           :: iloc, itot, i1, i2, n, n1, n2, n3, i, nkpts, my_nkpts
    integer           :: scriptunit, dataunit, loop_kpt
    real(kind=dp)     :: avec_2d(3, 3), bvec(3, 3), yvec(3), zvec(3), &
                         b1mod, b2mod, ymod, cosb1b2, kcorner_cart(3), &
                         areab1b2, cosyb2, kpt(3), kpt_x, kpt_y, k1, k2, &
                         imf_k_list(3, 3, fermi_n), img_k_list(3, 3, fermi_n), &
                         imh_k_list(3, 3, fermi_n), Morb_k(3, 3), curv(3), morb(3), &
                         spn_k(num_wann), del_eig(num_wann, 3), Delta_k, Delta_E, &
                         zhat(3), vdum(3), rdum, shc_k_fermi(fermi_n)
    logical           :: plot_fermi_lines, plot_curv, plot_morb, &
                         fermi_lines_color, heatmap, plot_shc
    character(len=120) :: filename, square

    complex(kind=dp), allocatable :: HH(:, :)
    complex(kind=dp), allocatable :: delHH(:, :, :)
    complex(kind=dp), allocatable :: UU(:, :)
    real(kind=dp), allocatable :: eig(:)

    ! Output data buffers
    real(kind=dp), allocatable    :: coords(:, :), my_coords(:, :), &
                                     spndata(:, :), my_spndata(:, :), &
                                     bandsdata(:, :), my_bandsdata(:, :), &
                                     zdata(:, :), my_zdata(:, :)
    logical, allocatable          :: spnmask(:, :), my_spnmask(:, :)

    integer, allocatable :: counts(:), displs(:)
    logical :: on_root = .false.
    integer :: my_node_id, num_nodes

    my_node_id = mpirank(comm)
    num_nodes = mpisize(comm)
    allocate (counts(0:num_nodes - 1))
    allocate (displs(0:num_nodes - 1))
    if (my_node_id == 0) on_root = .true.

    plot_fermi_lines = index(pw90_kslice%task, 'fermi_lines') > 0
    plot_curv = index(pw90_kslice%task, 'curv') > 0
    plot_morb = index(pw90_kslice%task, 'morb') > 0
    plot_shc = index(pw90_kslice%task, 'shc') > 0
    fermi_lines_color = pw90_kslice%fermi_lines_colour /= 'none'
    heatmap = plot_curv .or. plot_morb .or. plot_shc
    if (plot_fermi_lines .and. fermi_lines_color .and. heatmap) then
      call set_error_input(error, 'Error: spin-colored Fermi lines not allowed in ' &
                           //'curv/morb/shc heatmap plots', comm)
      return
    end if
    if (plot_shc) then
      if (pw90_berry%kubo_smearing%use_adaptive) then
        call set_error_input(error, 'Error: Must use fixed smearing when plotting spin Hall conductivity', comm)
        return
      end if
      if (fermi_n == 0) then
        call set_error_input(error, 'Error: must specify Fermi energy', comm)
        return
      else if (fermi_n /= 1) then
        call set_error_input(error, 'Error: kpath plot only accept one Fermi energy, ' &
                             //'use fermi_energy instead of fermi_energy_min', comm)
        return
      end if
    end if

    if (on_root) then
      call kslice_print_info(plot_fermi_lines, fermi_lines_color, plot_curv, plot_morb, plot_shc, &
                             stdout, pw90_berry, fermi_energy_list, error, comm)
      if (allocated(error)) return
    end if

    call get_HH_R(dis_manifold, kpt_latt, print_output, wigner_seitz, HH_R, u_matrix, v_matrix, &
                  eigval, real_lattice, scissors_shift, num_bands, num_kpts, num_wann, &
                  num_Valence_bands, effective_model, have_disentangled, seedname, ws_distance, ws_region, &
                  stdout, timer, error, comm)
    if (allocated(error)) return

    if (plot_curv .or. plot_morb) then
      if (effective_model) then
        call get_AA_R_effective(print_output, AA_R, HH_R, wigner_seitz%nrpts, num_wann, seedname, &
                                stdout, timer, error, comm)
      else
        call get_AA_R(pw90_berry, dis_manifold, kmesh_info, kpt_latt, print_output, wannier_data, AA_R, &
                      v_matrix, eigval, wigner_seitz, ws_distance, ws_region, num_bands, num_kpts, &
                      num_wann, have_disentangled, seedname, stdout, timer, error, comm)
      endif
      if (allocated(error)) return

    endif
    if (plot_morb) then
      call get_BB_R(pw90_berry, dis_manifold, kmesh_info, kpt_latt, print_output, HH_R, BB_R, v_matrix, &
                    eigval, scissors_shift, wigner_seitz, ws_distance, ws_region, num_bands, num_kpts, &
                    num_wann, have_disentangled, seedname, stdout, timer, error, comm)
      if (allocated(error)) return

      call get_CC_R(pw90_berry, dis_manifold, kmesh_info, kpt_latt, print_output, pw90_oper_read, &
                    HH_R, BB_R, CC_R, v_matrix, eigval, scissors_shift, wigner_seitz, ws_distance, &
                    ws_region, num_bands, num_kpts, num_wann, have_disentangled, seedname, stdout, &
                    timer, error, comm)
      if (allocated(error)) return

    endif

    if (plot_shc) then
      if (effective_model) then
        call get_AA_R_effective(print_output, AA_R, HH_R, wigner_seitz%nrpts, num_wann, seedname, &
                                stdout, timer, error, comm)
      else
        call get_AA_R(pw90_berry, dis_manifold, kmesh_info, kpt_latt, print_output, wannier_data, AA_R, &
                      v_matrix, eigval, wigner_seitz, ws_distance, ws_region, num_bands, num_kpts, &
                      num_wann, have_disentangled, seedname, stdout, timer, error, comm)
      endif
      if (allocated(error)) return

      call get_SS_R(dis_manifold, kpt_latt, print_output, pw90_oper_read, SS_R, v_matrix, eigval, &
                    wigner_seitz, ws_distance, ws_region, num_bands, num_kpts, num_wann, &
                    have_disentangled, seedname, stdout, timer, error, comm)
      if (allocated(error)) return

      call get_SHC_R(dis_manifold, kmesh_info, kpt_latt, print_output, pw90_oper_read, &
                     pw90_spin_hall, SH_R, SHR_R, SR_R, v_matrix, eigval, scissors_shift, &
                     wigner_seitz, ws_distance, ws_region, num_bands, num_kpts, num_wann, &
                     num_valence_bands, have_disentangled, seedname, stdout, timer, error, comm)
      if (allocated(error)) return

    end if

    if (fermi_lines_color) then
      call get_SS_R(dis_manifold, kpt_latt, print_output, pw90_oper_read, SS_R, v_matrix, eigval, &
                    wigner_seitz, ws_distance, ws_region, num_bands, num_kpts, num_wann, &
                    have_disentangled, seedname, stdout, timer, error, comm)
      if (allocated(error)) return

    endif
    call utility_recip_lattice_base(real_lattice, recip_lattice, volume)
    ! Set Cartesian components of the vectors (b1,b2) spanning the slice

    bvec(1, :) = matmul(pw90_kslice%b1(:), recip_lattice(:, :))
    bvec(2, :) = matmul(pw90_kslice%b2(:), recip_lattice(:, :))
    ! z_vec (orthogonal to the slice)
    zvec(1) = bvec(1, 2)*bvec(2, 3) - bvec(1, 3)*bvec(2, 2)
    zvec(2) = bvec(1, 3)*bvec(2, 1) - bvec(1, 1)*bvec(2, 3)
    zvec(3) = bvec(1, 1)*bvec(2, 2) - bvec(1, 2)*bvec(2, 1)
    ! y_vec (orthogonal to b1=x_vec)
    yvec(1) = zvec(2)*bvec(1, 3) - zvec(3)*bvec(1, 2)
    yvec(2) = zvec(3)*bvec(1, 1) - zvec(1)*bvec(1, 3)
    yvec(3) = zvec(1)*bvec(1, 2) - zvec(2)*bvec(1, 1)
    ! Area (modulus b1 x b2 = z_vec)
    areab1b2 = sqrt(zvec(1)**2 + zvec(2)**2 + zvec(3)**2)
    if (areab1b2 < eps8) then
      call set_error_fatal(error, 'Error in kslice: Vectors pw90_kslice_b1 and pw90_kslice_b2 ' &
                           //'not linearly independent', comm)
      return
    endif
    ! This is the unit vector zvec/|zvec| which completes the triad
    ! in the 2D case
    bvec(3, :) = zvec(:)/areab1b2
    ! Now that we have bvec(3,:), we can compute the dual vectors
    ! avec_2d as in the 3D case
    call utility_recip_lattice(bvec, avec_2d, rdum, error, comm)
    if (allocated(error)) return

    ! Moduli b1,b2,y_vec
    b1mod = sqrt(bvec(1, 1)**2 + bvec(1, 2)**2 + bvec(1, 3)**2)
    b2mod = sqrt(bvec(2, 1)**2 + bvec(2, 2)**2 + bvec(2, 3)**2)
    ymod = sqrt(yvec(1)**2 + yvec(2)**2 + yvec(3)**2)
    ! Cosine of the angle between y_vec and b2
    cosyb2 = yvec(1)*bvec(2, 1) + yvec(2)*bvec(2, 2) + yvec(3)*bvec(2, 3)
    cosyb2 = cosyb2/(ymod*b2mod)
    ! Cosine of the angle between b1=x_vec and b2
    cosb1b2 = bvec(1, 1)*bvec(2, 1) + bvec(1, 2)*bvec(2, 2) + bvec(1, 3)*bvec(2, 3)
    cosb1b2 = cosb1b2/(b1mod*b2mod)
    if (abs(cosb1b2) < eps8 .and. abs(b1mod - b2mod) < eps8) then
      square = 'True'
    else
      square = 'False'
    end if

    nkpts = (pw90_kslice%kmesh2d(1) + 1)*(pw90_kslice%kmesh2d(2) + 1)

    ! Partition set of k-points into junks
    call comms_array_split(nkpts, counts, displs, comm)
    my_nkpts = counts(my_node_id)

    allocate (my_coords(2, my_nkpts))
    if (heatmap) allocate (my_zdata(3, my_nkpts))

    if (plot_fermi_lines) then
      allocate (HH(num_wann, num_wann))
      allocate (UU(num_wann, num_wann))
      allocate (eig(num_wann))
      if (fermi_lines_color) then
        allocate (delHH(num_wann, num_wann, 3))
        allocate (my_spndata(num_wann, my_nkpts))
        allocate (my_spnmask(num_wann, my_nkpts))
        my_spnmask = .false.
      else
        allocate (my_bandsdata(num_wann, my_nkpts))
      endif
    endif

    ! Loop over local portion of uniform mesh of k-points covering the slice,
    ! including all four borders

    do iloc = 1, my_nkpts
      itot = iloc - 1 + displs(my_node_id)
      i2 = itot/(pw90_kslice%kmesh2d(1) + 1) ! slow
      i1 = itot - i2*(pw90_kslice%kmesh2d(1) + 1) !fast
      ! k1 and k2 are the coefficients of the k-point in the basis
      ! (pw90_kslice_b1,pw90_kslice_b2)
      k1 = i1/real(pw90_kslice%kmesh2d(1), dp)
      k2 = i2/real(pw90_kslice%kmesh2d(2), dp)
      kpt = pw90_kslice%corner + k1*pw90_kslice%b1 + k2*pw90_kslice%b2
      ! Add to (k1,k2) the projection of pw90_kslice_corner on the
      ! (pw90_kslice_b1,pw90_kslice_b2) plane, expressed as a linear
      ! combination of pw90_kslice_b1 and pw90_kslice_b2
      kcorner_cart(:) = matmul(pw90_kslice%corner(:), recip_lattice(:, :))
      k1 = k1 + dot_product(kcorner_cart, avec_2d(1, :))/twopi
      k2 = k2 + dot_product(kcorner_cart, avec_2d(2, :))/twopi
      ! Convert to (kpt_x,kpt_y), the 2D Cartesian coordinates
      ! with x along x_vec=b1 and y along y_vec
      kpt_x = k1*b1mod + k2*b2mod*cosb1b2
      kpt_y = k2*b2mod*cosyb2

      my_coords(:, iloc) = [kpt_x, kpt_y]

      if (plot_fermi_lines) then
        if (fermi_lines_color) then
          call spin_get_nk(ws_region, pw90_spin, wannier_data, ws_distance, wigner_seitz, HH_R, &
                           SS_R, kpt, real_lattice, spn_k, mp_grid, num_wann, error, comm)
          if (allocated(error)) return

          do n = 1, num_wann
            if (spn_k(n) > 1.0_dp - eps8) then
              spn_k(n) = 1.0_dp - eps8
            elseif (spn_k(n) < -1.0_dp + eps8) then
              spn_k(n) = -1.0_dp + eps8
            endif
          enddo

          call wham_get_eig_deleig(dis_manifold, kpt_latt, pw90_band_deriv_degen, ws_region, &
                                   print_output, wannier_data, ws_distance, wigner_seitz, delHH, &
                                   HH, HH_R, u_matrix, UU, v_matrix, del_eig, eig, eigval, kpt, &
                                   real_lattice, scissors_shift, mp_grid, num_bands, num_kpts, &
                                   num_wann, num_valence_bands, effective_model, &
                                   have_disentangled, seedname, stdout, timer, error, comm)
          if (allocated(error)) return

          Delta_k = max(b1mod/pw90_kslice%kmesh2d(1), b2mod/pw90_kslice%kmesh2d(2))
        else
          call pw90common_fourier_R_to_k(ws_region, wannier_data, ws_distance, wigner_seitz, HH, &
                                         HH_R, kpt, real_lattice, mp_grid, 0, num_wann, error, comm)
          if (allocated(error)) return

          call utility_diagonalize(HH, num_wann, eig, UU, error, comm)
          if (allocated(error)) return

        endif

        if (allocated(my_bandsdata)) then
          my_bandsdata(:, iloc) = eig(:)
        else
          my_spndata(:, iloc) = spn_k(:)
          do n = 1, num_wann
            ! vdum = dE/dk projected on the k-slice
            zhat = zvec/sqrt(dot_product(zvec, zvec))
            vdum(:) = del_eig(n, :) - dot_product(del_eig(n, :), zhat)*zhat(:)
            Delta_E = sqrt(dot_product(vdum, vdum))*Delta_k
!                   Delta_E=Delta_E*sqrt(2.0_dp) ! optimize this factor
            my_spnmask(n, iloc) = abs(eig(n) - fermi_energy_list(1)) < Delta_E
          end do
        end if
      end if

      if (plot_curv) then

        call berry_get_imf_klist(dis_manifold, fermi_energy_list, kpt_latt, ws_region, print_output, &
                                 wannier_data, ws_distance, wigner_seitz, AA_R, BB_R, CC_R, HH_R, &
                                 u_matrix, v_matrix, eigval, kpt, real_lattice, imf_k_list, &
                                 scissors_shift, mp_grid, num_bands, num_kpts, num_wann, &
                                 num_valence_bands, effective_model, have_disentangled, seedname, &
                                 stdout, timer, error, comm)
        if (allocated(error)) return

        curv(1) = sum(imf_k_list(:, 1, 1))
        curv(2) = sum(imf_k_list(:, 2, 1))
        curv(3) = sum(imf_k_list(:, 3, 1))
        if (pw90_berry%curv_unit == 'bohr2') curv = curv/bohr**2
        ! Print _minus_ the Berry curvature
        my_zdata(:, iloc) = -curv(:)
      else if (plot_morb) then
        call berry_get_imfgh_klist(dis_manifold, fermi_energy_list, kpt_latt, ws_region, &
                                   print_output, wannier_data, ws_distance, wigner_seitz, AA_R, &
                                   BB_R, CC_R, HH_R, u_matrix, v_matrix, eigval, kpt, &
                                   real_lattice, scissors_shift, mp_grid, fermi_n, num_bands, &
                                   num_kpts, num_wann, num_valence_bands, effective_model, &
                                   have_disentangled, seedname, stdout, timer, error, comm, &
                                   imf_k_list, img_k_list, imh_k_list)
        if (allocated(error)) return

        Morb_k = img_k_list(:, :, 1) + imh_k_list(:, :, 1) &
                 - 2.0_dp*fermi_energy_list(1)*imf_k_list(:, :, 1)
        Morb_k = -Morb_k/2.0_dp ! differs by -1/2 from Eq.97 LVTS12
        morb(1) = sum(Morb_k(:, 1))
        morb(2) = sum(Morb_k(:, 2))
        morb(3) = sum(Morb_k(:, 3))
        my_zdata(:, iloc) = morb(:)
      else if (plot_shc) then
        call berry_get_shc_klist(pw90_berry, dis_manifold, fermi_energy_list, kpt_latt, &
                                 pw90_band_deriv_degen, ws_region, pw90_spin_hall, print_output, &
                                 wannier_data, ws_distance, wigner_seitz, AA_R, HH_R, SH_R, SHR_R, &
                                 SR_R, SS_R, SAA_R, SBB_R, u_matrix, v_matrix, eigval, kpt, real_lattice, &
                                 scissors_shift, mp_grid, fermi_n, num_bands, num_kpts, num_wann, &
                                 num_valence_bands, effective_model, have_disentangled, seedname, &
                                 stdout, timer, error, comm, shc_k_fermi=shc_k_fermi)
        if (allocated(error)) return

        my_zdata(1, iloc) = shc_k_fermi(1)
      end if

    end do !iloc

    ! Send results to root process
    if (on_root) then
      allocate (coords(2, nkpts))
    else
      allocate (coords(1, 1))
    end if
    call comms_gatherv(my_coords, 2*my_nkpts, coords, 2*counts, 2*displs, error, comm)
    if (allocated(error)) return

    if (allocated(my_spndata)) then
      if (on_root) then
        allocate (spndata(num_wann, nkpts))
      else
        allocate (spndata(1, 1))
      end if
      call comms_gatherv(my_spndata, num_wann*my_nkpts, &
                         spndata, num_wann*counts, num_wann*displs, error, comm)
      if (allocated(error)) return
    end if

    if (allocated(my_spnmask)) then
      if (on_root) then
        allocate (spnmask(num_wann, nkpts))
      else
        allocate (spnmask(1, 1))
      end if
      call comms_gatherv(my_spnmask(1, 1), num_wann*my_nkpts, &
                         spnmask(1, 1), num_wann*counts, num_wann*displs, error, comm)
      if (allocated(error)) return
    end if

    if (allocated(my_bandsdata)) then
      if (on_root) then
        allocate (bandsdata(num_wann, nkpts))
      else
        allocate (bandsdata(1, 1))
      end if
      call comms_gatherv(my_bandsdata, num_wann*my_nkpts, &
                         bandsdata, num_wann*counts, num_wann*displs, error, comm)
      if (allocated(error)) return
    end if

    ! This holds either -curv or morb
    if (allocated(my_zdata)) then
      if (on_root) then
        allocate (zdata(3, nkpts))
      else
        allocate (zdata(1, 1))
      end if
      call comms_gatherv(my_zdata, 3*my_nkpts, &
                         zdata, 3*counts, 3*displs, error, comm)
      if (allocated(error)) return
    end if

    ! Write output files
    if (on_root) then
      ! set kpt_x and kpt_y to last evaluated point
      kpt_x = coords(1, nkpts)
      kpt_y = coords(2, nkpts)

      write (stdout, '(/,/,1x,a)') 'Output files:'

      if (.not. fermi_lines_color) then
        filename = trim(seedname)//'-kslice-coord.dat'
        call write_data_file(stdout, filename, '(2E16.8)', coords)
      end if

      if (allocated(bandsdata)) then
        ! For python
        filename = trim(seedname)//'-kslice-bands.dat'
        call write_data_file(stdout, filename, '(E16.8)', &
                             reshape(bandsdata, [1, nkpts*num_wann]))

        ! For gnuplot, using 'grid data' format
        if (.not. heatmap) then
          do n = 1, num_wann
            n1 = n/100
            n2 = (n - n1*100)/10
            n3 = n - n1*100 - n2*10
            filename = trim(seedname)//'-bnd_' &
                       //achar(48 + n1)//achar(48 + n2)//achar(48 + n3)//'.dat'

            call write_coords_file(stdout, filename, '(3E16.8)', coords, &
                                   reshape(bandsdata(n, :), [1, 1, nkpts]), &
                                   blocklen=pw90_kslice%kmesh2d(1) + 1)
          enddo
        endif
      end if

      if (allocated(spndata)) then
        filename = trim(seedname)//'-kslice-fermi-spn.dat'
        call write_coords_file(stdout, filename, '(3E16.8)', coords, &
                               reshape(spndata, [1, num_wann, nkpts]), &
                               spnmask)
      end if

      if (allocated(my_zdata)) then
        if (plot_curv .or. plot_morb .or. plot_shc) then
          if (plot_morb) then ! ugly. But to keep the logic the same as other places
            filename = trim(seedname)//'-kslice-morb.dat'
          elseif (plot_curv) then
            filename = trim(seedname)//'-kslice-curv.dat'
          elseif (plot_shc) then
            filename = trim(seedname)//'-kslice-shc.dat'
          endif
          write (stdout, '(/,3x,a)') filename
          open (newunit=dataunit, file=filename, form='formatted')
          if (plot_shc) then
            if (pw90_berry%curv_unit == 'bohr2') zdata = zdata/bohr**2
            do loop_kpt = 1, nkpts
              write (dataunit, '(1E16.8)') zdata(1, loop_kpt)
            end do
          else
            do loop_kpt = 1, nkpts
              write (dataunit, '(4E16.8)') zdata(:, loop_kpt)
            end do
          end if
          write (dataunit, *) ' '
          close (dataunit)
        end if
      endif

      if (plot_fermi_lines .and. .not. fermi_lines_color .and. .not. heatmap) then
        !
        ! gnuplot script for black Fermi lines
        !
        filename = trim(seedname)//'-kslice-fermi_lines.gnu'
        write (stdout, '(/,3x,a)') filename
        open (newunit=scriptunit, file=filename, form='formatted')
        write (scriptunit, '(a)') "unset surface"
        write (scriptunit, '(a)') "set contour"
        write (scriptunit, '(a)') "set view map"
        write (scriptunit, '(a,f9.5)') "set cntrparam levels discrete ", &
          fermi_energy_list(1)
        write (scriptunit, '(a)') "set cntrparam bspline"
        do n = 1, num_wann
          n1 = n/100
          n2 = (n - n1*100)/10
          n3 = n - n1*100 - n2*10
          write (scriptunit, '(a)') "set table 'bnd_" &
            //achar(48 + n1)//achar(48 + n2)//achar(48 + n3)//".dat'"
          write (scriptunit, '(a)') "splot '"//trim(seedname)//"-bnd_" &
            //achar(48 + n1)//achar(48 + n2)//achar(48 + n3)//".dat'"
          write (scriptunit, '(a)') "unset table"
        enddo
        write (scriptunit, '(a)') &
          "#Uncomment next two lines to create postscript"
        write (scriptunit, '(a)') "#set term post eps enh"
        write (scriptunit, '(a)') &
          "#set output '"//trim(seedname)//"-kslice-fermi_lines.eps'"
        write (scriptunit, '(a)') "set size ratio -1"
        write (scriptunit, '(a)') "unset tics"
        write (scriptunit, '(a)') "unset key"
        write (scriptunit, '(a)') &
          "#For postscript try changing lw 1 --> lw 2 in the next line"
        write (scriptunit, '(a)') "set style line 1 lt 1 lw 1"
        if (num_wann == 1) then
          write (scriptunit, '(a)') &
            "plot 'bnd_001.dat' using 1:2 w lines ls 1"
        else
          write (scriptunit, '(a)') &
            "plot 'bnd_001.dat' using 1:2 w lines ls 1,"//achar(92)
        endif
        do n = 2, num_wann - 1
          n1 = n/100
          n2 = (n - n1*100)/10
          n3 = n - n1*100 - n2*10
          write (scriptunit, '(a)') "     'bnd_" &
            //achar(48 + n1)//achar(48 + n2)//achar(48 + n3) &
            //".dat' using 1:2 w lines ls 1,"//achar(92)
        enddo
        n = num_wann
        n1 = n/100
        n2 = (n - n1*100)/10
        n3 = n - n1*100 - n2*10
        write (scriptunit, '(a)') "     'bnd_" &
          //achar(48 + n1)//achar(48 + n2)//achar(48 + n3) &
          //".dat' using 1:2 w lines ls 1"
        close (scriptunit)
        !
        ! Python script for black Fermi lines
        !
        filename = trim(seedname)//'-kslice-fermi_lines.py'
        write (stdout, '(/,3x,a)') filename
        open (newunit=scriptunit, file=filename, form='formatted')
        call script_common(scriptunit, areab1b2, square, seedname)
        call script_fermi_lines(scriptunit, seedname, fermi_energy_list)
        write (scriptunit, '(a)') " "
        write (scriptunit, '(a)') "# Remove the axes"
        write (scriptunit, '(a)') "ax = pl.gca()"
        write (scriptunit, '(a)') "ax.xaxis.set_visible(False)"
        write (scriptunit, '(a)') "ax.yaxis.set_visible(False)"
        write (scriptunit, '(a)') " "
        write (scriptunit, '(a)') "pl.axes().set_aspect('equal')"
        write (scriptunit, '(a)') " "
        write (scriptunit, '(a)') "outfile = '"//trim(seedname)// &
          "-fermi_lines.pdf'"
        write (scriptunit, '(a)') " "
        write (scriptunit, '(a)') " "
        write (scriptunit, '(a)') "pl.savefig(outfile,bbox_inches='tight')"
        write (scriptunit, '(a)') "pl.show()"
        close (scriptunit)
      endif !plot_fermi_lines .and. .not.fermi_lines_color .and. .not.heatmap

      if (plot_fermi_lines .and. fermi_lines_color .and. .not. heatmap) then
        !
        ! gnuplot script for spin-colored Fermi lines
        !
        filename = trim(seedname)//'-kslice-fermi_lines.gnu'
        write (stdout, '(/,3x,a)') filename
        open (newunit=scriptunit, file=filename, form='formatted')
        write (scriptunit, '(a)') "unset key"
        write (scriptunit, '(a)') "unset tics"
        write (scriptunit, '(a)') "set cbtics"
        write (scriptunit, '(a)') &
          "set palette defined (-1 'blue', 0 'green', 1 'red')"
        write (scriptunit, '(a)') "set pm3d map"
        write (scriptunit, '(a)') "set zrange [-1:1]"
        write (scriptunit, '(a)') "set size ratio -1"
        write (scriptunit, '(a)') &
          "#Uncomment next two lines to create postscript"
        write (scriptunit, '(a)') "#set term post eps enh"
        write (scriptunit, '(a)') "#set output '" &
          //trim(seedname)//"-kslice-fermi_lines.eps'"
        write (scriptunit, '(a)') "splot '" &
          //trim(seedname)//"-kslice-fermi-spn.dat' with dots palette"
        !
        ! python script for spin-colored Fermi lines
        !
        filename = trim(seedname)//'-kslice-fermi_lines.py'
        write (stdout, '(/,3x,a)') filename
        open (newunit=scriptunit, file=filename, form='formatted')
        write (scriptunit, '(a)') "import pylab as pl"
        write (scriptunit, '(a)') "import numpy as np"
        write (scriptunit, '(a)') "data = np.loadtxt('"//trim(seedname)// &
          "-kslice-fermi-spn.dat')"
        write (scriptunit, '(a)') "x=data[:,0]"
        write (scriptunit, '(a)') "y=data[:,1]"
        write (scriptunit, '(a)') "z=data[:,2]"
        write (scriptunit, '(a)') &
          "pl.scatter(x,y,c=z,marker='+',s=2,cmap=pl.cm.jet)"
        write (scriptunit, '(a,F12.6,a)') &
          "pl.plot([0,", kpt_x, "],[0,0],color='black',linestyle='-'," &
          //"linewidth=0.5)"
        write (scriptunit, '(a,F12.6,a,F12.6,a,F12.6,a)') &
          "pl.plot([", kpt_x, ",", kpt_x, "],[0,", kpt_y, "],color='black'," &
          //"linestyle='-',linewidth=0.5)"
        write (scriptunit, '(a,F12.6,a,F12.6,a,F12.6,a)') &
          "pl.plot([0,", kpt_x, "],[", kpt_y, ",", kpt_y, &
          "],color='black',linestyle='-',linewidth=0.5)"
        write (scriptunit, '(a,F12.6,a)') "pl.plot([0,0],[0,", kpt_y, &
          "],color='black',linestyle='-',linewidth=0.5)"
        write (scriptunit, '(a,F12.6,a)') "pl.xlim([0,", kpt_x, "])"
        write (scriptunit, '(a,F12.6,a)') "pl.ylim([0,", kpt_y, "])"
        write (scriptunit, '(a)') "cbar=pl.colorbar()"
        write (scriptunit, '(a)') "ax = pl.gca()"
        write (scriptunit, '(a)') "ax.xaxis.set_visible(False)"
        write (scriptunit, '(a)') "ax.yaxis.set_visible(False)"
        write (scriptunit, '(a)') "pl.savefig('"//trim(seedname)// &
          "-kslice-fermi_lines.pdf',bbox_inches='tight')"
        write (scriptunit, '(a)') "pl.show()"
        close (scriptunit)
      endif ! plot_fermi_lines .and. fermi_lines_color .and. .not.heatmap

      if (heatmap .and. (.not. plot_shc)) then
        !
        ! python script for curvature/Morb/SHC heatmaps [+ black Fermi lines]
        !
        do i = 1, 3

          if (plot_curv .and. .not. plot_fermi_lines) then
            filename = trim(seedname)//'-kslice-curv_'//achar(119 + i)//'.py'
            write (stdout, '(/,3x,a)') filename
            open (newunit=scriptunit, file=filename, form='formatted')
          elseif (plot_curv .and. plot_fermi_lines) then
            filename = trim(seedname)//'-kslice-curv_'//achar(119 + i)// &
                       '+fermi_lines.py'
            write (stdout, '(/,3x,a)') filename
            open (newunit=scriptunit, file=filename, form='formatted')
          elseif (plot_morb .and. .not. plot_fermi_lines) then
            filename = trim(seedname)//'-kslice-morb_'//achar(119 + i)//'.py'
            write (stdout, '(/,3x,a)') filename
            open (newunit=scriptunit, file=filename, form='formatted')
          elseif (plot_morb .and. plot_fermi_lines) then
            filename = trim(seedname)//'-kslice-morb_'//achar(119 + i)// &
                       '+fermi_lines.py'
            write (stdout, '(/,3x,a)') filename
            open (newunit=scriptunit, file=filename, form='formatted')
          endif
          call script_common(scriptunit, areab1b2, square, seedname)
          if (plot_fermi_lines) call script_fermi_lines(scriptunit, seedname, fermi_energy_list)

          if (plot_curv) then
            write (scriptunit, '(a)') " "
            write (scriptunit, '(a)') "outfile = '"//trim(seedname)// &
              "-kslice-curv_"//achar(119 + i)//".pdf'"
            write (scriptunit, '(a)') " "
            write (scriptunit, '(a)') &
              "val = np.loadtxt('"//trim(seedname)// &
              "-kslice-curv.dat', usecols=("//achar(47 + i)//",))"
            write (scriptunit, '(a)') " "
            write (scriptunit, '(a)') &
                 "val_log=np.array([np.log10(abs(elem))*np.sign(elem) &
                 &if abs(elem)>10 else elem/10.0 for elem in val])"
            write (scriptunit, '(a)') " "
            write (scriptunit, '(a)') "if square: "
            write (scriptunit, '(a)') "  Z=val_log.reshape(dimy,dimx)"
            write (scriptunit, '(a)') "  mn=int(np.floor(Z.min()))"
            write (scriptunit, '(a)') "  mx=int(np.ceil(Z.max()))"
            write (scriptunit, '(a)') "  ticks=range(mn,mx+1)"
            write (scriptunit, '(a)') "  pl.contourf(x_coord,y_coord,Z," &
              //"ticks,origin='lower')"
            write (scriptunit, '(a)') "  #pl.imshow(Z,origin='lower'," &
              //"extent=(min(x_coord),max(x_coord),min(y_coord)," &
              //"max(y_coord)))"
            write (scriptunit, '(a)') "else: "
            write (scriptunit, '(a)') "  valint = interpolate.griddata((points_x," &
              //"points_y), val_log, (grid_x,grid_y), method='nearest') # or 'cubic' or 'nearest'"
            write (scriptunit, '(a)') "  mn=int(np.floor(valint.min()))"
            write (scriptunit, '(a)') "  mx=int(np.ceil(valint.max()))"
            write (scriptunit, '(a)') "  ticks=range(mn,mx+1)"
            write (scriptunit, '(a)') "  pl.contourf(xint,yint,valint,ticks)"
            write (scriptunit, '(a)') "  #pl.imshow(valint,origin='lower'," &
              //"extent=(min(xint),max(xint),min(yint),max(yint)))"
            write (scriptunit, '(a)') " "
            write (scriptunit, '(a)') "ticklabels=[]"
            write (scriptunit, '(a)') "for n in ticks:"
            write (scriptunit, '(a)') " if n<0: "
            write (scriptunit, '(a)') &
              "  ticklabels.append('-$10^{%d}$' % abs(n))"
            write (scriptunit, '(a)') " elif n==0:"
            write (scriptunit, '(a)') "  ticklabels.append(' $%d$' %  n)"
            write (scriptunit, '(a)') " else:"
            write (scriptunit, '(a)') "  ticklabels.append(' $10^{%d}$' % n)"
            write (scriptunit, '(a)') " "
            write (scriptunit, '(a)') "cbar=pl.colorbar()"
            write (scriptunit, '(a)') "cbar.set_ticks(ticks)"
            write (scriptunit, '(a)') "cbar.set_ticklabels(ticklabels)"

          elseif (plot_morb) then

            write (scriptunit, '(a)') " "
            write (scriptunit, '(a)') "outfile = '"//trim(seedname)// &
              "-kslice-morb_"//achar(119 + i)//".pdf'"
            write (scriptunit, '(a)') " "
            write (scriptunit, '(a)') &
              "val = np.loadtxt('"//trim(seedname)// &
              "-kslice-morb.dat', usecols=("//achar(47 + i)//",))"
            write (scriptunit, '(a)') " "
            write (scriptunit, '(a)') "if square: "
            write (scriptunit, '(a)') "  Z=val.reshape(dimy,dimx)"
            write (scriptunit, '(a)') "  pl.imshow(Z,origin='lower'," &
              //"extent=(min(x_coord),max(x_coord),min(y_coord)," &
              //"max(y_coord)))"
            write (scriptunit, '(a)') "else: "
            write (scriptunit, '(a)') "  valint = interpolate.griddata((points_x," &
              //"points_y), val_log, (grid_x,grid_y), method='nearest') # or 'cubic' or 'nearest'"
            write (scriptunit, '(a)') "  pl.imshow(valint,origin='lower'," &
              //"extent=(min(xint),max(xint),min(yint),max(yint)))"
            write (scriptunit, '(a)') "cbar=pl.colorbar()"

          endif

          write (scriptunit, '(a)') " "
          write (scriptunit, '(a)') "ax = pl.gca()"
          write (scriptunit, '(a)') "ax.xaxis.set_visible(False)"
          write (scriptunit, '(a)') "ax.yaxis.set_visible(False)"
          write (scriptunit, '(a)') " "
          write (scriptunit, '(a)') "pl.savefig(outfile,bbox_inches='tight')"
          write (scriptunit, '(a)') "pl.show()"

          close (scriptunit)

        enddo !i

      endif !heatmap

      if (heatmap .and. plot_shc) then
        if (.not. plot_fermi_lines) then
          filename = trim(seedname)//'-kslice-shc'//'.py'
          write (stdout, '(/,3x,a)') filename
          open (newunit=scriptunit, file=filename, form='formatted')
        elseif (plot_fermi_lines) then
          filename = trim(seedname)//'-kslice-shc'//'+fermi_lines.py'
          write (stdout, '(/,3x,a)') filename
          open (scriptunit, file=filename, form='formatted')
        endif
        write (scriptunit, '(a)') "# uncomment these two lines if you are " &
          //"running in non-GUI environment"
        write (scriptunit, '(a)') "#import matplotlib"
        write (scriptunit, '(a)') "#matplotlib.use('Agg')"
        write (scriptunit, '(a)') "import matplotlib.pyplot as plt"
        call script_common(scriptunit, areab1b2, square, seedname)
        if (plot_fermi_lines) call script_fermi_lines(scriptunit, seedname, fermi_energy_list)

        write (scriptunit, '(a)') " "
        write (scriptunit, '(a)') "def shiftedColorMap(cmap, start=0, " &
          //"midpoint=0.5, stop=1.0, name='shiftedcmap'):"
        write (scriptunit, '(a)') "  '''"
        write (scriptunit, '(a)') '  Function to offset the "center" ' &
          //'of a colormap. Useful for'
        write (scriptunit, '(a)') '  data with a negative min and ' &
          //'positive max and you want the'
        write (scriptunit, '(a)') "  middle of the colormap's dynamic " &
          //"range to be at zero."
        write (scriptunit, '(a)') '  '
        write (scriptunit, '(a)') '  Input'
        write (scriptunit, '(a)') '  -----'
        write (scriptunit, '(a)') '  cmap : The matplotlib colormap to ' &
          //'be altered'
        write (scriptunit, '(a)') "  start : Offset from lowest point in " &
          //"the colormap's range."
        write (scriptunit, '(a)') '    Defaults to 0.0 (no lower offset). ' &
          //'Should be between'
        write (scriptunit, '(a)') '    0.0 and `midpoint`.'
        write (scriptunit, '(a)') '  midpoint : The new center of the ' &
          //'colormap. Defaults to '
        write (scriptunit, '(a)') '    0.5 (no shift). Should be between ' &
          //'0.0 and 1.0. In'
        write (scriptunit, '(a)') '    general, this should be  1 - ' &
          //'vmax / (vmax + abs(vmin))'
        write (scriptunit, '(a)') '    For example if your data range from ' &
          //'-15.0 to +5.0 and'
        write (scriptunit, '(a)') '    you want the center of the colormap ' &
          //'at 0.0, `midpoint`'
        write (scriptunit, '(a)') '    should be set to  1 - 5/(5 + 15)) ' &
          //'or 0.75'
        write (scriptunit, '(a)') "  stop : Offset from highest point in " &
          //"the colormap's range."
        write (scriptunit, '(a)') '    Defaults to 1.0 (no upper offset). ' &
          //'Should be between'
        write (scriptunit, '(a)') '    `midpoint` and 1.0.'
        write (scriptunit, '(a)') "  '''"
        write (scriptunit, '(a)') "  cdict = {'red': [],'green': []," &
          //"'blue': [],'alpha': []}"
        write (scriptunit, '(a)') '  # regular index to compute the colors'
        write (scriptunit, '(a)') '  reg_index = np.linspace(start, stop, 257)'
        write (scriptunit, '(a)') '  # shifted index to match the data'
        write (scriptunit, '(a)') '  shift_index = np.hstack(['
        write (scriptunit, '(a)') '    np.linspace(0.0, midpoint, 128, ' &
          //'endpoint=False),'
        write (scriptunit, '(a)') '    np.linspace(midpoint, 1.0, 129, ' &
          //'endpoint=True)'
        write (scriptunit, '(a)') '  ])'
        write (scriptunit, '(a)') '  for ri, si in zip(reg_index, shift_index):'
        write (scriptunit, '(a)') '    r, g, b, a = cmap(ri)'
        write (scriptunit, '(a)') "    cdict['red'].append((si, r, r))"
        write (scriptunit, '(a)') "    cdict['green'].append((si, g, g))"
        write (scriptunit, '(a)') "    cdict['blue'].append((si, b, b))"
        write (scriptunit, '(a)') "    cdict['alpha'].append((si, a, a))"
        write (scriptunit, '(a)') '  newcmap = matplotlib.colors' &
          //'.LinearSegmentedColormap(name, cdict)'
        write (scriptunit, '(a)') '  plt.register_cmap(cmap=newcmap)'
        write (scriptunit, '(a)') '  return newcmap'
        write (scriptunit, '(a)') " "
        write (scriptunit, '(a)') "outfile = '"//trim(seedname)//"-kslice-shc.pdf'"
        write (scriptunit, '(a)') " "
        write (scriptunit, '(a)') "val = np.loadtxt('"//trim(seedname) &
          //"-kslice-shc.dat', usecols=(0,))"
        write (scriptunit, '(a)') " "
        write (scriptunit, '(a)') "val_log=np.array([np.log10(abs(elem))*np.sign(elem)" &
          //"if abs(elem)>10 else elem/10.0 for elem in val])"
        write (scriptunit, '(a)') "#val_log = val"
        write (scriptunit, '(a)') "valmax=max(val_log)"
        write (scriptunit, '(a)') "valmin=min(val_log)"
        write (scriptunit, '(a)') "#cmnew=shiftedColorMap(matplotlib.cm.bwr," &
          //"0,1-valmax/(valmax+abs(valmin)),1)"
        write (scriptunit, '(a)') " "
        write (scriptunit, '(a)') "if square: "
        write (scriptunit, '(a)') "  Z=val_log.reshape(dimy,dimx)"
        write (scriptunit, '(a)') "  mn=int(np.floor(Z.min()))"
        write (scriptunit, '(a)') "  mx=int(np.ceil(Z.max()))"
        write (scriptunit, '(a)') "  ticks=range(mn,mx+1)"
        write (scriptunit, '(a)') "  #pl.contourf(x_coord,y_coord,Z," &
          //"ticks,origin='lower')"
        write (scriptunit, '(a)') "  pl.imshow(Z,origin='lower'," &
          //"extent=(min(x_coord),max(x_coord),min(y_coord)," &
          //"max(y_coord)))#,cmap=cmnew)"
        write (scriptunit, '(a)') "else: "
        write (scriptunit, '(a)') "  grid_x, grid_y = np.meshgrid(xint,yint)"
        write (scriptunit, '(a)') "  valint = interpolate.griddata((points_x," &
          //"points_y), val_log, (grid_x,grid_y), method='nearest')"
        write (scriptunit, '(a)') "  mn=int(np.floor(valint.min()))"
        write (scriptunit, '(a)') "  mx=int(np.ceil(valint.max()))"
        write (scriptunit, '(a)') "  ticks=range(mn,mx+1)"
        write (scriptunit, '(a)') "  #pl.contourf(xint,yint,valint,ticks)"
        write (scriptunit, '(a)') "  pl.imshow(valint,origin='lower'," &
          //"extent=(min(xint),max(xint),min(yint),max(yint)))#,cmap=cmnew)"
        write (scriptunit, '(a)') " "
        write (scriptunit, '(a)') "ticklabels=[]"
        write (scriptunit, '(a)') "for n in ticks:"
        write (scriptunit, '(a)') " if n<0: "
        write (scriptunit, '(a)') "  ticklabels.append('-$10^{%d}$' % abs(n))"
        write (scriptunit, '(a)') " elif n==0:"
        write (scriptunit, '(a)') "  ticklabels.append(' $%d$' %  n)"
        write (scriptunit, '(a)') " else:"
        write (scriptunit, '(a)') "  ticklabels.append(' $10^{%d}$' % n)"
        write (scriptunit, '(a)') " "
        write (scriptunit, '(a)') "cbar=pl.colorbar()"
        write (scriptunit, '(a)') "#cbar.set_ticks(ticks)"
        write (scriptunit, '(a)') "#cbar.set_ticklabels(ticklabels)"
        write (scriptunit, '(a)') " "
        write (scriptunit, '(a)') "ax = pl.gca()"
        write (scriptunit, '(a)') "ax.xaxis.set_visible(False)"
        write (scriptunit, '(a)') "ax.yaxis.set_visible(False)"
        write (scriptunit, '(a)') " "
        write (scriptunit, '(a)') "pl.savefig(outfile,bbox_inches='tight')"
        write (scriptunit, '(a)') "pl.show()"
      end if

      write (stdout, *) ' '

    end if ! on_root

  end subroutine k_slice

  !================================================!
  !                   PRIVATE PROCEDURES
  !================================================!

  subroutine kslice_print_info(plot_fermi_lines, fermi_lines_color, plot_curv, plot_morb, &
                               plot_shc, stdout, pw90_berry, fermi_energy_list, error, comm)
    !================================================!

    use w90_constants, only: dp
    use w90_postw90_types, only: pw90_berry_mod_type
    use w90_comms, only: w90_comm_type

    type(pw90_berry_mod_type), intent(in) :: pw90_berry
    real(kind=dp), allocatable, intent(in) :: fermi_energy_list(:)
    type(w90_error_type), allocatable, intent(out) :: error
    type(w90_comm_type), intent(in) :: comm
    integer, intent(in) :: stdout
    logical, intent(in) :: plot_fermi_lines, fermi_lines_color, plot_curv, plot_morb, plot_shc

    integer :: fermi_n

    write (stdout, '(/,/,1x,a)') &
      'Properties calculated in module  k s l i c e'
    write (stdout, '(1x,a)') &
      '--------------------------------------------'

    fermi_n = 0
    if (allocated(fermi_energy_list)) fermi_n = size(fermi_energy_list)
    if (plot_fermi_lines) then
      if (fermi_n /= 1) then
        call set_error_input(error, 'Must specify one Fermi level when kslice_task=fermi_lines', &
                             comm)
        return
      endif
      select case (fermi_lines_color)
      case (.false.)
        write (stdout, '(/,3x,a)') '* Fermi lines'
      case (.true.)
        write (stdout, '(/,3x,a)') '* Fermi lines coloured by spin'
      end select
      write (stdout, '(/,7x,a,f10.4,1x,a)') &
        '(Fermi level: ', fermi_energy_list(1), 'eV)'
    endif

    if (plot_curv) then
      if (pw90_berry%curv_unit == 'ang2') then
        write (stdout, '(/,3x,a)') '* Negative Berry curvature in Ang^2'
      elseif (pw90_berry%curv_unit == 'bohr2') then
        write (stdout, '(/,3x,a)') '* Negative Berry curvature in Bohr^2'
      endif
      if (fermi_n /= 1) then
        call set_error_input(error, 'Must specify one Fermi level when kslice_task=curv', comm)
        return
      endif
    elseif (plot_morb) then
      write (stdout, '(/,3x,a)') '* Orbital magnetization k-space integrand in eV.Ang^2'
      if (fermi_n /= 1) then
        call set_error_input(error, 'Must specify one Fermi level when kslice_task=morb', comm)
        return
      endif
    elseif (plot_shc) then
      if (pw90_berry%curv_unit == 'ang2') then
        write (stdout, '(/,3x,a)') '* Berry curvature-like term ' &
          //'of spin Hall conductivity in Ang^2'
      elseif (pw90_berry%curv_unit == 'bohr2') then
        write (stdout, '(/,3x,a)') '* Berry curvature-like term ' &
          //'of spin Hall conductivity in Bohr^2'
      endif
      if (fermi_n /= 1) then
        call set_error_input(error, 'Must specify one Fermi level when kslice_task=shc', comm)
        return
      endif
    endif

  end subroutine kslice_print_info

  subroutine write_data_file(stdout, filename, fmt, data)
    use w90_constants, only: dp

    integer, intent(in) :: stdout
    character(len=*), intent(in) :: filename, fmt
    real(kind=dp), intent(in) :: data(:, :)

    integer :: n, i, fileunit

    write (stdout, '(/,3x,a)') filename
    open (newunit=fileunit, file=filename, form='formatted')

    n = size(data, 2)
    do i = 1, n
      write (fileunit, fmt) data(:, i)
    end do

    write (fileunit, *) ''
    close (fileunit)
  end subroutine

  subroutine write_coords_file(stdout, filename, fmt, coords, vals, mask, blocklen)
    !================================================!

    use w90_constants, only: dp

    integer, intent(in) :: stdout
    character(len=*), intent(in) :: filename, fmt
    real(kind=dp), intent(in)  :: coords(:, :), vals(:, :, :)
    logical, intent(in), optional :: mask(:, :)
    integer, intent(in), optional :: blocklen

    integer :: n, m, i, j, fileunit, bl

    write (stdout, '(/,3x,a)') filename
    open (newunit=fileunit, file=filename, form='formatted')

    n = size(vals, 3)
    m = size(vals, 2)

    if (present(mask)) then
      do i = 1, n
        do j = 1, m
          if (mask(j, i)) then
            write (fileunit, fmt) coords(:, i), vals(:, j, i)
          end if
        end do
      end do
    else
      if (present(blocklen)) then
        bl = blocklen
      else
        bl = n
      end if

      do i = 1, n
        do j = 1, m
          write (fileunit, fmt) coords(:, i), vals(:, j, i)
        end do
        if (mod(i, bl) == 0) then
          write (fileunit, *) ''
        end if
      end do
    end if
    write (fileunit, *) ''
    close (fileunit)
  end subroutine

  subroutine script_common(scriptunit, areab1b2, square, seedname)
    !================================================!

    use w90_constants, only: dp

    integer, intent(in) :: scriptunit
    real(kind=dp), intent(in) :: areab1b2
    character(len=25), intent(in) :: square
    character(len=50), intent(in) :: seedname

    write (scriptunit, '(a)') "import pylab as pl"
    write (scriptunit, '(a)') "import numpy as np"
    write (scriptunit, '(a)') "from scipy import interpolate"
    write (scriptunit, '(a)') "from collections import OrderedDict"
    write (scriptunit, '(a)') " "
    write (scriptunit, '(a)') "points = np.loadtxt('"//trim(seedname)// &
      "-kslice-coord.dat')"
    write (scriptunit, '(a)') "# Avoid numerical noise"
    write (scriptunit, '(a)') "points_x=np.around(points[:,0],decimals=10)"
    write (scriptunit, '(a)') "points_y=np.around(points[:,1],decimals=10)"
    write (scriptunit, '(a)') "num_pt=len(points)"
    write (scriptunit, '(a)') " "
    write (scriptunit, '(a,f12.6)') "area=", areab1b2
    write (scriptunit, '(a)') " "
    write (scriptunit, '(a)') "square= "//square
    write (scriptunit, '(a)') " "
    write (scriptunit, '(a)') "if square:"
    write (scriptunit, '(a)') &
      "  x_coord=list(OrderedDict.fromkeys(points_x))"
    write (scriptunit, '(a)') &
      "  y_coord=list(OrderedDict.fromkeys(points_y))"
    write (scriptunit, '(a)') "  dimx=len(x_coord)"
    write (scriptunit, '(a)') "  dimy=len(y_coord)"
    write (scriptunit, '(a)') "else:"
    write (scriptunit, '(a)') "  xmin=np.min(points_x)"
    write (scriptunit, '(a)') "  ymin=np.min(points_y)"
    write (scriptunit, '(a)') "  xmax=np.max(points_x)"
    write (scriptunit, '(a)') "  ymax=np.max(points_y)"
    write (scriptunit, '(a)') &
      "  a=np.max(np.array([xmax-xmin,ymax-ymin]))"
    write (scriptunit, '(a)') &
      "  num_int=int(round(np.sqrt(num_pt*a**2/area)))"
    write (scriptunit, '(a)') "  xint = np.linspace(xmin,xmin+a,num_int)"
    write (scriptunit, '(a)') "  yint = np.linspace(ymin,ymin+a,num_int)"
    write (scriptunit, '(a)') " "

  end subroutine script_common

  subroutine script_fermi_lines(scriptunit, seedname, fermi_energy_list)
    !================================================!

    use w90_constants, only: dp
    implicit none
    real(kind=dp), allocatable, intent(in) :: fermi_energy_list(:)
    integer, intent(in) :: scriptunit
    character(len=50), intent(in)  :: seedname

    write (scriptunit, '(a)') &
      "# Energy level for isocontours (typically the Fermi level)"
    write (scriptunit, '(a,f12.6)') "ef=", fermi_energy_list(1)
    write (scriptunit, '(a)') " "
    write (scriptunit, '(a)') &
      "bands=np.loadtxt('"//trim(seedname)//"-kslice-bands.dat')"
    write (scriptunit, '(a)') "numbands=bands.size//num_pt"
    write (scriptunit, '(a)') "if square:"
    write (scriptunit, '(a)') &
      "  bbands=bands.reshape((dimy,dimx,numbands))"
    write (scriptunit, '(a)') "  for i in range(numbands):"
    write (scriptunit, '(a)') "    Z=bbands[:,:,i]"
    write (scriptunit, '(a)') "    pl.contour(x_coord,y_coord,Z," &
      //"[ef],colors='black')"
    write (scriptunit, '(a)') "else:"
    write (scriptunit, '(a)') "  bbands=bands.reshape((num_pt," &
      //"numbands))"
    write (scriptunit, '(a)') "  bandint=[]"
    write (scriptunit, '(a)') "  grid_x, grid_y = np.meshgrid(xint,yint)"
    write (scriptunit, '(a)') "  for i in range(numbands):"
    write (scriptunit, '(a)') "    bandint.append(interpolate.griddata" &
      //"((points_x,points_y), bbands[:,i], (grid_x,grid_y), " &
      //"method='nearest'))"
    write (scriptunit, '(a)') "    pl.contour(grid_x,grid_y," &
      //"bandint[i],[ef],colors='black')"

  end subroutine script_fermi_lines

end module w90_kslice
