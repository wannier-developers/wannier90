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
!  w90_tetrahedron: functions for tetrahedron method         !
!                                                            !
!------------------------------------------------------------!

module w90_tetrahedron

  use w90_constants, only: dp

  implicit none

  private
  !tetrahedron_fermidirac
  !tetrahedron_sort
  !tetrahedron_jacobian
  !tetrahedron_log1p

  public :: tetrahedron_spinhall
  ! To do: extend to other quantities(anomalous, ...)
  public :: tetrahedron_integral
  public :: tetrahedron_P_matrix_init
  public :: tetrahedron_array_init
  public :: tetrahedron_array_small_init

contains

  subroutine tetrahedron_P_matrix_init(P_matrix)

    use w90_constants, only: dp

    implicit none

    real(kind=dp), dimension(4, 20), intent(inout) :: P_matrix

    ! correction, Kawamura PRB 89 094515
    P_matrix(1, 1:4) = real((/1440, 0, 30, 0/), dp)
    P_matrix(2, 1:4) = real((/0, 1440, 0, 30/), dp)
    P_matrix(3, 1:4) = real((/30, 0, 1440, 0/), dp)
    P_matrix(4, 1:4) = real((/0, 30, 0, 1440/), dp)
    !
    P_matrix(1, 5:8) = real((/-38, 7, 17, -28/), dp)
    P_matrix(2, 5:8) = real((/-28, -38, 7, 17/), dp)
    P_matrix(3, 5:8) = real((/17, -28, -38, 7/), dp)
    P_matrix(4, 5:8) = real((/7, 17, -28, -38/), dp)
    !
    P_matrix(1, 9:12) = real((/-56, 9, -46, 9/), dp)
    P_matrix(2, 9:12) = real((/9, -56, 9, -46/), dp)
    P_matrix(3, 9:12) = real((/-46, 9, -56, 9/), dp)
    P_matrix(4, 9:12) = real((/9, -46, 9, -56/), dp)
    !
    P_matrix(1, 13:16) = real((/-38, -28, 17, 7/), dp)
    P_matrix(2, 13:16) = real((/7, -38, -28, 17/), dp)
    P_matrix(3, 13:16) = real((/17, 7, -38, -28/), dp)
    P_matrix(4, 13:16) = real((/-28, 17, 7, -38/), dp)
    !
    P_matrix(1, 17:20) = real((/-18, -18, 12, -18/), dp)
    P_matrix(2, 17:20) = real((/-18, -18, -18, 12/), dp)
    P_matrix(3, 17:20) = real((/12, -18, -18, -18/), dp)
    P_matrix(4, 17:20) = real((/-18, 12, -18, -18/), dp)
    !
    P_matrix(1:4, 1:20) = P_matrix(1:4, 1:20)/1260.0_dp

  end subroutine tetrahedron_P_matrix_init

  subroutine tetrahedron_array_small_init(tet_array_small)

    implicit none

    integer, dimension(6, 4), intent(inout) :: tet_array_small

    ! without correction
    !                                               z                     !
    !                                                 7----------8        ! tetrahedron #1: 1-5-6-8
    tet_array_small(1, :) = (/1, 5, 6, 8/)  !        /|         /| layer 2! tetrahedron #2: 1-7-5-8
    tet_array_small(2, :) = (/1, 7, 5, 8/)  !       5----------6 |        ! tetrahedron #3: 1-3-7-8
    tet_array_small(3, :) = (/1, 3, 7, 8/)  !       | | y      | |        ! tetrahedron #4: 1-4-3-8
    tet_array_small(4, :) = (/1, 4, 3, 8/)  !       | 3--------|-4        ! tetrahedron #5: 1-2-4-8
    tet_array_small(5, :) = (/1, 2, 4, 8/)  !       |/         |/  layer 1! tetrahedron #6: 1-6-2-8
    tet_array_small(6, :) = (/1, 6, 2, 8/)  !       1----------2 x        !
    ! 1-8 : main diagonal         !
    !=============================!

  end subroutine tetrahedron_array_small_init

  subroutine tetrahedron_array_init(tet_array)

    implicit none

    integer, dimension(6, 20), intent(inout) :: tet_array
    ! with correction, Kawamura PRB 89 094515
    !                                               z                     !
    !                                       !         61--62--63-64       ! tetrahedron #1: 22-38-39-43
    tet_array(1, 1:4) = (/22, 38, 39, 43/)  !        /|         /| layer 4! tetrahedron #2: 22-42-38-43
    tet_array(2, 1:4) = (/22, 42, 38, 43/)  !       49--50--51-52|        ! tetrahedron #3: 22-26-42-43
    tet_array(3, 1:4) = (/22, 26, 42, 43/)  !       | | y      | |        ! tetrahedron #4: 22-27-26-43
    tet_array(4, 1:4) = (/22, 27, 26, 43/)  !       | 13--14--1|-16       ! tetrahedron #5: 22-23-27-43
    tet_array(5, 1:4) = (/22, 23, 27, 43/)  !       |/         |/  layer 1! tetrahedron #6: 22-39-23-43
    tet_array(6, 1:4) = (/22, 39, 23, 43/)  !       1---2---3--4 x        !
    ! 22-43 : main diagonal       !
    !=============================!
    tet_array(1, 5:20) = (/6, 37, 35, 64, 5, 33, 56, 48, 1, 54, 40, 47, 59, 23, 42, 18/)
    tet_array(2, 5:20) = (/2, 46, 33, 64, 6, 41, 54, 44, 1, 62, 34, 48, 63, 18, 47, 17/)
    tet_array(3, 5:20) = (/18, 10, 41, 64, 2, 9, 62, 60, 1, 30, 58, 44, 47, 38, 27, 21/)
    tet_array(4, 5:20) = (/17, 28, 9, 64, 18, 11, 30, 59, 1, 32, 25, 60, 48, 21, 44, 5/)
    tet_array(5, 5:20) = (/21, 19, 11, 64, 17, 3, 32, 63, 1, 24, 31, 59, 44, 26, 39, 6/)
    tet_array(6, 5:20) = (/5, 55, 3, 64, 21, 35, 24, 47, 1, 56, 7, 63, 60, 6, 59, 2/)
    !note: 18 points(4, 8, 12, 13, 14, 15, 16, 20, 29, 36, 45, 49, 50, 51, 52, 53, 57, 61) are not used

  end subroutine tetrahedron_array_init

  !================================================!
  function tetrahedron_spinhall(F, E1, E2, t, hw, Ef, type, tet_cutoff)
    !=============================================================================
    ! Calculates contribution from a single tetrahedron, for the Kubo formula.   !
    ! "integral d3k (f_{nk} - f_{mk} * (...)
    ! Ref: Minsu Ghim and Cheol-Hwan Park, PRB ()                                !
    !=============================================================================

    use w90_constants, only: dp

    implicit none

    real(kind=dp) :: tetrahedron_spinhall
    real(kind=dp), dimension(4), intent(in) :: F, E1, E2
    real(kind=dp), dimension(3, 3), intent(in) :: t
    real(kind=dp), intent(in) :: hw, Ef
    integer, intent(in) :: type
    real(kind=dp), intent(in) :: tet_cutoff
    !intermediate variables, s:sorted
    real(kind=dp) :: t_s(3, 3), t_small(3, 3), x(3), y
    real(kind=dp), dimension(4) :: D, occ1, occ2, F_s, E1_s, E2_s, &
                                   F_small, D_small
    integer       :: i
    logical       :: flag1, flag2

    !fnk-fmk = 0 then quickly returns zero
    occ1 = 0.0_dp; occ2 = 0.0_dp; flag1 = .true.; flag2 = .true.
    do i = 1, 4
      if (E1(i) < Ef) occ1(i) = 1.0_dp
      if (E2(i) < Ef) occ2(i) = 1.0_dp
      if (occ1(i) /= 1.0 .or. occ2(i) /= 1.0) flag1 = flag1 .and. .false.
      if (occ1(i) /= 0.0 .or. occ2(i) /= 0.0) flag2 = flag2 .and. .false.
    enddo
    if (flag1 .or. flag2) then
      tetrahedron_spinhall = 0.0_dp
    else
      tetrahedron_spinhall = tetrahedron_fermidirac(F, E1, E2, t, hw, Ef, type, tet_cutoff) &
                             - tetrahedron_fermidirac(F, E2, E1, t, hw, Ef, type, tet_cutoff)
    endif

    return

  end function tetrahedron_spinhall

  function tetrahedron_fermidirac(F, E_ref, E2, t, hw, Ef, type, tet_cutoff)
    !=============================================================================
    ! Calculates contribution from a single tetrahedron, for the Kubo formula.   !
    ! "integral d3k f_{nk} * (...)
    ! Ref: Minsu Ghim and Cheol-Hwan Park, PRB ()                                !
    !=============================================================================

    use w90_constants, only: dp

    implicit none

    real(kind=dp) :: tetrahedron_fermidirac
    real(kind=dp), dimension(4), intent(in) :: F, E_ref, E2
    real(kind=dp), dimension(3, 3), intent(in) :: t
    real(kind=dp), intent(in) :: hw, Ef
    integer, intent(in) :: type
    real(kind=dp), intent(in) :: tet_cutoff
    !intermediate variables, s:sorted
    real(kind=dp) :: t_s(3, 3), t_small(3, 3), x(3), y
    real(kind=dp), dimension(4) :: D, F_s, E1_s, E2_s, F_small, D_small
    !result
    real(kind=dp) :: Ans
    integer :: i

    Ans = 0.0_dp
    !sorting vertices according to E_ref
    F_s = F; E1_s = E_ref; E2_s = E2; t_s = t !_s: sorted
    call tetrahedron_sort(E1_s, E2_s, F_s, t_s)
    D = E1_s - E2_s
    !case 1,2,3,4,5
    if (Ef < E1_s(1)) then       ! case 1: zero
      Ans = Ans + 0.0_dp
    else if (Ef < E1_s(2)) then  ! case 2: a small tet.

      x(1) = (Ef - E1_s(1))/(E1_s(2) - E1_s(1))
      x(2) = (Ef - E1_s(1))/(E1_s(3) - E1_s(1))
      x(3) = (Ef - E1_s(1))/(E1_s(4) - E1_s(1))

      F_small(1) = F_s(1); D_small(1) = D(1)
      do i = 1, 3
        F_small(i + 1) = F_s(1) + (F_s(i + 1) - F_s(1))*x(i)
        D_small(i + 1) = D(1) + (D(i + 1) - D(1))*x(i)
        t_small(:, i) = t_s(:, i)*x(i)
      enddo
      Ans = Ans + tetrahedron_integral(F_small, D_small, t_small, hw, type, tet_cutoff)

    else if (Ef < E1_s(3)) then  ! case 3: two tet.'s with cases 2 and 4

      x(1) = (Ef - E1_s(4))/(E1_s(2) - E1_s(4))
      x(2) = (Ef - E1_s(1))/(E1_s(3) - E1_s(1))
      x(3) = (Ef - E1_s(1))/(E1_s(4) - E1_s(1))
      y = (Ef - E1_s(3))/(E1_s(2) - E1_s(3))

      F_small = F_s; D_small = D; t_small = t_s
      F_small(4) = F_s(1) + (F_s(4) - F_s(1))*x(3)
      D_small(4) = D(1) + (D(4) - D(1))*x(3)
      t_small(:, 3) = t_s(:, 3)*x(3)
      Ans = Ans + tetrahedron_integral(F_small, D_small, t_small, hw, type, tet_cutoff)

      F_small(1) = F_s(1) + (F_s(3) - F_s(1))*x(2); F_small(2) = F_s(3) + (F_s(2) - F_s(3))*y
      D_small(1) = D(1) + (D(3) - D(1))*x(2); D_small(2) = D(3) + (D(2) - D(3))*y
      t_small(:, 1) = t_s(:, 1)*y + t_s(:, 2)*(1 - y - x(2))
      t_small(:, 2) = t_s(:, 2)*(1 - x(2))
      t_small(:, 3) = t_s(:, 3)*x(3) - t_s(:, 2)*x(2)
      Ans = Ans - tetrahedron_integral(F_small, D_small, t_small, hw, type, tet_cutoff)

      F_small(1) = F_s(4) + (F_s(2) - F_s(4))*x(1)
      F_small(3) = F_small(2); F_small(2) = F_s(2)
      D_small(1) = D(4) + (D(2) - D(4))*x(1)
      D_small(3) = D_small(2); D_small(2) = D(2)
      t_small(:, 1) = (t_s(:, 1) - t_s(:, 3))*(1 - x(1))
      t_small(:, 2) = t_s(:, 1)*(y - x(1)) + t_s(:, 2)*(1 - y) + t_s(:, 3)*(x(1) - 1)
      t_small(:, 3) = -t_s(:, 1)*x(1) + t_s(:, 3)*(x(1) + x(3) - 1)
      Ans = Ans + tetrahedron_integral(F_small, D_small, t_small, hw, type, tet_cutoff)

    else if (Ef < E1_s(4)) then  ! case 4: a large tet. - a small tet.

      x(1) = (Ef - E1_s(4))/(E1_s(2) - E1_s(4))
      x(2) = (Ef - E1_s(4))/(E1_s(3) - E1_s(4))
      x(3) = (Ef - E1_s(1))/(E1_s(4) - E1_s(1))

      F_small(1) = F_s(4); F_small(2) = F_s(1) + (F_s(4) - F_s(1))*x(3)
      F_small(3) = F_s(4) + (F_s(2) - F_s(4))*x(1)
      F_small(4) = F_s(4) + (F_s(3) - F_s(4))*x(2)
      D_small(1) = D(4); D_small(2) = D(1) + (D(4) - D(1))*x(3)
      D_small(3) = D(4) + (D(2) - D(4))*x(1)
      D_small(4) = D(4) + (D(3) - D(4))*x(2)
      t_small(:, 1) = -t_s(:, 3)*(1 - x(3))
      t_small(:, 2) = (t_s(:, 1) - t_s(:, 3))*x(1)
      t_small(:, 3) = (t_s(:, 2) - t_s(:, 3))*x(2)
      Ans = Ans + tetrahedron_integral(F_s, D, t_s, hw, type, tet_cutoff) &
            - tetrahedron_integral(F_small, D_small, t_small, hw, type, tet_cutoff)

    else                      ! case 5: a large tet.
      Ans = Ans + tetrahedron_integral(F_s, D, t_s, hw, type, tet_cutoff)
    endif
    tetrahedron_fermidirac = Ans

  end function tetrahedron_fermidirac

  !=======================================================================
  function tetrahedron_integral(F_in, D_in, t_in, hw, type, tet_cutoff)
    !=============================================================================
    ! Calculates contribution from a single tetrahedron, for the Kubo formula.   !
    ! To do: extend to complex numbers                                           !
    ! type 1: nondissipative, linear/linear                                      !
    ! type 2: dissipative                                                        !
    ! type 3: nondissipative, linear/(linear)**2                                 !
    !                                                                            !
    ! Ref: Minsu Ghim and Cheol-Hwan Park, PRB ()                                !
    !=============================================================================

    use w90_constants, only: dp
    use w90_utility, only: utility_inv3

    implicit none

    real(kind=dp) :: tetrahedron_integral
    real(kind=dp), dimension(4), intent(in) :: F_in, D_in
    real(kind=dp), dimension(3, 3), intent(in) :: t_in
    real(kind=dp), intent(in) :: hw
    integer, intent(in) :: type
    real(kind=dp), intent(in) :: tet_cutoff
    !result
    real(kind=dp) :: Ans, Det_t

    !intermediate variables
    integer :: i, j, k, l
    real(kind=dp) :: D(4), F(4), t(3, 3), DAV, D_small_prev, D_large_prev
    real(kind=dp) :: dd(3), ll(3), ff, bb(4), cc(4, 3)
    integer :: a, b, c
    REAL(kind=dp) :: t_inverse(3, 3), dummy(4)
    REAL(kind=dp) :: GradD, Jac, y
    REAL(kind=dp) :: x(3)
    REAL(kind=dp), DIMENSION(0:2) :: F_uv
    !small parameter to avoid a problem of degenearcy
    real(kind=dp), parameter :: avoid_deg = 1.e-4_dp

    D = D_in; F = F_in; t = t_in
    call tetrahedron_sort(D, F, dummy, t)
    Ans = 0.0_dp

    !case 1 and 3: nondissipative part, case 2: dissipative part
    if (type == 1 .or. type == 3) then
      if (type == 3) then
        !treatment for accidental small band spliting(but degenerate actually) cases
        do j = 1, 4
          if (abs(D(j)) < avoid_deg) then
            D(j) = avoid_deg*(abs(D(j))/D(j))
            F = 0
          endif
        enddo
      endif

      !cutoff treatment, hw == 0.0 for case 3
      DAV = (D(2) + D(3))/2.0_dp
      if (abs((D(2) - D(3))/(DAV + hw)) < tet_cutoff) then
        D_small_prev = D(2); D_large_prev = D(3)
        D(3) = DAV + 0.5_dp*abs(DAV + hw)*tet_cutoff
        D(2) = DAV - 0.5_dp*abs(DAV + hw)*tet_cutoff
        if (D(1) > D(2)) D(1) = D(1) + (D(2) - D_small_prev)
        if (D(3) > D(4)) D(4) = D(4) + (D(3) - D_large_prev)
      endif
      DAV = (D(1) + D(2))/2.0_dp
      if (abs((D(1) - D(2))/(DAV + hw)) < tet_cutoff) then
        if (D(2) > 0) then
          D(1) = D(2)*(2.0_dp - tet_cutoff)/(2.0_dp + tet_cutoff)
        else
          D(1) = D(2)*(2.0_dp + tet_cutoff)/(2.0_dp - tet_cutoff)
        endif
      endif
      DAV = (D(3) + D(4))/2.0_dp
      if (abs((D(3) - D(4))/(DAV + hw)) < tet_cutoff) then
        if (D(3) > 0) then
          D(4) = D(3)*(2.0_dp + tet_cutoff)/(2.0_dp - tet_cutoff)
        else
          D(4) = D(3)*(2.0_dp - tet_cutoff)/(2.0_dp + tet_cutoff)
        endif
      endif

      !intermediate variables for case 1 and 3
      do i = 1, 3
        dd(i) = (D(4) - D(i))/(D(i) + hw)
        ll(i) = tetrahedron_log1p(dd(i))
      enddo
      ff = 1.0_dp

      !determinant factor from parametrisation(tetrahedron volume)
      Det_t = ABS(t(1, 1)*t(2, 2)*t(3, 3) + t(1, 2)*t(2, 3)*t(3, 1) + t(1, 3)*t(2, 1)*t(3, 2) &
                  - t(1, 1)*t(2, 3)*t(3, 2) - t(1, 3)*t(2, 2)*t(3, 1) - t(1, 2)*t(2, 1)*t(3, 3))
    endif

    select case (type)
    case (1)
      do i = 1, 3
        a = i
        b = mod(i, 3) + 1
        c = mod(i + 1, 3) + 1
        cc(a, a) = -(1.0_dp + dd(a))*(3.0_dp*dd(a)**2 - 2.0_dp*(dd(b) + dd(c))*dd(a) + dd(b)*dd(c)) &
                   *((dd(b) - dd(c))*dd(b)*dd(c))**2
        cc(b, a) = -dd(a)*(1.0_dp + dd(b))*(dd(c) - dd(a))*((dd(b) - dd(c))*dd(b)*dd(c))**2
        cc(c, a) = dd(a)*(1.0_dp + dd(c))*(dd(a) - dd(b))*((dd(b) - dd(c))*dd(b)*dd(c))**2
        cc(4, a) = -(dd(a) - dd(b))*(dd(c) - dd(a))*((dd(b) - dd(c))*dd(b)*dd(c))**2
        bb(a) = cc(4, a)*dd(a)
        ff = ff*(1.0_dp + dd(a))/(dd(a)*(dd(a) - dd(b)))**2
      enddo
      bb(4) = -dd(1)*dd(2)*dd(3)*((dd(1) - dd(2))*(dd(2) - dd(3))*(dd(3) - dd(1)))**2
      ff = -ff/6.0_dp

      do i = 1, 4
        Ans = Ans + F(i)*(cc(i, 1)*ll(1) + cc(i, 2)*ll(2) + cc(i, 3)*ll(3) + bb(i))
      enddo
      Ans = Ans*ff/(D(4) + hw)

      tetrahedron_integral = Ans*Det_t

    case (2)
      ! calculate integ d2k / |grad D| * F(k)

      ! obtaining |grad D|
      call utility_inv3(t, t_inverse, Det_t)
      t_inverse = t_inverse/Det_t

      GradD = 0.0_dp
      DO i = 1, 3
        DO j = 1, 3
          DO k = 1, 3
            GradD = GradD + t_inverse(i, k)*t_inverse(j, k)*(D(i + 1) - D(1))*(D(j + 1) - D(1))
          END DO
        END DO
      END DO
      GradD = SQRT(ABS(GradD))

      IF (hw < D(1)) THEN
        Ans = 0.0_dp
      ELSE IF (hw < D(2)) THEN
        ! parametrization
        x(1) = (hw - D(1))/(D(2) - D(1))
        x(2) = (hw - D(1))/(D(3) - D(1))
        x(3) = (hw - D(1))/(D(4) - D(1))
        ! Jacobian factor
        Jac = tetrahedron_jacobian(t, x, 1)
        ! integration formula
        F_uv(0) = F(1) + (F(2) - F(1))*x(1)
        F_uv(1) = (F(3) - F(1))*x(2) - (F(2) - F(1))*x(1)
        F_uv(2) = (F(4) - F(1))*x(3) - (F(2) - F(1))*x(1)
        Ans = Jac*(F_uv(0)/2.0_dp + (F_uv(1) + F_uv(2))/6.0_dp)/GradD
      ELSE IF (hw < D(3)) THEN
        ! parametrization
        x(1) = (hw - D(4))/(D(2) - D(4))
        x(2) = (hw - D(1))/(D(3) - D(1))
        x(3) = (hw - D(1))/(D(4) - D(1))

        !! triangle 1
        ! Jacobian factor
        Jac = tetrahedron_jacobian(t, x, 2)
        ! integration formula
        F_uv(0) = F(1) + (F(4) - F(1))*x(3)
        F_uv(1) = (F(2) - F(1))*x(1) + (F(4) - F(1))*(1.0_dp - x(1) - x(3))
        F_uv(2) = (F(3) - F(1))*x(2) - (F(4) - F(1))*x(3)
        Ans = Jac*(F_uv(0)/2.0_dp + (F_uv(1) + F_uv(2))/6.0_dp)

        !! triangle 2
        y = (hw - D(3))/(D(2) - D(3))
        ! Jacobian factor
        Jac = tetrahedron_jacobian(t, x, 3)
        ! integration formula
        F_uv(0) = F(1) + (F(2) - F(1))*y + (F(3) - F(1))*(1.0_dp - y)
        F_uv(1) = (F(2) - F(1))*(x(1) - y) + &
                  (F(3) - F(1))*(y - 1.0_dp) + (F(4) - F(1))*(1.0_dp - x(1))
        F_uv(2) = -(F(2) - F(1))*y + (F(3) - F(1))*(y - 1.0_dp + x(2))
        Ans = Ans + Jac*(F_uv(0)/2.0_dp + (F_uv(1) + F_uv(2))/6.0_dp)

        Ans = Ans/GradD
      ELSE IF (hw < D(4)) THEN
        ! parametrization
        x(1) = (hw - D(4))/(D(2) - D(4))
        x(2) = (hw - D(4))/(D(3) - D(4))
        x(3) = (hw - D(1))/(D(4) - D(1))
        ! Jacobian factor
        Jac = tetrahedron_jacobian(t, x, 4)
        ! integration formula
        F_uv(0) = F(1) + (F(4) - F(1))*x(3)
        F_uv(1) = (F(2) - F(1))*x(1) + (F(4) - F(1))*(1.0_dp - x(1) - x(3))
        F_uv(2) = (F(3) - F(1))*x(2) + (F(4) - F(1))*(1.0_dp - x(2) - x(3))
        Ans = Jac*(F_uv(0)/2.0_dp + (F_uv(1) + F_uv(2))/6.0_dp)/GradD
      ELSE
        Ans = 0.0_dp
      END IF
      tetrahedron_integral = Ans

    case (3)
      do i = 1, 3
        a = i
        b = mod(i, 3) + 1
        c = mod(i + 1, 3) + 1
        cc(a, a) = -(1.0_dp + dd(a))*(2.0_dp*dd(a)**3 + (3.0_dp - dd(b) - dd(c))*dd(a)**2 &
                                      - 2.0_dp*(dd(b) + dd(c))*dd(a) + dd(b)*dd(c))*((dd(b) - dd(c))*dd(b)*dd(c))**2
        cc(b, a) = -(1.0_dp + dd(a))*dd(a)*(1.0_dp + dd(b))*(dd(c) - dd(a))*((dd(b) - dd(c))*dd(b)*dd(c))**2
        cc(c, a) = (1.0_dp + dd(a))*dd(a)*(1.0_dp + dd(c))*(dd(a) - dd(b))*((dd(b) - dd(c))*dd(b)*dd(c))**2
        cc(4, a) = -(1.0_dp + dd(a))*(dd(a) - dd(b))*(dd(c) - dd(a))*((dd(b) - dd(c))*dd(b)*dd(c))**2
        bb(a) = cc(4, a)*dd(a)
        ff = ff*(1.0_dp + dd(a))/(dd(a)*(dd(a) - dd(b)))**2
      enddo
      bb(4) = -dd(1)*dd(2)*dd(3)*((dd(1) - dd(2))*(dd(2) - dd(3))*(dd(3) - dd(1)))**2
      ff = ff/2.0_dp
      do i = 1, 4
        Ans = Ans + F(i)*(cc(i, 1)*ll(1) + cc(i, 2)*ll(2) + cc(i, 3)*ll(3) + bb(i))
      enddo
      Ans = Ans*ff/(D(4) + hw)**2

      tetrahedron_integral = Ans*Det_t
    case default
      tetrahedron_integral = 0.0_dp
    end select

    return
  end function tetrahedron_integral

  !=========================================================
  subroutine tetrahedron_sort(a, b1, b2, t)
    !===========================================!
    !A simple bubble sort subroutine for size 4 !
    !Assuming the size of a is four             !
    !===========================================!

    real(kind=dp), DIMENSION(4), INTENT(INOUT) :: a, b1, b2
    real(kind=dp), DIMENSION(3, 3), INTENT(INOUT) :: t
    real(kind=dp), DIMENSION(4) :: b1_temp, b2_temp
    real(kind=dp), DIMENSION(3, 4) :: t_temp
    INTEGER :: i, j
    INTEGER, DIMENSION(4) :: reference
    real(kind=dp) :: temp
    integer :: temp2

    DO i = 1, 4
      reference(i) = i
    ENDDO
    b1_temp = b1; b2_temp = b2
    DO j = 1, 3
      t_temp(j, 1) = 0
    ENDDO
    DO i = 2, 4
      DO j = 1, 3
        t_temp(j, i) = t(j, i - 1)
      ENDDO
    ENDDO
    DO i = 4, 1, -1
      DO j = 1, i - 1, +1
        IF (a(j) > a(j + 1)) THEN
          temp = a(j)
          a(j) = a(j + 1)
          a(j + 1) = temp
          temp2 = reference(j)
          reference(j) = reference(j + 1)
          reference(j + 1) = temp2
        END IF
      END DO
    END DO
    ! rearrange t(j,i)
    DO i = 2, 4
      DO j = 1, 3
        t(j, i - 1) = t_temp(j, reference(i)) - t_temp(j, reference(1))
      ENDDO
    ENDDO
    ! rearrange b(j,i)
    DO i = 1, 4
      b1(i) = b1_temp(reference(i))
      b2(i) = b2_temp(reference(i))
    ENDDO

  end subroutine tetrahedron_sort

!=========================================================
  function tetrahedron_jacobian(t, x, type)
    !=======================================!
    !                                       !
    ! Jacobian part of surface integrations !
    !                                       !
    !=======================================!
    use w90_constants, only: dp
    implicit none
    real(kind=dp) :: tetrahedron_jacobian
    real(kind=dp), DIMENSION(3, 3), INTENT(IN) :: t
    real(kind=dp), DIMENSION(3), INTENT(IN) :: x
    real(kind=dp), DIMENSION(3, 2) :: J
    real(kind=dp) :: y
    real(kind=dp) :: Ans
    INTEGER, INTENT(IN) :: type
    INTEGER :: j_, k, a, b, c, d

    IF (type == 1) THEN
      J(1, 1) = -x(1); J(1, 2) = -x(1)
      J(2, 1) = x(2); J(2, 2) = 0.0_dp
      J(3, 1) = 0.0_dp; J(3, 2) = x(3)
    ELSE IF (type == 2) THEN
      J(1, 1) = x(1); J(1, 2) = 0.0_dp
      J(2, 1) = 0.0_dp; J(2, 2) = x(2)
      J(3, 1) = 1.0_dp - x(1) - x(3); J(3, 2) = -x(3)
    ELSE IF (type == 3) THEN
      y = x(1)*(x(2) - 1.0_dp)*x(3)/(-x(2) + x(1)*x(2) + x(2)*x(3) - x(1)*x(3))
      J(1, 1) = x(1) - y; J(1, 2) = -y
      J(2, 1) = y - 1.0_dp; J(2, 2) = y - 1.0_dp + x(2)
      J(3, 1) = 1.0_dp - x(1); J(3, 2) = 0.0_dp
    ELSE
      J(1, 1) = x(1); J(1, 2) = 0.0_dp
      J(2, 1) = 0.0_dp; J(2, 2) = x(2)
      J(3, 1) = 1.0_dp - x(1) - x(3); J(3, 2) = 1.0_dp - x(2) - x(3)
    END IF
    Ans = 0.0_dp
    DO j_ = 1, 3
      DO k = 1, 3
        DO a = 1, 3
          DO b = 1, 3
            DO c = 1, 3
              DO d = 1, 3
                Ans = Ans + t(j_, a)*t(j_, b)*t(k, c)*t(k, d)*J(a, 1)*J(c, 2)*(J(b, 1)*J(d, 2) - J(b, 2)*J(d, 1))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    tetrahedron_jacobian = SQRT(ABS(Ans))

  end function tetrahedron_jacobian

  function tetrahedron_log1p(x)

    use w90_constants, only: dp

    implicit none

    real(kind=dp) :: tetrahedron_log1p
    real(kind=dp), intent(in) :: x
    real(kind=dp) :: y, z

    if (ABS(x) > 0.5_dp) then
      tetrahedron_log1p = LOG(ABS(1.0_dp + x))
    else
      y = 1.0_dp + x
      z = y - 1.0_dp
      if (z == 0) then
        tetrahedron_log1p = x
      else
        tetrahedron_log1p = x*LOG(y)/z
      endif
    endif

  end function tetrahedron_log1p

end module w90_tetrahedron
