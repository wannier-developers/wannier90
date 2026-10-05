program proj_overflow
  use mpi
  use w90_library

  implicit none

  type(lib_common_type) :: w90
  integer :: ierr, stdout, stderr
  integer :: num_proj
  integer :: mp_grid(3)
  integer :: l(8), m(8), s(8), rad(8)
  real(kind=8) :: kpoints(3, 1), unit_cell(3, 3)
  real(kind=8) :: atoms_frac(3, 1)
  real(kind=8) :: site(3, 8), x(3, 8), z(3, 8), sqa(3, 8), zona(8)
  character(len=2) :: symbols(1)

  mp_grid = (/1, 1, 1/)
  kpoints(:, 1) = (/0.0d0, 0.0d0, 0.0d0/)
  unit_cell = 0.0d0
  unit_cell(1, 1) = 10.0d0
  unit_cell(2, 2) = 10.0d0
  unit_cell(3, 3) = 10.0d0
  atoms_frac(:, 1) = (/0.25d0, 0.25d0, 0.25d0/)
  symbols(1) = 'Bi'

  call MPI_Init(ierr)
  call w90_get_fortran_stdout(stdout)
  call w90_get_fortran_stderr(stderr)

  call w90_set_comm(w90, MPI_COMM_WORLD)
  call w90_set_option(w90, 'num_bands', 4)
  call w90_set_option(w90, 'num_wann', 2)
  call w90_set_option(w90, 'mp_grid', mp_grid)
  call w90_set_option(w90, 'kpoints', kpoints)
  call w90_set_option(w90, 'unit_cell_cart', unit_cell)
  call w90_set_option(w90, 'symbols', symbols)
  call w90_set_option(w90, 'atoms_frac', atoms_frac)
  call w90_set_option(w90, 'projections', 'Bi:s,p')

  call w90_input_setopt(w90, 'proj_overflow', stdout, stderr, ierr)
  if (ierr /= 0) then
    write (stderr, '(a,i0)') 'w90_input_setopt failed with ierr=', ierr
    call MPI_Finalize(ierr)
    stop 1
  end if

  num_proj = size(l)
  call w90_get_proj(w90, num_proj, site, l, m, s, rad, x, z, sqa, zona, stdout, stderr, ierr)
  if (ierr /= 0) then
    write (stderr, '(a,i0)') 'w90_get_proj failed with ierr=', ierr
    call MPI_Finalize(ierr)
    stop 1
  end if

  if (num_proj /= 4) then
    write (stderr, '(a,i0,a)') 'Expected 4 projections, got ', num_proj, '.'
    call MPI_Finalize(ierr)
    stop 1
  end if

  call MPI_Finalize(ierr)
end program proj_overflow
