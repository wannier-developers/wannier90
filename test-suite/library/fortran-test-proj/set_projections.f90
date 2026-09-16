program set_projections
  use w90_library
  use mpi
  implicit none

  type(lib_common_type) :: w90
  integer :: ierr, stdout, stderr
  integer :: mp_grid(3) = [1, 1, 1]
  real(kind=8) :: unit_cell(3, 3), kpoints(3, 1), atoms_frac(3, 2)
  character(len=2) :: symbols(2)
  integer :: num_wann

  call w90_get_fortran_stdout(stdout)
  call w90_get_fortran_stderr(stderr)

  call MPI_Init(ierr)

  num_wann = 2

  unit_cell = reshape([3.0d0, 0.0d0, 0.0d0, &
                       0.0d0, 3.0d0, 0.0d0, &
                       0.0d0, 0.0d0, 3.0d0], [3, 3])
  kpoints(:, 1) = [0.0d0, 0.0d0, 0.0d0]
  atoms_frac(:, 1) = [0.0d0, 0.0d0, 0.0d0]
  atoms_frac(:, 2) = [0.5d0, 0.0d0, 0.0d0]
  symbols = ['H ', 'H ']

  call w90_set_comm(w90, MPI_COMM_WORLD)
  call w90_set_option(w90, 'num_wann', num_wann)
  call w90_set_option(w90, 'mp_grid', mp_grid)
  call w90_set_option(w90, 'kpoints', kpoints)
  call w90_set_option(w90, 'unit_cell_cart', unit_cell)
  call w90_set_option(w90, 'total_bands', num_wann)
  call w90_set_option(w90, 'atoms_frac', atoms_frac)
  call w90_set_option(w90, 'symbols', symbols)
  call w90_set_option(w90, 'projections', 'f=0.0,0.0,0.0:s')
  call w90_set_option(w90, 'projections', 'f=0.5,0.0,0.0:s')

  call w90_input_setopt(w90, 'wannier90', stdout, stderr, ierr)

  if (ierr /= 0) then
    write (*, '(a,i0)') 'FAIL: w90_input_setopt returned ierr = ', ierr
    call MPI_Finalize(ierr)
    stop 1
  end if

  call MPI_Finalize(ierr)

end program set_projections
