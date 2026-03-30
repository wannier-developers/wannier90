! Reproducer for bug in w90_readwrite_get_projections (library mode):
! the last projection entry passed via w90_set_option is always skipped.
!
! Root cause: in library mode line_s=0 and line_e=settings%num_entries,
! but the projection loop runs "do line = line_s+1, line_e-1", so the
! entry at index line_e is never visited.
!
! Fix: set line_e = settings%num_entries + 1 in the library branch of
! w90_readwrite_get_projections (readwrite.F90).
!
! To compile (from the wannier90 develop directory):
!   mpif90-openmpi-gcc14 -I src/obj bug_reproduce.F90 libwannier90_mpi.a -o bug_reproduce
! To run:
!   mpirun -n 1 ./bug_reproduce
!
program bug_reproduce
  use w90_library
  use mpi
  implicit none

  ! Minimal H2 system: 2 s orbitals → num_wann = 2.
  ! Projections block has 2 entries: "H:s" passed twice (one per atom).
  ! In library mode these become settings entries 3 and 4 (after num_wann
  ! and mp_grid are already stored). With the bug, line_e = 4 and the loop
  ! runs 1..3, skipping entry 4 ("H:s" for the second atom). Counter = 1
  ! < num_wann = 2 → "too few projection functions defined".

  type(lib_common_type) :: w90
  integer :: ierr, stdout, stderr
  integer :: mp_grid(3) = [1, 1, 1]
  real(kind=8) :: unit_cell(3,3), kpoints(3,1), atoms_frac(3,2)
  character(len=2) :: symbols(2)
  integer :: num_wann

  call w90_get_fortran_stdout(stdout)
  call w90_get_fortran_stderr(stderr)

  call MPI_Init(ierr)

  num_wann = 2

  unit_cell = reshape([3.0d0, 0.0d0, 0.0d0, &
                       0.0d0, 3.0d0, 0.0d0, &
                       0.0d0, 0.0d0, 3.0d0], [3,3])
  kpoints(:,1) = [0.0d0, 0.0d0, 0.0d0]
  atoms_frac(:,1) = [0.0d0, 0.0d0, 0.0d0]
  atoms_frac(:,2) = [0.5d0, 0.0d0, 0.0d0]
  symbols = ['H ', 'H ']

  call w90_set_comm(w90, MPI_COMM_WORLD)
  call w90_set_option(w90, 'num_wann',      num_wann)
  call w90_set_option(w90, 'mp_grid',       mp_grid)
  call w90_set_option(w90, 'kpoints',       kpoints)
  call w90_set_option(w90, 'unit_cell_cart', unit_cell)
  call w90_set_option(w90, 'total_bands',   num_wann)
  call w90_set_option(w90, 'atoms_frac',    atoms_frac)
  call w90_set_option(w90, 'symbols',       symbols)
  ! Pass two projection lines via w90_set_option (library interface):
  call w90_set_option(w90, 'projections', 'f=0.0,0.0,0.0:s')   ! entry N-1
  call w90_set_option(w90, 'projections', 'f=0.5,0.0,0.0:s')   ! entry N  ← skipped by bug

  call w90_input_setopt(w90, 'wannier90', stdout, stderr, ierr)

  if (ierr /= 0) then
    write(*, '(a,i0)') 'FAIL: w90_input_setopt returned ierr = ', ierr
    write(*, '(a)') 'Expected success. Bug: last projection entry skipped (line_e off-by-one).'
    call MPI_Finalize(ierr)
    stop 1
  end if

  write(*, '(a)') 'PASS: both projection entries were parsed correctly.'

  call MPI_Finalize(ierr)

end program bug_reproduce
