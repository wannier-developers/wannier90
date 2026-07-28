# `wannier90` as a library

Starting with Wannier90 version 4, all functionality of the standalone
executable is available via a Fortran library interface.

The library interface consists of a mechanism to set options and perform
initialisation.  After processing options, the library can be probed for the
${\mathbf{k}}$,${\mathbf{k'}}$ combinations required for calculation of the
$M_{mn}^{(\mathbf{k,k'})}$ matrix. Execution of the library then proceeds using
a set of functions which take pointers to the main matrix quantities and a set
of high level functions which may be called to carry out disentanglement, MLWF
optimisation, plotting and those transport features implemented in Wannier90.
Taking pointers to matrices allocated by the calling code prevents duplication
of potentially large matrices.

Variables set using the library interface take the same name as the
corresponding tokens in the Wannier90 input file.  A minimal set of variables
must be set using the function interface; all others may be set using the
function interface or may be read from the Wannier90 (.win) input file.  An
optional library function is provided to read the input file and collect any
such additional variables.

The library makes use of a data object to record state; this object must be
passed to all Wannier90 library calls.  Data members within the object may be
accessed directly by the calling code, but it is not envisioned that this will
be necessary in normal use.

The library functions in parallel using MPI, with the same performance as the
standalone code.  The largest matrix input, the overlaps matrix
$M_{mn}^{(\mathbf{k,b})}$, is supplied to the library in decomposed form.  The
$U$ and $U^{opt}$ matrices are replicated across ranks.  When the parallel
version of the library is used, a valid and initialised MPI communicator must
be passed to the library.  Compilation of the parallel MPI library must be done
using the same MPI system and interface level as the calling program.

Error handling in the library has been implemented to ensure that errors
encountered by the library need not cause termination of the calling program.
Each of the library functions returns (as an argument) an integer error variable
which is non-zero in case of error.  This allows the calling code to ignore
errors in the Wannier90 library or to exit using the calling program's own error
handling mechanism.  Synchronisation of error state across MPI ranks is
performed prior to all MPI operations.  Error handling internally adopts the
strategy developed by [Bálint Aradi](https://github.com/aradi/errorfx).

Section [Using the Library](#using-the-library) provides detailed instructions
on using the library interface.

Section [Compiling and Linking](#compiling-and-linking) describes how to compile
and link to the Wannier90 library, including the names of the library files.

Section [Examples](#examples) documents some minimal Fortran examples that use
the library in serial and parallel.  These examples are distributed in the
directory `test-suite/library-mode-test`.

Section [C-interface](#c-interface) lists the C interface functions and shows
how to use them.  A test is available in `test-suite/library-mode-test-C-interface`.

$M_{mn}^{(\mathbf{k,b})}$ overlaps (Ref. [@marzari-prb97], Eq. (25)) and
$A_{mn}^{(\mathbf{k})}=\left\langle \psi_{m\mathbf{k}}|g_{n}\right\rangle$
projections (Ref. [@marzari-prb97], Eq. (62); Ref. [@souza-prb01], Eq. (22)).

## Library Specific Keywords

<!-- markdownlint-disable MD013 -->
{{ read_csv('docs/parameters/w90-library-parameters.csv', colalign=('left', 'center', 'left')) }}
<!-- markdownlint-enable MD013 -->

- the `distk` array communicates to Wannier90 how k-points are distributed over
MPI ranks.  In particular this affects the decomposition of the $M$-matrix.

- `total_bands` is provided for convenience when working with `exclude_bands`
when these determine the number of bands in the Wannier90 calculation instead
of vice-versa.

- `dump_inputs` writes a minimal .win file and .mmn, .amn and .eig files for
use by standalone executable.

## Using the Library

```fortran title="Fortran"
  USE w90_library, ONLY : lib_common_type, w90_set_comm, w90_disentangle, &
                          w90_project_overlap, w90_wannierise, w90_set_option, &
                          w90_input_setopt, w90_get_nn, w90_get_nnkp, &
                          w90_get_gkpb, w90_get_proj, w90_get_centres, &
                          w90_get_spreads, w90_plot, w90_set_eigval, &
                          w90_set_u_opt, w90_set_m_local, w90_set_u_matrix, &
                          w90_input_reader, w90_transport
```

The library exposes a number of functions via a Fortran module that should be
*use*d.  Use follows these steps:

1. create an instance of the [library data structure](#lib_common_type)
2. set necessary control parameters/flags/options using
   [w90_set_option](#w90_set_option)
3. (optionally) set other parameters/flags/options using
   [w90_set_option](#w90_set_option)
4. process parameters/flags/options using [w90_input_setopt](#w90_input_setopt)
5. (optionally) set other parameters/flags/options by reading the (.win) input
   file using [w90_input_reader](#w90_input_reader)
6. obtain the number of finite difference neighbours using
   [w90_get_nn](#w90_get_nn)
7. obtain the finite difference neighbour indices using
   [w90_get_nnkp](#w90_get_nnkp)
8. obtain the finite difference neighbour BZ offsets using
   [w90_get_gkpb](#w90_get_gkpb)
9. get projector definition corresponding to input string
   [w90_get_proj](#w90_get_proj)  (your code may have its own mechanism for
choosing projectors; this is a convenience function that uses Wannier90 to
parse the .win projector string into a site,l,m,zaxis configuration)
10. calculate projections and overlap matrices
11. pass pointers to $U$, $M$, $U^{opt}$ and (if disentangling) eigenvalue
    matrices using [w90_set_eigval](#w90_set_eigval),
[w90_set_u_matrix](#w90_set_u_matrix), etc
12. initial projections, $A$, should be stored in the $U^{opt}$ matrix
13. if disentangling is required, call [w90_disentangle](#w90_disentangle)
14. prepare for MLWF algorithm by projecting $M$, $U^{opt}$ onto subspace by
    calling [w90_project_overlap](#w90_project_overlap)
15. call [w90_wannierise](#w90_wannierise) (this may be called multiple times
    in order to test convergence)
16. obtain centres and spreads with [w90_get_centres](#w90_get_centres) and
    [w90_get_spreads](#w90_get_spreads)

### lib_common_type

An instance of this type must be declared and passed to the following library
subroutines.

```fortran title="Fortran"
  type(lib_common_type) :: wannier_data
```

### w90_set_comm

Pass the MPI communicator to the library.

Arguments are an instance of the Wannier90 library datatype and an initialised
MPI communicator.

```fortran title="Fortran"
  subroutine w90_set_comm(common_data, comm)

    type(lib_common_type), intent(inout) :: common_data
#ifdef MPI08
    type(mpi_comm), intent(in) :: comm
#
```

!!! warning
    The MPI inclusion style (Fortran 77 header file, Fortran 90 module or
    Fortran 2008 module) used to compile Wannier90 in parallel must match that
    used to compile the calling code (because it determines the type of the
    communicator variable).

### w90_disentangle

Perform disentanglement.  Arguments are the Wannier90 library object, integer
Fortran unit numbers for standard error and output streams and an integer status
value.  Successful disentanglement returns ierr zero.

```fortran title="Fortran"
  subroutine w90_disentangle(common_data, istdout, istderr, ierr)

    type(lib_common_type), intent(inout) :: common_data
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
```

### w90_project_overlap

### w90_wannierise

Perform MLWF optimisation.  Arguments are the Wannier90 library object, integer
Fortran unit numbers for standard error and output streams and an integer status
value.  Successful optimisation returns ierr zero; a negative value indicates
that the maximum number of iterations was reached without convergence
(suggesting that w90_wannierise may usefully be called again).

```fortran title="Fortran"
  subroutine w90_wannierise(common_data, istdout, istderr, ierr)

    type(lib_common_type), intent(inout) :: common_data
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
```

### w90_set_option

Call w90_set_option to specify options guiding the operation of the library.
Variable names and acceptable data types are the same as for the wannier90.x
executable's input (.win) file (with some exceptions associated with
projections, see below).  Arguments are the Wannier90 library object, a string
specifying the option and the corresponding data.  The type of the data must
match that expected for the specified option.

w90_set_option can be called any number of times, in any order (it accumulates
data in an internal structure).

w90_input_setopt must be called after the last call to w90_set_option.

```fortran title="Fortran"
    CALL w90_set_option(w90main, 'kpoints', kpt_latt)   ! rank-2 real array
    CALL w90_set_option(w90main, 'mp_grid', mp_grid)    ! rank-1 integer array
    CALL w90_set_option(w90main, 'num_iter', num_iter)  ! integer scalar
    CALL w90_set_option(w90main, 'spinors', noncolin)   ! logical
    CALL w90_set_option(w90main, 'dis_froz_max', dis_froz_max) ! real scalar
```

### w90_input_setopt

Call w90_input_setopt to initialise the library after setting all options using
w90_set_option.

Arguments are the Wannier90 library object, a string specifying the base of
filenames to be created (or read) by the library, integer Fortran unit numbers
for standard error and output streams and an integer status value.  Successful
processing of input returns ierr zero.

!!! warning
    w90_input_setopt must follow all w90_set_option calls.

```fortran title="Fortran"
  subroutine w90_input_setopt(common_data, seedname, istdout, istderr, ierr)

    character(len=*), intent(in) :: seedname
    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(lib_common_type), intent(inout) :: common_data
```

### w90_get_nn

Returns number of points used in finite difference scheme (k' for each k).  This
should be used to allocate arrays.

Must follow w90_input_setopt

```fortran title="Fortran"
  subroutine w90_get_nn(common_data, nn, istdout, istderr, ierr)

    integer, intent(out) :: nn, ierr
    integer, intent(in) :: istdout, istderr
    type(lib_common_type), intent(inout) :: common_data
```

### w90_get_nnkp

Returns tabulation of finite difference scheme k' indices.

On entry, nnkp must be allocated with dimension(num_kpts,nn) where nn is the
number of points in the FD scheme (returned by w90_get_nn) and num_kpts is the
total number of k-points (FBZ).

Must follow w90_input_setopt

```fortran title="Fortran"
  subroutine w90_get_nnkp(common_data, nnkp, istdout, istderr, ierr)

    integer, intent(out) :: nnkp(:, :), ierr
    integer, intent(in) :: istdout, istderr
    type(lib_common_type), intent(inout) :: common_data
```

### w90_get_gkpb

Get BZ offsets corresponding to each k,k' pair described in nnkp.

On entry, gkpb must be allocated with dimension(3,num_kpts,nn) where nn is the
number of points in the FD scheme (returned by w90_get_nn) and num_kpts is the
total number of k-points (FBZ).

Must follow w90_input_setopt

```fortran title="Fortran"
  subroutine w90_get_gkpb(common_data, gkpb, istdout, istderr, ierr)

    integer, intent(out) :: gkpb(:, :, :), ierr
    integer, intent(in) :: istdout, istderr
    type(lib_common_type), intent(inout) :: common_data
```

### w90_get_proj

The Wannier90 library can be used to translate a set of strings describing
projections into a set of arrays containing the position, angular momentum,
spin projection, etc, that can be directly used for setting up the projection
matrix $A$.

This is a convenience function for plane-wave codes that do not have a natural
basis for projection, or to take advantage of Wannier90's flexible projection
description mechanism.

Parsing the projection tokens requires the specification of both the atom
positions and their symbols, as well as the projection strings.

Unlike for other array inputs, when specifying projections it is necessary to
call w90_set_option successively for each projection string.

w90_get_proj takes as arguments a series of pre-allocated arrays, which must be
adequately sized

```fortran title="Fortran"
  subroutine w90_get_proj(common_data, n, site, l, m, s, rad, x, z, sqa, zona, &
                          istdout, istderr, ierr)

    integer, intent(in) :: istdout, istderr
    integer, intent(inout) :: n, l(:), m(:), s(:)
    integer, intent(inout) :: rad(:)
    integer, intent(out) :: ierr
    real(kind=dp), intent(inout) :: site(:, :)
    real(kind=dp), intent(inout) :: sqa(:, :), z(:, :), x(:, :), zona(:)
    type(lib_common_type), intent(in), target :: common_data
```

| argument | type                             | purpose                        |
|----------|----------------------------------|--------------------------------|
| n        | integer                          | number of projections          |
| site     | double precision, dimension(3,n) | origin of projections          |
| l        | integer, dimension(n)            | angular quantum number         |
| m        | integer, dimension(n)            | orbital quantum number         |
| s        | integer, dimension(n)            | spin quantum number            |
| rad      | integer, dimension(n)            | model for radial function      |
| x        | double precision, dimension(3,n) | vector defining $x$            |
| z        | double precision, dimension(3,n) | vector defining $z$            |
| sqa      | double precision, dimension(3,n) | spin quantisation vector       |
| zona     | double precision, dimension(n)   |                                |
| istdout  | integer                          | Fortran unit number for stdout |
| istderr  | integer                          | Fortran unit number for stderr |
| ierr     | integer                          | error status, zero on success  |

Must follow w90_input_setopt.

### w90_get_centres

Probe library for position of WF centres.

On entry, centres must be a double precision array allocated with
dimension(3,num_wannier)

Must follow w90_input_setopt

(No error code.)

```fortran title="Fortran"
  subroutine w90_get_centres(common_data, centres)

    real(kind=dp), intent(out) :: centres(:, :)
    type(lib_common_type), intent(in) :: common_data
```

### w90_get_spreads

Probe library for position of WF spreads.

On entry, spreads must be a double precision array allocated with
dimension(num_wannier).

Must follow w90_input_setopt

(No error code.)

```fortran title="Fortran"
  subroutine w90_get_spreads(common_data, spreads)

    real(kind=dp), intent(out) :: spreads(:)
    type(lib_common_type), intent(in) :: common_data
```

### w90_plot

Perform plotting functions.  Arguments are the Wannier90 library object, integer
Fortran unit numbers for standard error and output streams and an integer status
value.  Successful optimisation returns ierr zero.

```fortran title="Fortran"
  subroutine w90_plot(common_data, istdout, istderr, ierr)

    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(lib_common_type), intent(inout) :: common_data
```

### w90_transport

Perform those transport property calculations that are implemented in the
Wannier90 main executable (note that this is a much smaller functionality than
that afforded by postw90.x).  Arguments are the Wannier90 library object,
integer Fortran unit numbers for standard error and output streams and an
integer status value.  Successful optimisation returns ierr zero.

```fortran title="Fortran"
  subroutine w90_transport(common_data, istdout, istderr, ierr)

    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(lib_common_type), intent(inout) :: common_data
```

### w90_set_eigval

Pass a pointer to a preexisting double precision array of dimension(nbands,
num_kpts).

This must be accomplished before calling w90_disentangle.  Contents of the array
(the eigenvalues of each band at each k-point) need not be valid at call to
set_eigval, but must be valid at call to w90_disentangle.

This must follow w90_input_setopt.

```fortran title="Fortran"
  subroutine w90_set_eigval(common_data, eigval)

    type(lib_common_type), intent(inout) :: common_data
    real(kind=dp), intent(in), target :: eigval(:, :)
```

### w90_set_u_opt

Pass a pointer to a preexisting double precision complex array of
dimension(nbands, num_wannier, num_kpts).

This must be accomplished before calling w90_disentangle.

This must follow w90_input_setopt.

```fortran title="Fortran"
  subroutine w90_set_u_opt(common_data, u_matrix_opt)

    type(lib_common_type), intent(inout) :: common_data
    complex(kind=dp), intent(inout), target :: u_matrix_opt(:, :, :)
```

### w90_set_m_local

Pass a pointer to a preexisting double precision complex array of
dimension(nbands, nbands, nn, nklocal), where nklocal is the number of kpoints
associated with this rank and nn is the number of finite difference k-point
neighbours.  In serial, nklocal equals num_kpoints.  The distribution of M
across k-points is described by the `distk` array.

This must be done before calling w90_disentangle or w90_wannierise and must
follow w90_input_setopt.

```fortran title="Fortran"
  subroutine w90_set_m_local(common_data, m_matrix_local)

    type(lib_common_type), intent(inout) :: common_data
    complex(kind=dp), intent(inout), target :: m_matrix_local(:, :, :, :)
```

### w90_set_u_matrix

Pass a pointer to a preexisting double precision complex array of
dimension(num_wannier, num_wannier, num_kpoints).

The full U matrix is duplicated on all ranks.

This must be accomplished before calling w90_disentangle or w90_wannierise.

This must follow w90_input_setopt.

```fortran title="Fortran"
  subroutine w90_set_u_matrix(common_data, u_matrix)

    type(lib_common_type), intent(inout) :: common_data
    complex(kind=dp), intent(inout), target :: u_matrix(:, :, :)
```

### w90_input_reader

w90_input_reader provides an optional mechanism for passing additional flags to
the library using the input (.win) file.  All valid input tokens of the main
program may be specified in this way, except for variables listed in
[w90_set_option](#w90_set_option), i.e. the most important variables defining the
calculation must be specified by w90_set_option.

w90_input_reader must be called after w90_input_setopt.

If optional variables are specified using w90_input_setopt and then subsequently
found in the input file during a call to w90_input_reader, then the input file
value overwrites the former value.

The name of the input file that is read is 'seedname.win' where seedname is the
string passed to w90_input_setopt.

Arguments are the Wannier90 library object, integer Fortran unit numbers for
standard error and output streams and an integer status value.  Successful
optimisation returns ierr zero.

```fortran title="Fortran"
  subroutine w90_input_reader(common_data, istdout, istderr, ierr)

    integer, intent(in) :: istdout, istderr
    integer, intent(out) :: ierr
    type(lib_common_type), intent(inout) :: common_data
```

## Compiling and Linking

Depending on whether a serial or MPI compilation has happened, different
library files are produced:

| filename              | description               |
|-----------------------|---------------------------|
| libwannier90.a        | static library, serial    |
| libwannier90_mpi.a    | static library, parallel  |
| libwannier90.so.4     | dynamic library, serial   |
| libwannier90_mpi.so.4 | dynamic library, parallel |
| w90_library.mod       | fortran module            |
| wannier90.h           | C header (compile with WANNIER90_WITH_C=ON) |

MPI operations are only supported by "libwannier90_mpi".

To generate dynamic libraries using GNU make, you need to build the target
"dynlib"

```bash
    make dynlib
```

The directory containing the libraries must be added to your link line
and library search path.

You need to "use" the library fortran module to have access to the type and
function definitions in your code: the directory containing the module file
must be added to your "include" path.

The library should be compiled with the same compiler as the calling code.

If you link to the parallel library, you must initialise the MPI system
(mpi_init) and must pass a valid communicator to the library (failure to do
this results in an error:  "Error: parallel Wannier90 library invoked with
invalid communicator, exiting.  Use w90_set_comm()!" )

## Examples

See directory: test-suite/library-mode-test/

## C Interface

A C interface is provided using fortran 2003's iso_c_binding; the functions
available are described in the header file "wannier90.h" and follow the
fortran interface except for the passing of multi-dimensional arrays as options,
where different functions must be called for 1-d and 2-d data.

An example that re-implements the main wannier90 executable is given in
directory test-suite/library-mode-test-C-interface/

## Python interface

The directory `wrap/` in the sources contains an Python wrapping of the Fortran
library.  It is constructed using the [f90wrap
package](https://github.com/jameskermode/f90wrap) (See also DOI
10.1088/1361-648X/ab82d2).  Specifiy `F90WRAP` in the build configuration.

### Build instructions

- Make sure f90wrap is installed
- cd wrap
- make -f Makefile

Edit the makefiles as appropriate.

### Use

- if wannier90 is installed in e.g. W90DIR then export
  PYTHONPATH=$(W90DIR)/wrap

- since the wannier90 wrapper is built as a separate shared lib, also export
  LD_LIBRARY_PATH=$(W90DIR)/wrap

- run python3 and type commands or 'python3 script.py' (the mpi version is
  something like 'mpirun -np n python3 script.py')

The serial example is for diamond; just edit the seedname in the read call and
alter what you want to run after wannierise. MPI the same with the extra
imports.

example-dos.py tests one of the DOS examples with the draft postw90 interface.

## Frequently Asked Questions

### Can symmetry adapted WF be calculated using the library?

Not yet.

### C-interface is not built

The library needs to be built in serial for c_interface and you must set the
environment variable WANNIER90_WITH_C before invoking make.

### Error messages "Fortran runtime error: File cannot be deleted"

Likely you are using the incorrect mpirun command (causing several serial
instances to be executed concurrently); the mpirun command must correspond to
the library used during compilation.

### Inconsistency in library compilation and compilation of code calling library

Fortran modules are compiler specific.  Errors such as "Fatal Error: Reading
module ‘w90_library.mod’ at line 1 column 2: Unexpected EOF" suggest that
different compilers (or compiler versions) are being used.

### Inconsistency in MPI flavour used

MPI libraries for Fortran may support both old style ('use mpi') and more
modern interfaces ('use mpi_f08'); depending on the interface used to build the
library, communicator objects may be specially typed or treated as integer and
this must be consistent with what is used in the calling code.  This behaviour
can be adjusted by compiling Wannier90 with use of COMMS=MPI90 , MPI08 or MPIH,
as needed.

```bash
696 | call input_setopt(w90main, filename, stdout, stderr, ierr, comm)
    |
Error: Type mismatch in argument ‘comm’ at (1); passed INTEGER(4) to TYPE(mpi_comm)
```

### Build requirements

1. Do not use 'no_realloc_lhs' option

2. [libasan](https://github.com/google/sanitizers/wiki/AddressSanitizer) may fail
for gcc/openmpi builds (at pmpi_init() call); similarly ifx/impi builds with
"-check" ("-check all")

### Parallel libraries cannot be called from serial code

Parallel libraries cannot be called from serial code because mpi_init() must be
called before any call to the Wannier90 library when it is built with MPI.
Recomple a serial version if necessary.
