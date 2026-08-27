# Need mpi library etc in paths as well as python requirements of serial example
# mpiexec -n 4 python mpi-example.py
import numpy
from mpi4py import MPI

# Maybe should have a common name...
import wan90mpi as wan90

ftn_output = wan90.w90_library.w90_get_fortran_stdout()
ftn_error = wan90.w90_library.w90_get_fortran_stderr()

data = wan90.w90_library.lib_common_type()
#w90data = wan90.w90_library.lib_wannier_type()

#wan90.w90_library.w90_set_comm(data, MPI.COMM_WORLD.py2f())
data.comm.comm = MPI.COMM_WORLD.py2f()

status = wan90.w90_library_extra.input_reader_special(data, "diamond", ftn_output, ftn_error)
status = wan90.w90_library.w90_input_reader(data, ftn_output, ftn_error)

if not data.kmesh_info.explicit_nnkpts :
    status = wan90.w90_library.w90_create_kmesh(data, ftn_output, ftn_error)

# attempt to duplicate comms_array_split packing
kpts = numpy.zeros(data.num_kpts, dtype=numpy.int32)
nproc = MPI.COMM_WORLD.Get_size()
my_proc = MPI.COMM_WORLD.Get_rank()
pts_per_rank = data.num_kpts // nproc
extra_pts = data.num_kpts - nproc * pts_per_rank
if extra_pts > 0:
    pts_per_rank = pts_per_rank + 1

k = 0
cur_proc = 0
cnt = 0
while k < data.num_kpts:
    kpts[k] = cur_proc
    k = k + 1
    cnt = cnt + 1
    if cnt == pts_per_rank:
        cnt = 0
        cur_proc = cur_proc + 1
        extra_pts = extra_pts - 1
        if extra_pts == 0:
            pts_per_rank = pts_per_rank - 1

wan90.w90_library_extra.set_kpoint_distribution(data, kpts, ftn_output, ftn_error)

counts = numpy.zeros(nproc, dtype=numpy.int32)
#displs = numpy.zeros(nproc, dtype=numpy.int32)
cur_proc = 0
k = 0
cnt = 0
#displs[0] = 0
while k < data.num_kpts:
    if kpts[k] != cur_proc:
        counts[cur_proc] = cnt
        cur_proc = cur_proc + 1
#        if cur_proc <= nproc:
#            displs[cur_proc] = k
        cnt = 0
    cnt = cnt + 1
    k = k + 1
counts[nproc-1] = cnt
#wan90.w90_library.w90_set_kpoint_block(data, counts, displs)

#m_matrix = numpy.zeros((data.num_wann, data.num_wann, data.kmesh_info.nntot, data.num_kpts), dtype=numpy.cdouble, order='F')
u_matrix = numpy.zeros((data.num_wann, data.num_wann, data.num_kpts), dtype=numpy.cdouble, order='F')

#wan90.w90_library.w90_set_m_matrix(m_matrix)
wan90.w90_library.w90_set_u_matrix(data, u_matrix)

m_matrix = numpy.zeros((data.num_bands, data.num_bands, data.kmesh_info.nntot, counts[my_proc]), dtype=numpy.cdouble, order='F')
wan90.w90_library.w90_set_m_local(data, m_matrix)
u_opt = numpy.zeros((data.num_bands, data.num_wann, data.num_kpts), dtype=numpy.cdouble, order='F')
wan90.w90_library.w90_set_u_opt(data, u_opt)
status = wan90.w90_library_extra.overlaps(data, ftn_output, ftn_error)

if data.num_wann == data.num_bands:
    status = wan90.w90_library.w90_project_overlap(data, ftn_output, ftn_error)
else:
    eigval = numpy.zeros((data.num_bands, data.num_kpts), dtype=numpy.double, order='F')
    status = wan90.w90_library_extra.read_eigvals(data, eigval, ftn_output, ftn_error)
    wan90.w90_library.w90_set_eigval(data, eigval)
    status = wan90.w90_library.w90_disentangle(data, ftn_output, ftn_error)
    if status == 1:
        exit

status = wan90.w90_library.w90_wannierise(data, ftn_output, ftn_error)

#wan90.w90_library.w90_checkpoint(data, "postwann", ftn_output, ftn_error)

#print (data.num_wann)

#wan90.w90_library.w90_plot(data, ftn_output, ftn_error)

#wan90.w90_library.w90_transport(data, ftn_output, ftn_error)

if my_proc == 0:
    wan90.w90_library_extra.print_times(data, ftn_output)
