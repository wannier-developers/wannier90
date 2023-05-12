# Need mpi library etc in paths as well as python requirements of serial example
# mpiexec -n 4 python mpi-example.py
import numpy
from mpi4py import MPI

# Maybe should have a common name...
import wan90mpi as wan90

ftn_output = wan90.w90_library.get_fortran_stdout()
ftn_error = wan90.w90_library.get_fortran_stderr()

data = wan90.w90_library.lib_common_type()
w90data = wan90.w90_library.lib_wannier_type()
#comm = wan90.w90_comms.w90_comm_type()

wan90.w90_library.set_parallel_comms(data, MPI.COMM_WORLD.py2f())

status = wan90.w90_library.input_reader(data, w90data, "diamond", ftn_output, ftn_error)
#status = wan90.w90_library.input_reader(data, w90data, "cnt55", ftn_output, ftn_error)

exit

if not data.kmesh_info.explicit_nnkpts :
    status = wan90.w90_library.create_kmesh(data, ftn_output, ftn_error)

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

wan90.w90_library.set_kpoint_distribution(data, kpts)

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
#wan90.w90_library.set_kpoint_block(data, counts, displs)

#m_matrix = numpy.zeros((data.num_wann, data.num_wann, data.kmesh_info.nntot, data.num_kpts), dtype=numpy.cdouble, order='F')
u_matrix = numpy.zeros((data.num_wann, data.num_wann, data.num_kpts), dtype=numpy.cdouble, order='F')

#wan90.w90_library.set_m_matrix(w90data, m_matrix)
wan90.w90_library.set_u_matrix(data, u_matrix)

m_matrix_loc = numpy.zeros((data.num_wann, data.num_wann, data.kmesh_info.nntot, counts[my_proc]), dtype=numpy.cdouble, order='F')
wan90.w90_library.set_m_matrix_local(w90data, m_matrix_loc)

if data.num_wann == data.num_bands:
    u_opt = numpy.zeros((1, 1, 1), dtype=numpy.cdouble, order='F')
    wan90.w90_library.set_u_opt(data, u_opt)
    status = wan90.w90_library.overlaps(data, w90data, ftn_output, ftn_error)
else:
    a_matrix = numpy.zeros((data.num_bands, data.num_wann, data.num_kpts), dtype=numpy.cdouble, order='F')
    wan90.w90_library.set_a_matrix(w90data, a_matrix)
    u_opt = numpy.zeros((data.num_bands, data.num_wann, data.num_kpts), dtype=numpy.cdouble, order='F')
    wan90.w90_library.set_u_opt(data, u_opt)
    m_orig = numpy.zeros((data.num_bands, data.num_bands, data.kmesh_info.nntot, data.num_kpts), dtype=numpy.cdouble, order='F')
    wan90.w90_library.set_m_orig(w90data, m_orig)
    status = wan90.w90_library.overlaps(data, w90data, ftn_output, ftn_error)
    status = wan90.w90_library.disentangle(data, w90data, ftn_output, ftn_error)
    if status == 1:
        exit

status = wan90.w90_library.wannierise(data, w90data, ftn_output, ftn_error)

#wan90.w90_library.checkpoint(data, w90data, "postwann", ftn_output, ftn_error)

#wan90.w90_library.transport(helper, w90data, ftn_output, ftn_error, status)

#print (data.num_wann)

#wan90.w90_library.plot_files(data, w90data, ftn_output, ftn_error, status)

if my_proc == 0:
    wan90.w90_library.print_times(data, ftn_output)
