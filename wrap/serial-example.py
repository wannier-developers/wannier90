# Need PYTHONPATH to contain wrapper directory
# Also need LD_LIBRARY_PATH to include it for libwannier_serial.so
import numpy
import wan90

ftn_output = wan90.w90_helper_types.get_fortran_stdout()
ftn_error = wan90.w90_helper_types.get_fortran_stderr()

data = wan90.w90_helper_types.lib_global_type()
w90data = wan90.w90_helper_types.lib_w90_type()
comm = wan90.w90_comms.w90_comm_type()
status = wan90.w90_helper_types.input_reader(data, w90data, "diamond", ftn_output, ftn_error, comm)

if not data.kmesh_info.explicit_nnkpts :
    status = wan90.w90_helper_types.create_kmesh(data, ftn_output, ftn_error, comm)

# create dummy distribution, all done on proc 0
kpts = numpy.zeros(data.num_kpts, dtype=numpy.int32)
wan90.w90_helper_types.set_kpoint_distribution(data, kpts, ftn_output, ftn_error, comm)

#counts = numpy.zeros(1, dtype=numpy.int32)
#counts[0]=data.num_kpts
#displs = numpy.zeros(1, dtype=numpy.int32)
#wan90.w90_helper_types.set_kpoint_block(data, counts, displs)

#m_matrix = numpy.zeros((data.num_wann, data.num_wann, data.kmesh_info.nntot, data.num_kpts), dtype=numpy.cdouble, order='F')
u_matrix = numpy.zeros((data.num_wann, data.num_wann, data.num_kpts), dtype=numpy.cdouble, order='F')
#wan90.w90_helper_types.set_m_matrix(w90data, m_matrix)
wan90.w90_helper_types.set_u_matrix(data, u_matrix)
m_matrix_loc = numpy.zeros((data.num_wann, data.num_wann, data.kmesh_info.nntot, data.num_kpts), dtype=numpy.cdouble, order='F')
wan90.w90_helper_types.set_m_matrix_local(w90data, m_matrix_loc)

#m_matrix.flags.f_contiguous should be true

if data.num_wann == data.num_bands:
    u_opt = numpy.zeros((1, 1, 1), dtype=numpy.cdouble, order='F')
    wan90.w90_helper_types.set_u_opt(data, u_opt)
    status = wan90.w90_helper_types.overlaps(data, w90data, ftn_output, ftn_error, comm)
    status = wan90.w90_helper_types.projovlp(data, w90data, ftn_output, ftn_error, comm)
else:
    a_matrix = numpy.zeros((data.num_bands, data.num_wann, data.num_kpts), dtype=numpy.cdouble, order='F')
    wan90.w90_helper_types.set_a_matrix(w90data, a_matrix)
    u_opt = numpy.zeros((data.num_bands, data.num_wann, data.num_kpts), dtype=numpy.cdouble, order='F')
    wan90.w90_helper_types.set_u_opt(data, u_opt)
    m_orig = numpy.zeros((data.num_bands, data.num_bands, data.kmesh_info.nntot, data.num_kpts), dtype=numpy.cdouble, order='F')
    wan90.w90_helper_types.set_m_orig(w90data, m_orig)
    status = wan90.w90_helper_types.overlaps(data, w90data, ftn_output, ftn_error, comm)
    # allocate problem here (and seedname, and ierr not out) + set_eigval() "cnt55"
    eigval = numpy.zeros((data.num_bands, data.num_kpts), dtype=numpy.double, order='F')
    status = wan90.w90_helper_types.read_eigvals(data, eigval, ftn_output, ftn_error, comm)
    wan90.w90_helper_types.set_eigval(data, eigval)
    status = wan90.w90_helper_types.disentangle(data, w90data, ftn_output, ftn_error, comm)

status = wan90.w90_helper_types.wannierise(data, w90data, ftn_output, ftn_error, comm)


#status = wan90.w90_helper_types.checkpoint(data, w90data, "postwann", ftn_output, ftn_error, comm)

status = wan90.w90_helper_types.plot_files(data, w90data, ftn_output, ftn_error, comm)

#status = wan90.w90_helper_types.transport(data, w90data, ftn_output, ftn_error, comm)

wan90.w90_helper_types.print_times(data, ftn_output)
