# Need PYTHONPATH to contain wrapper directory
# Also need LD_LIBRARY_PATH to include it for libwannier_serial.so
import numpy
import wan90

ftn_output = wan90.w90_library.w90_get_fortran_stdout()
ftn_error = wan90.w90_library.w90_get_fortran_stderr()

data = wan90.w90_library.lib_common_type()
status = wan90.w90_library_extra.input_reader_special(data, "diamond", ftn_output, ftn_error)
status = wan90.w90_library.w90_input_reader(data, ftn_output, ftn_error)

if not data.kmesh_info.explicit_nnkpts :
    status = wan90.w90_library.w90_create_kmesh(data, ftn_output, ftn_error)

# create dummy distribution, all done on proc 0
kpts = numpy.zeros(data.num_kpts, dtype=numpy.int32)
wan90.w90_library_extra.set_kpoint_distribution(data, kpts, ftn_output, ftn_error)

#counts = numpy.zeros(1, dtype=numpy.int32)
#counts[0]=data.num_kpts
#displs = numpy.zeros(1, dtype=numpy.int32)
#wan90.w90_library.w90_set_kpoint_block(data, counts, displs)

m_matrix = numpy.zeros((data.num_bands, data.num_bands, data.kmesh_info.nntot, data.num_kpts), dtype=numpy.cdouble, order='F')
wan90.w90_library.w90_set_m_local(data, m_matrix)
u_opt = numpy.zeros((data.num_bands, data.num_wann, data.num_kpts), dtype=numpy.cdouble, order='F')
wan90.w90_library.w90_set_u_opt(data, u_opt)
u_matrix = numpy.zeros((data.num_wann, data.num_wann, data.num_kpts), dtype=numpy.cdouble, order='F')
wan90.w90_library.w90_set_u_matrix(data, u_matrix)
status = wan90.w90_library_extra.overlaps(data, ftn_output, ftn_error)

if data.num_wann == data.num_bands:
    status = wan90.w90_library.w90_project_overlap(data, ftn_output, ftn_error)
else:
    # allocate problem here (and seedname, and ierr not out) + set_eigval() "cnt55"
    eigval = numpy.zeros((data.num_bands, data.num_kpts), dtype=numpy.double, order='F')
    status = wan90.w90_library_extra.read_eigvals(data, eigval, ftn_output, ftn_error)
    wan90.w90_library.w90_set_eigval(data, eigval)
    status = wan90.w90_library.w90_disentangle(data, ftn_output, ftn_error)

status = wan90.w90_library.w90_wannierise(data, ftn_output, ftn_error)

#status = wan90.w90_library.w90_checkpoint(data, "postwann", ftn_output, ftn_error)

status = wan90.w90_library.w90_plot(data, ftn_output, ftn_error)

#status = wan90.w90_library.w90_transport(data, ftn_output, ftn_error)

wan90.w90_library_extra.print_times(data, ftn_output)
