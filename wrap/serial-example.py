# Need PYTHONPATH to contain wrapper directory
# Also need LD_LIBRARY_PATH to include it for libwannier_serial.so
import wan90

ftn_output = wan90.w90_helper_types.get_fortran_stdout()

data = wan90.w90_helper_types.lib_global_type()
plot = wan90.w90_helper_types.lib_plot_type()
tran = wan90.w90_helper_types.lib_transport_type()
comm = wan90.w90_comms.w90comm_type()
wan90.w90_helper_types.input_reader(data, plot, tran, "diamond", ftn_output, comm)

if not data.kmesh_info.explicit_nnkpts :
    wan90.w90_helper_types.create_kmesh(data, ftn_output, comm)

import numpy

m_matrix = numpy.zeros((data.num_wann, data.num_wann, data.kmesh_info.nntot, data.num_kpts), dtype=numpy.cdouble, order='F')
u_matrix = numpy.zeros((data.num_wann, data.num_wann, data.num_kpts), dtype=numpy.cdouble, order='F')
a_matrix = numpy.zeros((data.num_bands, data.num_wann, data.num_kpts), dtype=numpy.cdouble, order='F')

#m_matrix.flags.f_contiguous should be true

if data.num_wann == data.num_bands:
    m_orig = numpy.zeros((1, 1, 1, 1), dtype=numpy.cdouble, order='F')
    u_opt = numpy.zeros((1, 1, 1), dtype=numpy.cdouble, order='F')
    wan90.w90_helper_types.overlaps(data, a_matrix, m_matrix, m_orig, u_matrix, ftn_output, comm)
else:
    u_opt = numpy.zeros((data.num_bands, data.num_wann, data.num_kpts), dtype=numpy.cdouble, order='F')
    m_orig = numpy.zeros((data.num_bands, data.num_bands, data.kmesh_info.nntot, data.num_kpts), dtype=numpy.cdouble, order='F')
    wan90.w90_helper_types.overlaps(data, a_matrix, m_matrix, m_orig, u_matrix, ftn_output, comm)
    wan90.w90_helper_types.disentangle(data, a_matrix, m_matrix, m_orig, u_matrix, u_opt, ftn_output, comm)
    
wan90.w90_helper_types.wannierise(data, plot, tran, m_matrix, u_matrix, u_opt, ftn_output, comm)


#wan90.w90_helper_types.checkpoint(data, "postwann", m_matrix, u_matrix, u_opt, ftn_output, comm)

#wan90.w90_helper_types.transport(helper, plot, tran, u_matrix, u_opt, ftn_output, comm)

print (data.num_wann)

wan90.w90_helper_types.plot_files(data, plot, tran, m_matrix, u_matrix, u_opt, ftn_output, comm)

#wan90.w90_helper_types.transport(data, plot, tran, u_matrix, u_opt, ftn_output, comm)

wan90.w90_helper_types.print_times(data, ftn_output)
