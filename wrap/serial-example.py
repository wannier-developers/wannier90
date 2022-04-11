# Need PYTHONPATH to contain wrapper directory
# Also need LD_LIBRARY_PATH to include it for libwannier_serial.so
import wan90

data = wan90.w90_helper_types.lib_global_type()
plot = wan90.w90_helper_types.lib_plot_type()
tran = wan90.w90_helper_types.lib_transport_type()
comm = wan90.w90_comms.w90comm_type()
wan90.w90_helper_types.input_reader(data, plot, tran, "diamond", 6, comm)

if not data.kmesh_info.explicit_nnkpts :
    wan90.w90_helper_types.create_kmesh(data, 6, comm)

import numpy

m_matrix = numpy.zeros((data.num_wann, data.num_wann, data.kmesh_info.nntot, data.num_kpts), dtype=numpy.cdouble, order='F')
u_matrix = numpy.zeros((data.num_wann, data.num_wann, data.num_kpts), dtype=numpy.cdouble, order='F')
a_matrix = numpy.zeros((data.num_bands, data.num_wann, data.num_kpts), dtype=numpy.cdouble, order='F')

u_opt = numpy.zeros((data.num_bands, data.num_wann, data.num_kpts), dtype=numpy.cdouble, order='F')

#m_matrix.flags.f_contiguous should be true

wan90.w90_helper_types.overlaps(data, a_matrix, m_matrix, u_matrix, 6, comm)
wan90.w90_helper_types.wannierise(data, plot, tran, m_matrix, u_matrix, u_opt, 6, comm)

#wan90.w90_helper_types.checkpoint(data, "postwann", m_matrix, u_matrix, u_opt, 6, comm)

#wan90.w90_helper_types.transport(helper, plot, tran, u_matrix, u_opt, 6, comm)

print (data.num_wann)

wan90.w90_helper_types.plot_files(data, plot, tran, m_matrix, u_matrix, u_opt, 6, comm)

wan90.w90_helper_types.print_times(data, 6)
