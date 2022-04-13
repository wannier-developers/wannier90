# Need PYTHONPATH to contain wrapper directory
# Also need LD_LIBRARY_PATH to include it for libwannier_serial.so
import wan90

data = wan90.w90_helper_types.lib_global_type()
pw90 = wan90.w90_lib_all.lib_postw90_type()
plot = wan90.w90_helper_types.lib_plot_type()
tran = wan90.w90_helper_types.lib_transport_type()
comm = wan90.w90_comms.w90comm_type()
wan90.w90_lib_all.read_all_input(data, pw90, plot, tran, "Fe", 6, comm)

if not pw90.effective_model :
    if not data.kmesh_info.explicit_nnkpts :
        wan90.w90_helper_types.create_kmesh(data, 6, comm)

import numpy

m_matrix = numpy.zeros((data.num_wann, data.num_wann, data.kmesh_info.nntot, data.num_kpts), dtype=numpy.cdouble, order='F')
u_matrix = numpy.zeros((data.num_wann, data.num_wann, data.num_kpts), dtype=numpy.cdouble, order='F')
#a_matrix = numpy.zeros((data.num_bands, data.num_wann, data.num_kpts), dtype=numpy.cdouble, order='F')

# have_disentangled is in the checkpoint file so unset!!!
if wann90.have_disentangled :
    u_opt = numpy.zeros((data.num_bands, data.num_wann, data.num_kpts), dtype=numpy.cdouble, order='F')

wan90.w90_lib_all.read_checkpoint(data, pw90, m_matrix, u_matrix, u_opt, 6, comm)

v_matrix = numpy.empty((data.num_bands, data.num_wann, data.num_kpts), dtype=numpy.cdouble, order='F')
wan90.w90_lib_all.calc_v_matrix(data, u_matrix, u_opt, v_matrix)

# should check pw90.dos.index for 'dos_plot'
if pw90.calculation.dos :
    wan90.w90_lib_all.calc_dos(data, pw90, u_matrix, v_matrix, 6, comm)
