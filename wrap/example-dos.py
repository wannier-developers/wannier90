# Need PYTHONPATH to contain wrapper directory
# Also need LD_LIBRARY_PATH to include it for libwannier_serial.so
import wan90

ftn_output = wan90.w90_helper_types.get_fortran_stdout()
ftn_error = wan90.w90_helper_types.get_fortran_stderr()

data = wan90.w90_helper_types.lib_global_type()
pw90 = wan90.w90_lib_all.lib_postw90_type()
w90data = wan90.w90_helper_types.lib_w90_type()
comm = wan90.w90_comms.w90_comm_type()
status = wan90.w90_lib_all.read_all_input(data, w90data, pw90, "Fe", ftn_output, ftn_error, comm)

if not pw90.effective_model :
    if not data.kmesh_info.explicit_nnkpts :
        status = wan90.w90_helper_types.create_kmesh(data, ftn_output, ftn_error, comm)

import numpy

# create dummy distribution, all done on proc 0
kpts = numpy.zeros(data.num_kpts, dtype=numpy.int32)
wan90.w90_helper_types.set_kpoint_distribution(data, kpts)

#m_matrix = numpy.zeros((data.num_wann, data.num_wann, data.kmesh_info.nntot, data.num_kpts), dtype=numpy.cdouble, order='F')
u_matrix = numpy.zeros((data.num_wann, data.num_wann, data.num_kpts), dtype=numpy.cdouble, order='F')
#wan90.w90_helper_types.set_m_matrix(w90data, m_matrix)
wan90.w90_helper_types.set_u_matrix(data, u_matrix)
#a_matrix = numpy.zeros((data.num_bands, data.num_wann, data.num_kpts), dtype=numpy.cdouble, order='F')

# have_disentangled is in the checkpoint file so unset!!!
if wann90.have_disentangled :
    u_opt = numpy.zeros((data.num_bands, data.num_wann, data.num_kpts), dtype=numpy.cdouble, order='F')
    wan90.w90_helper_types.set_u_opt(data, u_opt)

status = wan90.w90_lib_all.read_checkpoint(data, pw90, ftn_output, ftn_error, comm)

v_matrix = numpy.empty((data.num_bands, data.num_wann, data.num_kpts), dtype=numpy.cdouble, order='F')
wan90.w90_lib_all.calc_v_matrix(data, pw90, v_matrix)

# should check pw90.dos.index for 'dos_plot'
if pw90.calculation.dos :
    status = wan90.w90_lib_all.calc_dos(data, pw90, ftn_output, ftn_error, comm)
