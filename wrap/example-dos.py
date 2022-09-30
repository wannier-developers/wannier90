# Need PYTHONPATH to contain wrapper directory
# Also need LD_LIBRARY_PATH to include it for libwannier_serial.so
import wan90

ftn_output = wan90.w90_library.w90_get_fortran_stdout()
ftn_error = wan90.w90_library.w90_get_fortran_stderr()

data = wan90.w90_library.lib_common_type()
pw90 = wan90.w90_lib_all.lib_postw90_type()
#comm = wan90.w90_comms.w90_comm_type()
# if you know the eigvals then allocate here and set, otherwise let the reader get them
#eigval = numpy.zeros((data.num_bands, data.num_kpts), dtype=numpy.double, order='F')
#status = wan90.w90_lib_all.read_all_input_has_eigs(data, pw90, eigval, ftn_output, ftn_error)
#status = wan90.w90_lib_all.read_all_input_and_eigs(data, pw90, "Fe", ftn_output, ftn_error)
status = wan90.w90_lib_all.read_pw90_input(data, pw90, "Fe", ftn_output, ftn_error)

if not pw90.effective_model :
    if not data.kmesh_info.explicit_nnkpts :
        status = wan90.w90_library.w90_create_kmesh(data, ftn_output, ftn_error)

import numpy

# create dummy distribution, all done on proc 0
kpts = numpy.zeros(data.num_kpts, dtype=numpy.int32)
wan90.w90_library_extra.set_kpoint_distribution(data, kpts, ftn_output, ftn_error)

m_matrix = numpy.zeros((data.num_bands, data.num_bands, data.kmesh_info.nntot, data.num_kpts), dtype=numpy.cdouble, order='F')
u_opt = numpy.zeros((data.num_bands, data.num_wann, data.num_kpts), dtype=numpy.cdouble, order='F')
wan90.w90_library.w90_set_m_local(data, m_matrix)
wan90.w90_library.w90_set_u_opt(data, u_opt)
u_matrix = numpy.zeros((data.num_wann, data.num_wann, data.num_kpts), dtype=numpy.cdouble, order='F')
wan90.w90_library.w90_set_u_matrix(data, u_matrix)

# have_disentangled is in the checkpoint file so unset!!!
#if data.have_disentangled :
if data.num_bands > data.num_wann :
    #u_opt = numpy.zeros((data.num_bands, data.num_wann, data.num_kpts), dtype=numpy.cdouble, order='F')
    #wan90.w90_library.w90_set_u_opt(data, u_opt)

status = wan90.w90_lib_all.read_checkpoint(data, pw90, ftn_output, ftn_error)
status = wan90.w90_lib_all.pw_setup(data, pw90, ftn_output, ftn_error)

v_matrix = numpy.empty((data.num_bands, data.num_wann, data.num_kpts), dtype=numpy.cdouble, order='F')
wan90.w90_lib_all.calc_v_matrix(data, pw90, v_matrix)

# should check pw90.dos.index for 'dos_plot'
if pw90.calculation.dos :
    status = wan90.w90_lib_all.calc_dos(data, pw90, ftn_output, ftn_error)
