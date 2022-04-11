# Need mpi library etc in paths as well as python requirements of serial example
# mpiexec -n 4 python mpi-example.py
from mpi4py import MPI
fcomm = MPI.COMM_WORLD.py2f()

# Maybe should have a common name...
import wan90mpi as wan90

data = wan90.w90_helper_types.lib_global_type()
plot = wan90.w90_helper_types.lib_plot_type()
tran = wan90.w90_helper_types.lib_transport_type()
comm = wan90.w90_comms.w90comm_type()

# probably should just do
#comm.comm = MPI.COMM_WORLD.py2f()
comm.comm = fcomm

wan90.w90_helper_types.input_reader(data, plot, tran, "diamond", 6, comm)

wan90.w90_helper_types.create_kmesh(data, 6, comm)

import numpy

m_matrix = numpy.zeros((data.num_wann, data.num_wann, data.kmesh_info.nntot, data.num_kpts), dtype=numpy.cdouble, order='F')
u_matrix = numpy.zeros((data.num_wann, data.num_wann, data.num_kpts), dtype=numpy.cdouble, order='F')
a_matrix = numpy.zeros((data.num_bands, data.num_wann, data.num_kpts), dtype=numpy.cdouble, order='F')

u_opt = numpy.zeros((data.num_bands, data.num_wann, data.num_kpts), dtype=numpy.cdouble, order='F')

#m_matrix = numpy.asfortranarray(m_matrix)
#m_matrix.flags.f_contiguous

wan90.w90_helper_types.overlaps(data, a_matrix, m_matrix, u_matrix, 6, comm)
wan90.w90_helper_types.wannierise(data, plot, tran, m_matrix, u_matrix, u_opt, 6, comm)

#wan90.w90_helper_types.checkpoint(data, "postwann", m_matrix, u_matrix, u_opt, 6, comm)

#wan90.w90_helper_types.transport(helper, plot, tran, u_matrix, u_opt, 6, comm)

#print (data.num_wann)

#wan90.w90_helper_types.plot_files(data, plot, tran, m_matrix, u_matrix, u_opt, 6, comm)

#wan90.w90_helper_types.print_times(data)
