#! /usr/bin/env python3
#
# This file is distributed as part of the Wannier90 code and
# under the terms of the GNU General Public License. See the
# file `LICENSE' in the root directory of the Wannier90
# distribution, or http://www.gnu.org/copyleft/gpl.txt
#
# The webpage of the Wannier90 code is www.wannier.org
#
# The Wannier90 code is hosted on GitHub:
#
# https://github.com/wannier-developers/wannier90
#
# Python3 script to find the indexes of a coarse mesh into a finer mesh
# provided they are commensurate.
#
# Written by Antimo Marrazzo (EPFL)
# Last update September 13th, 2016
#
import sys
import numpy as np
import os


def prepare_mesh(coarse_grid, nscf_output_file):
    import subprocess
    with open(nscf_output_file, 'r') as f:
        read_data = f.readlines()
    f.close()
    read_kpts = False
    k_fine_list = []
    for line in read_data:
        if "number of k points=" in line:
            numk_line = line.strip('\n').split()
            num_kpoints = int(numk_line[4])
            print(
                'Number of kpoints provided to Yambo through a NSCF calculation',
                num_kpoints)
        if read_kpts == True:
            kline = line.strip('\n').split()
            if 'wk' in kline:
                a = kline[4:6]
                b = kline[6].split(')')[0]
                k_vec = [float(a[0]), float(a[1]), float(b)]
                k_fine_list.append(k_vec)
            else:
                read_kpts = False
        if "cryst. coord." in line and 'site' not in line:
            read_kpts = True
    coarse_text = [str(i) + ' ' for i in coarse_grid]
    k_coarse_mesh = subprocess.check_output([
        os.path.join(os.path.dirname(os.path.realpath(__file__)), 'kmesh.pl'),
        coarse_text[0], coarse_text[1], coarse_text[2], 'wan'
    ],
                                            universal_newlines=True)
    k_coarse_mesh = k_coarse_mesh.split('\n')
    k_coarse_list = []
    for i in range(coarse_grid[0] * coarse_grid[1] * coarse_grid[2]):
        line = k_coarse_mesh[i].split()
        k_coarse_list.append([float(j) for j in line])
    return (k_fine_list, k_coarse_list)


def indexes_list(fine_mesh, coarse_mesh):
    import numpy as np
    opt = np.array([0, 1, -1])
    k_list = []
    for i in coarse_mesh:
        count = 1
        for j in fine_mesh:
            q = i - j
            q = np.around(q, decimals=5)
            if (q[0] in opt and q[1] in opt and q[2] in opt):
                k_list.append(count)
            count = count + 1
    return k_list


print('####################')
print('####Mesh mapper#####')
print('####################')
coarse_grid = sys.argv[1:4]
coarse_grid = [int(i) for i in coarse_grid]
print('Input coarse mesh:', coarse_grid)
nscf_output_file = sys.argv[4]
print('Path of the QE NSCF output file', nscf_output_file)

(k_fine_list, k_coarse_list) = prepare_mesh(coarse_grid, nscf_output_file)

ind_list = indexes_list(np.array(k_fine_list), np.array(k_coarse_list))
print('List of k-indexes to pass to Yambo', ind_list)
for i in ind_list:
    print(str(i) + '|' + str(i) + '|' + 'first band' + '|' + 'last band' + '|')
