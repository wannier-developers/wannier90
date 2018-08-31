#!/usr/bin/env python
from __future__ import print_function
import numpy
import sys

results_fname = 'results.dat'
ref_fname = 'ref/results_ref.dat'

err_threshold = 1.e-4 # Do not decrease too much, due to the precision in output

results = numpy.loadtxt(results_fname)
ref = numpy.loadtxt(ref_fname)

abs_err = numpy.abs(results - ref)
err_positions = abs_err > err_threshold

if err_positions.any():
    print(">>> Error, some values are different from the reference")
    print(">>> Error positions:")
    print(err_positions)
    print(">>> RESULT")
    with open(results_fname) as f:
        print(f.read())
    print(">>> REFERENCE")
    with open(ref_fname) as f:
        print(f.read())
    sys.exit(1)
else:
    print("Centres and spreads are OK")
