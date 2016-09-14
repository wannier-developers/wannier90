#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:  Dominik Gresch <greschd@gmx.ch>
# Date:    14.09.2016 09:54:19 CEST
# File:    run_test.py

import os
import re
import subprocess
import contextlib

import numpy as np

W90 = os.path.abspath('../../../wannier90.x')

TEST_MAPPING = {
    'test1': 'Input parameter nnkpts_block is allowed only if postproc_setup = .true.',
    'test2': 'The number of rows in nnkpts must be a multiple of num_kpts',
}

def clean(path):
    to_clean = ['wannier.werr']
    for f in to_clean:
        try:
            os.remove(os.path.join(path, f))
        except OSError:
            pass

def check(path, string):
    with open(os.path.join(path, 'wannier.werr'), 'r') as f:
        res = f.read()
    assert string in res

def run_test(path):
    with open(os.devnull, 'w') as null:
        subprocess.call(W90, cwd=path, stderr=null, stdout=null)
    check(path, TEST_MAPPING[path])

if __name__ == '__main__':
    test_cases = (t for t in os.listdir('.') if t.startswith('test'))
    for t in test_cases:
        try:
            run_test(t)
            print(t, 'passed')
        except (AssertionError, OSError):
            print(t, 'failed')
