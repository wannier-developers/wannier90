#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:  Dominik Gresch <greschd@gmx.ch>
# Date:    14.09.2016 09:54:19 CEST
# File:    run_test.py

import os
import subprocess

import numpy as np

W90 = os.path.abspath('../../../wannier90.x')

def clean(path):
    to_clean = ['wannier.wout', 'wannier.nnkp']
    for f in to_clean:
        try:
            os.remove(os.path.join(path, f))
        except OSError:
            pass

#~ def find_block(tag):


#~ def check(path):
    

def run_test(path):
    clean(path)
    subprocess.call([W90, '-pp'], cwd=path)
    

if __name__ == '__main__':
    test_cases = (t for t in os.listdir('.') if t.startswith('test'))
    for t in test_cases:
        run_test(t)
