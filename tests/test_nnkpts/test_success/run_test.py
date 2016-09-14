#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:  Dominik Gresch <greschd@gmx.ch>
# Date:    14.09.2016 09:54:19 CEST
# File:    run_test.py

import os
import re
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

def get_block(lines, tag):
    start, end = None, None
    for i, l in enumerate(lines):
        if re.match('[\s]*begin[\s]+' + tag + '[\s]*', l):
            start = i
        elif re.match('[\s]*end[\s]+' + tag + '[\s]*', l):
            end = i
    if start is None or end is None:
        return []
    else:
        return lines[start + 1:end]

def parse_block(lines, **kwargs):
    return [list(np.fromstring(l, sep=' ', **kwargs)) for l in lines]
    
def compare_sorted(l1, l2):
    assert list(sorted(l1)) == list(sorted(l2))
    
def check(path):
    with open(os.path.join(path, 'wannier.win'), 'r') as f:
        win = f.readlines()
    with open(os.path.join(path, 'wannier.nnkp'), 'r') as f:
        out = f.readlines()
        
    in_nnkpts = parse_block(get_block(win, 'nnkpts'), dtype=int)
    out_nnkpts = parse_block(get_block(out, 'nnkpts')[1:], dtype=int) # skip nntot
    compare_sorted(in_nnkpts, out_nnkpts)

    in_kpts = parse_block(get_block(win, 'kpoints'))
    out_kpts = parse_block(get_block(out, 'kpoints')[1:])
    assert np.isclose(np.array(in_kpts), np.array(out_kpts)).all()
    print(path, 'passed')

def run_test(path):
    clean(path)
    subprocess.call([W90, '-pp'], cwd=path)
    check(path)
    

if __name__ == '__main__':
    test_cases = (t for t in os.listdir('.') if t.startswith('test'))
    for t in test_cases:
        run_test(t)
