# -*- coding: utf-8 -*-
"""
Tests for functions calculating lamination parameters
"""
__version__ = '2.0'
__author__ = 'Noemie Fedon'

import sys
import pytest
import numpy as np

sys.path.append(r'C:\LAYLA')
from src.CLA.lp_functions_2 import calc_lampamA, calc_lampamD
from src.LAYLA_V02.constraints import Constraints

# test 1
constraints1 = Constraints(
    sym=False,
    set_of_angles=np.array([0, 45, 90, -45]))
stack1 = np.array([0, 0, 0, 0, 0, 0, 0, 0], int)
lampam1A = np.array([1, 1, 0, 0], float)
lampam1D = np.array([1, 1, 0, 0], float)

# test 2
constraints2 = Constraints(
    sym=True,
    set_of_angles=np.array([0, 45, 90, -45]))
stack2 = np.array([90, 0, 0, 0, 90], int)
lampam2A = np.array([0.2, 1, 0, 0], float)
lampam2D = np.array([-0.568, 1, 0, 0], float)

@pytest.mark.parametrize(
    """stack, constraints, expect""", [
        (stack1, constraints1, lampam1A),
        (stack2, constraints2, lampam2A)
        ])

def test_calc_lampamA(stack, constraints, expect):
    "tests function calc_lampamA"
    output = calc_lampamA(stack, constraints)
    assert np.allclose(output, expect)

@pytest.mark.parametrize(
    """stack, constraints, expect""", [
        (stack1, constraints1, lampam1D),
        (stack2, constraints2, lampam2D)
        ])

def test_calc_lampamD(stack, constraints, expect):
    "tests function calc_lampamD"
    output = calc_lampamD(stack, constraints)
    assert np.allclose(output, expect)
