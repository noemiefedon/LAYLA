#  - * -  coding: utf - 8  - * -
"""
This module test the functions for the 10% rule.
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import pytest
import numpy as np

sys.path.append(r'C:\LAYLA')
from src.LAYLA_V02.constraints import Constraints
from src.guidelines.ten_percent_rule import is_ten_percent_rule

@pytest.mark.parametrize(
    """constraints, stack, ply_queue, n_plies_per_angle, equality_45_135,
equality_0_90, expect""", [
        (Constraints(rule_10_percent=True, percent_0=50),
         np.array([0, 45, 90], int), [], None, False, False, False),
        (Constraints(rule_10_percent=True, percent_0=50),
         np.array([0, 45, 90], int), None, None, False, False, False),
        (Constraints(rule_10_percent=True, percent_0=50),
         np.array([0, 666, 666], int), [0, 45], None, False, False, True),
        (Constraints(rule_10_percent=True, percent_0=50),
         None, None, np.array([3, 0, 3, 0]), False, False, False),
        (Constraints(rule_10_percent=True, percent_0=50),
         None, None, np.array([0, 3, 0, 3]), False, False, True)
        ])

def test_is_ten_percent_rule(
        constraints, stack, ply_queue, n_plies_per_angle, equality_45_135,
        equality_0_90, expect):
    output = is_ten_percent_rule(
        constraints, stack, ply_queue, n_plies_per_angle, equality_45_135,
        equality_0_90)
    assert output == expect
