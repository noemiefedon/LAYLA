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
from src.CLA.lp_functions_1 import calc_lampam_from_delta_lp_matrix

from src.LAYLA_V02.parameters import Parameters
from src.LAYLA_V02.constraints import Constraints
from src.LAYLA_V02.divide_laminate_sym import divide_laminate_sym
from src.LAYLA_V02.divide_laminate_asym import divide_laminate_asym
from src.LAYLA_V02.ply_order import calc_ply_order, calc_levels
from src.LAYLA_V02.moment_of_areas import calc_mom_of_areas
from src.LAYLA_V02.targets import Targets

# test 1
constraints1 = Constraints(
    sym=False,
    set_of_angles=np.array([0, 45, 90, -45]))
parameters1 = Parameters(
    constraints=constraints1,
    group_size_min=1,
    group_size_max=100000)
stack1 = np.array([0, 0, 0, 0, 0, 0, 0, 0], int)
targets1 = Targets(n_plies=stack1.size, lampam=np.array(()))
ply_order1 = calc_ply_order(constraints1, targets1)
n_plies_in_groups1, pos_first_ply_groups1, n_groups1 = \
divide_laminate_asym(parameters1, targets1, step=0)
levels_in_groups1 = calc_levels(ply_order1, n_plies_in_groups1, n_groups1)
mom_areas1, cummul_mom_areas1, group_mom_areas1 = calc_mom_of_areas(
    constraints1, targets1, ply_order1, n_plies_in_groups1)
delta_lampams1 = np.empty((stack1.size,
                           constraints1.n_set_of_angles, 12), float)
for ind in range(delta_lampams1.shape[0]):
    delta_lampams1[ind, :, 0:4] = mom_areas1[ind, 0] * constraints1.cos_sin
    delta_lampams1[ind, :, 4:8] = mom_areas1[ind, 1] * constraints1.cos_sin
    delta_lampams1[ind, :, 8:12] = mom_areas1[ind, 2] * constraints1.cos_sin
lampam1 = np.array([1, 1, 0, 0,
                    0, 0, 0, 0,
                    1, 1, 0, 0], float)

# test 2
constraints2 = Constraints(
    sym=True,
    set_of_angles=np.array([0, 45, 90, -45]))
parameters2 = Parameters(
    constraints=constraints2,
    group_size_min=1,
    group_size_max=100000)
stack2 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0], int)
targets2 = Targets(n_plies=stack2.size, lampam=np.array(()))
ply_order2 = calc_ply_order(constraints2, targets2)
n_plies_in_groups2, pos_first_ply_groups2, n_groups2 = \
divide_laminate_sym(parameters2, targets2, step=0)
levels_in_groups2 = calc_levels(ply_order2, n_plies_in_groups2, n_groups2)
mom_areas2, cummul_mom_areas2, group_mom_areas2 = calc_mom_of_areas(
    constraints2, targets2, ply_order2, n_plies_in_groups2)
delta_lampams2 = np.empty((stack2.size // 2 + stack2.size % 2,
                           constraints2.n_set_of_angles, 12), float)
for ind in range(delta_lampams2.shape[0]):
    delta_lampams2[ind, :, 0:4] = mom_areas2[ind, 0] * constraints2.cos_sin
    delta_lampams2[ind, :, 4:8] = 0
    delta_lampams2[ind, :, 8:12] = mom_areas2[ind, 2] * constraints2.cos_sin
lampam2 = np.array([1, 1, 0, 0,
                    0, 0, 0, 0,
                    1, 1, 0, 0], float)

@pytest.mark.parametrize(
    """stack, constraints, delta_lampams, expect""", [
        (stack1, constraints1, delta_lampams1, lampam1),
        (stack2, constraints2, delta_lampams2, lampam2)
        ])

def test_calc_lampam_from_delta_lp_matrix(
        stack, constraints, delta_lampams, expect):
    "tests function calc_lampam_from_delta_lp_matrix"
    output = calc_lampam_from_delta_lp_matrix(
        stack, constraints, delta_lampams)
    assert np.allclose(output, expect)
