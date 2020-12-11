# -*- coding: utf-8 -*-
"""
script to generate symmetric stacking sequence populations
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import itertools
import pandas as pd
import numpy as np

sys.path.append(r'C:\LAYLA')

from src.LAYLA_V02.constraints import Constraints
from src.LAYLA_V02.materials import Material
from src.CLA.lampam_functions import calc_lampam
from src.CLA.ABD import A_from_lampam, B_from_lampam, D_from_lampam
from src.LAYLA_V02.save_set_up import save_constraints_LAYLA
from src.LAYLA_V02.save_set_up import save_materials
from src.divers.excel import append_df_to_excel, autofit_column_widths, delete_file
from src.guidelines.contiguity import is_contig
from src.guidelines.disorientation import is_diso_ss
from src.guidelines.balance import is_balanced
from src.guidelines.dam_tol import is_dam_tol
from src.guidelines.ten_percent_rule import is_ten_percent_rule

#==============================================================================
# output set up
#==============================================================================
# number of plies in a laminate
n_plies = 9
# file name
filename = 'trad-sym-' + str(n_plies) + '-plies.xlsx'
#==============================================================================
# Material properties
#==============================================================================
# Elastic modulus in the fibre direction (Pa)
E11 = 130e9
# Elastic modulus in the transverse direction (Pa)
E22 = 9e9
# Poisson's ratio relating transverse deformation and axial loading (-)
nu12 = 0.3
# In-plane shear modulus (Pa)
G12 = 4e9
mat_prop = Material(E11 = E11, E22 = E22, G12 = G12, nu12 = nu12)
#==============================================================================
# Design guidelines
#==============================================================================
# set of admissible fibre orientations
set_of_angles = np.array([-45, 0, 45, 90], dtype=int)
#set_of_angles = np.array([-45, 0, 45, 90, +30, -30, +60, -60], dtype=int)

# symmetry
sym = True

# balance and in-plane orthotropy requirements
bal = True
ipo = True

# out-of-plane orthotropy requirements
oopo = False

# damage tolerance
dam_tol = True
# rule 1: one outer ply at + or -45 deg at laminate surfaces
# rule 2: [+45, -45] or [-45, +45] plies at laminate surfaces
dam_tol_rule = 2

# 10% rule
rule_10_percent = True
percent_0 = 10 # percentage used in the 10% rule for 0 deg plies
percent_45 =  0 # percentage used in the 10% rule for +45 deg plies
percent_90 = 10 # percentage used in the 10% rule for 90 deg plies
percent_135 = 0 # percentage used in the 10% rule for -45 deg plies
percent_45_135 = 10 # percentage used in the 10% rule for +-45 deg plies

# disorientation
diso = True

# Upper bound of the variation of fibre orientation between two
# contiguous plies if the disorientation constraint is active
delta_angle = 45

# contiguity
contig = True

n_contig = 4
# No more that constraints.n_contig plies with same fibre orientation should be
# next to each other if the contiguity constraint is active. The value taken
# can only be 2, 3, 4 or 5, otherwise test functions should be modified

constraints = Constraints(
    sym=sym,
    bal=bal,
    ipo=ipo,
    oopo=oopo,
    dam_tol=dam_tol,
    dam_tol_rule=dam_tol_rule,
    rule_10_percent=rule_10_percent,
    percent_0=percent_0,
    percent_45=percent_45,
    percent_90=percent_90,
    percent_135=percent_135,
    percent_45_135=percent_45_135,
    diso=diso,
    contig=contig,
    n_contig=n_contig,
    delta_angle=delta_angle,
    set_of_angles=set_of_angles)
#==============================================================================
#
#==============================================================================

delete_file(filename)
popu = pd.DataFrame()

for ind, stack in enumerate(itertools.product(
        constraints.set_of_angles, repeat=n_plies)):

    if not ind % 100:
        print('ind', ind)

    popu.loc[ind, 'ply_count'] = 2*n_plies

    ss_flatten = np.array(stack, dtype=str)
    ss_flatten = ' '.join(ss_flatten)
    popu.loc[ind, 'half_stack'] = ss_flatten

    stack = np.hstack((stack, np.flip(stack)))

    lampam = calc_lampam(stack, constraints)

    popu.loc[ind, 'lampam1'] = lampam[0]
    popu.loc[ind, 'lampam2'] = lampam[1]
    popu.loc[ind, 'lampam3'] = lampam[2]
    popu.loc[ind, 'lampam4'] = lampam[3]
    popu.loc[ind, 'lampam5'] = lampam[4]
    popu.loc[ind, 'lampam6'] = lampam[5]
    popu.loc[ind, 'lampam7'] = lampam[6]
    popu.loc[ind, 'lampam8'] = lampam[7]
    popu.loc[ind, 'lampam9'] = lampam[8]
    popu.loc[ind, 'lampam10'] = lampam[9]
    popu.loc[ind, 'lampam11'] = lampam[10]
    popu.loc[ind, 'lampam12'] = lampam[11]

    A = A_from_lampam(lampam, mat_prop)
    B = B_from_lampam(lampam, mat_prop, constraints.sym)
    D = D_from_lampam(lampam, mat_prop)

#    filter_ABD(A=A, B=B, D=D, sym=constraints.sym, ply_t=mat_prop.ply_t)

    popu.loc[ind, 'A11'] = A[0, 0]
    popu.loc[ind, 'A22'] = A[1, 1]
    popu.loc[ind, 'A12'] = A[0, 1]
    popu.loc[ind, 'A66'] = A[2, 2]
    popu.loc[ind, 'A16'] = A[0, 2]
    popu.loc[ind, 'A26'] = A[1, 2]

    popu.loc[ind, 'B11'] = B[0, 0]
    popu.loc[ind, 'B22'] = B[1, 1]
    popu.loc[ind, 'B12'] = B[0, 1]
    popu.loc[ind, 'B66'] = B[2, 2]
    popu.loc[ind, 'B16'] = B[0, 2]
    popu.loc[ind, 'B26'] = B[1, 2]

    popu.loc[ind, 'D11'] = D[0, 0]
    popu.loc[ind, 'D22'] = D[1, 1]
    popu.loc[ind, 'D12'] = D[0, 1]
    popu.loc[ind, 'D66'] = D[2, 2]
    popu.loc[ind, 'D16'] = D[0, 2]
    popu.loc[ind, 'D26'] = D[1, 2]

#    popu.loc[ind, 'N0'] = N0
#    popu.loc[ind, 'N45'] = N45
#    popu.loc[ind, 'N90'] = N90
#    popu.loc[ind, 'N-45'] = N135

    popu.loc[ind, 'dam_tol'] = is_dam_tol(stack, constraints)

    popu.loc[ind, 'diso'] = is_diso_ss(
        stack,
        delta_angle=constraints.delta_angle,
        dam_tol=constraints.dam_tol,
        dam_tol_rule=constraints.dam_tol_rule)

    if (abs(lampam[2:4]) > 1e-10).any():
        popu.loc[ind, 'ipo'] = False
    else:
        popu.loc[ind, 'ipo'] = True

    if is_balanced(stack, constraints):
        popu.loc[ind, 'balance'] = True
    else:
        popu.loc[ind, 'balance'] = False

    if is_contig(stack, constraints.n_contig):
        popu.loc[ind, 'contig'] = True
    else:
        popu.loc[ind, 'contig'] = False

    if is_ten_percent_rule(constraints, stack=stack):
        popu.loc[ind, 'rule_10_percent'] = True
    else:
        popu.loc[ind, 'rule_10_percent'] = False

save_constraints_LAYLA(filename, constraints)
save_materials(filename, mat_prop)
append_df_to_excel(filename, popu, 'stacks', index=True, header=True)
autofit_column_widths(filename)
