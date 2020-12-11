# -*- coding: utf-8 -*-
"""
Function to create a population of 200 symmetric laminates of
40 plies satisfying the design and manufacturing guidelines C0
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
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

n_plies_in_panels = 40
n_pop = 57
filename = 'pop_sym_C0_40plies.xlsx'
delete_file(filename)
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
mat = Material(E11=E11, E22=E22, G12=G12, nu12=nu12)
#==============================================================================
# Design guidelines
#==============================================================================
### Set of design and manufacturing constraints:
constraints_set = 'C0'

# Set of admissible fibre orientations
set_of_angles = np.array([-45, 0, 45, 90], dtype=int)
#set_of_angles= np.array([-45, 0, 45, 90, +30, -30, +60, -60], dtype=int)

# symmetry
sym = True

# balance and in-plane orthotropy requirements
if constraints_set == 'C0':
    bal = False
    ipo = False
else:
    bal = True
    ipo = True

# out-of-plane orthotropy requirements
oopo = False

# damage tolerance
dam_tol = False

# 10% rule
if constraints_set == 'C0':
    rule_10_percent = False
else:
    rule_10_percent = True

percent_0 = 10 # percentage used in the 10% rule for 0 deg plies
percent_45 = 10 # percentage used in the 10% rule for +45 deg plies
percent_90 = 10 # percentage used in the 10% rule for 90 deg plies
percent_135 =10 # percentage used in the 10% rule for -45 deg plies

# disorientation
if constraints_set == 'C0':
    diso = False
else:
    diso = True

# Upper bound of the variation of fibre orientation between two
# contiguous plies if the disorientation constraint is active
delta_angle = 45
# delta_angle must be at least 45 to ensure convergence when the 10% is active

# contiguity
if constraints_set == 'C0':
    contig = False
else:
    contig = True

n_contig = 5
# No more that constraints.n_contig plies with same fibre orientation should be
# next to each other if the contiguity constraint is active. The value taken
# can only be 2, 3, 4 or 5, otherwise test functions should be modified

constraints = Constraints(
    sym=sym,
    bal=bal,
    ipo=ipo,
    oopo=oopo,
    dam_tol=dam_tol,
    rule_10_percent=rule_10_percent,
    percent_0=percent_0,
    percent_45=percent_45,
    percent_90=percent_90,
    percent_135=percent_135,
    diso=diso,
    contig=contig,
    n_contig=n_contig,
    delta_angle=delta_angle,
    set_of_angles=set_of_angles)
#print(constraints)

Pop = pd.DataFrame()

# Initialisation
print(' Creating Pop ... ')
ipop = 1

# Loop until Pop is complete
while ipop < n_pop + 1:

    if ipop == 1:
        ss=0*np.ones(20,)
    if ipop == 2:
        ss=90*np.ones(20,)
    if ipop == 3:
        ss=45*np.ones(20,)
    if ipop == 4:
        ss=-45*np.ones(20,)
    if ipop == 5:
        ss = np.array([0 ,0 ,0 ,0 , 0, 0 ,0 ,0 ,0 , 0, -45, -45, -45, 45, -45, -45, 45, -45, -45, -45])
    if ipop == 6:
        ss = np.array([0 ,0 ,0 ,0 , 0, 0 ,0 ,0 ,0 , 0, -45, -45, -45, 45, -45, -45, 45, -45, -45, -45])
    if ipop == 7:
        ss = np.array([45 ,45 ,45 ,0 , 0, 0 ,0 ,0 ,0 , 0, -45, -45, -45, 45, -45, -45, 45, -45, -45, -45])
    if ipop == 8:
        ss = np.array([0 ,0 ,0 ,0 , 0, 0 ,0 ,90 ,90 , 90, -45, -45, -45, 45, -45, -45, 45, -45, 90, 90])
    if ipop == 9:
        ss = np.array([0 ,0 ,0 ,0 , 0, 90 ,90 ,90 ,90 , 90, 45, -45, 45, 45, -45, 45 -45, 45, -45, -45, -45])
    if ipop == 10:
        ss = np.array([90 ,0 ,0 ,90 , 90, 0 ,90 ,0 ,0 , 90, 0 ,90 ,90 ,90 , 0, 90 ,90 ,90 ,0 , 0,])
    if ipop == 11:
        ss = np.array([90 ,45 ,45 ,90 , 90, 0 ,90 ,0 ,0 , 90, 0 ,90 ,90 ,90 , 0, 90 ,90 ,90 ,0 , 0,])
    if ipop == 12:
        ss = np.array([90 ,45 ,-45 ,90 , 90, 0 ,90 ,90 ,90 , 90, 90 ,90 ,90 ,90 , 0, 90 ,90 ,90 ,90 , 90,])
    if ipop == 13:
        ss = np.array([90 ,45 ,90 ,90 , 90, 0 ,90 ,90 ,90 , 90, 90 ,90 ,90 ,90 , 0, 90 ,90 ,90 ,90 , 90,])
    if ipop == 14:
        ss = np.array([90 ,45 ,90 ,90 , 90, 45 ,90 ,90 ,90 , 45, 90 ,-45 ,90 ,-45 , 0, 90 ,90 ,90 ,90 , 90,])
    if ipop == 15:
        ss = np.array([90 ,45 ,90 ,90 , 0, 45 ,90 ,0 ,0 , 45, 90 ,-45 ,90 ,45 , 0, 90 ,0 ,90 ,90 , 90,])
    if ipop == 16:
        ss = np.array([90 ,45 ,90 ,90 , 0, 0 ,90 ,0 ,0 , 45, 90 ,-45 ,90 ,45 , 0, 90 ,0 ,90 ,90 , 90,])
    if ipop == 17:
        ss = np.array([0 ,45 ,90 ,90 , 0, 0 ,90 ,0 ,0 , 45, 90 ,-45 ,90 ,45 , 0, 90 ,0 ,90 ,90 , 0,])
    if ipop == 18:
        ss = np.array([90 ,45 ,90 ,90 , 90, 90 ,90 ,0 ,0 , 45, 90 , 0,0 ,45 , 0, 0 ,0 ,90 ,0 , 45,])
    if ipop == 19:
        ss = np.array([-45 ,45 ,-45 ,-45 , 90, 90 ,90 ,0 ,45 , 45, 90 , 45,0 ,45 , 0, 0 ,0 ,90 ,0 , 45,])
    if ipop == 20:
        ss = np.array([-45 ,45 ,-45 ,-45 , 45, 90 ,45 ,0 ,45 , 45, -45 , 45,0 ,45 , 0, -45 ,45 ,-45 ,-45 , 45,])
    if ipop == 21:
        ss = np.array([-45 ,45 ,-45 ,-45 , 45, 45 ,45 ,0 ,45 , 45, -45 , 45,45 ,45 , 0, -45 ,45 ,-45 ,-45 , 45,])
    if ipop == 22:
        ss = np.array([-45 ,45 ,0 ,-45 , 45, 45 ,45 ,0 ,45 , 45, -45 , 45,45 ,45 , 0, -45 ,45 ,-45 ,0 , 45,])
    if ipop == 23:
        ss = np.array([45 ,45 ,0 ,-45 , 45, 45 ,0 ,0 ,0 , 0, -45 , 45,45 ,45 , 0, -45 ,45 ,-45 ,0 , -45,])
    if ipop == 24:
        ss = np.array([45 ,45 ,90 ,-45 , 45, 45 ,90 ,90 ,90 , 0, -45 , 45,45 ,45 , 0, -45 ,45 ,-45 ,0 , -45,])
    if ipop == 25:
        ss = np.array([45 ,45 ,45 ,45 , 45, 45 ,90 ,90 ,90 , 0, -45 , 45,45 ,45 , 0, 0 ,45 ,0 ,0 , 0,])
    if ipop == 26:
        ss = np.array([45 ,45 ,90 ,45 , 45, 45 ,90 ,90 ,90 , 0, 0 , 45,45 ,45 , 0, 0 ,45 ,0 ,0 , 0,])
    if ipop == 27:
        ss = np.array([45 ,0 ,90 ,-45 , -45, 45 ,90 ,0 ,90 , 0, 90 , 45,45 ,45 , 0, 0 ,45 ,-45 ,90 , 0,])
    if ipop == 28:
        ss = np.array([45 ,0 ,90 ,-45 , 45, 45 ,90 ,0 ,90 , 0, 90 , -45,45 ,45 , 0, 0 ,45 ,-45 ,90 , 0,])
    if ipop == 29:
        ss = np.array([45 ,0 ,0 ,-45 , 45, 45 ,90 ,0 ,90 , 0, 90 , -45,45 ,0 , 0, 0 ,45 ,-45 ,90 , 0,])
    if ipop == 30:
        ss = np.array([45 ,90 ,90 ,-45 , -45, 45 ,90 ,0 ,90 , 0, 90 , -45,45 ,0 , 0, 0 ,45 ,-45 ,90 , 0,])
    if ipop == 31:
        ss = np.array([45 ,45, 45 ,-45 , -45, 45 ,90 ,0 ,90 , 0, 90 , -45,45 ,-45 , 0, 0 ,90, 90 ,90 , 0,])
    if ipop == 32:
        ss = np.array([90 ,45 ,90 ,90 , 90, 90 ,90 ,0 ,0 , 45, 90 , -45,45 ,-45 , 0, 0 ,90, 90 ,90 , 0,])
    if ipop == 33:
        ss = np.array([90 ,45 ,90 ,90 , 90, 90 ,90 ,0 ,0 , 45, 90 ,45 ,0 ,45 , 45, -45 , 45,0 ,45 , 0,])
    if ipop == 34:
        ss = np.array([45 ,45 ,45 ,0 , 0, 0 ,0 ,0 ,0 , 0, -45, -45, -45, 45, 90, 90, 90, -45, -45, -45])
    if ipop == 35:
        ss = np.array([45 ,90 ,45 ,90 , 90, 90 ,0 ,0 ,0 , 0, -45, -45, -45, 45, 90, 90, 90, -45, -45, -45])
    if ipop == 36:
        ss = np.array([0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, ])
    if ipop == 37:
        ss = np.array([0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, 0, 45, -45, 90, ])
    if ipop == 38:
        ss = np.array([0, 45, -45, 45, 0, 45, -45, 90, 0, 45, -45, 90, -45, 45, -45, 90, 0, 45, -45, 90, ])
    if ipop == 39:
        ss = np.array([0, 45, 90, 0, 0, 0, 0, 90, 0, 45, -45, 90, 90, 90, -45, 90, 0, 45, -45, 90, ])
    if ipop == 40:
        ss = np.array([0, 45, 90, 0, 0, 0, 0, 0, 0, 45, -45, 90, 90, 0, -45, 90, 0, 45, -45, 90, ])
    if ipop == 41:
        ss = np.array([0, 45, 90, 0, 0, 0, 0, 0, 90, 45, -45, 0, 0, 0, -45, 0, 0, 45, -45, 0, ])
    if ipop == 42:
        ss = np.array([0, 45, 90, 0, 0, 0, 0, 0, 90, 45, -45, 0, 0, 0, -45, 0, 0, 45, -45, 0, ])
    if ipop == 43:
        ss = np.array([0, 45, 90, 0, 0, 0, 0, 0, 90, 45, 90, 0, 0, 0, -45, 0, 0, 45, -45, 0, ])
    if ipop == 44:
        ss = np.array([0 ,45, 45 ,0 , 0, 0, 0, -45, 0 , 0, -45, -45, -45, 45, -45, -45, 45, -45, -45, -45])
    if ipop == 45:
        ss = np.array([0 ,-45, -45 ,0 , 0, 0, 0, -45, 0 , 0, -45, -45, 45, 45, 45, 45, 45, -45, 45, 45])
    if ipop == 46:
        ss = np.array([90 ,90 ,90 ,0 , 0, 0 ,0 ,0 ,0 , 90, 45, 45, 45, 45, -45, -45, 45, -45, -45, -45])
    if ipop == 47:
        ss = np.array([0 ,0 ,0 ,0 , 90, 0 ,0 ,90 ,90 , 90, 45, 45, 45, 45, -45, -45, 45, -45, -45, -45])
    if ipop == 48:
        ss = np.array([-45, -45, -45 ,0 , 0, 0 ,0 ,0 ,0 , 0, -45, -45, 90, 90, 90, 90, 45, 90, 90, 90])
    if ipop == 49:
        ss = np.array([-45, -45, -45 ,0 , 0, 0 ,0 ,0 ,0 , 0, -45, -45, 0, 0, 0, 90, 45, 90, 0, 0])
    if ipop == 50:
        ss = np.array([45, 45, 45 ,0 , 0, 0 ,0 ,0 ,0 , 45, -45, -45, 0, 0, 0, 90, 45, 0, 0, 0])
    if ipop == 51:
        ss = np.array([ 0, 45, -45, 90, 90, 0, -45, 90, 0 , 45, -45, -45, 0, 0, 0, 90, 45, 0, 0, 0])
    if ipop == 52:
        ss = np.array([0 ,0 ,0 ,0 , 0, 0, 0 ,0 ,0 ,0 , 0, 90, 90, 90, 45, -45, 45, -45, -45, -45])
    if ipop == 53:
        ss = np.array([0 ,0 ,0 ,0 , 0, 45, 45 ,0 ,0 , 0, 90, 90, 90,45, 90, -45, 45, -45, -45, -45])
    if ipop == 54:
        ss = np.array([ -45, 45, -45, -45, -45, 0, 0, 0, 0, 45, -45, 90, 90, 0, -45, 90, 0, 45, 0, 90, ])
    if ipop == 55:
        ss = np.array([0, 45, 0, 0, -45, 45, -45, -45, -45, 45, -45, 90, 90, 0, -45, 90, 0, 45, 0, 90, ])
    if ipop == 56:
        ss = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 90, 90, 90, 90, 90 ,90 , ])
    if ipop == 57:
        ss = np.array([0, 45, 0, 0, 0, 0, 0, 0, 0, 45, -45, 90, 90, 0, -45, 90, 0, 45, 0, 90, ]) # useless from here
    if ipop == 58:
        ss = np.hstack((0*np.ones(15,), 90*np.ones(45,), 0*np.ones(25,), -45*np.ones(15,)))
    if ipop == 59:
        ss = np.hstack((0*np.ones(35,), -45*np.ones(35,), 45*np.ones(25,), 0*np.ones(5,)))
    if ipop == 60:
        ss = np.hstack((45*np.ones(42,), 0*np.ones(58,)))
    if ipop == 61:
        ss = np.hstack((0*np.ones(33,), 45*np.ones(67,)))
    if ipop == 62:
        ss = np.hstack((0*np.ones(17,), 90*np.ones(83,)))
    if ipop == 63:
        ss = np.hstack((90*np.ones(33,), 45*np.ones(33,), 0*np.ones(34,)))
    if ipop == 64:
        ss = np.hstack((90*np.ones(27,), 45*np.ones(45,), 0*np.ones(28,)))
    if ipop == 65:
        ss = np.hstack((45*np.ones(21,), -45*np.ones(79,)))

    # Complete the stacking sequence
    #print(ipop)
    ss.reshape((20,))
#    print(ipop, ss.ndim, ss)
    if ss.ndim == 2:
        ss = np.hstack((ss, np.flip(ss, axis=1)))
    else:
        ss = np.hstack((ss, np.flip(ss, axis=0)))
#    print(ss)

    # For the correct total number of plies
    if ss.size != n_plies_in_panels:
        raise Exception('Laminate with incorrect ply count at line ', ipop)

     # Storage of the stacking sequence
    ss = np.ravel(ss)
    N0 = sum(ss == 0)
    N90 = sum(ss==90)
    N45 = sum(ss==45)
    N135 = sum(ss==-45)
    ss = ss.astype(int)

    #print(ipop, N0, N90, N45, N135)
    lampam = calc_lampam(ss, constraints)

    Pop.loc[ipop, 'ply_counts'] = n_plies_in_panels

    Pop.loc[ipop, 'lampam[1]'] = lampam[0]
    Pop.loc[ipop, 'lampam[2]'] = lampam[1]
    Pop.loc[ipop, 'lampam[3]'] = lampam[2]
    Pop.loc[ipop, 'lampam[4]'] = lampam[3]
    Pop.loc[ipop, 'lampam[5]'] = lampam[4]
    Pop.loc[ipop, 'lampam[6]'] = lampam[5]
    Pop.loc[ipop, 'lampam[7]'] = lampam[6]
    Pop.loc[ipop, 'lampam[8]'] = lampam[7]
    Pop.loc[ipop, 'lampam[9]'] = lampam[8]
    Pop.loc[ipop, 'lampam[10]'] = lampam[9]
    Pop.loc[ipop, 'lampam[11]'] = lampam[10]
    Pop.loc[ipop, 'lampam[12]'] = lampam[11]

    A = A_from_lampam(lampam, mat)

    Pop.loc[ipop, 'A11'] = A[0, 0]
    Pop.loc[ipop, 'A22'] = A[1, 1]
    Pop.loc[ipop, 'A12'] = A[0, 1]
    Pop.loc[ipop, 'A66'] = A[2, 2]
    Pop.loc[ipop, 'A16'] = A[0, 2]
    Pop.loc[ipop, 'A26'] = A[1, 2]

    B = B_from_lampam(lampam, mat, constraints.sym)

    Pop.loc[ipop, 'B11'] = B[0, 0]
    Pop.loc[ipop, 'B22'] = B[1, 1]
    Pop.loc[ipop, 'B12'] = B[0, 1]
    Pop.loc[ipop, 'B66'] = B[2, 2]
    Pop.loc[ipop, 'B16'] = B[0, 2]
    Pop.loc[ipop, 'B26'] = B[1, 2]

    D = D_from_lampam(lampam, mat)

    Pop.loc[ipop, 'D11'] = D[0, 0]
    Pop.loc[ipop, 'D22'] = D[1, 1]
    Pop.loc[ipop, 'D12'] = D[0, 1]
    Pop.loc[ipop, 'D66'] = D[2, 2]
    Pop.loc[ipop, 'D16'] = D[0, 2]
    Pop.loc[ipop, 'D26'] = D[1, 2]

    Pop.loc[ipop, 'N0'] = N0
    Pop.loc[ipop, 'N45'] = N45
    Pop.loc[ipop, 'N90'] = N90
    Pop.loc[ipop, 'N-45'] = N135

    ss_flatten = np.array(ss, dtype=str)
    ss_flatten = '  '.join(ss_flatten)
    Pop.loc[ipop, 'ss'] = ss_flatten

    ipop += 1

print(' Pop Created ... ')

for i in range(1, len(Pop.index) +1):

    ss = Pop.loc[i, 'ss']
    # unflatenning
    ss = ss.split('  ')
    ss_ini = np.array(ss, dtype=int)
    ss_ini = ss_ini.astype(int)

    # Creating more stacking sequences with:
    # 0->90
    # 90-> 0
    ss = np.copy(ss_ini)
    ss[ss == 0] = 80
    ss[ss == 90] = 0
    ss[ss == 80] = 90

    lampam = calc_lampam(ss, constraints)
    N0 = sum(ss == 0)
    N90 = sum(ss==90)
    N45 = sum(ss==45)
    N135 = sum(ss==-45)

    Pop.loc[ipop, 'ply_counts'] = n_plies_in_panels

    Pop.loc[ipop, 'lampam[1]'] = lampam[0]
    Pop.loc[ipop, 'lampam[2]'] = lampam[1]
    Pop.loc[ipop, 'lampam[3]'] = lampam[2]
    Pop.loc[ipop, 'lampam[4]'] = lampam[3]
    Pop.loc[ipop, 'lampam[5]'] = lampam[4]
    Pop.loc[ipop, 'lampam[6]'] = lampam[5]
    Pop.loc[ipop, 'lampam[7]'] = lampam[6]
    Pop.loc[ipop, 'lampam[8]'] = lampam[7]
    Pop.loc[ipop, 'lampam[9]'] = lampam[8]
    Pop.loc[ipop, 'lampam[10]'] = lampam[9]
    Pop.loc[ipop, 'lampam[11]'] = lampam[10]
    Pop.loc[ipop, 'lampam[12]'] = lampam[11]

    A = A_from_lampam(lampam, mat)

    Pop.loc[ipop, 'A11'] = A[0, 0]
    Pop.loc[ipop, 'A22'] = A[1, 1]
    Pop.loc[ipop, 'A12'] = A[0, 1]
    Pop.loc[ipop, 'A66'] = A[2, 2]
    Pop.loc[ipop, 'A16'] = A[0, 2]
    Pop.loc[ipop, 'A26'] = A[1, 2]

    B = B_from_lampam(lampam, mat, constraints.sym)

    Pop.loc[ipop, 'B11'] = B[0, 0]
    Pop.loc[ipop, 'B22'] = B[1, 1]
    Pop.loc[ipop, 'B12'] = B[0, 1]
    Pop.loc[ipop, 'B66'] = B[2, 2]
    Pop.loc[ipop, 'B16'] = B[0, 2]
    Pop.loc[ipop, 'B26'] = B[1, 2]

    D = D_from_lampam(lampam, mat)

    Pop.loc[ipop, 'D11'] = D[0, 0]
    Pop.loc[ipop, 'D22'] = D[1, 1]
    Pop.loc[ipop, 'D12'] = D[0, 1]
    Pop.loc[ipop, 'D66'] = D[2, 2]
    Pop.loc[ipop, 'D16'] = D[0, 2]
    Pop.loc[ipop, 'D26'] = D[1, 2]

    Pop.loc[ipop, 'N0'] = N0
    Pop.loc[ipop, 'N45'] = N45
    Pop.loc[ipop, 'N90'] = N90
    Pop.loc[ipop, 'N-45'] = N135
    ss = np.array(ss, dtype=str)
    ss = '  '.join(ss)
    Pop.loc[ipop, 'ss'] = ss

    ipop +=1

    # Creating more stacking sequences with:
    # 0->90
    # 90-> 0
    # 45->-45
    # -45-> 45
    ss = np.copy(ss_ini)
    ss[ss == 45] = -10
    ss[ss == -45] = 45
    ss[ss == -10] = -45
    ss[ss == 0] = 80
    ss[ss == 90] = 0
    ss[ss == 80] = 90
    lampam = calc_lampam(ss, constraints)
    N0 = sum(ss == 0)
    N90 = sum(ss==90)
    N45 = sum(ss==45)
    N135 = sum(ss==-45)

    Pop.loc[ipop, 'ply_counts'] = n_plies_in_panels

    Pop.loc[ipop, 'lampam[1]'] = lampam[0]
    Pop.loc[ipop, 'lampam[2]'] = lampam[1]
    Pop.loc[ipop, 'lampam[3]'] = lampam[2]
    Pop.loc[ipop, 'lampam[4]'] = lampam[3]
    Pop.loc[ipop, 'lampam[5]'] = lampam[4]
    Pop.loc[ipop, 'lampam[6]'] = lampam[5]
    Pop.loc[ipop, 'lampam[7]'] = lampam[6]
    Pop.loc[ipop, 'lampam[8]'] = lampam[7]
    Pop.loc[ipop, 'lampam[9]'] = lampam[8]
    Pop.loc[ipop, 'lampam[10]'] = lampam[9]
    Pop.loc[ipop, 'lampam[11]'] = lampam[10]
    Pop.loc[ipop, 'lampam[12]'] = lampam[11]

    A = A_from_lampam(lampam, mat)

    Pop.loc[ipop, 'A11'] = A[0, 0]
    Pop.loc[ipop, 'A22'] = A[1, 1]
    Pop.loc[ipop, 'A12'] = A[0, 1]
    Pop.loc[ipop, 'A66'] = A[2, 2]
    Pop.loc[ipop, 'A16'] = A[0, 2]
    Pop.loc[ipop, 'A26'] = A[1, 2]

    B = B_from_lampam(lampam, mat, constraints.sym)

    Pop.loc[ipop, 'B11'] = B[0, 0]
    Pop.loc[ipop, 'B22'] = B[1, 1]
    Pop.loc[ipop, 'B12'] = B[0, 1]
    Pop.loc[ipop, 'B66'] = B[2, 2]
    Pop.loc[ipop, 'B16'] = B[0, 2]
    Pop.loc[ipop, 'B26'] = B[1, 2]

    D = D_from_lampam(lampam, mat)

    Pop.loc[ipop, 'D11'] = D[0, 0]
    Pop.loc[ipop, 'D22'] = D[1, 1]
    Pop.loc[ipop, 'D12'] = D[0, 1]
    Pop.loc[ipop, 'D66'] = D[2, 2]
    Pop.loc[ipop, 'D16'] = D[0, 2]
    Pop.loc[ipop, 'D26'] = D[1, 2]

    Pop.loc[ipop, 'N0'] = N0
    Pop.loc[ipop, 'N45'] = N45
    Pop.loc[ipop, 'N90'] = N90
    Pop.loc[ipop, 'N-45'] = N135
    ss = np.array(ss, dtype=str)
    ss = '  '.join(ss)
    Pop.loc[ipop, 'ss'] = ss

    ipop +=1

    # Creating more stacking sequences with:
    # 45->-45
    # -45-> 45
    ss = np.copy(ss_ini)
    ss[ss == 45] = -10
    ss[ss == -45] = 45
    ss[ss == -10] = -45
    lampam = calc_lampam(ss, constraints)
    N0 = sum(ss == 0)
    N90 = sum(ss==90)
    N45 = sum(ss==45)
    N135 = sum(ss==-45)

    Pop.loc[ipop, 'ply_counts'] = n_plies_in_panels

    Pop.loc[ipop, 'lampam[1]'] = lampam[0]
    Pop.loc[ipop, 'lampam[2]'] = lampam[1]
    Pop.loc[ipop, 'lampam[3]'] = lampam[2]
    Pop.loc[ipop, 'lampam[4]'] = lampam[3]
    Pop.loc[ipop, 'lampam[5]'] = lampam[4]
    Pop.loc[ipop, 'lampam[6]'] = lampam[5]
    Pop.loc[ipop, 'lampam[7]'] = lampam[6]
    Pop.loc[ipop, 'lampam[8]'] = lampam[7]
    Pop.loc[ipop, 'lampam[9]'] = lampam[8]
    Pop.loc[ipop, 'lampam[10]'] = lampam[9]
    Pop.loc[ipop, 'lampam[11]'] = lampam[10]
    Pop.loc[ipop, 'lampam[12]'] = lampam[11]

    A = A_from_lampam(lampam, mat)

    Pop.loc[ipop, 'A11'] = A[0, 0]
    Pop.loc[ipop, 'A22'] = A[1, 1]
    Pop.loc[ipop, 'A12'] = A[0, 1]
    Pop.loc[ipop, 'A66'] = A[2, 2]
    Pop.loc[ipop, 'A16'] = A[0, 2]
    Pop.loc[ipop, 'A26'] = A[1, 2]

    B = B_from_lampam(lampam, mat, constraints.sym)

    Pop.loc[ipop, 'B11'] = B[0, 0]
    Pop.loc[ipop, 'B22'] = B[1, 1]
    Pop.loc[ipop, 'B12'] = B[0, 1]
    Pop.loc[ipop, 'B66'] = B[2, 2]
    Pop.loc[ipop, 'B16'] = B[0, 2]
    Pop.loc[ipop, 'B26'] = B[1, 2]

    D = D_from_lampam(lampam, mat)

    Pop.loc[ipop, 'D11'] = D[0, 0]
    Pop.loc[ipop, 'D22'] = D[1, 1]
    Pop.loc[ipop, 'D12'] = D[0, 1]
    Pop.loc[ipop, 'D66'] = D[2, 2]
    Pop.loc[ipop, 'D16'] = D[0, 2]
    Pop.loc[ipop, 'D26'] = D[1, 2]

    Pop.loc[ipop, 'N0'] = N0
    Pop.loc[ipop, 'N45'] = N45
    Pop.loc[ipop, 'N90'] = N90
    Pop.loc[ipop, 'N-45'] = N135

    ss = np.array(ss, dtype=str)
    ss = '  '.join(ss)
    Pop.loc[ipop, 'ss'] = ss

    ipop +=1

Pop.drop_duplicates(keep='first', inplace=True)
Pop.reset_index(inplace=True)

print(f'The population consist of {len(Pop.index)} individuals')


save_constraints_LAYLA(filename, constraints)
save_materials(filename, mat)
append_df_to_excel(filename, Pop, 'stacks', index=True, header=True)
autofit_column_widths(filename)
