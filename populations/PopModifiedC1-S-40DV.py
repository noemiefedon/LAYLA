# -*- coding: utf-8 -*-
"""
Function to check that the modified layerwise approach produce varied
solution stacking sequences with different initial stacking sequences

@author: Noemie Fedon
"""

import numpy as np
import pandas as pd
import sys
sys.path.append(r'C:\LAYLA')

from src.LAYLA_V02.constraints import Constraints
from src.CLA.lampam_functions import calc_lampam
from src.guidelines.ten_percent_rule import ten_percent_rule
from src.guidelines.internal_diso_contig import internal_diso_contig

n_plies_in_panels = 80
n_pop = 7
filename = 'SS_C1_40DV_Modified.xlsx'

############################################################
#
# --- Design guidelines ---
#
############################################################

# The various composite design guidelines are set in this section with the helampam of the constraints structures:
# 1. symmetry. Stacking sequence is mirrored about the mid-plane.
# 2. balance. All fibre angles, except 0deg and 90deg, occur in +- pairs.
# 3. dam_tol. Damage Tolerance, a block of two plies [45, -45] or [-45, 45] are covering the upper and lower surface of the laminate.
# 4. Rule10percent.
#       A rule with ply counts enforces a  minimum of of plies in each of the 0deg, 45deg, -45deg and 90deg.
#       A restriction of the in-plane lamination parameters 1 and 3 is expected to approach the 10# rule.
# 5. disorientation. The change of angles between two consecutive plies should not exceed delta_angle (set by user).
# 6. contiguity. No more that 'X' plies with same fibre orientation should be next to each other (set by user).
# 7. set_of_angles. The set of fibre orientations can be [-45, 0, 45, 90] or [-60, -45, -30, 0, 30, 45, 60, 90].


sym = True # symmetry
bal = False # balance
dam_tol = False # Damage tolerance

# -----     -----     -----     -----
#             10% rule
# -----     -----     -----     -----
rule_10_percent = False

# constraints.percent_0 = percentage used in the 10# rule for plies in the 0 deg fibre direction
# constraints.percent_45 = percentage used in the 10# rule for plies in the 45 deg fibre direction
# constraints.percent_90 = percentage used in the 10# rule for plies in the 90 deg fibre direction
# constraints.percent_135 =percentage used in the 10# rule for plies in the -45 deg fibre direction
# constraints.percent_45_135 =percentage used in the 10# rule for sum	of  the plies in the45 and -45 deg fibre direction
 percent_0 = 10
 percent_45 = 10
 percent_90 = 10
 percent_135 =10
 percent_45_135 = 0
# Set the precentages to 0 if they do not matter

# Type of 10# rule
ply_counts = True
# parameters.penalty_10_pc_switch = True 10# rule with counts for the ply orientations
# otherwise 10# rule with restriction of the lampam design space with hope of final result satisfying the 10# rule
# with counts for the ply orientations, otherwise LAYLA restarts
# considering the 10# rule with ply counts.

# -----     -----     -----     -----
#        disorientation
# -----     -----     -----     -----
diso = False

# Upper bound of the variation of fibre orientation between two
# contiguous plies if the disorientation constraint is active
delta_angle = 45
# delta_angle must be at least 45

# -----     -----     -----     -----
#           contiguity
# -----     -----     -----     -----
contig = False

# No more that constraints.n_contiguity plies with same fibre orientation
# should be next to each other if the contiguity constraint is active
n_contig = 5
# Optional (only needed if the contiguity constraint is active) The value
# taken can only be 2, 3, 4 or 5, otherwise modify
# anglesCombinationsFunction1 and  anglesCombinationsFunction2

# -----     -----     -----     -----     -----     -----     -----
#
#                           FIXED VALUES
#                         do not change it
#
# -----     -----     -----     -----     -----     -----     -----

set_of_angles= np.array([-45, 0, 45, 90], dtype=int)

constraints = Constraints(
    sym=sym,
    bal=bal,
    dam_tol = dam_tol,
    rule_10_percent = rule_10_percent,
    ply_counts = ply_counts,
    percent_0 = percent_0,
    percent_45 = percent_45,
    percent_90 = percent_90,
    percent_135 = percent_135,
    percent_45_135 = percent_45_135,
    diso = diso,
    contig = contig,
    n_contig = n_contig ,
    delta_angle = delta_angle,
    set_of_angles = set_of_angles)
#print(constraints)

Pop = pd.DataFrame()

# Initialisation
print(' Creating Pop ... ')
ipop = 1

# Loop until Pop is complete
while ipop < n_pop + 1:

    if ipop == 1:
        ss = 0*np.ones(40,)
    if ipop == 2:
        ss = 45*np.ones(40,)
    if ipop == 3:
        ss = np.matlib.repmat(np.array([0, 45, 90, -45]), 1, 10)
    if ipop == 4:
        ss = np.hstack((0*np.ones(10,), 45*np.ones(10,), 0*np.ones(10,), -45*np.ones(10,)))
    if ipop == 5:
        ss = np.hstack((0*np.ones(28,), 45*np.ones(4,), 0*np.ones(4,), -45*np.ones(4,)))
    if ipop == 6:
        ss = np.array([0, 0, 0, 0, 0, 45 , 90, 90, 90, 90, 90, -45,
                0, 0, 0, 0, 0, 45 , 90, 90, 90, 90, 90, -45,
                0, 0, 0, 0, 0, 45 , 90, 90, 90, 90, 90, -45,
                0, 45, 90, -45])
    if ipop == 7:
        ss = np.matlib.repmat(np.array([0, 0, 45, 0, 0, -45, 90, 45, 90, -45]), 1, 4)
#    if ipop == 8:
#        ss = np.matlib.repmat(np.array([90, 45, 0, -45]), 1, 10)
#    if ipop == 9:
#        ss = np.matlib.repmat(np.array([0 ,0 ,0 ,0 , -45, -45, -45, -45]), 1, 5)
#    if ipop == 10:
#        ss = np.matlib.repmat(np.array([90, 90, 90, 90, 0 ,0 ,0 ,0 ]), 1, 5)
#    if ipop == 11:
#        ss = np.matlib.repmat(np.array([45, 45, 45, 45, 90 ,90 ,90 ,90 ]), 1, 5)
#    if ipop == 12:
#        ss = np.matlib.repmat(np.array([-45, -45, -45, -45, 45 ,45 ,45 ,45 ]), 1, 5)
#    if ipop == 13:
#        ss = np.matlib.repmat(np.array([90, 90, 90, 90, 0 ,0 ,0 ,0 ]), 1, 5)
#    if ipop == 14:
#        ss = np.matlib.repmat(np.array([-45, 0, 0, 0 ]), 1, 10)
#    if ipop == 15:
#        ss = np.matlib.repmat(np.array([0, 90, 90, 90 ]), 1, 10)
#    if ipop == 16:
#        ss = np.matlib.repmat(np.array([90, 45, 45, 45]), 1, 10)
#    if ipop == 17:
#        ss = np.matlib.repmat(np.array([45, -45, -45, -45]), 1, 10)
#    if ipop == 18:
#        ss = np.matlib.repmat(np.array([45, -45]), 1, 20)
#    if ipop == 19:
#        ss = np.matlib.repmat(np.array([90, 45, 45, 45, 45, 45, 45, 45]), 1, 5)
#    if ipop == 20:
#        ss = np.matlib.repmat(np.array([0, 90, 90, 90 , 90, 90, 90, 90 ]), 1, 5)
#    if ipop == 21:
#        ss = np.matlib.repmat(np.array([-45, 0, 0, 0 , 0, 0, 0 , 0]), 1, 5)
#    if ipop == 22:
#        ss = np.matlib.repmat(np.array([45, -45, -45, -45, -45, -45, 0, 0]), 1, 5)
#    if ipop == 23:
#        ss = np.matlib.repmat(np.array([45, -45, -45, -45, -45, -45, 90, 90]), 1, 5)
#    if ipop == 24:
#        ss = np.matlib.repmat(np.array([45, -45, -45, -45, -45, -45, 0, 90]), 1, 5)
#    if ipop == 25:
#        ss = np.matlib.repmat(np.array([90, 45, 45, 45, 45, 45, 0, 0]), 1, 5)
#    if ipop == 26:
#        ss = np.matlib.repmat(np.array([90, 45, 45, 45, 45, 45, 90, 90]), 1, 5)
#    if ipop == 27:
#        ss = np.matlib.repmat(np.array([90, 45, 45, 45, 45, 45, 0, 90]), 1, 5)
#    if ipop == 28:
#        ss = np.matlib.repmat(np.array([0, 90, 90, 90 , 90, 90, 0, 0 ]), 1, 5)
#    if ipop == 29:
#        ss = np.matlib.repmat(np.array([0, 90, 90, 90 , 90, 90, 45, 45 ]), 1, 5)
#    if ipop == 30:
#        ss = np.matlib.repmat(np.array([0, 90, 90, 90 , 90, 90, 0, 45 ]), 1, 5)
#    if ipop == 31:
#        ss = np.matlib.repmat(np.array([-45, 0, 0, 0 , 0, 0, 90 , 90]), 1, 5)
#    if ipop == 32:
#        ss = np.matlib.repmat(np.array([-45, 0, 0, 0 , 0, 0, 45 , 45]), 1, 5)
#    if ipop == 33:
#        ss = np.hstack((0*np.ones(15,), 45*np.ones(10,), 0*np.ones(15,)))
#    if ipop == 34:
#        ss = np.matlib.repmat(np.array([0, 0, 0 ,0 ,0 ,0 ,0 , 90 , 90 ,90]), 1, 4)
#    if ipop == 35:
#        ss = np.matlib.repmat(np.array([90, 45]), 1, 20)
#    if ipop == 36:
#        ss = np.matlib.repmat(np.array([45, 90]), 1, 20)
#    if ipop == 37:
#        ss = np.hstack((0*np.ones(15,), 45*np.ones(10,), 45*np.ones(15,)))
#    if ipop == 38:
#        ss = np.matlib.repmat(np.array([45, -45]), 1, 20)
#    if ipop == 39:
#        ss = np.hstack((45*np.ones(15,), 45*np.ones(10,), 0*np.ones(15,)))
#    if ipop == 40:
#        ss = np.hstack((45*np.ones(15,), 0*np.ones(10,), 45*np.ones(15,)))
#    if ipop == 41:
#        ss = np.hstack((-45*np.ones(15,), 90*np.ones(10,), 90*np.ones(15,)))
#    if ipop == 42:
#        ss = np.matlib.repmat(np.array([45, 45, 0, 0 , 0, 90, 90 , -45, 0, 90]), 1, 4)
#    if ipop == 43:
#        ss = np.matlib.repmat(np.array([-45, -45, 90, 0 , 0, 0, 90 , -45, 0, 90]), 1, 4)
#    if ipop == 44:
#        ss = np.matlib.repmat(np.array([-45, -45, 0, 0 , 0, 45, 90 , -45, 0, 90]), 1, 4)
#    if ipop == 45:
#        ss = np.matlib.repmat(np.array([-45, -45, 45, 0 , 0, 0, 90 , -45, 0, 90]), 1, 4)
#    if ipop == 46:
#        ss = np.matlib.repmat(np.array([-45, 90, 90, 90 , 0, 90, 90 , +45, 0, 90]), 1, 4)
#    if ipop == 47:
#        ss = np.matlib.repmat(np.array([-45, -45, 90, 90 , 90, 0, 90 , -45, 0, 90]), 1, 4)
#    if ipop == 48:
#        ss = np.matlib.repmat(np.array([-45,90, 0, 0 , 0, 0, 90 , 90, 0, 90]), 1, 4)
#    if ipop == 49:
#        ss = np.matlib.repmat(np.array([90, -45, 0, 0 , 0, 0, 90 , 0, 0, 90]), 1, 4)
#    if ipop == 50:
#        ss = np.matlib.repmat(np.array([-45, -45, 0, 0 , 0, 90, 90 , 90, 0, 90]), 1, 4)
#    if ipop == 51:
#        ss = np.matlib.repmat(np.array([0, -45, 0, 0 , 0, 0, 90 , -45, 0, 90]), 1, 4)
#    if ipop == 52:
#        ss = np.hstack((0*np.ones(5,), 45*np.ones(20,), 0*np.ones(15,)))
#    if ipop == 53:
#        ss = np.hstack((0*np.ones(15,), 45*np.ones(20,), 90*np.ones(5,)))
#    if ipop == 54:
#        ss = np.hstack((0*np.ones(15,), 45*np.ones(15,), 0*np.ones(5,), 45*np.ones(5,)))
#    if ipop == 55:
#        ss = np.hstack((0*np.ones(15,), 45*np.ones(15,), 0*np.ones(5,), -45*np.ones(5,)))
#    if ipop == 56:
#        ss = np.hstack((90*np.ones(15,), 45*np.ones(15,), 0*np.ones(5,), -45*np.ones(5,)))
#    if ipop == 57:
#        ss = np.hstack((0*np.ones(15,), -45*np.ones(15,), 0*np.ones(5,), -45*np.ones(5,)))
#    if ipop == 58:
#        ss = np.hstack((0*np.ones(15,), 90*np.ones(15,), 0*np.ones(5,), -45*np.ones(5,)))
#    if ipop == 59:
#        ss = np.hstack((0*np.ones(15,), -45*np.ones(15,), 45*np.ones(5,), 0*np.ones(5,)))
#    if ipop == 60:
#        ss = np.hstack((45*np.ones(12,), 0*np.ones(28,)))
#    if ipop == 61:
#        ss = np.hstack((0*np.ones(12,), 45*np.ones(28,)))
#    if ipop == 62:
#        ss = np.hstack((0*np.ones(12,), 90*np.ones(28,)))
#    if ipop == 63:
#        ss = np.hstack((90*np.ones(6,), 45*np.ones(6,), 0*np.ones(28,)))
#    if ipop == 64:
#        ss = np.hstack((90*np.ones(8,), 45*np.ones(4,), 0*np.ones(28,)))
#    if ipop == 65:
#        ss = np.hstack((45*np.ones(12,), -45*np.ones(28,)))

    # Complete the stacking sequence
#    print(ss)
    ss.reshape((40,))
#    print(ipop, ss.ndim, ss)
    if ss.ndim == 2:
        ss = np.hstack((ss, np.flip(ss, axis=1)))
    else:
        ss = np.hstack((ss, np.flip(ss, axis=0)))


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

    Pop.loc[ipop, 'N0'] = N0
    Pop.loc[ipop, 'N45'] = N45
    Pop.loc[ipop, 'N90'] = N90
    Pop.loc[ipop, 'N-45'] = N135

    ss_flatten = np.array(ss, dtype=str)
    ss_flatten = '  '.join(ss_flatten)
    Pop.loc[ipop, 'ss'] = ss_flatten

    ipop += 1

print(' Pop Created ... ')

#for i in range(1, len(Pop.index) +1):
#
#    ss = Pop.loc[i, 'ss']
#    # unflatenning
#    ss = ss.split('  ')
#    ss_ini = np.array(ss, dtype=int)
#    ss_ini = ss_ini.astype(int)
#
#    # Creating more stacking sequences with:
#    # 0->90
#    # 90-> 0
#    ss = np.copy(ss_ini)
#    ss[ss == 0] = 80
#    ss[ss == 90] = 0
#    ss[ss == 80] = 90
#
#    lampam = calc_lampam(ss, constraints)
#    N0 = sum(ss == 0)
#    N90 = sum(ss==90)
#    N45 = sum(ss==45)
#    N135 = sum(ss==-45)
#
#    Pop.loc[ipop, 'ply_counts'] = n_plies_in_panels
#
#    Pop.loc[ipop, 'lampam[1]'] = lampam[0]
#    Pop.loc[ipop, 'lampam[2]'] = lampam[1]
#    Pop.loc[ipop, 'lampam[3]'] = lampam[2]
#    Pop.loc[ipop, 'lampam[4]'] = lampam[3]
#    Pop.loc[ipop, 'lampam[5]'] = lampam[4]
#    Pop.loc[ipop, 'lampam[6]'] = lampam[5]
#    Pop.loc[ipop, 'lampam[7]'] = lampam[6]
#    Pop.loc[ipop, 'lampam[8]'] = lampam[7]
#    Pop.loc[ipop, 'lampam[9]'] = lampam[8]
#    Pop.loc[ipop, 'lampam[10]'] = lampam[9]
#    Pop.loc[ipop, 'lampam[11]'] = lampam[10]
#    Pop.loc[ipop, 'lampam[12]'] = lampam[11]
#    Pop.loc[ipop, 'N0'] = N0
#    Pop.loc[ipop, 'N45'] = N45
#    Pop.loc[ipop, 'N90'] = N90
#    Pop.loc[ipop, 'N-45'] = N135
#    ss = np.array(ss, dtype=str)
#    ss = '  '.join(ss)
#    Pop.loc[ipop, 'ss'] = ss
#
#    ipop +=1
#
#    # Creating more stacking sequences with:
#    # 0->90
#    # 90-> 0
#    # 45->-45
#    # -45-> 45
#    ss = np.copy(ss_ini)
#    ss[ss == 45] = -10
#    ss[ss == -45] = 45
#    ss[ss == -10] = -45
#    ss[ss == 0] = 80
#    ss[ss == 90] = 0
#    ss[ss == 80] = 90
#    lampam = calc_lampam(ss, constraints)
#    N0 = sum(ss == 0)
#    N90 = sum(ss==90)
#    N45 = sum(ss==45)
#    N135 = sum(ss==-45)
#
#    Pop.loc[ipop, 'ply_counts'] = n_plies_in_panels
#
#    Pop.loc[ipop, 'lampam[1]'] = lampam[0]
#    Pop.loc[ipop, 'lampam[2]'] = lampam[1]
#    Pop.loc[ipop, 'lampam[3]'] = lampam[2]
#    Pop.loc[ipop, 'lampam[4]'] = lampam[3]
#    Pop.loc[ipop, 'lampam[5]'] = lampam[4]
#    Pop.loc[ipop, 'lampam[6]'] = lampam[5]
#    Pop.loc[ipop, 'lampam[7]'] = lampam[6]
#    Pop.loc[ipop, 'lampam[8]'] = lampam[7]
#    Pop.loc[ipop, 'lampam[9]'] = lampam[8]
#    Pop.loc[ipop, 'lampam[10]'] = lampam[9]
#    Pop.loc[ipop, 'lampam[11]'] = lampam[10]
#    Pop.loc[ipop, 'lampam[12]'] = lampam[11]
#    Pop.loc[ipop, 'N0'] = N0
#    Pop.loc[ipop, 'N45'] = N45
#    Pop.loc[ipop, 'N90'] = N90
#    Pop.loc[ipop, 'N-45'] = N135
#    ss = np.array(ss, dtype=str)
#    ss = '  '.join(ss)
#    Pop.loc[ipop, 'ss'] = ss
#
#    ipop +=1
#
#    # Creating more stacking sequences with:
#    # 45->-45
#    # -45-> 45
#    ss = np.copy(ss_ini)
#    ss[ss == 45] = -10
#    ss[ss == -45] = 45
#    ss[ss == -10] = -45
#    lampam = calc_lampam(ss, constraints)
#    N0 = sum(ss == 0)
#    N90 = sum(ss==90)
#    N45 = sum(ss==45)
#    N135 = sum(ss==-45)
#
#    Pop.loc[ipop, 'ply_counts'] = n_plies_in_panels
#
#    Pop.loc[ipop, 'lampam[1]'] = lampam[0]
#    Pop.loc[ipop, 'lampam[2]'] = lampam[1]
#    Pop.loc[ipop, 'lampam[3]'] = lampam[2]
#    Pop.loc[ipop, 'lampam[4]'] = lampam[3]
#    Pop.loc[ipop, 'lampam[5]'] = lampam[4]
#    Pop.loc[ipop, 'lampam[6]'] = lampam[5]
#    Pop.loc[ipop, 'lampam[7]'] = lampam[6]
#    Pop.loc[ipop, 'lampam[8]'] = lampam[7]
#    Pop.loc[ipop, 'lampam[9]'] = lampam[8]
#    Pop.loc[ipop, 'lampam[10]'] = lampam[9]
#    Pop.loc[ipop, 'lampam[11]'] = lampam[10]
#    Pop.loc[ipop, 'lampam[12]'] = lampam[11]
#    Pop.loc[ipop, 'N0'] = N0
#    Pop.loc[ipop, 'N45'] = N45
#    Pop.loc[ipop, 'N90'] = N90
#    Pop.loc[ipop, 'N-45'] = N135
#    ss = np.array(ss, dtype=str)
#    ss = '  '.join(ss)
#    Pop.loc[ipop, 'ss'] = ss
#
#    ipop +=1

Pop.drop_duplicates(keep='first', inplace=True)
Pop.reset_index(inplace=True)

print(f'The population consist of {len(Pop.index)} individuals')

## Write results in a Excell sheet
writer = pd.ExcelWriter(filename)
Pop.to_excel(writer)
writer.save()
