# -*- coding: utf-8 -*-
"""
Calculation of inhomogeneity factor for a population of stacking sequence 

@author: Noemie Fedon
"""

import sys
sys.path.append(r'C:\LAYLA')
import numpy as np
import pandas as pd
from src.CLA.lampam_functions import calc_lampam

# Creation of a table of stacking sequences
ss = np.array([0, 45, 90, -45, 0, 45, 90, -45])
ss = np.hstack((ss, np.flip(ss, axis=0)))
ss = np.array(ss, dtype=str)
ss = '  '.join(ss)
sst = ss

ss = np.array([0, 45, 90, -45, 0, 45, 90, -45])
ss = np.matlib.repmat(ss, 1, 2)
ss = np.ravel(ss)
ss = np.hstack((ss, np.flip(ss, axis=0)))
ss = np.array(ss, dtype=str)
ss = '  '.join(ss)
sst = np.vstack((sst, ss))

ss = np.array([0, 45, 90, -45, 0, 45, 90, -45])
ss = np.matlib.repmat(ss, 1, 3)
ss = np.ravel(ss)
ss = np.hstack((ss, np.flip(ss, axis=0)))
ss = np.array(ss, dtype=str)
ss = '  '.join(ss)
sst = np.vstack((sst, ss))

ss = np.array([0, 45, 90, -45, 0, 45, 90, -45])
ss = np.matlib.repmat(ss, 1, 4)
ss = np.ravel(ss)
ss = np.hstack((ss, np.flip(ss, axis=0)))
ss = np.array(ss, dtype=str)
ss = '  '.join(ss)
sst = np.vstack((sst, ss))

ss = np.array([0, 0, 0, 0, 90, 90, 90, 90])
ss = np.hstack((ss, np.flip(ss, axis=0)))
ss = np.array(ss, dtype=str)
ss = '  '.join(ss)
sst = np.vstack((sst, ss))

ss = np.array([90, 90, 90, 90, 0, 0, 0, 0])
ss = np.hstack((ss, np.flip(ss, axis=0)))
ss = np.array(ss, dtype=str)
ss = '  '.join(ss)
sst = np.vstack((sst, ss))

ss = np.array([0, 0, 0, 0, 0, 90, 90, 90])
ss = np.hstack((ss, np.flip(ss, axis=0)))
ss = np.array(ss, dtype=str)
ss = '  '.join(ss)
sst = np.vstack((sst, ss))

ss = np.array([90, 90, 90, 0, 0, 0, 0, 0])
ss = np.hstack((ss, np.flip(ss, axis=0)))
ss = np.array(ss, dtype=str)
ss = '  '.join(ss)
sst = np.vstack((sst, ss))

ss = np.array([0, 0, 0, 0, 0, 0, 90, 90])
ss = np.hstack((ss, np.flip(ss, axis=0)))
ss = np.array(ss, dtype=str)
ss = '  '.join(ss)
sst = np.vstack((sst, ss))

ss = np.array([90, 90, 0, 0, 0, 0, 0, 0])
ss = np.hstack((ss, np.flip(ss, axis=0)))
ss = np.array(ss, dtype=str)
ss = '  '.join(ss)
sst = np.vstack((sst, ss))

ss = np.array([90, 0, 0, 0, 0, 0, 0, 0])
ss = np.hstack((ss, np.flip(ss, axis=0)))
ss = np.array(ss, dtype=str)
ss = '  '.join(ss)
sst = np.vstack((sst, ss))

ss = np.array([0, 0, 0, 0, 0, 0, 0, 90])
ss = np.hstack((ss, np.flip(ss, axis=0)))
ss = np.array(ss, dtype=str)
ss = '  '.join(ss)
sst = np.vstack((sst, ss))

ss = np.array([0, 0, 0, 0, 0, 0, 0, 0])
ss = np.hstack((ss, np.flip(ss, axis=0)))
ss = np.array(ss, dtype=str)
ss = '  '.join(ss)
sst = np.vstack((sst, ss))

ss = np.array([90, 90, 90, 90, 90, 0, 0, 0])
ss = np.hstack((ss, np.flip(ss, axis=0)))
ss = np.array(ss, dtype=str)
ss = '  '.join(ss)
sst = np.vstack((sst, ss))

ss = np.array([0, 0, 0, 90, 90, 90, 90, 90])
ss = np.hstack((ss, np.flip(ss, axis=0)))
ss = np.array(ss, dtype=str)
ss = '  '.join(ss)
sst = np.vstack((sst, ss))

ss = np.array([90, 90, 90, 90, 90, 90, 0, 0])
ss = np.hstack((ss, np.flip(ss, axis=0)))
ss = np.array(ss, dtype=str)
ss = '  '.join(ss)
sst = np.vstack((sst, ss))

ss = np.array([0, 0, 90, 90, 90, 90, 90, 90])
ss = np.hstack((ss, np.flip(ss, axis=0)))
ss = np.array(ss, dtype=str)
ss = '  '.join(ss)
sst = np.vstack((sst, ss))

ss = np.array([0, 90, 90, 90, 90, 90, 90, 90])
ss = np.hstack((ss, np.flip(ss, axis=0)))
ss = np.array(ss, dtype=str)
ss = '  '.join(ss)
sst = np.vstack((sst, ss))

ss = np.array([90, 90, 90, 90, 90, 90, 90, 0])
ss = np.hstack((ss, np.flip(ss, axis=0)))
ss = np.array(ss, dtype=str)
ss = '  '.join(ss)
sst = np.vstack((sst, ss))

ss = np.array([90, 90, 90, 90, 90, 90, 90, 0])
ss = np.hstack((ss, np.flip(ss, axis=0)))
ss = np.array(ss, dtype=str)
ss = '  '.join(ss)
sst = np.vstack((sst, ss))

ss = np.array([45, 45, 0, 0, 0, 0, 0, 0])
ss = np.hstack((ss, np.flip(ss, axis=0)))
ss = np.array(ss, dtype=str)
ss = '  '.join(ss)
sst = np.vstack((sst, ss))

ss = np.array([-45, -45, 0, 0, 0, 0, 0, 0])
ss = np.hstack((ss, np.flip(ss, axis=0)))
ss = np.array(ss, dtype=str)
ss = '  '.join(ss)
sst = np.vstack((sst, ss))

ss = np.array([0, 0, 0, 0, 0, 0, 45, 45])
ss = np.hstack((ss, np.flip(ss, axis=0)))
ss = np.array(ss, dtype=str)
ss = '  '.join(ss)
sst = np.vstack((sst, ss))

ss = np.array([0, 0, 0, 0, 0, 0, -45, -45])
ss = np.hstack((ss, np.flip(ss, axis=0)))
ss = np.array(ss, dtype=str)
ss = '  '.join(ss)
sst = np.vstack((sst, ss))

Pop = pd.DataFrame()

for ipop in range(sst.shape[0]):

    # stacking sequence
    ss = sst[ipop][0]
    Pop.loc[ipop, 'ss'] = ss
    ss = ss.split('  ')
    ss = np.array(ss)
    ss = ss.astype(int)
    
    # lamination parameter
    lampam = calc_lampam(ss, constraints)
    
    # inhomogeneity factor
    inh = np.linalg.norm(lampam[0:4] - lampam[8:12]) 
    Pop.loc[ipop, 'inh'] = inh
    
    Pop.loc[ipop, 'lampam[9]- lampam[1]'] = lampam[8] - lampam[0] 
    Pop.loc[ipop, 'lampam[11]- lampam[3]'] = lampam[10] - lampam[2] 
    
    Pop.loc[ipop, 'lampam[1]'] = lampam[0]
    Pop.loc[ipop, 'lampam[2]'] = lampam[1] 
    Pop.loc[ipop, 'lampam[3]'] = lampam[2]    
    
    Pop.loc[ipop, 'lampam[9]'] = lampam[8] 
    Pop.loc[ipop, 'lampam[10]'] = lampam[9] 
    Pop.loc[ipop, 'lampam[11]'] = lampam[10]  
    Pop.loc[ipop, 'lampam[12]'] = lampam[11] 

    ipop += 1
    
print(f'The population consist of {len(Pop.index)} individuals')

## Write results in a Excell sheet
writer = pd.ExcelWriter('Inhomogeneity factors.xlsx')
Pop.to_excel(writer, 'Inhomogeneity factors')
writer.save()