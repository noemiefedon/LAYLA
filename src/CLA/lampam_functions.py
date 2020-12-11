# -*- coding: utf-8 -*-
"""
Functions to calculate lamination parameters

- filter_lampam
    returns lamination parameters where a few numerical approximations have
    been filtered regarding the constraints for laminate symmetry and balance,
    and the fibre orientations used

- calc_lampam, calc_lampam_2
    calculates the lamination parameters of one or more laminates from their
    stacking sequence

- calc_lampam_sym
    returns the lamination parameters of one or more symmetric laminates
    with even ply counts

- calc_lampam_10_sym
    returns the 10th lamination parameters of a symmetric laminates
    with even ply counts from the half stacking sequences

- calc_lampam_11_sym
    returns the 11th lamination parameters of a symmetric laminates
    with even ply counts from the half stacking sequences

- test_lampam
    returns the lamination parameter components associated to a symmetric
    multi-panel structure with the orientations of its outer plies known, the
    rest of the plies assumed as quasi-isotopic

- calc_lampam_mp
    returns the lamination parameters of a multipanel structure

- calc_delta_lampam_1
    returns ply partial lamination parameters
    (considers the two symmetric parts for a symmetric laminate)

- calc_delta_lampam
    returns the lamination parameters of ply groups plies
    (considers the two symmetric parts for a symmetric laminate)

- calc_delta_lampamA
    returns the in-plane lamination parameters of ply groups plies
    (considers the two symmetric parts for a symmetric laminate)

- calc_delta_lampamD
    returns the out-of-plane lamination parameters of ply groups plies
    (considers the two symmetric parts for a symmetric laminate)

- calc_delta_lampam_mp
    returns lamination parameters associated with the sublaminate
    corresponding to a group of  plies in a multi-panel structure

- calc_delta_lampam_mp_2
    returns the partial lamination parameters associated to the outer plies
    used to improve the damage tolerance of a multi-panel structure

- calc_delta_lampam_mp_3
    returns the lamination parameters associated to a single ply in a
    multi-panel structure
    (when the ply does not cover a panel, the lamination parameter are zeros)

- calc_delta_lampam_mp_3A
    returns the in-plane lamination parameters associated to a single ply in a
    multi-panel structure
    (when the ply does not cover a panel, the lamination parameter are zeros)

- calc_delta_lampam_mp_3D
    returns the out-of-plane lamination parameters associated to a single ply
    in a multi-panel structure
    (when the ply does not cover a panel, the lamination parameter are zeros)

- calc_delta_lampam_tab
    returns lamination parameter components associated with one group of plies
    with uniform thickness
    (considers the two symmetric parts for a symmetric laminate)

- calc_delta_lampam_tab_t
    returns lamination parameter components associated with groups of plies
    with varying thickness
    (considers the two symmetric parts for a symmetric laminate)

- calc_delta_lampam_tab_t_1
    returns lamination parameter components associated with one group of plies
    with varying thickness
    (considers the two symmetric parts for a symmetric laminate)
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import numpy as np

sys.path.append(r'C:\LAYLA')
from src.divers.pretty_print import print_lampam, print_ss
from src.LAYLA_V02.parameters import Parameters
from src.LAYLA_V02.constraints import Constraints

def filter_lampam(lampam, constraints):
    """
    returns lamination parameters where a few numerical approximations have
    been filtered regarding the constraints for laminate symmetry and balance,
    and the fibre orientations used

    OUTPUTS

    - lampam: array storing the filtered lamination parameters

    INPUTS

    - lampam: array storing the laminate lamination parameters
    - constraints: design and manufacturing guidelines
    """
    if lampam.ndim == 1:
        if constraints.sym:
            lampam[4:8] = 0
#        if constraints.bal:
#            lampam[2:4] = 0
        if constraints.n_set_of_angles:
            sett = set([0, 45, -45, 90, -90, 135, -135])
            if np.all([el in sett  for el in constraints.set_of_angles]):
                lampam[3] = 0
                lampam[7] = 0
                lampam[11] = 0
    elif lampam.ndim == 2:
        if constraints.sym:
            lampam[:, 4:8] = 0
#        if constraints.bal:
#            lampam[:, 2:4] = 0
        if constraints.n_set_of_angles:
            sett = set([0, 45, -45, 90, -90, 135, -135])
            if np.all([el in sett  for el in constraints.set_of_angles]):
                lampam[:, 3] = 0
                lampam[:, 7] = 0
                lampam[:, 11] = 0
    else:
        raise Exception('This should not happen.')
    return lampam


def calc_lampam_from_delta_lp_matrix(stack, constraints, delta_lampams):
    """
    returns the lamination parameters of a laminate

    INPUTS

    - ss: laminate stacking sequences
    - constraints: design and manufacturing guidelines
    - delta_lampams: ply partial lamination parameters
    """
    lampam = np.zeros((12,), float)
    for ind_ply in range(delta_lampams.shape[0]):
        lampam += delta_lampams[
            ind_ply, constraints.ind_angles_dict[stack[ind_ply]]]
    return lampam


def calc_lampam_2(ss):
    """
    returns the lamination parameters of one or more laminates

    OUTPUTS

    - lampam: laminate lamination parameters

    INPUTS

    - ss: laminate stacking sequences
    - constraints: design and manufacturing guidelines
    """
    if isinstance(ss, list):
        lampam = np.zeros((len(ss), 12), float)
        for index in range(len(ss)):
            lampam[index] = calc_lampam_2(ss[index])
        return lampam
    if ss.ndim == 2 and ss.shape[0] > 1:
        lampam = np.zeros((ss.shape[0], 12), float)
        for index in range(ss.shape[0]):
            lampam[index] = calc_lampam_2(ss[index])
        return lampam

    n_plies_in_panels = np.size(ss) # laminate ply count

    theta2 = np.deg2rad(2*ss.astype(float))
    theta4 = 2*theta2
    cos_sin = np.concatenate((
        np.cos(theta2),
        np.cos(theta4),
        np.sin(theta2),
        np.sin(theta4))).reshape((4, n_plies_in_panels))

    for_the_top = np.arange(n_plies_in_panels)
    z_0 = np.ones(n_plies_in_panels)
    z_2 = ((1-n_plies_in_panels/2)*z_0+for_the_top)**3 \
        - ((1-n_plies_in_panels/2)*z_0+for_the_top - 1)**3
    z_1 = ((1-n_plies_in_panels/2)*z_0+for_the_top)**2 \
        - ((1-n_plies_in_panels/2)*z_0+for_the_top - 1)**2

    return np.array([
        (1/n_plies_in_panels)*np.matmul(cos_sin, z_0),
        (2/n_plies_in_panels**2)*np.matmul(cos_sin, z_1),
        (4/n_plies_in_panels**3)*np.matmul(cos_sin, z_2)]).reshape(12)

def calc_lampam(ss, constraints=None):
    """
    returns the lamination parameters of one or more laminates

    OUTPUTS

    - lampam: laminate lamination parameters

    INPUTS

    - ss: laminate stacking sequences
    - constraints: design and manufacturing guidelines
    """
    if constraints is None:
        return calc_lampam_2(ss)

    if isinstance(ss, list):
        lampam = np.zeros((len(ss), 12), float)
        for index in range(len(ss)):
            lampam[index] = calc_lampam(ss[index], constraints)
        return lampam
    if ss.ndim == 2 and ss.shape[0] > 1:
        lampam = np.zeros((ss.shape[0], 12), float)
        for index in range(ss.shape[0]):
            lampam[index] = calc_lampam(ss[index], constraints)
        return lampam
    n_plies_in_panels = np.size(ss) # laminate ply count

    if not constraints.sym:
        cos_sin = np.empty((4, n_plies_in_panels), float)
        for ind in range(n_plies_in_panels):
            cos_sin[:, ind] = np.copy(constraints.cos_sin[
                constraints.ind_angles_dict[ss[ind]]].reshape((4, )))

        for_the_top = np.arange(n_plies_in_panels)
        z_0 = np.ones(n_plies_in_panels)
        z_2 = ((1-n_plies_in_panels/2)*z_0+for_the_top)**3 \
            - ((1-n_plies_in_panels/2)*z_0+for_the_top - 1)**3
        z_1 = ((1-n_plies_in_panels/2)*z_0+for_the_top)**2 \
            - ((1-n_plies_in_panels/2)*z_0+for_the_top - 1)**2
        return np.array([
            (1/n_plies_in_panels)*np.matmul(cos_sin, z_0),
            (2/n_plies_in_panels**2)*np.matmul(cos_sin, z_1),
            (4/n_plies_in_panels**3)*np.matmul(cos_sin, z_2)]).reshape(12)

    cos_sin = np.empty((4, np.size(ss) // 2), float)
    for ind in range(np.size(ss) // 2):
        cos_sin[:, ind] = constraints.cos_sin[
            constraints.ind_angles_dict[ss[ind]]].reshape((4,))

    for_the_top = np.arange(np.size(ss) // 2)
    z_0 = np.ones(np.size(ss) // 2)
    z_2 = ((1 - n_plies_in_panels / 2) * z_0 + for_the_top) ** 3 \
    - ((1 - n_plies_in_panels / 2) * z_0 + for_the_top - 1) ** 3
    lampam = np.array([
        (2/n_plies_in_panels)*np.matmul(cos_sin, z_0),
        np.array([0, 0, 0, 0]),
        (8/n_plies_in_panels**3)*np.matmul(cos_sin, z_2)]).reshape(12)

    if np.size(ss) % 2:
        cos_sin_mid = constraints.cos_sin[
            constraints.ind_angles_dict[ss[n_plies_in_panels // 2]]]
        lampam += np.array([
            (1/n_plies_in_panels)*cos_sin_mid,
            np.zeros((4,), dtype=float),
            (1/n_plies_in_panels**3)*cos_sin_mid]).reshape(12)
    return lampam

def calc_lampam_sym(ss, constraints):
    """
    returns the lamination parameters of one or more symmetric laminates
    with even ply counts from the half stacking sequences

    OUTPUTS

    - lampam: laminate lamination parameters

    INPUTS

    - ss: half laminate stacking sequences
    - constraints: design and manufacturing guidelines
    """
    if isinstance(ss, list):
        lampam = np.zeros((len(ss), 12), float)
        for index in range(len(ss)):
            lampam[index] = calc_lampam_sym(ss[index], constraints)
        return lampam
    if ss.ndim == 2 and ss.shape[0] > 1:
        lampam = np.zeros((ss.shape[0], 12), float)
        for index in range(ss.shape[0]):
            lampam[index] = calc_lampam_sym(ss[index], constraints)
        return lampam

    n_plies_in_panels = 2 * np.size(ss) # laminate ply count

    cos_sin = np.empty((4, n_plies_in_panels // 2), float)
    for ind in range(n_plies_in_panels // 2):
        cos_sin[:, ind] = constraints.cos_sin[
            constraints.ind_angles_dict[ss[ind]]].reshape((4, ))

    for_the_top = np.arange(n_plies_in_panels // 2)
    z_0 = np.ones(n_plies_in_panels // 2)
    z_2 = ((1 - n_plies_in_panels / 2) * z_0 + for_the_top) ** 3 \
    - ((1 - n_plies_in_panels / 2) * z_0 + for_the_top - 1) ** 3
    lampam = np.array([
        (2 / n_plies_in_panels)*np.matmul(cos_sin, z_0),
        np.array([0, 0, 0, 0]),
        (8 / n_plies_in_panels**3)*np.matmul(cos_sin, z_2)]).reshape(12)
    return lampam


def calc_lampam_10_sym(ss, constraints):
    """
    returns the 10th lamination parameters of a symmetric laminates
    with even ply counts from the half stacking sequences

    INPUTS

    - ss: half laminate stacking sequences
    - constraints: design and manufacturing guidelines
    """
    n_plies_in_panels = 2 * np.size(ss) # laminate ply count

    cos_sin = np.empty((n_plies_in_panels // 2), float)
    for ind in range(n_plies_in_panels // 2):
        cos_sin[ind] = constraints.cos_sin[
            constraints.ind_angles_dict[ss[ind]]][1]

    for_the_top = np.arange(n_plies_in_panels // 2)
    z_0 = np.ones(n_plies_in_panels // 2)
    z_2 = ((1 - n_plies_in_panels / 2) * z_0 + for_the_top) ** 3 \
    - ((1 - n_plies_in_panels / 2) * z_0 + for_the_top - 1) ** 3

    return (8 / n_plies_in_panels**3)*np.matmul(cos_sin, z_2)

def calc_lampam_11_sym(ss, constraints):
    """
    returns the 11th lamination parameters of a symmetric laminates
    with even ply counts from the half stacking sequences

    INPUTS

    - ss: half laminate stacking sequences
    - constraints: design and manufacturing guidelines
    """
    n_plies_in_panels = 2 * np.size(ss) # laminate ply count

    cos_sin = np.empty((n_plies_in_panels // 2), float)
    for ind in range(n_plies_in_panels // 2):
        cos_sin[ind] = constraints.cos_sin[
            constraints.ind_angles_dict[ss[ind]]][2]

    for_the_top = np.arange(n_plies_in_panels // 2)
    z_0 = np.ones(n_plies_in_panels // 2)
    z_2 = ((1 - n_plies_in_panels / 2) * z_0 + for_the_top) ** 3 \
    - ((1 - n_plies_in_panels / 2) * z_0 + for_the_top - 1) ** 3

    return (8 / n_plies_in_panels**3)*np.matmul(cos_sin, z_2)



def calc_lampam_mp(sslist, constraints):
    """
    returns lamination parameter components associated to an entire
    multipanel structure

    OUTPUTS

     - lampam: array storing the lamination parameters of the laminates (array)

    INPUTS

     - sslist: list of the stacking sequences for each panel
     - constraints: design and manufacturing guidelines
    """
    n_panels = len(sslist)
    lampam = np.zeros((n_panels, 12), dtype=float)
    for ind_panel in range(n_panels):
        lampam[ind_panel] = calc_lampam(sslist[ind_panel], constraints)
    return lampam


def test_lampam(ss_top, n_plies_per_panel):
    """
    returns the lamination parameter components associated to a symmetric
    multi-panel structure with the orientations of its outer plies known
    (in ss_top), the rest of the plies assumed as quasi-isotopic (lampam = 0)

    INPUTS

    - n_plies_per_panel: list of the number of plies per panels
    - ss_top: list of the panel partial half stacking sequences
    """
    n_panels = len(ss_top)
    if n_panels != len(n_plies_per_panel):
        raise Exception('This should not happen!')
    lampam = np.zeros((n_panels, 12), dtype=float)
    for ind_panel in range(n_panels):
        if n_plies_per_panel[ind_panel] // 2 + 1 == len(ss_top[ind_panel]):
            middle_ply = len(ss_top[ind_panel])
            n_plies_group = len(ss_top[ind_panel]) - 1
        else:
            middle_ply = 0
            n_plies_group = len(ss_top[ind_panel])
        lampam[ind_panel] = calc_delta_lampam(
            ss_top[ind_panel],
            n_first_ply=1,
            n_plies_group=n_plies_group,
            n_plies_in_panels=n_plies_per_panel[ind_panel],
            constraints=constraints,
            middle_ply=middle_ply)
    return lampam


def calc_delta_lampam(ss, n_first_ply, n_plies_group, n_plies_in_panels,
                      constraints, middle_ply=0):
    """
    returns the lamination parameters of ply groups plies taking into account
    the two symmetric part for a symmetric sublaminate.
    Attention: if a middle ply + X plies are accounted for, enter X for the
    number of plies to consider and not X + 1/2

    OUTPUTS

    - delta_lampam: array storing the sublaminate partial lamination parameters

    INPUTS

    - ss: array storing the sublaminate stacking sequence
    - n_first_ply is the number of the top ply in the sublaminate with a
    numbering starting from the bottom to the top of the laminate (int)
    - n_plies_group: ply count of the sublaminate (int),
      BEWARE: n_plies_group does not account for any middle ply!!!
    - n_plies_in_panels: ply count of the laminate (int)
    - constraints: design and manufacturing guidelines
    - middle_ply = 0 if there is no ply overlapping the mid-surface,
    otherwise middle_ply is equal to the number of this ply
    """
#    print('n_first_ply', n_first_ply)
#    print('n_plies_group', n_plies_group)
#    print('n_plies_in_panels', n_plies_in_panels)
#    print('ss', ss)
    if n_plies_group > ss.size:
        raise Exception("""
The stacking sequence of the sublaminate does not have enough plies.""")
    if n_plies_group + n_first_ply - 1 > n_plies_in_panels:
        raise Exception("""
The sublaminate is not defined as to be within the laminate.""")

    cos_sin = np.empty((4, n_plies_group), float)
    for ind in range(n_plies_group):
        cos_sin[:, ind] = constraints.cos_sin[
            constraints.ind_angles_dict[ss[ind]]].reshape((4, ))

    for_the_top = np.arange(n_plies_group)
    z_0 = np.ones(n_plies_group)
    z_2 = ((n_first_ply-n_plies_in_panels/2)*z_0+for_the_top)**3 \
        - ((n_first_ply-n_plies_in_panels/2)*z_0+for_the_top - 1)**3
    if not constraints.sym:
        z_1 = ((n_first_ply - n_plies_in_panels/2)*z_0+for_the_top)**2 \
            - ((n_first_ply - n_plies_in_panels/2)*z_0+for_the_top - 1)**2
        delta_lampam = np.array([
            (1/n_plies_in_panels)*np.matmul(cos_sin, z_0),
            (2/n_plies_in_panels**2)*np.matmul(cos_sin, z_1),
            (4/n_plies_in_panels**3)*np.matmul(cos_sin, z_2)]).reshape(12)
    else:
        delta_lampam = np.array([
            (2/n_plies_in_panels)*np.matmul(cos_sin, z_0),
            np.zeros((4,), dtype=float),
            (8/n_plies_in_panels**3)*np.matmul(cos_sin, z_2)]).reshape(12)

        # Add the contribution of a ply overlapping the middle surface
        if n_first_ply + n_plies_group == middle_ply:
            cos_sin_mid = constraints.cos_sin[
                constraints.ind_angles_dict[ss[n_plies_group]]].reshape(4)
            delta_lampam += np.array([
                (1/n_plies_in_panels)*cos_sin_mid,
                np.zeros((4,), dtype=float),
                (1/n_plies_in_panels**3)*cos_sin_mid]).reshape(12)
    return delta_lampam

def calc_delta_lampamA(ss, n_first_ply, n_plies_group, n_plies_in_panels,
                       constraints, middle_ply=0):
    """
    returns the in-pnae lamination parameters of ply groups plies taking into
    account the two symmetric part for a symmetric sublaminate.
    Attention: if a middle ply + X plies are accounted for, enter X for the
    number of plies to consider and not X + 1/2

    OUTPUTS

    - delta_lampam: array storing the sublaminate partial lamination parameters

    INPUTS

    - ss: array storing the sublaminate stacking sequence
    - n_first_ply is the number of the top ply in the sublaminate with a
    numbering starting from the bottom to the top of the laminate (int)
    - n_plies_group: ply count of the sublaminate (int),
      BEWARE: n_plies_group does not account for any middle ply!!!
    - n_plies_in_panels: ply count of the laminate (int)
    - constraints: design and manufacturing guidelines
    - middle_ply = 0 if there is no ply overlapping the mid-surface,
    otherwise middle_ply is equal to the number of this ply
    """
    if n_plies_group > ss.size:
        raise Exception("""
The stacking sequence of the sublaminate does not have enough plies.""")
    if n_plies_group + n_first_ply - 1 > n_plies_in_panels:
        raise Exception("""
The sublaminate is not defined as to be within the laminate.""")

    cos_sin = np.zeros((4,), float)
    for ind in range(n_plies_group):
        cos_sin += constraints.cos_sin[
            constraints.ind_angles_dict[ss[ind]]].reshape(4)

    if not constraints.sym:
        return (1 / n_plies_in_panels) * cos_sin

    # Add the contribution of a ply overlapping the middle surface
    if n_first_ply + n_plies_group == middle_ply:
        cos_sin += 0.5 * constraints.cos_sin[
            constraints.ind_angles_dict[ss[n_plies_group]]].reshape(4)
    return (2/n_plies_in_panels) * cos_sin

def calc_delta_lampamD(ss, n_first_ply, n_plies_group, n_plies_in_panels,
                       constraints, middle_ply=0):
    """
    returns the out-of-plane lamination parameters of ply groups plies taking
    into account the two symmetric part for a symmetric sublaminate.
    Attention: if a middle ply + X plies are accounted for, enter X for the
    number of plies to consider and not X + 1/2

    OUTPUTS

    - delta_lampam: array storing the sublaminate partial lamination parameters

    INPUTS

    - ss: array storing the sublaminate stacking sequence
    - n_first_ply is the number of the top ply in the sublaminate with a
    numbering starting from the bottom to the top of the laminate (int)
    - n_plies_group: ply count of the sublaminate (int),
      BEWARE: n_plies does not account for any middle ply!!!
    - n_plies_in_panels: ply count of the laminate (int)
    - constraints: design and manufacturing guidelines
    - middle_ply = 0 if there is no ply overlapping the mid-surface,
    otherwise middle_ply is equal to the number of this ply
    """
    if n_plies_group > ss.size:
        raise Exception("""
The stacking sequence of the sublaminate does not have enough plies.""")
    if n_plies_group + n_first_ply - 1 > n_plies_in_panels:
        raise Exception("""
The sublaminate is not defined as to be within the laminate.""")

    cos_sin = np.empty((4, n_plies_group), float)
    for ind in range(n_plies_group):
        cos_sin[:, ind] = constraints.cos_sin[
            constraints.ind_angles_dict[ss[ind]]].reshape((4, ))

    for_the_top = np.arange(n_plies_group)
    z_0 = np.ones(n_plies_group)
    z_2 = ((n_first_ply-n_plies_in_panels/2)*z_0+for_the_top)**3 \
        - ((n_first_ply-n_plies_in_panels/2)*z_0+for_the_top - 1)**3
    if not constraints.sym:
        delta_lampam = np.array([
            (4/n_plies_in_panels**3)*np.matmul(cos_sin, z_2)]).reshape(4)
    else:
        delta_lampam = np.array([
            (8/n_plies_in_panels**3)*np.matmul(cos_sin, z_2)]).reshape(4)
        # Add the contribution of a ply overlapping the middle surface
        if n_first_ply + n_plies_group == middle_ply:
            cos_sin_mid = constraints.cos_sin[
                constraints.ind_angles_dict[ss[n_plies_group]]]
            delta_lampam += (1/n_plies_in_panels**3) * cos_sin_mid.reshape((4))
    return delta_lampam

def calc_delta_lampam_mp(ss, multipanel, constraints, inner_step=-1):
    """
    returns the lamination parameters associated with a group of plies of a
    multi-panel structure

    OUTPUTS

    - delta_lampam: array storing the sublaminate partial lamination parameters

    INPUTS

    - ss: array storing the sublaminate stacking sequence
    - multipanel: multi-panel structure
    - constraints: design and manufacturing guidelines
    - inner_step: index of the step for the inner loop
    """
    # incorrect number of panels ?
    if len(ss) != multipanel.n_panels:
        raise Exception("""
Incorrect number of input stacking sequences for the multipanel structure.""")
    # incorrect ply counts for the partial stacking sequences ?
    for ind_panel, panel in enumerate(multipanel.panels):
        if not constraints.sym:
            pass
        elif (ss[ind_panel].size == panel.n_plies_per_group[inner_step] \
            and (panel.middle_ply == 0 \
                 or inner_step != multipanel.reduced.n_groups - 1)):
            pass
        elif (ss[ind_panel].size == panel.n_plies_per_group[inner_step] + 1\
            and (panel.middle_ply != 0 \
                 and inner_step == multipanel.reduced.n_groups - 1)):
            pass
        elif ss[ind_panel].size == panel.n_plies \
        and inner_step == multipanel.reduced.n_groups - 1:
            pass
        else:
#            print('ind_panel', ind_panel)
#            print('ss[ind_panel].size,  panel.n_plies',
#                  ss[ind_panel].size, panel.n_plies)
#            print('panel.n_plies_per_group[inner_step]',
#                  panel.n_plies_per_group[inner_step])
            raise Exception("""
The input stacking sequences do not have the correct number of plies.""")

    delta_lampam = np.zeros((multipanel.n_panels, 12), dtype=float)
    for ind_panel, panel in enumerate(multipanel.panels):
        if panel.middle_ply != 0 and inner_step == multipanel.reduced.n_groups - 1:
            delta_lampam[ind_panel] = calc_delta_lampam(
                ss=ss[ind_panel],
                n_first_ply=panel.n_first_plies[inner_step],
                n_plies_group=panel.n_plies_per_group[inner_step],
                n_plies_in_panels=panel.n_plies, constraints=constraints,
                middle_ply=panel.middle_ply)
        else:
            delta_lampam[ind_panel] = calc_delta_lampam(
                ss=np.array(ss[ind_panel], int),
                n_first_ply=panel.n_first_plies[inner_step],
                n_plies_group=panel.n_plies_per_group[inner_step],
                n_plies_in_panels=panel.n_plies, constraints=constraints,
                middle_ply=0)
    return delta_lampam


def calc_delta_lampam_mp_2(ss, multipanel, constraints, ss2=None):
    """
    returns lamination parameter components associated to the outer plies used
    for damage tolerance for a multi-panel structure

    OUTPUTS

    - delta_lampam: array storing the sublaminate partial lamination parameters

    INPUTS

    - ss: list of arrays storing the fibre orientations of the two outer
    plies of each panel
    - multipanel: multi-panel structure
    - constraints: design and manufacturing guidelines
    """
    if len(ss) != multipanel.n_panels:
        raise Exception("""
Incorrect number of input stacking sequences for the multipanel structure.""")
    for ind_panel, panel in enumerate(multipanel.panels):
        if ss[ind_panel].size != 2:
            raise Exception("""
The input stacking sequences do not have a correct number of plies.""")
    delta_lampam = np.zeros((multipanel.n_panels, 12), dtype=float)
    ss = ss[0]
    if constraints.sym:
        for ind_panel, panel in enumerate(multipanel.panels):
            delta_lampam[ind_panel] = calc_delta_lampam(
                ss=ss,
                n_first_ply=1,
                n_plies_group=2,
                n_plies_in_panels=panel.n_plies,
                constraints=constraints,
                middle_ply=0)
    else:
        ss2 = ss2[0]
        for ind_panel, panel in enumerate(multipanel.panels):
            delta_lampam[ind_panel] = calc_delta_lampam(
                ss=ss,
                n_first_ply=1,
                n_plies_group=2,
                n_plies_in_panels=panel.n_plies,
                constraints=constraints,
                middle_ply=0)
            delta_lampam[ind_panel] += calc_delta_lampam(
                ss=ss2,
                n_first_ply=panel.n_plies - 1,
                n_plies_group=2,
                n_plies_in_panels=panel.n_plies,
                constraints=constraints,
                middle_ply=0)
    return delta_lampam


def calc_delta_lampam_mp_3(
        ss, n_first_ply, n_plies, constraints, middle_ply=0):
    """
    returns lamination parameter components associated to a single ply in a
    multi-panel structure. When the ply does not cover a panel, (indicated with
    n_first_ply[ind_panel <= 0) the associated lamination parameters are zeros

    OUTPUTS

    - delta_lampam: array storing the partial lamination parameters

    INPUTS

    - ss: fibre orientation of the ply added to the multi-panel structure
    - n_first_ply[ind_panel] = the number of the ply in each panel,
    otherwise 0 if the patch does not cover some panels.
    - n_plies: number of plies of the laminates of each panel
    - constraints: design and manufacturing guidelines
    - middle_ply[ind_panel] = 0 if there is no ply overlapping the mid-surface,
    otherwise, middle_ply is equal to the position number of this ply
    """
    n_panels = n_first_ply.size # number of panels
    ss = np.array([ss])
    delta_lampam = np.zeros((n_panels, 12), dtype=float)
    if middle_ply == 0:
        middle_ply = np.zeros((n_panels,))
    for ind_panel in range(n_panels):
#        print('ind_panel', ind_panel)
#        print('n_first_ply[ind_panel]', n_first_ply[ind_panel])
        if n_first_ply[ind_panel] > 0:
            if middle_ply[ind_panel] == n_first_ply[ind_panel]: # a middle ply
                delta_lampam[ind_panel] = calc_delta_lampam(
                    ss,
                    n_first_ply=n_first_ply[ind_panel],
                    n_plies_group=0,
                    n_plies_in_panels=n_plies[ind_panel],
                    constraints=constraints,
                    middle_ply=middle_ply[ind_panel])
            else: # not a middle ply
                delta_lampam[ind_panel] = calc_delta_lampam(
                    ss,
                    n_first_ply=n_first_ply[ind_panel],
                    n_plies_group=1,
                    n_plies_in_panels=n_plies[ind_panel],
                    constraints=constraints,
                    middle_ply=0)
    return delta_lampam

def calc_delta_lampam_mp_3A(
        ss, n_first_ply, n_plies, constraints, middle_ply=0):
    """
    returns the in-plane lamination parameters associated to a single ply in a
    multi-panel structure. When the ply does not cover a panel, (indicated with
    n_first_ply[ind_panel <= 0) the associated lamination parameters are zeros

    OUTPUTS

    - delta_lampam: array storing the partial lamination parameters

    INPUTS

    - ss: fibre orientation of the ply added to the multi-panel structure
    - n_first_ply[ind_panel] = the number of the ply in each panel,
    otherwise 0 if the patch does not cover some panels.
    - n_plies: number of plies of the laminates of each panel
    - constraints: design and manufacturing guidelines
    - middle_ply[ind_panel] = 0 if there is no ply overlapping the mid-surface,
    otherwise, middle_ply is equal to the position number of this ply
    """
    n_panels = n_first_ply.size # number of panels
    ss = np.array([ss])
    delta_lampam = np.zeros((n_panels, 4), dtype=float)
    if middle_ply == 0:
        middle_ply = np.zeros((n_panels,))
    for ind_panel in range(n_panels):
        if n_first_ply[ind_panel] > 0:
            if middle_ply[ind_panel] == n_first_ply[ind_panel]: # a middle ply
                delta_lampam[ind_panel] = calc_delta_lampamA(
                    ss,
                    n_first_ply=n_first_ply[ind_panel],
                    n_plies_group=0,
                    n_plies_in_panels=n_plies[ind_panel],
                    constraints=constraints,
                    middle_ply=middle_ply[ind_panel])
            else: # not a middle ply
                delta_lampam[ind_panel] = calc_delta_lampamA(
                    ss,
                    n_first_ply=n_first_ply[ind_panel],
                    n_plies_group=1,
                    n_plies_in_panels=n_plies[ind_panel],
                    constraints=constraints,
                    middle_ply=0)
    return delta_lampam

def calc_delta_lampam_mp_3D(
        ss, n_first_ply, n_plies, constraints, middle_ply=0):
    """
    returns the out-of-plane lamination parameters associated to a single ply
    in a multi-panel structure. When the ply does not cover a panel,
    (indicated with n_first_ply[ind_panel <= 0) the associated lamination
    parameters are zeros

    OUTPUTS

    - delta_lampam: array storing the partial lamination parameters

    INPUTS

    - ss: fibre orientation of the ply added to the multi-panel structure
    - n_first_ply[ind_panel] = the number of the ply in each panel,
    otherwise 0 if the patch does not cover some panels.
    - n_plies: number of plies of the laminates of each panel
    - constraints: design and manufacturing guidelines
    - middle_ply[ind_panel] = 0 if there is no ply overlapping the mid-surface,
    otherwise, middle_ply is equal to the position number of this ply
    """
    n_panels = n_first_ply.size # number of panels
    ss = np.array([ss])
    delta_lampam = np.zeros((n_panels, 4), dtype=float)
    if middle_ply == 0:
        middle_ply = np.zeros((n_panels,))
    for ind_panel in range(n_panels):
        if n_first_ply[ind_panel] > 0:
            if middle_ply[ind_panel] == n_first_ply[ind_panel]: # a middle ply
                delta_lampam[ind_panel] = calc_delta_lampamD(
                    ss,
                    n_first_ply=n_first_ply[ind_panel],
                    n_plies_group=0,
                    n_plies_in_panels=n_plies[ind_panel],
                    constraints=constraints,
                    middle_ply=middle_ply[ind_panel])
            else: # not a middle ply
                delta_lampam[ind_panel] = calc_delta_lampamD(
                    ss,
                    n_first_ply=n_first_ply[ind_panel],
                    n_plies_group=1,
                    n_plies_in_panels=n_plies[ind_panel],
                    constraints=constraints,
                    middle_ply=0)
    return delta_lampam

def calc_delta_lampam_tab(
        angle, n_first_ply, n_plies_group, N, constraints, middle_ply=0):
    '''
    returns the partial lamination parameters for groups of plies of uniform
    thickness taking into account the two symmetric parts for symmetric
    laminates

    OUTPUTS

    - delta_lampam_tab: partial lamination parameters

    INPUTS

    - angle: the sublaminate stacking sequences columnby column
    - n_first_ply is the phe position of the first ply of the sublaminate with
    a numbering starting from the bottom to the top of the laminate (scalar)
    - n_plies_group: number of plies consisting the sublaminate (scalar)
    - N: total number of plies for the laminate (scalar)
    - middle_ply: 0 if there is no ply overlapping the mid-surface, otherwise,
    middle_ply is equal to the position number of this ply
    '''
    if n_plies_group > angle.size:
        raise Exception("""
The input set of angles have fewer elements that what is asked to be checked
""")
    if n_plies_group + n_first_ply - 1 > N:
        raise Exception("""
The sublaminate is not properly defined as to be contained within the laminate
""")

    size_delta_lampam_tab = angle.shape[0]
    delta_lampam_tab = np.empty((size_delta_lampam_tab, 12), dtype=float)
    for ii in np.arange(size_delta_lampam_tab):
        delta_lampam_tab[ii] = calc_delta_lampam(
            angle[ii],
            n_first_ply,
            n_plies_group,
            N,
            constraints,
            middle_ply)
    return delta_lampam_tab

def calc_delta_lampam_tab_t(angle, position_top, thickness, n_plies_group,
                            constraints, middle_ply=0):
    '''
    returns the partial lamination parameters for groups of plies of varying
    thickness taking into account the two symmetric parts for symmetric
    laminates

    OUTPUTS

    - delta_lampam_tab: sublaminate lamination parameters (line by line)

    INPUTS

    - angle: sublaminate stacking sequence
    - position_top normalized position of the top of the sublaminate
    - thickness: thicknesses of the plies
    - n_plies_group: ply count of the sublaminate (scalar), does not account
    for middle_ply !
    - constraints: set of constraints
    - middle_ply: 0 if there is no ply overlapping the mid-surface, otherwise,
    middle_ply is equal to the position number of this ply
    '''
    if n_plies_group > angle.size:
        raise Exception("""
The input set of angles have fewer elements that what is asked to be checked
""")
    size_delta_lampam_tab = angle.shape[0]
    delta_lampam_tab = np.empty((size_delta_lampam_tab, 12), float)
    for ii in np.arange(size_delta_lampam_tab):
        delta_lampam_tab[ii] = calc_delta_lampam_tab_t_1(
            np.array(angle[ii]),
            position_top,
            thickness,
            n_plies_group,
            constraints,
            middle_ply)
    return delta_lampam_tab


def calc_delta_lampam_tab_t_1(angle, position_top, thickness, n_plies_group,
                              constraints, middle_ply=0):
    '''
    returns the partial lamination parameters for a group of plies of varying
    thickness taking into account the two symmetric parts for symmetric
    laminates

    OUTPUTS

    - delta_lampam: sublaminate partial lamination parameters

    INPUTS

    - angle: sublaminate stacking sequence
    - position_top: normalized position of the top of the sublaminate
    - thickness: thicknesses of the plies
    - n_plies_group: ply count of the sublaminate (scalar), does not account
    for middle_ply !
    - constraints: set of constraints
    - middle_ply: 0 if there is no ply overlapping the mid-surface, otherwise,
    middle_ply is equal to the position number of this ply
    '''
#    print('angle.size', angle.size)
#    print('angle.shape', angle.shape)
    if angle.size == 0:
        return np.zeros((12,), float)

    if constraints.sym:
        if position_top - sum(thickness) < -1e-14:
            raise Exception("""
The sublaminate is not properly defined as to be contained within the laminate
""")
    else:
        if position_top - sum(thickness) < -1 -1e-14:
            raise Exception("""
The sublaminate is not properly defined as to be contained within the laminate
""")
    for_the_top = np.array([position_top])
    for i in range(n_plies_group):
        for_the_top = np.hstack((for_the_top, for_the_top[-1] - thickness[i]))
    for_the_bot = for_the_top[1:]
    for_the_top = np.delete(for_the_top, np.s_[-1], axis=0)

    cos_sin = np.empty((4, n_plies_group), float)
    for ind in range(n_plies_group):
        cos_sin[:, ind] = constraints.cos_sin[
            constraints.ind_angles_dict[angle[ind]]].reshape((4, ))

    z_0 = for_the_top - for_the_bot
    z_2 = for_the_top**3 - for_the_bot**3
    if constraints.sym: #
        delta_lampam = np.array([
            np.matmul(cos_sin, z_0),
            np.zeros((4,), dtype=float),
            np.matmul(cos_sin, z_2)]).reshape(12)
        # Add the contribution of a ply overlapping the middle surface
        if middle_ply != 0:
            if angle.size == n_plies_group + 1:
                cos_sin_mid = constraints.cos_sin[
                    constraints.ind_angles_dict[angle[-1]]]
                delta_lampam += np.array([
                    (for_the_bot[-1])*cos_sin_mid,
                    np.zeros((4,), dtype=float),
                    (for_the_bot[-1]**3)*cos_sin_mid]).reshape(12)
            else:
                raise Exception("""
The ply orientation of the middle-ply is not given as input""")
        return delta_lampam

    z_1 = -(for_the_top**2 - for_the_bot**2)
    # - correction because the gradual approach uses a top-to-bottom convention
    delta_lampam = 0.5*np.array([
        np.matmul(cos_sin, z_0),
        np.matmul(cos_sin, z_1),
        np.matmul(cos_sin, z_2)]).reshape(12)
    return delta_lampam


if __name__ == "__main__":
    constraints = Constraints(
        sym=False,
        set_of_angles=np.array([0, 45, 90, -45]))
    parameters = Parameters(constraints=constraints,
                            group_size_min=10, group_size_max=20)

    print('\n*** Test for the function filter_lampam ***\n')
#    lampam = np.arange(1, 13)
#    lampam = np.arange(1, 25).reshape((2, 12))
#    print('Lamination parameters:\n')
#    lampam = filter_lampam(lampam, constraints)
#    print_lampam(lampam[1])

    print('\n*** Test for the function calc_lampam ***\n')
    print('Input stacking sequence:\n')
    ss = np.array([45, 90, 45, 45, 0, -45, -45, 0, 90, -45])
    print(f'{ss}\n')
    print('Lamination parameters:\n')
    lampam = calc_lampam(ss)
    print_lampam(lampam)

    print('\n*** Test for the function calc_lampam_mp ***\n')
#    ss_target1 = np.array([0, 0, 0, 0])
#    ss_target2 = np.array([0, 0, 90, 0, 0])
#    sslist = [ss_target1, ss_target2]
#    print(f'sslist: {sslist}')
#    print('Lamination parameter outputs:\n')
#    lampam = calc_lampam_mp(sslist, constraints)
#    print_lampam(lampam[0])
#    print_lampam(lampam[1])

    print('\n*** Test for the function test_lampam ***\n')
#    print('Input stacking sequence:\n')
#    ss = np.array([0, 0, 90, 0, 0])
#    ss_top = [ss, ss]
#    n_plies_per_panel = [20, 10]
#    print(f'{ss}\n')
#    print('Lamination parameters:\n')
#    lampam = test_lampam(ss_top, n_plies_per_panel)
#    print_lampam(lampam[0], lampam[1])

    print("""\n*** Test for the functions:
          calc_delta_lampam
          calc_delta_lampamA
          calc_delta_lampamD ***\n""")
#    print('Input stacking sequence:\n')
#    ss = np.array([0, 0, 0, 0, 90, 90, 90, 0, 0, 0, 0])
#    n_first_ply = 1
#    n_plies_group = 5
#    n_plies_in_panels = 10
#    print(f'{ss}\n')
#    print('Lamination parameters:\n')
#    lampam = calc_delta_lampam(
#        ss, n_first_ply, n_plies_group, n_plies_in_panels, constraints,
#        middle_ply=0)
#    print_lampam(lampam)
#    lampamA = calc_delta_lampamA(
#        ss, n_first_ply, n_plies_group, n_plies_in_panels, constraints,
#        middle_ply=0)
#    print(lampamA)
#    lampamD = calc_delta_lampamD(
#        ss, n_first_ply, n_plies_group, n_plies_in_panels, constraints,
#        middle_ply=0)
#    print(lampamD)

    print('\n*** Test for the function calc_delta_lampam_mp ***\n')
#    print('Inputs:\n')
#    group_size_min = 4
#    # Desired number of plies for the groups at each outer loop
#    group_size_max = 10
#    # Maximum number of ply drop layouts to test for each group search
#    n_pdl_max = 5
#    # Relative importance of the lamination parameters from the global level
#    global_sensitivities = np.array([1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1])
#    parameters = Parameters(
#        constraints=constraints, n_outer_step=1, group_size_min=4,
#        group_size_max=10, sensitivities=global_sensitivities,
#        n_pdl_max=2, n_panels=2)
#    ss_target1 = np.zeros((10,))
#    n_plies_target1 = ss_target1.size
#    ss_target2 = np.zeros((8,))
#    n_plies_target2 = ss_target2.size
#    lampam_target1 = calc_lampam(ss_target1, constraints)
#    lampam_target2 = calc_lampam(ss_target2, constraints)
#    boundaries = np.array([[1, 0]])
#    panel_1 = Panel(
#        lampam_target=lampam_target1, n_plies=n_plies_target1, area=1,
#        constraints=constraints)
#    panel_2 = Panel(
#        lampam_target=lampam_target2, n_plies=n_plies_target2, area=1,
#        constraints)
#    multipanel = MultiPanel(panels=[panel_1, panel_2])
#    constraints.sym = False
#    divide_panels_2(multipanel, parameters, constraints, 0)
#    ss = [np.array([0, 0, 0, 0, 90, 90, 0, 0, 0, 0]),
#          np.array([0, 0, -45, -45, 90, 0, 45, 0])]
#    inner_step = -1
#    print(f'ss: {ss}')
#    print(f'sym: {sym}\n')
#    print('Lamination parameter outputs:\n')
#    print(calc_delta_lampam_mp(
#        ss, multipanel, constraints=constraints, inner_step=inner_step))

#    print('\n*** Test for the function calc_delta_lampam_mp_2 ***\n')
#    print('Inputs:\n')
#    constraints.sym = True
#    ss_target1 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
#    n_plies_target1 = ss_target1.size
#    ss_target2 = np.array([0, 0, 0, 0, 0, 0, 0, 0])
#    n_plies_target2 = ss_target2.size
#    lampam_target1 = calc_lampam(ss_target1, constraints)
#    lampam_target2 = calc_lampam(ss_target2, constraints)
#    boundaries = np.array([[1, 0]])
#    panel_1 = Panel(
#        lampam_target=lampam_target1, n_plies=n_plies_target1, area=1,
#        constraints=constraints)
#    panel_2 = Panel(
#        lampam_target=lampam_target2, n_plies=n_plies_target2, area=1,
#        constraints=constraints)
#    multipanel = MultiPanel(panels=[panel_1, panel_2])
#    divide_panels_2(multipanel, parameters, constraints, 0)
#    a = np.array([45, -45])
#    ss = [np.copy(a) for ind_panel in range(multipanel.n_panels)]
#    sym = True
#    print(f'ss: {ss}')
#    print(f'sym: {sym}\n')
#    print('Lamination parameter outputs:\n')
#    print(calc_delta_lampam_mp_2(ss, multipanel, constraints))
#
#    print("""\n*** Test for the functions:
#          calc_delta_lampam_mp_3
#          calc_delta_lampam_mp_3A
#          calc_delta_lampam_mp_3D ***\n""")
#    print('Inputs:\n')
#    ss = 0
#    n_first_ply = np.array([1, 5, 0])
#    n_plies = np.array([10, 10, 10])
#    sym = True
#    print(f'ss: {ss}')
#    print(f'n_plies: {n_plies}')
#    print(f'n_first_ply: {n_first_ply}')
#    print(f'sym: {sym}\n')
#    print('Lamination parameter outputs:\n')
#    print_lampam(calc_delta_lampam_mp_3(
#        ss, n_first_ply, n_plies, constraints, middle_ply=0)[0])
#    print_lampam(calc_delta_lampam_mp_3(
#        ss, n_first_ply, n_plies, constraints, middle_ply=0)[1])
#    print_lampam(calc_delta_lampam_mp_3(
#        ss, n_first_ply, n_plies, constraints, middle_ply=0)[2])
#    ss = 0
#    n_first_ply = np.array([2, 2])
#    n_plies = np.array([10, 6])
#    sym = True
#    print_lampam(calc_delta_lampam_mp_3(
#        ss, n_first_ply, n_plies, constraints, middle_ply=0)[0])
#    print_lampam(calc_delta_lampam_mp_3(
#        ss, n_first_ply, n_plies, constraints, middle_ply=0)[1])
#    ss = 0
#    n_first_ply = np.array([1, 5, 0])
#    n_plies = np.array([10, 10, 10])
#    sym = True
#    print(calc_delta_lampam_mp_3A(
#        ss, n_first_ply, n_plies, constraints, middle_ply=0)[0])
#    print(calc_delta_lampam_mp_3A(
#        ss, n_first_ply, n_plies, constraints, middle_ply=0)[1])
#    print(calc_delta_lampam_mp_3A(
#        ss, n_first_ply, n_plies, constraints, middle_ply=0)[2])
#    ss = 0
#    n_first_ply = np.array([2, 2])
#    n_plies = np.array([10, 6])
#    sym = True
#    print(calc_delta_lampam_mp_3A(
#        ss, n_first_ply, n_plies, constraints, middle_ply=0)[0])
#    print(calc_delta_lampam_mp_3A(
#        ss, n_first_ply, n_plies, constraints, middle_ply=0)[1])
#    ss = 0
#    n_first_ply = np.array([1, 5, 0])
#    n_plies = np.array([10, 10, 10])
#    sym = True
#    print(calc_delta_lampam_mp_3D(
#        ss, n_first_ply, n_plies, constraints, middle_ply=0)[0])
#    print(calc_delta_lampam_mp_3D(
#        ss, n_first_ply, n_plies, constraints, middle_ply=0)[1])
#    print(calc_delta_lampam_mp_3D(
#        ss, n_first_ply, n_plies, constraints, middle_ply=0)[2])
#    ss = 0
#    n_first_ply = np.array([2, 2])
#    n_plies = np.array([10, 6])
#    sym = True
#    print(calc_delta_lampam_mp_3D(
#        ss, n_first_ply, n_plies, constraints, middle_ply=0)[0])
#    print(calc_delta_lampam_mp_3D(
#        ss, n_first_ply, n_plies, constraints, middle_ply=0)[1])
