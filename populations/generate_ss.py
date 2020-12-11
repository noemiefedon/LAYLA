# -*- coding: utf-8 -*-
"""
Function to randomly create stacking sequences
"""
import sys
import random
import math as ma
import numpy as np

sys.path.append(r'C:\LAYLA')
from src.guidelines.internal_diso_contig import internal_diso_contig
from src.guidelines.one_stack import check_lay_up_rules
from src.guidelines.contiguity import is_contig
from src.guidelines.disorientation import is_diso_ss
from src.guidelines.ten_percent_rule import is_ten_percent_rule
from src.divers.pretty_print import print_lampam, print_ss
from src.CLA.lampam_functions import calc_lampam
from src.LAYLA_V02.parameters import Parameters
from src.LAYLA_V02.optimisation import optimisation

def generate_ss_LAYLA(n_plies, constraints, not_constraints):
    """
    randomly generates lamination parameters and retrieves stacking sequences
    satisfying the design guidelines using LAYLA

    Args:
        n_plies (scalar): number of plies in the laminate
        constraints (instance of the class Constraints): set of constraints
            that must be satisfied
        not_constraints (instance of the class Constraints): set of constraints
            that must not be satisfied

    Returns:
        a manufacturable stacking sequence and its lamination parameters
        created by the optimiser LAYLA aiming at matching lamination parameter
        targets randomly generated
    """
    if constraints.diso and not_constraints.diso:
        raise Exception("""
Decide whether or not the stacking sequences must satisfy disorientation""")
    if constraints.contig and not_constraints.contig:
        raise Exception("""
Decide whether or not the stacking sequences must satisfy contiguity""")
    if constraints.ipo and not_constraints.ipo:
        raise Exception("""
Decide whether or not the stacking sequences must satisfy balance""")
    if constraints.rule_10_percent and not_constraints.rule_10_percent:
        raise Exception("""
Decide whether or not the stacking sequences must satisfy 10% rule""")

    #==========================================================================
    # LAYLA parameters
    #==========================================================================
    n_outer_step = 1

    first_outer_loop_assumption = 'layerwise_QI'

    method = 'beam_search' # with beam search
    global_branching_limit = 20
    local_branching_limit = 10
    global_branching_limit_p = 20

    pseudo = True

    repair_membrane_switch = True
    repair_flexural_switch = True

    penalty_10_pc_switch = False
    penalty_10_lampam_switch = False
    penalty_ipo_switch = True
    penalty_bal_switch = False
    balanced_scheme = False

    coeff_10 = 1
    coeff_bal_ipo = 1
    coeff_oopo = 1

    p_A = 100
    n_D1 = 6
    n_D2 = 10

    group_size_min = 5
    group_size_max = np.array([1000, 12, 12, 12, 12])

    if constraints.set_of_angles is np.array([-45, 0, 45, 90], int):
        lampam_to_be_optimised = np.array([1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0])
    else:
        lampam_to_be_optimised = np.array([1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1])

    first_level_sensitivities = np.ones((12,), float)

    while True:

        lampam_target = 2 * np.random.random_sample((12,)) - 1
#        lampam_target = np.array([
#            -0.25662112, -0.01727515, -0.73962959, 0.08359081,
#            0.43671968, -0.7901057, 0.98404481, 0.65070345,
#            0.97056517, 0.7023994, 0.32539113, -0.35851357])
#        print('lampam_target', lampam_target)

        parameters = Parameters(
            constraints=constraints,
            lampam_target=lampam_target,
            coeff_10=coeff_10,
            coeff_bal_ipo=coeff_bal_ipo,
            coeff_oopo=coeff_oopo,
            p_A\
            =p_A,
            n_D1=n_D1,
            n_D2=n_D2,
            n_outer_step=n_outer_step,
            group_size_min=group_size_min,
            group_size_max=group_size_max,
            first_level_sensitivities=first_level_sensitivities,
            lampam_to_be_optimised=lampam_to_be_optimised,
            method=method,
            first_outer_loop_assumption=first_outer_loop_assumption,
            global_branching_limit=global_branching_limit,
            local_branching_limit=local_branching_limit,
            global_branching_limit_p=global_branching_limit_p,
            pseudo=pseudo,
            repair_membrane_switch=repair_membrane_switch,
            repair_flexural_switch=repair_flexural_switch,
            penalty_10_lampam_switch=penalty_10_lampam_switch,
            penalty_10_pc_switch=penalty_10_pc_switch,
            penalty_ipo_switch=penalty_ipo_switch,
            penalty_bal_switch=penalty_bal_switch,
            balanced_scheme=balanced_scheme)

        try:
            result = optimisation(
                parameters, constraints, np.copy(lampam_target), n_plies,
                not_constraints=not_constraints)
        except SystemExit:
            break

        ss = result[0] # solution stacking sequence
        lampam = result[1] # solution lamination parameters

        if not ss.size == n_plies:
            raise Exception('Wrong laminate ply count')

        try:
            check_lay_up_rules(ss, constraints)
        except:
            continue


        if hasattr(constraints, 'dam_tol_rule'):
            if not_constraints.diso and is_diso_ss(
                    ss,
                    not_constraints.delta_angle,
                    dam_tol=constraints.dam_tol,
                    dam_tol_rule=constraints.dam_tol_rule):
                continue
        else:
            if not_constraints.diso and is_diso_ss(
                    ss,
                    not_constraints.delta_angle,
                    dam_tol=constraints.dam_tol,
                    n_plies_dam_tol=constraints.n_plies_dam_tol):
                continue

        if not_constraints.contig and is_contig(ss, not_constraints.n_contig):
            continue

        if not_constraints.rule_10_percent \
        and is_ten_percent_rule(not_constraints, stack=ss):
            continue

        if not_constraints.ipo \
        and abs(lampam[2]) < 1e-10 and abs(lampam[3]) < 1e-10:
            continue
#
#        if is_ten_percent_rule(not_constraints, stack=ss) \
#        and abs(lampam[2]) < 1e-10 and abs(lampam[3]) < 1e-10:
#            continue

        break

    return ss, lampam


def generate_ss_randomly(n_plies, constraints, not_constraints):
    """
    randomly generates stacking sequences satisfying the design guidelines

    INPUTS

    n_plies: number of plies in the laminate
    constraints: set of constraints that must be satisfied
    not_constraints: set of constraints that must not be satisfied
    """
    if constraints.diso and not_constraints.diso:
        raise Exception("""
Decide whether or not the stacking sequences must satisfy disorientation""")
    if constraints.contig and not_constraints.contig:
        raise Exception("""
Decide whether or not the stacking sequences must satisfy contiguity""")
    if constraints.ipo and not_constraints.ipo:
        raise Exception("""
Decide whether or not the stacking sequences must satisfy balance""")
    if constraints.rule_10_percent and not_constraints.rule_10_percent:
        raise Exception("""
Decide whether or not the stacking sequences must satisfy 10% rule""")


    if constraints.dam_tol and constraints.dam_tol_rule not in {1, 2}:
        raise Exception("""
Only coded for damage tolerance rules 0, 1 and 2.""")

    if not_constraints.rule_10_percent:
        n_plies_0_lim = ma.ceil(not_constraints.percent_0 * n_plies)
        n_plies_90_lim = ma.ceil(not_constraints.percent_90 * n_plies)
        n_plies_45_lim = ma.ceil(not_constraints.percent_45 * n_plies)
        n_plies_135_lim = ma.ceil(not_constraints.percent_135 * n_plies)
        n_plies_45_135_lim = ma.ceil(not_constraints.percent_45_135 * n_plies)
        if constraints.sym:
            n_plies_0_lim = ma.ceil(n_plies_0_lim / 2)
            n_plies_90_lim = ma.ceil(n_plies_90_lim / 2)
            n_plies_45_lim = ma.ceil(n_plies_45_lim / 2)
            n_plies_135_lim = ma.ceil(n_plies_135_lim / 2)
            n_plies_45_135_lim = ma.ceil(n_plies_45_135_lim / 2)
#    print('n_plies_0_lim', n_plies_0_lim)
#    print('n_plies_45_lim', n_plies_45_lim)

    while True:

        ss = np.array((), int)
        n0 = 0
        n90 = 0

        if not_constraints.rule_10_percent:
            random_for_10 = random.randint(0, 4)

        if constraints.dam_tol:

            if constraints.dam_tol_rule == 1:

                reduced_n_plies = n_plies - 2

                if random.randint(0, 1):
                    ss = np.array(([45]), int)
                    if constraints.sym:
                        n45 = 2
                        n135= 0
                    else:
                        n45 = 4
                        n135= 0
                else:
                    ss = np.array(([-45]), int)
                    if constraints.sym:
                        n45 = 0
                        n135= 2
                    else:
                        n45 = 0
                        n135= 4

            elif constraints.dam_tol_rule == 2:

                reduced_n_plies = n_plies - 4

                if random.randint(0, 1):
                    ss = np.array(([45, -45]), int)
                else:
                    ss = np.array(([-45, 45]), int)
                if constraints.sym:
                    n45 = 2
                    n135= 2
                else:
                    n45 = 4
                    n135= 4

        else:
            reduced_n_plies = n_plies
            n45 = 0
            n135= 0

        if constraints.sym:
            reduced_n_plies //= 2


        if constraints.dam_tol and constraints.dam_tol_rule == 2:
            ss = ss[1:]

        for ind in range(reduced_n_plies):

            angle_added = False

            for angle in np.random.permutation(constraints.set_of_angles):

                ss_test = np.hstack((ss, angle))
#                print_ss(ss_test)

                if internal_diso_contig(ss_test, constraints)[0].size == 0:
                    continue

                if not_constraints.rule_10_percent:
                    if random_for_10 == 0 and n0 + 1 >= n_plies_0_lim:
                        continue
                    if random_for_10 == 1 and n90 + 1 >= n_plies_90_lim:
                        continue
                    if random_for_10 == 2 and n45 + 1 >= n_plies_45_lim:
                        continue
                    if random_for_10 == 3 and n135 + 1 >= n_plies_135_lim:
                        continue
                    if random_for_10 == 4 \
                    and n45 + n135 + 1 >= n_plies_45_135_lim:
                        continue

                angle_added = True
                if angle == 0:
                    n0 += 1
                elif angle == 90:
                    n90 += 1
                elif angle == 45:
                    n45 += 1
                elif angle == -45:
                    n135 += 1
                ss = ss_test

#                print_ss(ss)
                break

            if not angle_added:
                break

        if not angle_added:
            continue

        if constraints.dam_tol and constraints.dam_tol_rule == 2:
            ss = ss = np.hstack((-np.array(ss[0]), ss))

        if constraints.sym:
            if n_plies % 2:
                if random.randint(0, 1):
                    ss = np.hstack((ss, 0, np.flip(ss)))
                    n0 += 1
                else:
                    ss = np.hstack((ss, 90, np.flip(ss)))
                    n90 += 1
            else:
                ss = np.hstack((ss, np.flip(ss)))
        else:
            if constraints.dam_tol:
                if random.randint(0, 1):
                    ss = np.hstack((ss, 45, -45))
                else:
                    ss = np.hstack((ss, -45, 45))


        if not ss.size == n_plies:
            raise Exception('Wrong laminate ply count')

        try:
            check_lay_up_rules(ss, constraints)
        except:
            continue

        if hasattr(constraints, 'dam_tol_rule'):
            if not_constraints.diso and is_diso_ss(
                    ss,
                    not_constraints.delta_angle,
                    dam_tol=constraints.dam_tol,
                    dam_tol_rule=constraints.dam_tol_rule):
                continue
        else:
            if not_constraints.diso and is_diso_ss(
                    ss,
                    not_constraints.delta_angle,
                    dam_tol=constraints.dam_tol,
                    n_plies_dam_tol=constraints.n_plies_dam_tol):
                continue

        if not_constraints.contig and is_contig(ss, not_constraints.n_contig):
            continue

        if not_constraints.rule_10_percent \
        and is_ten_percent_rule(not_constraints, stack=ss):
            continue

        lampam = calc_lampam(ss)

        if not_constraints.bal \
        and abs(lampam[2]) < 1e-10 and abs(lampam[3]) < 1e-10:
            continue

        break

    return ss, lampam
