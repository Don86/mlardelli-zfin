import numpy as np
import pandas as pd

import cobra
import cobra.test
from cobra.test import create_test_model
from cobra.util.solver import linear_reaction_coefficients
from cobra.flux_analysis import flux_variability_analysis

import os
import re
from functools import reduce
import itertools


def get_gpr_list(gpr_str):
    """given a single gene_reaction_rule string, do whatever to make it into a list of gene conditions.

    PARAMS
    ------
    gpr_str: str. The gene reaction rule.

    RETURNS
    -------
    gpr_ls: list of string. The input gene reaction rule represented as a list of string.
    """
    gpr_ls = list(set(re.split('and|or', gpr_str.replace("(", "").replace(")", "").replace(" ", ""))))
    return gpr_ls


def eval_gpr(gpr_str, affected_genes_ls):
    """Gets inputs for eval, and runs eval().
    """
    gpr_ncbi_ids = get_gpr_list(gpr_str)
    # Make dict that maps ncbi_ids to a unique cond_id, because eval is weird like that
    cond_id_dict = {}
    cond_num = 0
    for k in gpr_ncbi_ids:
        cond_id_dict[str(k)] = "c"+str(cond_num)
        cond_num += 1

    # Init a dict of input boolean values for eval()
    eval_input_dict = {}
    for k in gpr_ncbi_ids:
        eval_input_dict[cond_id_dict[k]] = False

    # Form readable eval_string for eval(), by replacing ncbi_ids with cond_ids
    eval_str = gpr_str
    for k in cond_id_dict.keys():
        eval_str = eval_str.replace(k, cond_id_dict[k])

    # Update eval_input_dict{} with True values as necessary
    true_cond_id_ls = []
    for gne in affected_genes_ls:
        if gne in gpr_ncbi_ids:
            true_cond_id_ls.append(cond_id_dict[gne])
    #true_cond_id_ls = [cond_id_dict[k] for k in affected_genes_ls] # PROBLEM HERE
    for cond_id in true_cond_id_ls:
        eval_input_dict[cond_id]= True

    eval_bool = eval(eval_str, {}, eval_input_dict)

    return eval_bool
