#
# License-Identifier: GPL
#
# Copyright (C) 2024 The Yambo Team
#
# Authors: FP
#
# This file is part of the yambopy project
#
"""
This function finds degeneracies in a list of energies

Input:
    Energies --> energy array
    thr      --> threshold below which states are considered degenerate
    max_deg  --> max size of degenerate subspace

Output:
    A list of indices such as [[0,1],[2,3],[4],[5],[6,7]..], where
    degenerate indices are grouped in len>1 sublists.
"""
