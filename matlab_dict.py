#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# matlab_dict.py
"""
Functions for translating python struct of pyphi mip into a struct that can be saved as a .mat file.
"""

import numpy as np
import pyphi
from pyphi import convert

def mice2dict(mice):
    """Convert a PyPhi Mice to a dictionary suitable for conversion to a
    Matlab structure with scipy."""
    return {
        'phi': mice.phi,
        'direction': mice.direction,
        'purview': convert.nodes2indices(mice.purview),
        'partition': (
            ((convert.nodes2indices(mice.mip.partition[0].mechanism),
              convert.nodes2indices(mice.mip.partition[0].purview)),
             (convert.nodes2indices(mice.mip.partition[1].mechanism),
              convert.nodes2indices(mice.mip.partition[1].purview)))
            if mice.mip.partition is not None else 'None'
        ),
        'repertoire': (
            mice.mip.unpartitioned_repertoire if
            mice.mip.unpartitioned_repertoire is not None else 'None'
        ),
        'partitioned_repertoire': (
            mice.mip.partitioned_repertoire if mice.mip.partitioned_repertoire
            is not None else 'None'
        )
    }


def concept2dict(c):
    """Convert a PyPhi Concept to a dictionary suitable for conversion to a
    Matlab structure with scipy."""
    return {
        'phi': c.phi,
        'mechanism': convert.nodes2indices(c.mechanism),
        'cause': mice2dict(c.cause) if c.cause is not None else 'None',
        'effect': mice2dict(c.effect) if c.effect is not None else 'None'
    }

def bigmip2dict(mip, time):
    """Convert a BigMip to a dictionary suitable for conversion to a Matlab
    structure with scipy."""
    if mip is None:
        return np.array([])
    matlab_data = {
        'PhiMip': mip.phi,
        'main_complex': convert.nodes2indices(mip.subsystem.nodes),
        'MIP1': mip.cut.intact,
        'MIP2': mip.cut.severed,
        'current_state': mip.subsystem.network.current_state,
        #'past_state': mip.subsystem.network.past_state,
        'num_concepts': len(mip.unpartitioned_constellation),
        'calculation_time': time,
        'sum_small_phi': sum(c.phi for c in mip.unpartitioned_constellation),
        'partition_only_concepts': [
            concept2dict(c) for c in mip.partitioned_constellation if
            convert.nodes2indices(c.mechanism) not in
            [
                convert.nodes2indices(upc.mechanism) for upc in
                mip.unpartitioned_constellation
            ]
        ],
        'concepts': [concept2dict(c) for c in
                     mip.unpartitioned_constellation]
    }
    return matlab_data

def constellation2dict(constellation):
    """Convert a constellation to a dictionary suitable for conversion to a Matlab
    structure with scipy."""
    if constellation is None:
        return np.array([])
    matlab_data = {
        'num_concepts': len(constellation),
        'sum_small_phi': sum(c.phi for c in constellation),
        'concepts': [concept2dict(c) for c in constellation]
    }
    return matlab_data