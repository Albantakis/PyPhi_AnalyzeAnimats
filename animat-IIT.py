#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# animat.py
"""
Main File to evaluate Phi over a set of animats, using .txt output files from the C++ Evolution code
This used PyPhi version 0.7.0
taking only the current state
It was used for the first sets of changing environments, presented at CNS 2015
"""

import sys
import logging
import os
import pickle
import scipy.io
import numpy as np
from joblib import Parallel, delayed
import pyphi
from pyphi import convert
from time import time
import animalysis
import matlab_dict

formatter = logging.Formatter(
    fmt='%(asctime)s [%(name)s] %(levelname)s: %(message)s')

# Global settings
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parallel computation settings
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PARALLEL_VERBOSITY = 20
# Use all but three processors.
NUMBER_OF_CORES = 6
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Flag to indicate whether we're debugging the script.
DEBUG = False

# Only compute rule 105 if debugging.
if DEBUG:
    animats_to_compute = [3]
else:
    animats_to_compute = [20, 21, 31, 39, 41, 46, 58, 64, 75, 76, 82, 93, 95, 97, 100, 102, 107, 118, 129, 142, 152, 161, 167, 168, 184, 193, 197, 198]
    #[0, 34, 40, 69, 82, 85, 89, 93, 123, 128, 150, 156, 159, 202, 220, 231, 235, 265, 276, 284, 310, 318, 336, 355, 380, 382]
    #[206, 213, 220, 240, 258, 285, 330, 344, 368, 372, 399]#range(0,50)
    #[0, 3, 4, 5, 6, 10, 18, 27, 28, 29, 31, 32, 35, 38, 47, 48, 49, 52, 57, 58, 61, 62, 64, 65, 66, 68, 70, 72, 74, 76, 79, 83, 85, 91, 92, 96, 100, 110, 111, 118, 125, 129, 132, 133, 134, 135, 141, 143, 150, 152, 153, 157, 158, 162, 164, 167, 171, 174, 180, 185, 186, 192, 199]
    #[0, 3, 4, 5, 6, 10, 18, 27, 28, 29, 31, 32, 35, 38, 47, 48, 49, 52, 57, 58, 61, 62, 64, 65, 66, 68, 70, 72, 74, 76, 79, 83, 85, 91, 92, 96]
#[100, 110, 111, 118, 125, 129, 132, 133, 134, 135, 141, 143, 150, 152, 153, 157, 158, 162, 164, 167, 171, 174, 180, 185, 186, 192, 199]
#[52, 57, 58, 61, 62, 64, 65, 66, 68, 70, 72, 74, 76, 79, 83, 85, 91, 92, 96]
#[0, 3, 4, 5, 6, 10, 18, 27, 28, 29, 31, 32, 35, 38, 47, 48, 49]
#[103, 104, 118, 130, 131, 132, 136, 145, 147, 157, 163, 170, 175, 182, 186, 187, 191, 193, 195, 197]#range(50)
# 52 57 58 61 62 64 65 66 68 70 72 74 76 79 83 85 91 92 96 100 110 111 118 125 129 132 133 134 135 141 143 150 152 153 157 158 162 164 167 171 174 180 185 186 192 199

# Specify Task 1-4
task_condition = 'c2a4_change_c36a45'
animat_generations = [29696, 30208, 59904]#range(0,60000,512)
# Flag to indicate whether we calculate all complexes and find the main
# complex, or only calculate the full system.
calculate_all_complexes_and_whole_system = True

#ANIMAT_PATH = '/Volumes/Macintosh HD 2/Simulations/Arend_XCodeAnimat/temporalSpatialIntegrationLite/work_' + task_condition + '/trial'
ANIMAT_PATH = '/Users/larissa_alb/dev/animats/results/work_' + task_condition + '/trial'
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def compute_phi_data(animat, state, log):
    log.info("\n[Animat] " + ("current state " + str(state)).center(60, " "))

    main_complex = None
    whole_system = None
    elapsed = 0
    
    animat.current_state = state
    #import pdb; pdb.set_trace()

    tic = time()
    main_complex = pyphi.compute.main_complex(animat.network)
    toc = time()

    log.info("\n[Animat] Found main_complex:")
    log.info("\n[Animat]\tNodes: " + str(main_complex.subsystem.nodes))
    log.info("\n[Animat]\tPhi: " + str(main_complex.phi))
    if calculate_all_complexes_and_whole_system:
            subsystem = pyphi.Subsystem(range(len(animat.used_nodes)),
                                        animat.network)
            whole_system = pyphi.compute.constellation(subsystem)

    #import pdb; pdb.set_trace()
    elapsed = toc - tic
    log.info('\n[Animat] Elapsed time:', elapsed)

    return main_complex, whole_system, elapsed

def evaluate_animat(animat, gen, log, animat_string, matlab_results_dir):
    log.info("\n[Animat] " + (" Generation " + str(gen) + " ").center(60, "-"))
    print("running Animat: " + str(animat) + " Generation: " + str(gen) + "\n")
    
    data_path = ANIMAT_PATH + str(animat) + "_" + str(gen) + "_"  

    #------- Make animat network for PyPhi analysis using animalysis ---------------
    # Animalysis takes care of all the issues with the raw data:
    # - Motors shouldn't connect back
    # - Sensors shoudn't receive feedback
    # - every 36 trials there is a transition which should not be counted as a causal transition
    try:
        animat_network = animalysis.Animat(data_path)
        if animat_network.used_nodes:
            # In DEBUG take only first 2 states (this is hacky, cause below l. 173 is must be more than 1 to work)
            if DEBUG:
                # TODO: not transitons but current state only now
                current_states = animat_network.ranked_states[0:2]
            else:
                current_states = animat_network.ranked_states
            
            # Compute results for this animat.
            tic = time()
            results = {
                'animat': animat_network,
                'state': tuple(filter(lambda x: x is not None, [
                    compute_phi_data(animat_network, state, log)
                    for state, count in current_states
                ]))
            }
            
            #import pdb; pdb.set_trace()
            toc = time()
            elapsed = toc - tic

            used_nodes = tuple(animat_network.used_nodes)
            # We have binary nodes only
            num_states = 2**len(used_nodes)
            number_of_connections = np.count_nonzero(animat_network.cm)
            # Get the results in a form suitable for saving in a matfile.
            matlab_results = {
                'tpm': animat_network.tpm.reshape([num_states]+[len(used_nodes)], order='F'),
                'connectivity_matrix': animat_network.cm,
                'number_of_connections': number_of_connections,
                'used_nodes': used_nodes,
                'LifeStates': animat_network.ranked_states,
                'mip_states': [
                    matlab_dict.bigmip2dict(mip, t) for mip, whole, t in results['state']
                ],
                'whole_states': [
                    matlab_dict.constellation2dict(whole) for mip, whole, t in results['state']
                ],
                'PhiMip': [mip.phi for mip, whole, t in results['state']],
                #Todo: check what the result is if main complex is empty
                'main_complex': [convert.nodes2indices(mip.subsystem.nodes) for mip, whole, t in results['state']],
                'ssphi_whole_system': [sum(c.phi for c in whole) for mip, whole, t in results['state']],
                'whole_num_concepts': [len(whole) for mip, whole, t in results['state']],
                'mip_num_concepts': [len(mip.unpartitioned_constellation) for mip, whole, t in results['state']]
            }

            log.info('\n[Animat] Total time elapsed: ' + str(elapsed))
            #import pdb; pdb.set_trace()

            # Save the matlab results in a matfile for analysis with Matlab.
            matfile_filename = os.path.join(matlab_results_dir, animat_string + "_" + str(gen))
            scipy.io.savemat(matfile_filename, matlab_results, do_compression=True)
  
    except ValueError:
        print("[Animat] Animat network cannot be build.")
        pass

# Iterate over generations
def evaluate_generations(animat):
    # This is the name of the folder (directory) where results will be saved.
    matlab_results_dir = 'matlab_results/' + task_condition + '_trial' + str(animat)
    # Make the results directory if it doesn't exist.
    os.makedirs(matlab_results_dir, exist_ok=True)

    animat_string = task_condition + '_animat_' + str(animat)
        
    log = logging.getLogger(animat_string)
    handler = logging.FileHandler('logs/' + animat_string + '.log')
    handler.setFormatter(formatter)
    log.addHandler(handler)
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(formatter)
    log.addHandler(handler)
    log.setLevel(logging.INFO)

    log.info("\n\n[Animat] " + (" Task " + task_condition + " Trial " + str(animat) + " ").center(60, '='))
    
    [evaluate_animat(animat, gen, log, animat_string, matlab_results_dir) for gen in animat_generations]

# Run everything if this file is being executed.
if __name__ == "__main__":
    tic = time()
    
    Parallel(n_jobs=(NUMBER_OF_CORES), verbose=PARALLEL_VERBOSITY)(
        delayed(evaluate_generations)(animat) for animat in animats_to_compute)
    
    toc = time()
    elapsed = toc-tic
    print("Finished in" + str(elapsed) + "\n")

