"""Build regionalized voxel model with excluded experiments.

This script rebuilds the Knox et al. regionalized voxel model after excluding
experiments identified during quality control. It uses hyperparameters from
previous cross-validation and outputs regionalized connectivity metrics.
"""
from __future__ import division
import os
import logging

import numpy as np
import pandas as pd
import sys

sys.path.insert(0, "./mouse_connectivity_models/paper/figures/model_comparison/")

import allensdk.core.json_utilities as ju

from mcmodels.core import VoxelModelCache, Mask
from mcmodels.models.voxel import RegionalizedModel
from mcmodels.regressors.nonparametric.kernels import Polynomial
from mcmodels.utils import padded_diagonal_fill

from helpers.model_data import ModelData
from helpers.error import VoxelModelError
from helpers.utils import get_structure_id, get_ordered_summary_structures

# Path definitions
EXP_INPUT_DIR = './mouse_connectivity_models/paper/'
INPUT_JSON = os.path.join(EXP_INPUT_DIR, 'input.json')
OUTPUT_DIR = './mouse_connectivity_models/paper/figures/model_comparison/output/'
TOP_DIR = EXP_INPUT_DIR
LOG = False  # Use standard model (not log-transformed)


def fit_structure(cache, structure_id, experiments_exclude, kernel_params,
                  model_option='standard'):
    """Fit a voxel connectivity model for a single brain structure.
    
    Args:
        cache: VoxelModelCache instance
        structure_id: Allen CCF structure identifier
        experiments_exclude: List of experiment IDs to exclude
        kernel_params: Dictionary of kernel hyperparameters
        model_option: Model type ('standard' or 'log')
        
    Returns:
        tuple: (ModelData instance, fitted RegionalizedModel)
    """
    data = ModelData(cache, structure_id).get_voxel_data(
        experiments_exclude=experiments_exclude)

    # Nested cross validation
    nw_kwargs = {}
    if 'shape' in kernel_params:
        nw_kwargs['kernel'] = Polynomial(**kernel_params)
    else:
        nw_kwargs['kernel'] = 'rbf'
        nw_kwargs['gamma'] = kernel_params.pop('gamma')

    error = VoxelModelError(cache, data)
    return data, error.fit(**nw_kwargs, option=model_option)


def main():
    """Main function to build regionalized voxel model with excluded experiments."""
    # Validate command line arguments
    if len(sys.argv) < 3:
        print("Usage: python build_model_new_excluded.py <experiments_exclude.json> <output_suffix>",
              file=sys.stderr)
        print("\nExample: python build_model_new_excluded.py experiments_exclude.json rebuilt_oh_211_regions",
              file=sys.stderr)
        sys.exit(1)
    
    input_data = ju.read(INPUT_JSON)

    structures = input_data.get('structures')
    manifest_file = input_data.get('manifest_file')
    manifest_file = os.path.join(TOP_DIR, manifest_file)

    log_level = input_data.get('log_level', logging.DEBUG)
    logging.getLogger().setLevel(log_level)

    exp_exc_path = sys.argv[1]
    EXPERIMENTS_EXCLUDE_JSON = os.path.join("./mouse_connectivity_models/paper/", exp_exc_path)
    
    # Validate experiments exclude file exists
    if not os.path.exists(EXPERIMENTS_EXCLUDE_JSON):
        logging.error(f"Experiments exclude file not found: {EXPERIMENTS_EXCLUDE_JSON}")
        sys.exit(1)
        
    experiments_exclude = ju.read(EXPERIMENTS_EXCLUDE_JSON)

    # Load hyperparameter dict
    suffix = 'log' if LOG else 'standard'
    hyperparameter_json = os.path.join(OUTPUT_DIR, f'hyperparameters-{suffix}.json')
    
    if not os.path.exists(hyperparameter_json):
        logging.error(f"Hyperparameters file not found: {hyperparameter_json}")
        sys.exit(1)
        
    hyperparameters = ju.read(hyperparameter_json)

    # Get caching object
    cache = VoxelModelCache(manifest_file=manifest_file)

    # Get structure ids
    structure_ids = [get_structure_id(cache, s) for s in structures]

    # Mask for reordering source
    annotation = cache.get_annotation_volume()[0]
    cumm_source_mask = np.zeros(annotation.shape, dtype=np.int)

    offset = 1  # Start @ 1 so that nonzero can be used
    weights, nodes = [], []
    for sid, sac in zip(structure_ids, structures):
        logging.debug("Building model for structure: %s", sac)

        data, reg = fit_structure(cache, sid, experiments_exclude,
                                  hyperparameters[sac], model_option=suffix)

        w = reg.get_weights(data.injection_mask.coordinates)

        # Assign ordering to full source
        ordering = np.arange(offset, w.shape[0] + offset, dtype=np.int)
        offset += w.shape[0]

        # Get source mask
        data.injection_mask.fill_volume_where_masked(cumm_source_mask, ordering)

        # Append to list
        weights.append(w)
        nodes.append(reg.nodes)

    # Stack
    weights = padded_diagonal_fill(weights)
    nodes = np.vstack(nodes)

    # Need to reorder weights (subtract 1 to get proper index)
    permutation = cumm_source_mask[cumm_source_mask.nonzero()] - 1
    weights = weights[permutation, :]

    # Regionalized
    logging.debug('Regionalizing voxel weights')
    ontological_order = get_ordered_summary_structures(cache)

    outfile_suffix = sys.argv[2]
    if outfile_suffix in ("original_oh_211_regions", "rebuilt_oh_211_regions"):
        oh_regions_file = "../preprocessed/allen_template_inputs/oh_connectome_rgn_numbers_ccfv3.txt"
        if not os.path.exists(oh_regions_file):
            logging.error(f"Oh connectome regions file not found: {oh_regions_file}")
            sys.exit(1)
        ontological_order = np.loadtxt(oh_regions_file)
    
    source_mask = Mask.from_cache(cache, structure_ids=structure_ids, hemisphere_id=2)
    source_key = source_mask.get_key(structure_ids=ontological_order)
    ipsi_key = data.projection_mask.get_key(structure_ids=ontological_order, hemisphere_id=2)
    contra_key = data.projection_mask.get_key(structure_ids=ontological_order, hemisphere_id=1)
    ipsi_model = RegionalizedModel(weights, nodes, source_key, ipsi_key,
                                   ordering=ontological_order, dataframe=True)
    contra_model = RegionalizedModel(weights, nodes, source_key, contra_key,
                                     ordering=ontological_order, dataframe=True)
    get_metric = lambda s: pd.concat((getattr(ipsi_model, s), getattr(contra_model, s)),
                                     keys=('ipsi', 'contra'), axis=1)

    # Write results
    output_dir = os.path.join(TOP_DIR, f'connectivity/voxel-{suffix}-model')
    os.makedirs(output_dir, exist_ok=True)

    # Regionalized
    logging.debug('Saving to directory: %s', output_dir)
    get_metric('connection_density').to_csv(
        os.path.join(output_dir, f'connection_density_{outfile_suffix}.csv'))
    get_metric('connection_strength').to_csv(
        os.path.join(output_dir, f'connection_strength_{outfile_suffix}.csv'))
    get_metric('normalized_connection_density').to_csv(
        os.path.join(output_dir, f'normalized_connection_density_{outfile_suffix}.csv'))
    get_metric('normalized_connection_strength').to_csv(
        os.path.join(output_dir, f'normalized_connection_strength_{outfile_suffix}.csv'))

    # Voxel
    ju.write(os.path.join(output_dir, 'target_mask_params.json'),
             dict(structure_ids=structure_ids, hemisphere_id=3))
    ju.write(os.path.join(output_dir, 'source_mask_params.json'),
             dict(structure_ids=structure_ids, hemisphere_id=2))
    np.savetxt(os.path.join(output_dir, f'weights_{outfile_suffix}.csv.gz'),
               weights.astype(np.float32), delimiter=',')
    np.savetxt(os.path.join(output_dir, f'nodes_{outfile_suffix}.csv.gz'),
               nodes.astype(np.float32), delimiter=',')
    
    logging.info("Model building completed successfully")


if __name__ == "__main__":
    main()
