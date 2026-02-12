"""Build homogeneous connectivity model with excluded experiments.

This script rebuilds the homogeneous connectivity model (Oh et al., 2014) after
excluding experiments identified during quality control. Unlike the regionalized
voxel model, this uses a simpler homogeneous connectivity assumption.
"""
from __future__ import division
import os
import logging

import numpy as np
import sys
import pandas as pd

sys.path.insert(0, "./mouse_connectivity_models/paper/figures/model_comparison/")
import allensdk.core.json_utilities as ju

from mcmodels.core import VoxelModelCache
from mcmodels.models import HomogeneousModel
from mcmodels.utils import nonzero_unique

# model_data_updated accounts for Oh regions
from model_data_updated import ModelData
from helpers.utils import get_structure_id, get_ordered_summary_structures

# Configuration constants
EXP_INPUT_DIR = './mouse_connectivity_models/paper/'
INPUT_JSON = os.path.join(EXP_INPUT_DIR, 'input.json')
TOP_DIR = EXP_INPUT_DIR
ROOT_ID = 997  # Root structure ID for whole brain
HIGH_RES = False
OUTPUT_DIR = './mouse_connectivity_models/paper/figures/model_comparison/output/'
THRESHOLD_INJECTION = True


def get_summary_structure_ids(rgn_list_path):
    """Load summary structure IDs from region list file.
    
    Args:
        rgn_list_path: Path to CSV file containing region IDs
        
    Returns:
        pd.Series: Region IDs consistent with Oh et al., 2014
    """
    if not os.path.exists(rgn_list_path):
        raise FileNotFoundError(f"Region list file not found: {rgn_list_path}")
    
    structures = pd.read_csv(rgn_list_path, header=None).loc[:, 0]
    return structures


def fit(cache, rgn_list_path, eid_set=None, experiments_exclude=None,
        high_res=False, threshold_injection=True):
    """Fit homogeneous connectivity model.
    
    Args:
        cache: VoxelModelCache instance
        rgn_list_path: Path to region list file
        eid_set: Optional set of experiment IDs to include
        experiments_exclude: List of experiment IDs to exclude
        high_res: Use high resolution data if True
        threshold_injection: Apply injection thresholding if True
        
    Returns:
        pd.DataFrame: Fitted connectivity weights (ipsi and contra)
    """
    if experiments_exclude is None:
        experiments_exclude = []
    
    logging.debug('Getting ipsilateral data')
    ipsi_data = ModelData(cache, ROOT_ID).get_regional_data(
        rgn_list_path,
        eid_set=eid_set,
        experiments_exclude=experiments_exclude,
        high_res=high_res,
        threshold_injection=threshold_injection,
        projection_hemisphere_id=2
    )

    logging.debug('Getting contralateral data')
    contra_data = ModelData(cache, ROOT_ID).get_regional_data(
        rgn_list_path,
        eid_set=eid_set,
        experiments_exclude=experiments_exclude,
        high_res=high_res,
        projection_hemisphere_id=1,
        threshold_injection=threshold_injection
    )

    X = ipsi_data.injections
    y = np.hstack((ipsi_data.projections, contra_data.projections))

    logging.debug('Fitting homogeneous model')
    reg = HomogeneousModel(kappa=np.inf)
    reg.fit(X, y)

    # Get region IDs
    ss_ids = get_summary_structure_ids(rgn_list_path)
    injection_key = ipsi_data.injection_mask.get_key(structure_ids=ss_ids, hemisphere_id=2)
    ipsi_key = ipsi_data.projection_mask.get_key(structure_ids=ss_ids, hemisphere_id=2)
    contra_key = contra_data.projection_mask.get_key(structure_ids=ss_ids, hemisphere_id=1)

    injection_regions = nonzero_unique(injection_key)
    ipsi_regions = nonzero_unique(ipsi_key)
    contra_regions = nonzero_unique(contra_key)

    ipsi_w = pd.DataFrame(
        data=reg.weights[:, :len(ipsi_regions)],
        index=injection_regions,
        columns=ipsi_regions
    )
    contra_w = pd.DataFrame(
        data=reg.weights[:, len(ipsi_regions):],
        index=injection_regions,
        columns=contra_regions
    )

    return pd.concat((ipsi_w, contra_w), keys=('ipsi', 'contra'), axis=1)


def main():
    """Main function to build homogeneous model with excluded experiments."""
    # Validate command line arguments
    if len(sys.argv) < 4:
        print("Usage: python build_homogeneous_model_new_excluded.py <experiments_exclude.json> "
              "<region_list.csv> <output_suffix>",
              file=sys.stderr)
        print("\nExample: python build_homogeneous_model_new_excluded.py "
              "experiments_exclude.json oh_regions.txt rebuilt",
              file=sys.stderr)
        sys.exit(1)
    
    input_data = ju.read(INPUT_JSON)

    manifest_file = input_data.get('manifest_file')
    manifest_file = os.path.join(TOP_DIR, manifest_file)

    log_level = input_data.get('log_level', logging.DEBUG)
    logging.getLogger().setLevel(log_level)

    # Load experiments to exclude
    exp_exc_path = sys.argv[1]
    EXPERIMENTS_EXCLUDE_JSON = os.path.join("./mouse_connectivity_models/paper/", exp_exc_path)
    
    if not os.path.exists(EXPERIMENTS_EXCLUDE_JSON):
        logging.error(f"Experiments exclude file not found: {EXPERIMENTS_EXCLUDE_JSON}")
        sys.exit(1)
    
    experiments_exclude = ju.read(EXPERIMENTS_EXCLUDE_JSON)

    # Load resolution suffix
    suffix = 'high_res' if HIGH_RES else 'standard'

    # Get caching object
    cache = VoxelModelCache(manifest_file=manifest_file)

    fit_kwargs = dict(
        high_res=HIGH_RES,
        threshold_injection=THRESHOLD_INJECTION,
        experiments_exclude=experiments_exclude
    )

    # Load region list and fit model
    rgn_list_path = sys.argv[2]
    
    if not os.path.exists(rgn_list_path):
        logging.error(f"Region list file not found: {rgn_list_path}")
        sys.exit(1)
    
    logging.info(f"Building homogeneous model with region list: {rgn_list_path}")
    model = fit(cache, rgn_list_path, **fit_kwargs)

    # Write results
    logging.debug('Saving model to disk')
    outfile_suffix = sys.argv[3]
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    output_file = os.path.join(
        OUTPUT_DIR,
        f'homogeneous-{suffix}-model_{outfile_suffix}.csv'
    )

    model.to_csv(output_file)
    logging.info(f"Model saved to: {output_file}")


if __name__ == "__main__":
    main()
