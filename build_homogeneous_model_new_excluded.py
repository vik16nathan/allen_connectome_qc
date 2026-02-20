"""Build homogeneous connectivity model with updated experiment exclusions.

Fits a homogeneous (Oh et al., 2014-style) connectivity model with custom
region lists and experiment exclusions, then saves results as CSV.
"""

from __future__ import division
import argparse
import os
import logging

import numpy as np
import sys
import pandas as pd
sys.path.insert(0,"./mouse_connectivity_models/paper/figures/model_comparison/")
import allensdk.core.json_utilities as ju

from mcmodels.core import VoxelModelCache
from mcmodels.models import HomogeneousModel
from mcmodels.utils import nonzero_unique

####model_data_updated accounts for Oh regions

from model_data_updated import ModelData
from helpers.utils import get_structure_id, get_ordered_summary_structures


exp_input_dir='./mouse_connectivity_models/paper/'
INPUT_JSON = os.path.join(exp_input_dir, 'input.json')
TOP_DIR=exp_input_dir
ROOT_ID = 997
HIGH_RES = False
OUTPUT_DIR='./mouse_connectivity_models/paper/figures/model_comparison/output/'
THRESHOLD_INJECTION = True

def get_summary_structure_ids(rgn_list_path):
    """Load region structure IDs from a text file.

    Args:
        rgn_list_path: Path to a headerless CSV file with region IDs.

    Returns:
        Series of structure IDs.
    """
    structures = pd.read_csv(rgn_list_path, header=None).loc[:,0]
    return structures


def fit(cache, rgn_list_path, eid_set=None, experiments_exclude=None, high_res=False, threshold_injection=True):
    """Fit the homogeneous connectivity model.

    Args:
        cache: VoxelModelCache instance for data access.
        rgn_list_path: Path to file listing region numbers.
        eid_set: Optional set of experiment IDs to restrict to.
        experiments_exclude: List of experiment IDs to exclude.
        high_res: Whether to use high-resolution data.
        threshold_injection: Whether to threshold injection volumes.

    Returns:
        DataFrame of ipsi/contra connectivity weights.
    """
    if experiments_exclude is None:
        experiments_exclude = []
    logging.debug('getting data')
    ipsi_data = ModelData(cache, ROOT_ID).get_regional_data(rgn_list_path,
        eid_set=eid_set, experiments_exclude=experiments_exclude, high_res=high_res,
        threshold_injection=threshold_injection, projection_hemisphere_id=2)

    contra_data = ModelData(cache, ROOT_ID).get_regional_data(rgn_list_path,
        eid_set=eid_set, experiments_exclude=experiments_exclude, high_res=high_res,
        projection_hemisphere_id=1, threshold_injection=threshold_injection)

    X = ipsi_data.injections
    y = np.hstack((ipsi_data.projections, contra_data.projections))

    logging.debug('fitting')
    reg = HomogeneousModel(kappa=np.inf)
    reg.fit(X, y)

    # get ids
    ss_ids = get_summary_structure_ids(rgn_list_path)
    injection_key = ipsi_data.injection_mask.get_key(structure_ids=ss_ids, hemisphere_id=2)
    ipsi_key = ipsi_data.projection_mask.get_key(structure_ids=ss_ids, hemisphere_id=2)
    contra_key = contra_data.projection_mask.get_key(structure_ids=ss_ids, hemisphere_id=1)

    injection_regions = nonzero_unique(injection_key)
    ipsi_regions = nonzero_unique(ipsi_key)
    contra_regions = nonzero_unique(contra_key)

    ipsi_w = pd.DataFrame(data=reg.weights[:, :len(ipsi_regions)],
                          index=injection_regions,
                          columns=ipsi_regions)
    contra_w = pd.DataFrame(data=reg.weights[:, len(ipsi_regions):],
                            index=injection_regions,
                            columns=contra_regions)

    return pd.concat((ipsi_w, contra_w), keys=('ipsi', 'contra'), axis=1)


def parse_args():
    """Parse command-line arguments.

    Returns:
        Parsed argument namespace.
    """
    parser = argparse.ArgumentParser(
        description='Build homogeneous connectivity model with updated experiment exclusions.')
    parser.add_argument('excluded_exps_file', type=str,
                        help='Path to JSON file listing experiments to exclude.')
    parser.add_argument('rgn_list_path', type=str,
                        help='Path to file listing region numbers.')
    parser.add_argument('suffix', type=str,
                        help='Output file suffix (e.g., "original", "rebuilt").')
    return parser.parse_args()


def main():
    """Build and save the homogeneous connectivity model."""
    args = parse_args()

    input_data = ju.read(INPUT_JSON)

    manifest_file = input_data.get('manifest_file')
    manifest_file = os.path.join(TOP_DIR, manifest_file)

    log_level = input_data.get('log_level', logging.DEBUG)
    logging.getLogger().setLevel(log_level)

    # experiments to exclude
    ###CHANGED FROM ORIGINAL LIST TO ACCOUNT FOR QC FAILURE MODES
    EXPERIMENTS_EXCLUDE_JSON = os.path.join("./mouse_connectivity_models/paper/", args.excluded_exps_file)
    experiments_exclude = ju.read(EXPERIMENTS_EXCLUDE_JSON)

    # load hyperparameter dict
    suffix = 'high_res' if HIGH_RES else 'standard'

    # get caching object
    cache = VoxelModelCache(manifest_file=manifest_file)

    fit_kwargs = dict(high_res=HIGH_RES, threshold_injection=THRESHOLD_INJECTION,
                      experiments_exclude=experiments_exclude)
    
    model = fit(cache, args.rgn_list_path, **fit_kwargs)

    # write results
    logging.debug('saving')
    output_file = os.path.join(
        OUTPUT_DIR,
        'homogeneous-%s-model_%s.csv' % (suffix, args.suffix)
     )


    model.to_csv(output_file)


if __name__ == "__main__":
    main()
