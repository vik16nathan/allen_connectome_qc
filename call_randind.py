"""Calculate mutual information stability for community detection analysis.

This script computes adjusted mutual information scores between community 
assignments to identify stable community structures across different resolution
parameters (gammas) in Louvain community detection.
"""
import numpy as np
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import adjusted_mutual_info_score
import argparse
import scipy.io
import pandas as pd
import itertools
import os
import sys


def call_rand(ar1, ar2):
    """Calculate adjusted mutual information score between two community assignments.
    
    Args:
        ar1: First community assignment array
        ar2: Second community assignment array
        
    Returns:
        float: Adjusted mutual information score
    """
    return adjusted_mutual_info_score(ar1, ar2)


def parse_args():
    """Parse command line arguments.
    
    Returns:
        argparse.Namespace: Parsed arguments containing matfile, flag, and suffix
    """
    parser = argparse.ArgumentParser(
        description='Calculate mutual information for community stability analysis'
    )
    parser.add_argument('--matfile', type=str, required=True,
                        help='Path to MATLAB file containing community data')
    parser.add_argument('--flag', type=int, required=True, choices=[0, 1],
                        help='Processing mode: 0 for citemp, 1 for ciu')
    parser.add_argument('--suffix', type=str,
                        help='Suffix for output filename (required when flag=1)')
    args = parser.parse_args()
    
    # Validate arguments
    if args.flag == 1 and args.suffix is None:
        parser.error('--suffix is required when --flag=1')
    
    return args

if __name__ == '__main__':
    output_dir = '../derivatives/community_louvain/'
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    inputs = parse_args()
    
    # Validate input file exists
    if not os.path.exists(inputs.matfile):
        print(f"Error: Input file '{inputs.matfile}' not found", file=sys.stderr)
        sys.exit(1)
    
    try:
        mat = scipy.io.loadmat(inputs.matfile)
    except Exception as e:
        print(f"Error loading MATLAB file: {e}", file=sys.stderr)
        sys.exit(1)
    
    if inputs.flag == 0:
        # Process citemp data
        if 'citemp' not in mat:
            print("Error: 'citemp' key not found in MATLAB file", file=sys.stderr)
            sys.exit(1)
            
        communities = mat['citemp']
        means = []
        stds = []
        
        # Fixed: range should be communities.shape[1], not +1
        for col in range(2, communities.shape[1]):
            print(f'Processing column {col}')
            r = []
            
            for subset in itertools.combinations(range(col), 2):
                r.append(call_rand(communities[subset[0]], communities[subset[1]]))
            
            means.append(np.average(r))
            stds.append(np.std(r))
            print(f'Processed {len(means)} columns')
        
        output_file = os.path.join(output_dir, 'rand_output.mat')
        scipy.io.savemat(output_file, mdict={'means': means, 'stds': stds})
        print(f'Results saved to {output_file}')

    elif inputs.flag == 1:
        # Process ciu data
        if 'ciu' not in mat:
            print("Error: 'ciu' key not found in MATLAB file", file=sys.stderr)
            sys.exit(1)
            
        communities = mat['ciu']
        rands = {}
        
        for col in range(communities.shape[1]):
            r = []
            ar1 = communities[:, col]
            for col2 in range(communities.shape[1]):
                ar2 = communities[:, col2]
                r.append(call_rand(ar1, ar2))
            rands[f'{col}'] = r
        
        rands_df = pd.DataFrame(rands)
        output_file = os.path.join(output_dir, f'rand_output_{inputs.suffix}.mat')
        scipy.io.savemat(output_file, mdict={'a': rands_df.to_numpy()})
        print(f'Results saved to {output_file}')