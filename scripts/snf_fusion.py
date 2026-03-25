import argparse
import pandas as pd
import numpy as np
from snf import compute, snf  # Added snf for the fusion step
import sys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputs', nargs='+', required=True, help='List of omics CSV files')
    parser.add_argument('--output', required=True, help='Path to save fused .npy matrix')
    args = parser.parse_args()

    dataframes = []
    sample_sets = []

    # 1. Load data and track index (patient IDs)
    for f in args.inputs:
        try:
            df = pd.read_csv(f, index_col=0)
            dataframes.append(df)
            sample_sets.append(set(df.index))
            print(f"Loaded {f}: {df.shape[0]} samples, {df.shape[1]} features.")
        except Exception as e:
            print(f"Error loading {f}: {e}")
            sys.exit(1)

    # 2. Find common samples across ALL layers
    common_samples = sorted(list(set.intersection(*sample_sets)))
    print(f"Aligning datasets: {len(common_samples)} common samples found.")

    if len(common_samples) == 0:
        print("Error: No overlapping samples found between input files.")
        sys.exit(1)

    # 3. Reindex and extract values (Ensures identical ordering)
    aligned_datasets = [df.loc[common_samples].values for df in dataframes]

    # 4. Run SNF
    try:
        # Step A: Convert each dataset into an affinity matrix
        # K = number of neighbors, mu = hyperparameter for affinity scaling
        print("Calculating affinity matrices...")
        affinity_matrices = [
            compute.make_affinity(ds, K=20, mu=0.5) 
            for ds in aligned_datasets
        ]
        
        # Step B: Fuse the affinity matrices into one similarity network
        # t = number of iterations for the diffusion process
        print("Fusing networks...")
        fused_network = snf(affinity_matrices, K=20, t=20)
        
        # 5. Save the fused similarity matrix
        np.save(args.output, fused_network)
        print(f"Successfully saved fused network to {args.output}")
        
    except Exception as e:
        print(f"SNF Fusion failed: {e}")
        # This will now give you a useful error if parameters are mismatched
        sys.exit(1)

if __name__ == "__main__":
    main()
