import argparse
import pandas as pd
import numpy as np
from snf import compute, snf
import sys
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputs', nargs='+', required=True, help='List of omics TSV files')
    parser.add_argument('--output', required=True, help='Path to save fused .npy matrix')
    parser.add_argument('--sample_ids', required=False, help='Path to save sample IDs (txt)')
    args = parser.parse_args()

    dataframes = []
    sample_sets = []

    # 1. Load data and track sample IDs
    for f in args.inputs:
        try:
            df = pd.read_csv(f, sep='\t', index_col=0)

            # Drop first column if non-numeric metadata
            if df.shape[1] > 1 and df.iloc[:, 0].dtype == object:
                df = df.iloc[:, 1:]

            # Transpose: rows = samples, columns = features
            df = df.T

            dataframes.append(df)
            sample_sets.append(set(df.index))
            print(f"Loaded {f}: {df.shape[0]} samples, {df.shape[1]} features.")

        except Exception as e:
            print(f"Error loading {f}: {e}")
            sys.exit(1)

    # 2. Find common samples
    common_samples = sorted(list(set.intersection(*sample_sets)))
    print(f"Aligning datasets: {len(common_samples)} common samples found.")

    if len(common_samples) == 0:
        print("Error: No overlapping samples found between input files.")
        sys.exit(1)

    # 3. Reindex datasets to common samples and extract values
    aligned_datasets = [df.loc[common_samples].values for df in dataframes]

    # 4. Run SNF
    try:
        # Convert each dataset to an affinity matrix
        print("Calculating affinity matrices...")
        affinity_matrices = [compute.make_affinity(ds, K=20, mu=0.5) for ds in aligned_datasets]

        # Fuse affinity matrices
        print("Fusing networks...")
        fused_network = snf(affinity_matrices, K=20, t=20)

        # 5. Save fused network
        os.makedirs(os.path.dirname(args.output), exist_ok=True)
        np.save(args.output, fused_network)
        print(f"Successfully saved fused network to {args.output}")

        # 6. Save sample IDs if requested
        if args.sample_ids:
            os.makedirs(os.path.dirname(args.sample_ids), exist_ok=True)
            with open(args.sample_ids, "w") as f:
                for s in common_samples:
                    f.write(s + "\n")
            print(f"Successfully saved sample IDs to {args.sample_ids}")

    except Exception as e:
        print(f"SNF Fusion failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()

