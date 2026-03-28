import argparse
import numpy as np
import pandas as pd
from sklearn.decomposition import NMF
from scipy.cluster.hierarchy import linkage, cophenet
from scipy.spatial.distance import pdist
import sys

def calculate_cophenetic(W):
    # Generate consensus-like distance from the Basis matrix W
    # High values in W mean strong membership
    z = linkage(W, method='average')
    c, coph_dists = cophenet(z, pdist(W))
    return c

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--rank', type=int, required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    try:
        data = np.load(args.input)
        
        # NMF Setup
        model = NMF(n_components=args.rank, init='nndsvd', random_state=42, max_iter=1000)
        W = model.fit_transform(data)
        
        # Calculate CCC
        ccc = calculate_cophenetic(W)
        print(f"Rank {args.rank} Cophenetic Coefficient: {ccc:.4f}")

        clusters = np.argmax(W, axis=1) + 1
        
        results = pd.DataFrame({
            'sample_index': range(len(clusters)),
            'cluster': clusters,
            'cophenetic_coeff': ccc,  # Store for Snakemake to read
            'rank': args.rank
        })
        
        results.to_csv(args.output, index=False)

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

