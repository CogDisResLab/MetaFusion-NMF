import argparse
import numpy as np
import pandas as pd
from sklearn.decomposition import NMF
import sys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--summary', required=True, help='Path to meta_best_rank_summary.txt')
    parser.add_argument('--fused', required=True, help='Path to cohort-specific fused_network.npy')
    parser.add_argument('--omics', nargs='+', required=True, help='List of cohort-specific omics files')
    parser.add_argument('--output', required=True, help='Path to save top features CSV')
    args = parser.parse_args()

    try:
        # 1. Get the best K from the Meta-Summary
        # Format in file: "Meta-Best Rank (k): 3"
        best_k = None
        with open(args.summary, 'r') as f:
            for line in f:
                if "Meta-Best Rank (k):" in line:
                    best_k = int(line.split(':')[-1].strip())
                    break
        
        if best_k is None:
            raise ValueError(f"Could not find 'Meta-Best Rank (k)' in {args.summary}")
            
        print(f"Projecting for Meta-Best k={best_k}")

        # 2. Re-run NMF on the cohort-specific fused matrix to get W (Samples x Clusters)
        fused_data = np.load(args.fused)
        model = NMF(n_components=best_k, init='nndsvd', random_state=42, max_iter=1000)
        W = model.fit_transform(fused_data) 

        # 3. Align Samples (Logic must match snf_fusion.py)
        # Re-identifying common samples to ensure indexing of W matches the CSVs
        sample_sets = []
        dataframes = []
        for file_path in args.omics:
            df = pd.read_csv(file_path, index_col=0)
            dataframes.append(df)
            sample_sets.append(set(df.index))
        
        common_samples = sorted(list(set.intersection(*sample_sets)))
        
        # 4. Calculate Correlations per feature
        all_results = []
        for i, df in enumerate(dataframes):
            file_path = args.omics[i]
            # Align this specific cohort's dataframe to its own common samples
            df_aligned = df.loc[common_samples]
            
            for cluster_idx in range(best_k):
                # Membership weights for this specific cluster from the fused matrix
                weights = W[:, cluster_idx]
                
                # Correlate every feature (gene/probe) with the cluster membership weights
                cluster_series = pd.Series(weights, index=common_samples)
                corrs = df_aligned.apply(lambda col: col.corr(cluster_series))
                
                # Extract top 50 biomarkers for this cluster in this cohort
                top_features = corrs.sort_values(ascending=False).head(50)
                for feature, score in top_features.items():
                    all_results.append({
                        'cluster': cluster_idx + 1,
                        'feature': feature,
                        'correlation': score,
                        'source': file_path
                    })

        pd.DataFrame(all_results).to_csv(args.output, index=False)
        print(f"Successfully saved top features to {args.output}")

    except Exception as e:
        print(f"Error in feature projection: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
