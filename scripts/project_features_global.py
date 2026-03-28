import argparse
import numpy as np
import pandas as pd
import sys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--summary', required=True,
                        help='Path to meta_best_rank_summary.txt')
    parser.add_argument('--global_w', required=True,
                        help='Path to global NMF W matrix (.npy)')
    parser.add_argument('--index_map', required=True,
                        help='Path to sample_index.csv')
    parser.add_argument('--cohort', required=True,
                        help='Cohort name (must match index_map)')
    parser.add_argument('--fused', required=True,
                        help='Path to cohort-specific fused_network.npy (for sanity check)')
    parser.add_argument('--omics', nargs='+', required=True,
                        help='List of cohort-specific omics TSV files')
    parser.add_argument('--output', required=True,
                        help='Path to save top features CSV')
    args = parser.parse_args()

    try:
        # ---------------------------
        # 1. Get best K
        # ---------------------------
        best_k = None
        with open(args.summary, 'r') as f:
            for line in f:
                if "Meta-Best Rank (k):" in line:
                    best_k = int(line.split(':')[-1].strip())
                    break

        if best_k is None:
            raise ValueError("Could not parse Meta-Best Rank (k)")

        print(f"Using global k={best_k}")

        # ---------------------------
        # 2. Load global W + index map
        # ---------------------------
        W_global = np.load(args.global_w)
        index_df = pd.read_csv(args.index_map)

        # Filter rows for this cohort
        cohort_df = index_df[index_df["cohort"] == args.cohort]

        if cohort_df.empty:
            raise ValueError(f"No samples found for cohort {args.cohort} in index map")

        # Extract indices + sample names
        global_indices = cohort_df["global_index"].values
        cohort_samples = cohort_df["sample"].values

        # Subset W
        W_cohort = W_global[global_indices, :]

        print(f"{args.cohort}: {W_cohort.shape[0]} samples, {W_cohort.shape[1]} clusters")

        # ---------------------------
        # 3. Load and align omics
        # ---------------------------
        sample_sets = []
        dataframes = []

        for file_path in args.omics:
            df = pd.read_csv(file_path, sep='\t', index_col=0)

            # Drop metadata column if present
            if df.shape[1] > 1 and df.iloc[:, 0].dtype == object:
                df = df.iloc[:, 1:]

            # Transpose → samples as rows
            df = df.T

            dataframes.append(df)
            sample_sets.append(set(df.index))

        # Find common samples across omics
        common_samples = sorted(list(set.intersection(*sample_sets)))

        if len(common_samples) == 0:
            raise ValueError("No overlapping samples across omics layers")

        # Intersect with cohort samples from global index
        common_samples = [s for s in common_samples if s in set(cohort_samples)]

        if len(common_samples) == 0:
            raise ValueError("No overlap between omics samples and global index map")

        print(f"Aligned samples: {len(common_samples)}")

        # ---------------------------
        # 4. Align W to omics sample order
        # ---------------------------
        cohort_df = cohort_df.set_index("sample")

        # reorder to match omics sample order
        aligned_indices = cohort_df.loc[common_samples]["global_index"].values
        W_aligned = W_global[aligned_indices, :]

        # ---------------------------
        # 5. Correlation analysis
        # ---------------------------
        all_results = []

        for i, df in enumerate(dataframes):
            file_path = args.omics[i]

            df_aligned = df.loc[common_samples]

            for cluster_idx in range(best_k):
                weights = W_aligned[:, cluster_idx]
                cluster_series = pd.Series(weights, index=common_samples)

                corrs = df_aligned.apply(lambda col: col.corr(cluster_series))

                top_features = corrs.sort_values(ascending=False).head(50)

                for feature, score in top_features.items():
                    all_results.append({
                        'cluster': cluster_idx + 1,
                        'feature': feature,
                        'correlation': score,
                        'source': file_path
                    })

        # ---------------------------
        # 6. Save
        # ---------------------------
        pd.DataFrame(all_results).to_csv(args.output, index=False)
        print(f"Saved: {args.output}")

    except Exception as e:
        print(f"Error in feature projection: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

