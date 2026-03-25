import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputs', nargs='+', required=True, help='List of top_features_best_k.csv files')
    parser.add_argument('--output', required=True, help='Path to save heatmap (e.g., results/meta_heatmap.png)')
    args = parser.parse_args()

    # Load data and pivot to create "Signatures"
    cohort_data = {}
    for f in args.inputs:
        cohort_name = f.split('/')[-2]
        df = pd.read_csv(f)
        # Create a matrix of [Features x Clusters] for this cohort
        pivot_df = df.pivot_table(index='feature', columns='cluster', values='correlation').fillna(0)
        cohort_data[cohort_name] = pivot_df

    cohort_names = list(cohort_data.keys())
    if len(cohort_names) < 2:
        print("Error: Need at least 2 cohorts to compare.")
        return

    # Find common features across both cohorts to compare apples to apples
    c1, c2 = cohort_names[0], cohort_names[1]
    common_features = cohort_data[c1].index.intersection(cohort_data[c2].index)
    
    if len(common_features) == 0:
        print("Error: No common features found between cohorts. Check gene naming (Human vs Mouse?)")
        return

    # Calculate Correlation between Cluster Signatures
    sig1 = cohort_data[c1].loc[common_features]
    sig2 = cohort_data[c2].loc[common_features]
    
    # Correlation matrix between columns of sig1 and sig2
    corr_matrix = sig1.corrwith(sig2, axis=0, drop=True) # Simple approach
    # More robust: Pearson correlation between every cluster pair
    final_corr = pd.DataFrame(index=sig1.columns, columns=sig2.columns)
    for col1 in sig1.columns:
        for col2 in sig2.columns:
            final_corr.loc[col1, col2] = sig1[col1].corr(sig2[col2])

    # Plotting
    plt.figure(figsize=(10, 8))
    sns.heatmap(final_corr.astype(float), annot=True, cmap='RdBu_r', center=0)
    plt.title(f'Cluster Similarity: {c1} (Rows) vs {c2} (Cols)')
    plt.xlabel(f'{c2} Clusters')
    plt.ylabel(f'{c1} Clusters')
    plt.tight_layout()
    plt.savefig(args.output)
    print(f"Heatmap saved to {args.output}")

if __name__ == "__main__":
    main()
