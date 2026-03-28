import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--nmf_w', required=True, help='Path to global NMF W matrix (.npy)')
    parser.add_argument('--sample_index', required=True, help='CSV mapping samples to cohort')
    parser.add_argument('--output', required=True, help='Path to save clustered heatmap')
    args = parser.parse_args()

    # 1. Load data
    W = np.load(args.nmf_w)  # shape: (n_samples, n_clusters)
    sample_index = pd.read_csv(args.sample_index, index_col=0)

    if 'cohort' not in sample_index.columns:
        raise ValueError("Sample index CSV must contain a 'cohort' column")

    if len(sample_index) != W.shape[0]:
        raise ValueError(f"W shape {W.shape[0]} doesn't match index length {len(sample_index)}")

    # 2. Create DataFrame directly using the original index to ensure alignment
    W_df = pd.DataFrame(
        W,
        index=sample_index.index,
        columns=[f"Cluster_{i+1}" for i in range(W.shape[1])]
    )

    # 3. Map Colors to the index
    # We create the color series based on the same index as W_df
    unique_cohorts = sample_index['cohort'].unique()
    palette = sns.color_palette("Set2", n_colors=len(unique_cohorts))
    cohort_colors_map = {c: palette[i] for i, c in enumerate(unique_cohorts)}
    
    # This series is indexed identically to W_df
    row_colors = sample_index['cohort'].map(cohort_colors_map)

    # 4. Plot clustered heatmap
    # row_cluster=True will mathematically group samples. 
    # Seaborn will automatically reorder 'row_colors' to stay aligned with the rows.
    g = sns.clustermap(
        W_df,
        row_cluster=True,
        col_cluster=True,
        row_colors=row_colors,
        cmap="viridis",
        figsize=(12, 10),
        linewidths=0, # Thinner lines often look cleaner with many samples
        method='average',
        metric='euclidean',
        # z_score=0  # Uncomment this to normalize clusters across rows
    )

    # 5. Clean up labels
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=8)

    # 6. Legend for cohorts
    handles = [Patch(facecolor=cohort_colors_map[c], edgecolor='k', label=c) for c in unique_cohorts]
    g.ax_heatmap.legend(
        handles=handles,
        title="Cohort",
        bbox_to_anchor=(1.2, 1),
        loc='upper left'
    )

    g.fig.suptitle("Global NMF: Samples × Clusters (Clustered)", y=1.02)

    # Save
    g.savefig(args.output, bbox_inches='tight')
    print(f"Clustered heatmap saved to {args.output}")

if __name__ == "__main__":
    main()
