import argparse
import pandas as pd
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputs', nargs='+', required=True, help='List of top_features_best_k.csv files')
    parser.add_argument('--output', required=True, help='Path to meta_comparison.csv')
    args = parser.parse_args()

    # Load all cohort results
    all_dfs = []
    for f in args.inputs:
        # Extract cohort name from path (e.g., results/Cohort_A/...)
        cohort_name = f.split('/')[-2] 
        df = pd.read_csv(f)
        df['cohort'] = cohort_name
        all_dfs.append(df)

    combined = pd.concat(all_dfs)

    # Find features that appear in the SAME cluster across MULTIPLE cohorts
    # We group by cluster and feature, then count how many cohorts they appear in
    summary = combined.groupby(['cluster', 'feature']).agg({
        'correlation': 'mean',
        'cohort': 'count',
        'source': lambda x: ', '.join(x.unique())
    }).rename(columns={'cohort': 'cohort_count', 'correlation': 'avg_correlation'})

    # Filter for features found in more than 1 cohort
    conserved_drivers = summary[summary['cohort_count'] > 1].sort_values(
        by=['cluster', 'avg_correlation'], ascending=[True, False]
    )

    conserved_drivers.to_csv(args.output)
    print(f"Meta-comparison complete. Found {len(conserved_drivers)} conserved drivers.")

if __name__ == "__main__":
    main()
