import argparse
import pandas as pd
import os

def sanitize_feature(name: str) -> str:
    """Replace or escape characters that break Markdown tables."""
    if pd.isna(name):
        return ""
    # Replace pipe with dash
    return f"`{name.replace('|', '-')}`"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--summary', required=True, help='Path to meta_best_rank_summary.txt')
    parser.add_argument('--drivers', required=True, help='Path to meta_conserved_drivers.csv')
    parser.add_argument('--output', required=True, help='Path to output Markdown report')
    args = parser.parse_args()

    # Read Meta-Best K summary
    with open(args.summary, 'r') as f:
        best_k_line = f.readline().strip()
        derived_line = f.readline().strip()

    # Read conserved drivers
    drivers_df = pd.read_csv(args.drivers)

    with open(args.output, 'w') as f:
        # Header
        f.write("# Meta-SNF-NMF Integration Report\n")
        f.write(f"Generated on: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}\n\n")
        
        # Section 1: Mathematical Stability
        f.write("## 1. Mathematical Stability\n")
        f.write(f"* **{best_k_line}**\n")
        f.write(f"* **{derived_line}**\n\n")
        
        # Section 2: Conserved Biological Drivers
        f.write("## 2. Conserved Biological Drivers (Top 10 per Cluster)\n")
        f.write("These features were found to be top-ranked in BOTH cohorts.\n\n")
        
        for cluster in sorted(drivers_df['cluster'].unique()):
            f.write(f"### Cluster {cluster}\n")
            
            # Sort by avg_correlation descending and take top 10
            cluster_top = drivers_df[drivers_df['cluster'] == cluster] \
                            .sort_values('avg_correlation', ascending=False) \
                            .head(10)
            
            if not cluster_top.empty:
                f.write("| Feature | Avg Correlation | Source Layers |\n")
                f.write("| :--- | :--- | :--- |\n")
                for _, row in cluster_top.iterrows():
                    feature = sanitize_feature(row['feature'])
                    avg_corr = f"{row['avg_correlation']:.4f}"
                    source = row['source']
                    f.write(f"| {feature} | {avg_corr} | {source} |\n")
            else:
                f.write("*No conserved drivers found for this cluster across all cohorts.*\n")
            f.write("\n")

    print(f"Executive summary generated at {args.output}")

if __name__ == "__main__":
    main()

