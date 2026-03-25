import argparse
import pandas as pd
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--summary', required=True)
    parser.add_argument('--drivers', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    # Read the Meta-Best K
    with open(args.summary, 'r') as f:
        best_k_line = f.readline().strip()
        ccc_line = f.readline().strip()

    # Read the Conserved Drivers
    drivers_df = pd.read_csv(args.drivers)

    with open(args.output, 'w') as f:
        f.write("# Meta-SNF-NMF Integration Report\n")
        f.write(f"Generated on: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}\n\n")
        
        f.write("## 1. Mathematical Stability\n")
        f.write(f"* **{best_k_line}**\n")
        f.write(f"* **{ccc_line}**\n\n")
        
        f.write("## 2. Conserved Biological Drivers (Top 10 per Cluster)\n")
        f.write("These features were found to be top-ranked in BOTH cohorts.\n\n")
        
        for cluster in sorted(drivers_df['cluster'].unique()):
            f.write(f"### Cluster {cluster}\n")
            cluster_top = drivers_df[drivers_df['cluster'] == cluster].head(10)
            if not cluster_top.empty:
                f.write("| Feature | Avg Correlation | Source Layers |\n")
                f.write("| :--- | :--- | :--- |\n")
                for _, row in cluster_top.iterrows():
                    f.write(f"| {row['feature']} | {row['avg_correlation']:.4f} | {row['source']} |\n")
            else:
                f.write("*No conserved drivers found for this cluster across all cohorts.*\n")
            f.write("\n")

    print(f"Executive summary generated at {args.output}")

if __name__ == "__main__":
    main()
