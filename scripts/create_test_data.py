import pandas as pd
import numpy as np
import os

# ---------------------------
# Config
# ---------------------------
cohorts = ["cohort_a", "cohort_b"]
n_samples = 50

# HGNC and miRNA pools
m_genes = ["TP53", "EGFR", "MYC", "BRCA1", "BRCA2", "PTEN", "MTOR", "STAT3", "VEGFA", "GAPDH"]

def generate_hgnc_list(prefix, count):
    return m_genes + [f"{prefix}{i}" for i in range(count - len(m_genes))]

omics_config = {
    "expression.tsv": 1000,
    "methylation.tsv": 500,
    "mirna.tsv": 200
}

# ---------------------------
# Generate synthetic data with signal
# ---------------------------
for cohort_idx, cohort in enumerate(cohorts, start=1):
    cohort_dir = f"data/{cohort}"
    os.makedirs(cohort_dir, exist_ok=True)

    # Unique Sample IDs per cohort
    patients = [f"{cohort.upper()}_Patient_{i:03d}" for i in range(1, n_samples+1)]
    
    for filename, n_features in omics_config.items():
        if filename == "mirna.tsv":
            features = [f"hsa-mir-{i}" for i in range(1, n_features+1)]
        else:
            features = generate_hgnc_list(f"{cohort.upper()}_GENE", n_features)
        
        # 1. Start with baseline noise
        data = np.random.rand(n_features, n_samples) * 2 

        # 2. Inject "Clusters" (Biological Signal)
        # We split patients into two sub-groups (e.g., Subtype 1 and Subtype 2)
        # Group 1 (first 25 patients): High expression in features 0-50
        data[0:50, 0:25] += np.random.normal(loc=8, scale=1.5, size=(50, 25))
        
        # Group 2 (remaining 25 patients): High expression in features 50-100
        data[50:100, 25:] += np.random.normal(loc=12, scale=2.0, size=(50, 25))

        # 3. Add a Cohort-specific bias (so A looks different from B)
        if cohort == "cohort_a":
            data[100:150, :] += 5
        else:
            data[150:200, :] += 7

        # Ensure no negative values (NMF requirement)
        data = np.clip(data, 0, None)

        df = pd.DataFrame(data, index=features, columns=patients)

        # Save as TSV
        df.to_csv(f"{cohort_dir}/{filename}", sep="\t")
        print(f"Created {cohort_dir}/{filename} with signal ({df.shape[0]}x{df.shape[1]})")

print("\n✅ Synthetic data with latent patterns ready! Your NMF dendrogram should now show clear branching.")
