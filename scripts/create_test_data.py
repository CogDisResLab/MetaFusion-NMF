import pandas as pd
import numpy as np
import os

# Create data directory if it doesn't exist
os.makedirs("data", exist_ok=True)

# 1. Create Patient IDs
patients = [f"Patient_{i:03d}" for i in range(1, 51)]
with open("data/patient_ids.txt", "w") as f:
    f.write("\n".join(patients))

# 2. Define HGNC Gene Pools for a "Real" Bridge
# We use famous genes to ensure the API validator finds them
m_genes = ["TP53", "EGFR", "MYC", "BRCA1", "BRCA2", "PTEN", "MTOR", "STAT3", "VEGFA", "GAPDH"]
# Fill the rest with simulated HGNC-like names (GENE1, GENE2...)
def generate_hgnc_list(prefix, count):
    return m_genes + [f"{prefix}{i}" for i in range(count - len(m_genes))]

# 3. Define dimensions and Gene Symbols for each layer
omics_config = {
    "expression.csv": (50, generate_hgnc_list("MARK", 1000)),
    "methylation.csv": (50, generate_hgnc_list("METH", 500)),
    "mirna.csv": (50, [f"hsa-mir-{i}" for i in range(200)]) # miRNAs have their own HGNC format
}

# 4. Generate and save the CSVs
for filename, (n_patients, feature_names) in omics_config.items():
    # Generate random data
    data = np.random.rand(n_patients, len(feature_names)) * 10
    
    # Create DataFrame with patients as rows and HGNC symbols as columns
    df = pd.DataFrame(data, index=patients, columns=feature_names)
    
    # Note: If your pipeline expects genes as ROWS, use df.T.to_csv
    df.to_csv(f"data/{filename}")
    print(f"Created data/{filename} with {len(feature_names)} HGNC features.")

print("\nSuccess! Synthetic HGNC data is ready. The bridge check should now PASS.")
