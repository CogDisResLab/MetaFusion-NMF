import pandas as pd
import numpy as np
import os

# Create data directory if it doesn't exist
os.makedirs("data", exist_ok=True)

# 1. Create Patient IDs (The "Source of Truth")
patients = [f"Patient_{i:03d}" for i in range(1, 51)]
with open("data/patient_ids.txt", "w") as f:
    f.write("\n".join(patients))

# 2. Define dimensions for our fake omics
# mrna: 1000 genes, methy: 500 sites, mirna: 200 miRNAs
omics_config = {
    "expression.csv": (50, 1000),
    "methylation.csv": (50, 500),
    "mirna.csv": (50, 200)
}

# 3. Generate and save the CSVs
for filename, shape in omics_config.items():
    # Generate random data between 0 and 10
    data = np.random.rand(shape[0], shape[1]) * 10
    
    # Create DataFrame with patients as rows
    df = pd.DataFrame(data, index=patients, columns=[f"Feature_{i}" for i in range(shape[1])])
    
    df.to_csv(f"data/{filename}")
    print(f"Created data/{filename} with shape {shape}")

print("\nSuccess! Your 'data/' folder is ready for the Snakemake pipeline.")
