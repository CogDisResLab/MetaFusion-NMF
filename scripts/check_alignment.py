# scripts/check_alignment_real.py
import pandas as pd
import sys
import mygene
import os

def validate_hgnc_real(gene_list, cohort_name):
    """Pings MyGene.info to verify if symbols are valid HGNC entities."""
    if not gene_list:
        return set(), set()

    mg = mygene.MyGeneInfo()
    print(f"[{cohort_name}] Verifying {len(gene_list)} symbols...")
    
    try:
        results = mg.querymany(gene_list, scopes='symbol', fields='symbol', species='human', verbose=False)
        found_genes = {item['query'] for item in results if 'notfound' not in item}
        invalid_genes = set(gene_list) - found_genes
        return found_genes, invalid_genes
    except Exception as e:
        print(f"Error connecting to MyGene.info for {cohort_name}: {e}")
        return set(), set(gene_list)

def main(files):
    cohort_data = {}
    
    # 1. Process each file dynamically
    for fpath in files:
        # Use filename as the cohort label
        label = os.path.basename(fpath).replace('.csv', '')
        df = pd.read_csv(fpath)
        
        # Identify gene columns
        genes = [str(c) for c in df.columns if not c.lower().startswith(('unnamed', 'sample', 'id'))]
        
        found, invalid = validate_hgnc_real(genes, label)
        cohort_data[label] = {
            "found": found,
            "invalid": invalid,
            "total_cols": len(genes)
        }

    # 2. Global Bridge Calculation (Intersection of ALL cohorts)
    all_valid_sets = [data["found"] for data in cohort_data.values()]
    global_bridge = set.intersection(*all_valid_sets) if all_valid_sets else set()
    
    # 3. Report Generation
    print(f"\n{'='*60}")
    print(f"   MULTI-COHORT HGNC VALIDATION & BRIDGE REPORT")
    print(f"{'='*60}")

    max_found_count = 0
    for label, data in cohort_data.items():
        count = len(data["found"])
        max_found_count = max(max_found_count, count)
        print(f"\nCohort: {label}")
        print(f"  - Total Features: {data['total_cols']}")
        print(f"  - Valid HGNC:    {count}")
        if data["invalid"] and count < 10:
             print(f"  - Invalid (5):   {list(data['invalid'])[:5]}")

    print(f"\n{'='*60}")
    print(f"FINAL BRIDGE ANALYSIS")
    print(f"  - Shared Genes (All Cohorts): {len(global_bridge)}")
    
    if max_found_count == 0:
        print("\nCRITICAL: No valid HGNC symbols found in any cohort.")
        sys.exit(1)
    
    # Coverage relative to the largest cohort
    coverage = (len(global_bridge) / max_found_count) * 100
    print(f"  - Global Bridge Coverage:    {coverage:.2f}%")
    
    if coverage < 20:
        print("\nCRITICAL: Bridge is too small (<20%). Integration will likely fail.")
        sys.exit(1) 
    else:
        print("\nSUCCESS: Biological bridge confirmed across all cohorts.")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python scripts/check_alignment.py data/c1.csv data/c2.csv ... data/cN.csv")
        sys.exit(1)
    
    # Pass all arguments from the command line starting from the first file
    main(sys.argv[1:])
