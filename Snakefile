configfile: "config.yaml"
import itertools

# Extract cohorts from the nested config structure
COHORTS = config["cohorts"]
RANKS = [2, 3, 4, 5]

# Generate all unique pairs of cohorts for validation
# e.g., if cohorts are [A, B, C], pairs are [(A,B), (A,C), (B,C)]
COHORT_PAIRS = list(itertools.combinations(COHORTS, 2))

rule all:
    input:
        # 1. New Gatekeeper: Validation reports for every pair
        expand("results/validation/{c1}_vs_{c2}_bridge.txt", 
               zip, 
               c1=[p[0] for p in COHORT_PAIRS], 
               c2=[p[1] for p in COHORT_PAIRS]),
        # 2. Fused networks for every cohort
        expand("results/{cohort}/fused_network.npy", cohort=COHORTS),
        # 3. NMF results for every cohort and every rank
        expand("results/{cohort}/nmf_results_k{k}.csv", cohort=COHORTS, k=RANKS),
        # 4. Final meta-summary and projections
        "results/meta_best_rank_summary.txt",
        expand("results/{cohort}/top_features_best_k.csv", cohort=COHORTS),
        "results/meta_conserved_drivers.csv",
        "results/meta_similarity_heatmap.png",
        "results/FINAL_META_REPORT.md"

# --- [NEW] Rule: Validates biological identity for every pair ---
rule check_biological_bridge:
    input:
        a = lambda wildcards: config["data"][wildcards.c1]["mrna"],
        b = lambda wildcards: config["data"][wildcards.c2]["mrna"]
    output:
        report = "results/validation/{c1}_vs_{c2}_bridge.txt"
    conda:
        "envs/environment.yaml"
    shell:
        "python scripts/check_alignment.py {input.a} {input.b} > {output.report}"

rule run_snf:
    input:
        # This forces the pipeline to validate EVERY bridge before starting math
        bridges = expand("results/validation/{c1}_vs_{c2}_bridge.txt", 
                         zip, 
                         c1=[p[0] for p in COHORT_PAIRS], 
                         c2=[p[1] for p in COHORT_PAIRS]),
        layers = lambda wildcards: config["data"][wildcards.cohort].values()
    output:
        fused = "results/{cohort}/fused_network.npy"
    log:
        "logs/{cohort}/snf_fusion.log"
    conda:
        "envs/environment.yaml"
    shell:
        "python scripts/snf_fusion.py "
        "--inputs {input.layers} "
        "--output {output.fused} > {log} 2>&1"

rule run_nmf:
    input:
        fused = "results/{cohort}/fused_network.npy"
    output:
        csv = "results/{cohort}/nmf_results_k{k}.csv"
    conda:
        "envs/environment.yaml"
    shell:
        "python scripts/nmf_decompose.py "
        "--input {input.fused} "
        "--rank {wildcards.k} "
        "--output {output.csv}"

rule select_meta_best_k:
    input:
        expand("results/{cohort}/nmf_results_k{k}.csv", cohort=COHORTS, k=RANKS)
    output:
        summary = "results/meta_best_rank_summary.txt"
    run:
        import pandas as pd
        all_data = []
        for f in input:
            df = pd.read_csv(f)
            all_data.append(df[['rank', 'cophenetic_coeff']])
        
        res_df = pd.concat(all_data)
        mean_ccc = res_df.groupby('rank')['cophenetic_coeff'].mean()
        best_k = mean_ccc.idxmax()
        
        with open(output.summary, "w") as f:
            f.write(f"Meta-Best Rank (k): {int(best_k)}\n")
            f.write(f"Average Cophenetic Correlation: {mean_ccc[best_k]:.4f}\n")

rule project_features:
    input:
        summary = "results/meta_best_rank_summary.txt",
        fused = "results/{cohort}/fused_network.npy",
        mrna = lambda wildcards: config["data"][wildcards.cohort]["mrna"],
        methy = lambda wildcards: config["data"][wildcards.cohort].get("methy", ""),
        mirna = lambda wildcards: config["data"][wildcards.cohort].get("mirna", "")
    output:
        top_features = "results/{cohort}/top_features_best_k.csv"
    log:
        "logs/{cohort}/project_features.log"
    conda:
        "envs/environment.yaml"
    shell:
        "python scripts/project_features.py "
        "--summary {input.summary} "
        "--fused {input.fused} "
        "--omics {input.mrna} {input.methy} {input.mirna} "
        "--output {output.top_features} > {log} 2>&1"

rule compare_meta_cohorts:
    input:
        expand("results/{cohort}/top_features_best_k.csv", cohort=COHORTS)
    output:
        csv = "results/meta_conserved_drivers.csv"
    conda:
        "envs/environment.yaml"
    shell:
        "python scripts/compare_meta_cohorts.py "
        "--inputs {input} "
        "--output {output.csv}"

rule plot_meta_similarity:
    input:
        expand("results/{cohort}/top_features_best_k.csv", cohort=COHORTS)
    output:
        png = "results/meta_similarity_heatmap.png"
    conda:
        "envs/environment.yaml"
    shell:
        "python scripts/plot_meta_correlation.py "
        "--inputs {input} "
        "--output {output.png}"

rule generate_meta_report:
    input:
        summary = "results/meta_best_rank_summary.txt",
        drivers = "results/meta_conserved_drivers.csv"
    output:
        report = "results/FINAL_META_REPORT.md"
    conda:
        "envs/environment.yaml"
    shell:
        "python scripts/generate_meta_report.py "
        "--summary {input.summary} "
        "--drivers {input.drivers} "
        "--output {output.report}"
