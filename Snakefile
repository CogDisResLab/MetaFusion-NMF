configfile: "config.yaml"

COHORTS = config["cohorts"]

rule all:
    input:
        expand("results/{cohort}/fused_network.npy", cohort=COHORTS),
        expand("results/{cohort}/sample_ids.txt", cohort=COHORTS),
        "results/global/merged_fused.npy",
        "results/global/sample_index.csv",
        "results/global/global_nmf_W.npy",
        "results/global/global_nmf_H.npy",
        "results/meta_best_rank_summary.txt",
        expand("results/{cohort}/top_features_best_k.csv", cohort=COHORTS),
        "results/meta_conserved_drivers.csv",
        "results/meta_cohort_separation.png",
        "results/FINAL_META_REPORT.md"


# ---------------------------
# 1. SNF per cohort
# ---------------------------
rule run_snf:
    input:
        layers = lambda wc: list(config["data"][wc.cohort].values())
    output:
        fused = "results/{cohort}/fused_network.npy",
        samples = "results/{cohort}/sample_ids.txt"
    log:
        "logs/{cohort}/snf_fusion.log"
    conda:
        "envs/environment.yaml"
    shell:
        "python scripts/snf_fusion.py "
        "--inputs {input.layers} "
        "--output {output.fused} "
        "--sample_ids {output.samples} "
        "> {log} 2>&1"


# ---------------------------
# 2. Merge fused networks
# ---------------------------
rule merge_fused_networks:
    input:
        fused = expand("results/{cohort}/fused_network.npy", cohort=COHORTS),
        samples = expand("results/{cohort}/sample_ids.txt", cohort=COHORTS)
    output:
        merged = "results/global/merged_fused.npy",
        index = "results/global/sample_index.csv"
    conda:
        "envs/environment.yaml"
    shell:
        "python scripts/merge_fused_networks.py "
        "--inputs {input.fused} "
        "--sample_ids {input.samples} "
        "--output {output.merged} "
        "--index_map {output.index}"


# ---------------------------
# 3. Global NMF
# ---------------------------
rule run_global_nmf:
    input:
        fused = "results/global/merged_fused.npy"
    output:
        W = "results/global/global_nmf_W.npy",
        H = "results/global/global_nmf_H.npy"
    conda:
        "envs/environment.yaml"
    shell:
        "python scripts/nmf_global.py "
        "--input {input.fused} "
        "--output_w {output.W} "
        "--output_h {output.H}"


# ---------------------------
# 4. Select meta best k
# ---------------------------
rule select_meta_best_k:
    input:
        "results/global/global_nmf_W.npy"
    output:
        summary = "results/meta_best_rank_summary.txt"
    run:
        import numpy as np

        W = np.load(input[0])
        k = W.shape[1]

        with open(output.summary, "w") as f:
            f.write(f"Meta-Best Rank (k): {k}\n")
            f.write("Derived from global NMF\n")


# ---------------------------
# 5. Project features
# ---------------------------
rule project_features:
    input:
        summary = "results/meta_best_rank_summary.txt",
        W = "results/global/global_nmf_W.npy",
        index = "results/global/sample_index.csv",
        fused = "results/{cohort}/fused_network.npy",
        omics_layers = lambda wc: list(config["data"][wc.cohort].values())
    output:
        top_features = "results/{cohort}/top_features_best_k.csv"
    log:
        "logs/{cohort}/project_features.log"
    conda:
        "envs/environment.yaml"
    params:
        omics = lambda wc, input: ' '.join(input.omics_layers)
    shell:
        "python scripts/project_features_global.py "
        "--summary {input.summary} "
        "--global_w {input.W} "
        "--index_map {input.index} "
        "--cohort {wildcards.cohort} "
        "--fused {input.fused} "
        "--omics {params.omics} "
        "--output {output.top_features} > {log} 2>&1"


# ---------------------------
# 6. Compare cohorts
# ---------------------------
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


# ---------------------------
# 7. Plot similarity / cohort separation
# ---------------------------
rule plot_meta_similarity:
    input:
        W = "results/global/global_nmf_W.npy",
        index = "results/global/sample_index.csv"
    output:
        png = "results/meta_cohort_separation.png"
    conda:
        "envs/environment.yaml"
    shell:
        "python scripts/plot_meta_correlation.py "
        "--nmf_w {input.W} "
        "--sample_index {input.index} "
        "--output {output.png}"


# ---------------------------
# 8. Final report
# ---------------------------
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

