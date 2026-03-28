import argparse
import numpy as np
import os
import pandas as pd
import sys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputs', nargs='+', required=True,
                        help='List of fused .npy files')
    parser.add_argument('--sample_ids', nargs='+', required=True,
                        help='List of sample ID .txt files (same order as inputs)')
    parser.add_argument('--output', required=True,
                        help='Output merged .npy file')
    parser.add_argument('--index_map', default="results/global/sample_index.csv",
                        help='Path to save sample index mapping (CSV)')
    args = parser.parse_args()

    if len(args.inputs) != len(args.sample_ids):
        print("Error: --inputs and --sample_ids must have same length")
        sys.exit(1)

    fused_matrices = []
    cohort_names = []
    sizes = []

    global_sample_records = []

    # 1. Load all fused matrices + sample IDs
    for f, sid_file in zip(args.inputs, args.sample_ids):
        try:
            mat = np.load(f)

            if mat.shape[0] != mat.shape[1]:
                raise ValueError(f"{f} is not square: {mat.shape}")

            # Load sample IDs
            sample_ids = pd.read_csv(sid_file, header=None)[0].tolist()

            if len(sample_ids) != mat.shape[0]:
                raise ValueError(
                    f"Sample ID count mismatch in {sid_file}: "
                    f"{len(sample_ids)} vs matrix {mat.shape[0]}"
                )

            fused_matrices.append(mat)

            # Extract cohort name from path
            cohort = os.path.basename(os.path.dirname(f))
            cohort_names.append(cohort)
            sizes.append(mat.shape[0])

            print(f"Loaded {f}: {mat.shape} ({len(sample_ids)} samples)")

            # Store temporarily
            global_sample_records.append((cohort, sample_ids))

        except Exception as e:
            print(f"Error loading {f}: {e}")
            sys.exit(1)

    # 2. Build block-diagonal matrix
    total_size = sum(sizes)
    merged = np.zeros((total_size, total_size))

    index_rows = []
    start = 0

    for (cohort, mat, size), (_, sample_ids) in zip(
        zip(cohort_names, fused_matrices, sizes),
        global_sample_records
    ):
        end = start + size
        merged[start:end, start:end] = mat

        print(f"{cohort}: rows {start}–{end-1}")

        # Record mapping: GLOBAL index → cohort + sample
        for i, sample in enumerate(sample_ids):
            index_rows.append({
                "global_index": start + i,
                "cohort": cohort,
                "sample": sample
            })

        start = end

    # 3. Save outputs
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    np.save(args.output, merged)

    index_df = pd.DataFrame(index_rows)
    os.makedirs(os.path.dirname(args.index_map), exist_ok=True)
    index_df.to_csv(args.index_map, index=False)

    print(f"\nSaved merged matrix: {args.output}")
    print(f"Saved sample index map: {args.index_map}")


if __name__ == "__main__":
    main()

