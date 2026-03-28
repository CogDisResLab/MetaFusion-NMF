import argparse
import numpy as np
from sklearn.decomposition import NMF
import os
import sys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help='Merged fused matrix (.npy)')
    parser.add_argument('--k', type=int, default=3, help='Number of clusters')
    parser.add_argument('--output_w', required=True, help='Output W matrix (.npy)')
    parser.add_argument('--output_h', required=True, help='Output H matrix (.npy)')
    args = parser.parse_args()

    try:
        # 1. Load merged matrix
        X = np.load(args.input)

        if X.shape[0] != X.shape[1]:
            raise ValueError("Input fused matrix must be square")

        print(f"Loaded merged matrix: {X.shape}")
        print(f"Running NMF with k={args.k}")

        # 2. Run NMF
        model = NMF(
            n_components=args.k,
            init='nndsvd',
            random_state=42,
            max_iter=1000
        )

        W = model.fit_transform(X)
        H = model.components_

        print(f"W shape: {W.shape}")
        print(f"H shape: {H.shape}")

        # 3. Save outputs
        os.makedirs(os.path.dirname(args.output_w), exist_ok=True)
        np.save(args.output_w, W)

        os.makedirs(os.path.dirname(args.output_h), exist_ok=True)
        np.save(args.output_h, H)

        print(f"Saved W → {args.output_w}")
        print(f"Saved H → {args.output_h}")

    except Exception as e:
        print(f"NMF failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

