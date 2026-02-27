# There must be a file in the directory which has columns: 1. Distance and 2. Median MP, which should be entered as input argument.
# Base on median MP of each distance from TSS, this script will calculate where whole promoter starts & ends.
# You must input size of core promoter (bp) of the species
# Usage: python3 speed_promoter_delineations_calculator_keep_arg.py <Median_file> <Core_size>

import pandas as pd
import numpy as np
import sys

median_file = sys.argv[1]

# Position 2: Core_size (optional; can be empty or absent)
if len(sys.argv) < 3 or sys.argv[2] == '':
    core_size = 0
else:
    try:
        core_size = int(sys.argv[2])
    except ValueError:
        core_size = 0   # or exit with error message

# Position 3: Suffix (optional)
suffix = sys.argv[3] if len(sys.argv) > 3 else ""

# Read data once: position (col1), median methylation (col2)
df = pd.read_csv(median_file,
                 sep="\t", header=None, names=["pos", "median"])

# Precompute prefix sums
pos_means = df.groupby("pos")["median"].mean().sort_index()
positions = pos_means.index.values
prefix_sum = np.cumsum(pos_means.values)
prefix_count = np.cumsum(np.ones(len(pos_means)))

total_sum = prefix_sum[-1]
total_count = prefix_count[-1]

# Find best windows [i,j]
max_authenticity = -100
best_windows = []  # keep all (i,j) with current max

print("Computing optimal promoter window...", file=sys.stderr, flush=True)

base_name = "promoter_delineations"
out_path = f"{base_name}_{suffix}.tsv"
with open(out_path, "w") as out:
    for i in range(-10000, -1*core_size+1):
        for j in range(core_size, 10001):
            # O(1) window sum using prefix sums
            left_idx = np.searchsorted(positions, i)
            right_idx = np.searchsorted(positions, j, side="right") - 1

            if left_idx > right_idx:
                continue

            sum_prom = prefix_sum[right_idx] - (prefix_sum[left_idx - 1] if left_idx > 0 else 0)
            n_prom = prefix_count[right_idx] - (prefix_count[left_idx - 1] if left_idx > 0 else 0)
            avg_prom = sum_prom / n_prom if n_prom > 0 else 0

            sum_non = total_sum - sum_prom
            n_non = total_count - n_prom
            avg_non = sum_non / n_non if n_non > 0 else 0

            authenticity = avg_non - avg_prom

            # continuously show what is being calculated and current best
            print(
                f"i={i}, j={j}, authenticity={authenticity:.6f}; "
                f"current best={max_authenticity:.6f}",
                file=sys.stderr,
                flush=True,
            )

            if authenticity > max_authenticity:
                # new strictly better max: reset list and file
                max_authenticity = authenticity
                best_windows = [(i, j)]
                #out.seek(0)
                #out.truncate()
                out.write(f"{i}\t{j}\t{authenticity}\n")
                out.flush()
                # also show best on stdout
                print(
                    f"NEW BEST: [{i}, {j}] authenticity={authenticity:.6f}",
                    flush=True,
                )
            elif authenticity == max_authenticity:
                # tie: append to list and file
                best_windows.append((i, j))
                out.write(f"{i}\t{j}\t{authenticity}\n")
                out.flush()
                print(
                    f"TIE BEST: [{i}, {j}] authenticity={authenticity:.6f}",
                    flush=True,
                )

print(f"Finished. Max authenticity = {max_authenticity:.6f}")
print(f"Number of best windows found: {len(best_windows)}")
