
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from cect.readers import White_whole

# ============================================================
# C. elegans connectome: coupling ablation experiment
# ============================================================

OUTPUT = Path("coupling_ablation_results.csv")

DT = 0.02
DURATION = 10.0
TAU = 0.5

CHEMICAL_STRENGTHS = [0.4, 0.8, 1.2]
GAP_STRENGTHS = [0.125, 0.25, 0.375]
DECAYS = [0.075, 0.15, 0.225]

STIMULUS_CELLS = ["ASHL", "ASHR"]
STIMULUS_DURATION = 1.0
STIMULUS_AMPLITUDE = 1.0

# Keep consistent with the previous experiments.
# Chemical edges are directed; gap junctions are undirected.
print("Loading White et al. 1986 connectome...")
ds = White_whole.get_instance()
cells, connections = ds.get_neuron_to_neuron_conns()

cell_names = list(cells)
index = {name: i for i, name in enumerate(cell_names)}
n = len(cell_names)

chemical = np.zeros((n, n), dtype=float)
gap = np.zeros((n, n), dtype=float)

for conn in connections:
    pre = conn.pre_cell
    post = conn.post_cell

    if pre not in index or post not in index:
        continue

    i, j = index[pre], index[post]
    weight = float(conn.number)

    if weight <= 0:
        continue

    if "Generic_CS" in str(conn.synclass):
        chemical[i, j] += weight
    elif "Generic_GJ" in str(conn.synclass):
        # Canonicalize reciprocal records to one undirected edge.
        # Use max rather than sum to avoid double-counting reciprocal entries.
        gap[i, j] = max(gap[i, j], weight)
        gap[j, i] = max(gap[j, i], weight)

def normalize_matrix(matrix):
    row_sums = matrix.sum(axis=1)
    result = matrix.copy()
    nonzero = row_sums > 0
    result[nonzero] /= row_sums[nonzero, None]
    return result

chemical_norm = normalize_matrix(chemical)
gap_norm = normalize_matrix(gap)

# Keep the same normalization convention as the previous model.
# Chemical input is postsynaptic: chemical_norm.T @ activity.
# Gap input is a neighbor-weighted average minus current activity.
def simulate(stimulated_cell, chemical_strength, gap_strength, decay,
             use_chemical=True, use_gap=True):
    if stimulated_cell not in index:
        raise ValueError(
            f"{stimulated_cell} is missing from the connectome."
        )

    steps = int(DURATION / DT)
    activity = np.zeros((steps, n), dtype=float)
    stim_idx = index[stimulated_cell]

    for t in range(1, steps):
        previous = activity[t - 1]

        chemical_input = (
            chemical_norm.T @ previous if use_chemical
            else np.zeros(n)
        )
        electrical_input = (
            gap_norm @ previous - previous if use_gap
            else np.zeros(n)
        )

        external_input = np.zeros(n)
        if t * DT < STIMULUS_DURATION:
            external_input[stim_idx] = STIMULUS_AMPLITUDE

        drive = (
            chemical_strength * chemical_input
            + gap_strength * electrical_input
            - decay * previous
            + external_input
        )

        target = np.maximum(drive, 0.0)
        activity[t] = previous + (DT / TAU) * (target - previous)
        activity[t] = np.clip(activity[t], 0.0, 1.0)

    return activity

# Ablation modes:
# full, chemical only, gap only, and neither recurrent pathway.
MODES = [
    ("full", True, True),
    ("chemical_only", True, False),
    ("gap_only", False, True),
    ("no_recurrent_coupling", False, False),
]

def integrated_profile(activity):
    # Integral of activity over time for each neuron.
    return activity.sum(axis=0) * DT

def cosine_similarity(a, b):
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    if denom == 0:
        return np.nan
    return float(np.dot(a, b) / denom)

rows = []
total_runs = (
    len(CHEMICAL_STRENGTHS)
    * len(GAP_STRENGTHS)
    * len(DECAYS)
    * len(MODES)
)
run_number = 0

for cs, gs, decay in itertools.product(
    CHEMICAL_STRENGTHS, GAP_STRENGTHS, DECAYS
):
    for mode, use_chemical, use_gap in MODES:
        run_number += 1
        print(
            f"[{run_number}/{total_runs}] {mode}: "
            f"chemical={cs}, gap={gs}, decay={decay}"
        )

        left = simulate(
            "ASHL", cs, gs, decay,
            use_chemical=use_chemical,
            use_gap=use_gap,
        )
        right = simulate(
            "ASHR", cs, gs, decay,
            use_chemical=use_chemical,
            use_gap=use_gap,
        )

        left_profile = integrated_profile(left)
        right_profile = integrated_profile(right)

        similarity = cosine_similarity(left_profile, right_profile)

        rows.append({
            "mode": mode,
            "chemical_strength": cs,
            "gap_strength": gs,
            "decay": decay,
            "cosine_similarity": similarity,
            "left_integrated_network_activity": float(left_profile.sum()),
            "right_integrated_network_activity": float(right_profile.sum()),
            "left_peak_activity": float(left.max()),
            "right_peak_activity": float(right.max()),
            "left_active_neurons": int(np.count_nonzero(left_profile > 1e-8)),
            "right_active_neurons": int(np.count_nonzero(right_profile > 1e-8)),
            "left_max_activity": float(left.max()),
            "right_max_activity": float(right.max()),
        })

results = pd.DataFrame(rows)
results.to_csv(OUTPUT, index=False)

print("\n" + "=" * 70)
print("COUPLING ABLATION SUMMARY")
print("=" * 70)

summary = (
    results.groupby("mode")["cosine_similarity"]
    .agg(["count", "mean", "min", "max"])
    .sort_values("mean")
)
print(summary.to_string(float_format=lambda x: f"{x:.6f}"))

print("\nMean integrated network activity by mode:")
print(
    results.groupby("mode")[
        [
            "left_integrated_network_activity",
            "right_integrated_network_activity",
        ]
    ].mean().to_string(float_format=lambda x: f"{x:.6f}")
)

print("\nFull model versus ablations at chemical strength 1.2:")
high = results[results["chemical_strength"] == 1.2]
print(
    high.groupby("mode")["cosine_similarity"]
    .agg(["mean", "min", "max"])
    .to_string(float_format=lambda x: f"{x:.6f}")
)

print(f"\nSaved: {OUTPUT.resolve()}")
print(
    "\nInterpretation caution: this identifies behavior of the current "
    "exploratory model, not experimentally validated worm neural dynamics."
)
