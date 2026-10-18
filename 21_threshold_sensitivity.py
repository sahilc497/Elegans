
import numpy as np
import pandas as pd
from pathlib import Path
from cect.readers import White_whole

# ============================================================
# Experiment 21: activity-threshold sensitivity
# ============================================================

DT = 0.02
DURATION = 10.0
TAU = 0.5
CHEMICAL_STRENGTH = 0.8
GAP_STRENGTH = 0.25
DECAY = 0.15
STIMULUS_DURATION = 1.0
STIMULUS_AMPLITUDE = 1.0

THRESHOLDS = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 5e-2]
OUTPUT = Path("threshold_sensitivity_results.csv")

print("Loading connectome...")
ds = White_whole.get_instance()
cells, connections = ds.get_neuron_to_neuron_conns()

names = list(cells)
index = {name: i for i, name in enumerate(names)}
n = len(names)

chemical = np.zeros((n, n), dtype=float)
gap = np.zeros((n, n), dtype=float)

for conn in connections:
    pre, post = conn.pre_cell, conn.post_cell
    if pre not in index or post not in index:
        continue

    i, j = index[pre], index[post]
    weight = float(conn.number)
    if weight <= 0:
        continue

    if "Generic_CS" in str(conn.synclass):
        chemical[i, j] += weight
    elif "Generic_GJ" in str(conn.synclass):
        gap[i, j] = max(gap[i, j], weight)
        gap[j, i] = max(gap[j, i], weight)

def normalize_rows(matrix):
    result = matrix.copy()
    sums = result.sum(axis=1)
    mask = sums > 0
    result[mask] /= sums[mask, None]
    return result

chemical_norm = normalize_rows(chemical)
gap_norm = normalize_rows(gap)

def simulate(stimulus, use_chemical=True, use_gap=True):
    steps = int(DURATION / DT)
    activity = np.zeros((steps, n), dtype=float)

    for t in range(1, steps):
        previous = activity[t - 1]

        chem_input = (
            chemical_norm.T @ previous
            if use_chemical else np.zeros(n)
        )
        gap_input = (
            gap_norm @ previous - previous
            if use_gap else np.zeros(n)
        )

        external = np.zeros(n)
        if t * DT < STIMULUS_DURATION:
            external[index[stimulus]] = STIMULUS_AMPLITUDE

        drive = (
            CHEMICAL_STRENGTH * chem_input
            + GAP_STRENGTH * gap_input
            - DECAY * previous
            + external
        )

        target = np.maximum(drive, 0.0)
        activity[t] = previous + (DT / TAU) * (
            target - previous
        )
        activity[t] = np.clip(activity[t], 0.0, 1.0)

    return activity

rows = []
modes = {
    "full": (True, True),
    "chemical_only": (True, False),
    "gap_only": (False, True),
}

for stimulus in ["ASHL", "ASHR"]:
    for mode, (use_chemical, use_gap) in modes.items():
        print(f"Simulating {stimulus}: {mode}")
        activity = simulate(stimulus, use_chemical, use_gap)

        peak_by_neuron = activity.max(axis=0)
        integrated_by_neuron = activity.sum(axis=0) * DT

        for threshold in THRESHOLDS:
            rows.append({
                "stimulus": stimulus,
                "mode": mode,
                "threshold": threshold,
                "neurons_ever_above_threshold": int(
                    np.count_nonzero(peak_by_neuron > threshold)
                ),
                "neurons_above_threshold_at_end": int(
                    np.count_nonzero(activity[-1] > threshold)
                ),
                "peak_network_activity": float(activity.max()),
                "integrated_network_activity": float(activity.sum() * DT),
                "neurons_with_integrated_activity_above_threshold": int(
                    np.count_nonzero(integrated_by_neuron > threshold)
                ),
            })

results = pd.DataFrame(rows)
results.to_csv(OUTPUT, index=False)

print("\n" + "=" * 75)
print("THRESHOLD SENSITIVITY RESULTS")
print("=" * 75)

for stimulus in ["ASHL", "ASHR"]:
    for mode in modes:
        subset = results[
            (results["stimulus"] == stimulus)
            & (results["mode"] == mode)
        ]
        print(f"\n{stimulus} / {mode}")
        print(
            subset[
                [
                    "threshold",
                    "neurons_ever_above_threshold",
                    "neurons_above_threshold_at_end",
                    "neurons_with_integrated_activity_above_threshold",
                ]
            ].to_string(
                index=False,
                float_format=lambda x: f"{x:.6g}"
            )
        )

print(f"\nSaved: {OUTPUT.resolve()}")
print(
    "Caution: these thresholds are numerical sensitivity checks, "
    "not biological firing thresholds."
)
