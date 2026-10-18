
import csv
import numpy as np
import matplotlib.pyplot as plt
from cect.readers import White_whole

# Load the same connectome used in the main simulation.
ds = White_whole.get_instance()
cells, connections = ds.get_neuron_to_neuron_conns()
cell_index = {name: i for i, name in enumerate(cells)}
n = len(cells)

chemical = np.zeros((n, n), dtype=float)
gap_weights = {}

for conn in connections:
    pre = conn.pre_cell
    post = conn.post_cell
    weight = float(conn.number)

    if pre not in cell_index or post not in cell_index:
        continue

    i, j = cell_index[pre], cell_index[post]

    if conn.syntype == "Chemical":
        chemical[i, j] += weight
    elif conn.syntype == "GapJunction":
        pair = tuple(sorted((i, j)))
        gap_weights[pair] = max(gap_weights.get(pair, 0.0), weight)

gap = np.zeros((n, n), dtype=float)
for (i, j), weight in gap_weights.items():
    gap[i, j] = weight
    gap[j, i] = weight


def row_normalize(matrix):
    totals = matrix.sum(axis=1, keepdims=True)
    return np.divide(
        matrix, totals,
        out=np.zeros_like(matrix),
        where=totals != 0
    )


chemical_norm = row_normalize(chemical)
gap_norm = row_normalize(gap)

# Each experiment stimulates one neuron for the first second.
# These are computational probes, not validated biological stimuli.
stimulus_groups = {
    "ASH_pair": ["ASHL", "ASHR"],
    "ADFL": ["ADFL"],
    "AFDL": ["AFDL"],
}

for label, names in stimulus_groups.items():
    missing = [name for name in names if name not in cell_index]
    if missing:
        raise ValueError(f"{label}: missing cells {missing}")

dt = 0.02
duration = 10.0
steps = int(duration / dt)
time = np.arange(steps) * dt

tau = 0.5
chemical_strength = 0.8
gap_strength = 0.25
decay = 0.15

summary = []
traces = {}

for label, names in stimulus_groups.items():
    stimulus = np.zeros(n)
    for name in names:
        stimulus[cell_index[name]] = 1.0

    activity = np.zeros((steps, n))

    for t in range(1, steps):
        previous = activity[t - 1]

        chemical_input = chemical_norm.T @ previous
        electrical_input = gap_norm @ previous - previous

        external_input = stimulus if time[t] < 1.0 else 0.0

        drive = (
            chemical_strength * chemical_input
            + gap_strength * electrical_input
            - decay * previous
            + external_input
        )

        target = np.maximum(drive, 0.0)
        activity[t] = previous + (dt / tau) * (target - previous)
        activity[t] = np.clip(activity[t], 0.0, 1.0)

    traces[label] = activity

    mean_trace = activity.mean(axis=1)
    peak_trace = activity.max(axis=1)

    summary.append({
        "experiment": label,
        "stimulated_cells": ";".join(names),
        "peak_network_activity": float(peak_trace.max()),
        "max_mean_network_activity": float(mean_trace.max()),
        "integrated_mean_activity": float(np.trapezoid(mean_trace, x=time)),
        "final_mean_activity": float(mean_trace[-1]),
    })

with open("stimulus_experiment_summary.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=summary[0].keys())
    writer.writeheader()
    writer.writerows(summary)

# Plot network mean activity for each experiment.
plt.figure(figsize=(10, 5))
for label, activity in traces.items():
    plt.plot(time, activity.mean(axis=1), label=label)

plt.xlabel("Time (seconds)")
plt.ylabel("Mean activity (arbitrary units)")
plt.title("Connectome Response to Different Input Neurons")
plt.legend()
plt.tight_layout()
plt.savefig("11_stimulus_comparison.png", dpi=160)
plt.close()

# Plot peak activity for each experiment.
plt.figure(figsize=(10, 5))
for label, activity in traces.items():
    plt.plot(time, activity.max(axis=1), label=label)

plt.xlabel("Time (seconds)")
plt.ylabel("Peak activity (arbitrary units)")
plt.title("Peak Network Response to Different Inputs")
plt.legend()
plt.tight_layout()
plt.savefig("11_peak_comparison.png", dpi=160)
plt.close()

print("=== STIMULUS EXPERIMENTS ===")
for row in summary:
    print(f"\n{row['experiment']}")
    for key, value in row.items():
        if key != "experiment":
            print(f"  {key}: {value}")

print("\nSaved:")
print("  stimulus_experiment_summary.csv")
print("  11_stimulus_comparison.png")
print("  11_peak_comparison.png")
