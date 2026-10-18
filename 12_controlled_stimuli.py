
import csv
import numpy as np
import matplotlib.pyplot as plt
from cect.readers import White_whole

# ============================================================
# 1. LOAD C. ELEGANS CONNECTOME
# ============================================================
ds = White_whole.get_instance()
cells, connections = ds.get_neuron_to_neuron_conns()

cell_index = {name: i for i, name in enumerate(cells)}
n = len(cells)

chemical = np.zeros((n, n), dtype=float)
gap_pairs = {}

for conn in connections:
    pre = conn.pre_cell
    post = conn.post_cell
    weight = float(conn.number)

    if pre not in cell_index or post not in cell_index:
        continue

    i = cell_index[pre]
    j = cell_index[post]

    if conn.syntype == "Chemical":
        # Directed chemical connection.
        chemical[i, j] += weight

    elif conn.syntype == "GapJunction":
        # Treat reciprocal gap-junction records as one pair.
        pair = tuple(sorted((i, j)))
        gap_pairs[pair] = max(gap_pairs.get(pair, 0.0), weight)

gap = np.zeros((n, n), dtype=float)

for (i, j), weight in gap_pairs.items():
    gap[i, j] = weight
    gap[j, i] = weight


def normalize_rows(matrix):
    """Normalize each row while safely handling empty rows."""
    totals = matrix.sum(axis=1, keepdims=True)

    return np.divide(
        matrix,
        totals,
        out=np.zeros_like(matrix),
        where=totals != 0,
    )


chemical_norm = normalize_rows(chemical)
gap_norm = normalize_rows(gap)

print("=" * 60)
print("C. ELEGANS CONTROLLED STIMULUS EXPERIMENTS")
print("=" * 60)
print("Cells:", n)
print("Chemical directed matrix entries:", np.count_nonzero(chemical))
print("Undirected gap-junction pairs:", len(gap_pairs))

# ============================================================
# 2. EXPERIMENT CONFIGURATION
# ============================================================
conditions = {
    "ASHL_full": {"ASHL": 1.0},
    "ASHR_full": {"ASHR": 1.0},
    "ASH_pair_half": {"ASHL": 0.5, "ASHR": 0.5},
    "ASH_pair_full": {"ASHL": 1.0, "ASHR": 1.0},
}

for label, inputs in conditions.items():
    missing = [name for name in inputs if name not in cell_index]
    if missing:
        raise ValueError(
            f"Missing cells in condition {label}: {missing}"
        )

dt = 0.02
duration = 10.0
steps = int(duration / dt)
time = np.arange(steps) * dt

# These are exploratory model parameters, not fitted
# biological constants.
tau = 0.5
chemical_strength = 0.8
gap_strength = 0.25
decay = 0.15

# Store experiment outputs.
traces = {}
summary_rows = []
ranking_rows = []

# ============================================================
# 3. RUN EACH STIMULUS EXPERIMENT
# ============================================================
for label, inputs in conditions.items():

    stimulus = np.zeros(n)

    for neuron, amplitude in inputs.items():
        stimulus[cell_index[neuron]] = amplitude

    activity = np.zeros((steps, n), dtype=float)

    for t in range(1, steps):
        previous = activity[t - 1]

        # Activity arriving through directed chemical connections.
        chemical_input = chemical_norm.T @ previous

        # Electrical coupling responds to differences in activity.
        electrical_input = gap_norm @ previous - previous

        # Stimulus is applied for the first second.
        external_input = stimulus if time[t] < 1.0 else 0.0

        drive = (
            chemical_strength * chemical_input
            + gap_strength * electrical_input
            - decay * previous
            + external_input
        )

        # Simplified nonnegative rate-based dynamics.
        target = np.maximum(drive, 0.0)

        activity[t] = previous + (dt / tau) * (
            target - previous
        )

        activity[t] = np.clip(activity[t], 0.0, 1.0)

    traces[label] = activity

    # --------------------------------------------------------
    # 3A. Network-wide measurements
    # --------------------------------------------------------
    mean_trace = activity.mean(axis=1)
    peak_trace = activity.max(axis=1)

    peak_network = float(peak_trace.max())
    maximum_mean = float(mean_trace.max())
    integrated_mean = float(np.trapezoid(mean_trace, x=time))
    final_mean = float(mean_trace[-1])

    summary_rows.append({
        "condition": label,
        "stimulated_cells": ";".join(inputs.keys()),
        "total_input_amplitude": sum(inputs.values()),
        "peak_network_activity": peak_network,
        "maximum_mean_activity": maximum_mean,
        "integrated_mean_activity": integrated_mean,
        "final_mean_activity": final_mean,
    })

    # --------------------------------------------------------
    # 3B. Rank individual neurons
    # --------------------------------------------------------
    peak_per_neuron = activity.max(axis=0)
    auc_per_neuron = np.trapezoid(
        activity, x=time, axis=0
    )

    # Rank by total integrated activity over the experiment.
    ranking = np.argsort(auc_per_neuron)[::-1]

    print("\n" + "=" * 60)
    print(f"EXPERIMENT: {label}")
    print("=" * 60)
    print("Stimulated cells:", ", ".join(inputs.keys()))
    print("Total input amplitude:", sum(inputs.values()))
    print(f"Peak network activity: {peak_network:.8f}")
    print(f"Maximum mean activity: {maximum_mean:.8f}")
    print(f"Integrated mean activity: {integrated_mean:.8f}")
    print(f"Final mean activity: {final_mean:.8f}")

    print("\nTOP 10 NEURONS BY INTEGRATED ACTIVITY")

    for rank, i in enumerate(ranking[:10], start=1):
        print(
            f"{rank:2}. {cells[i]:8s} "
            f"integrated={auc_per_neuron[i]:.8f} "
            f"peak={peak_per_neuron[i]:.8f}"
        )

    for rank, i in enumerate(ranking, start=1):
        ranking_rows.append({
            "condition": label,
            "rank": rank,
            "neuron": cells[i],
            "integrated_activity": float(auc_per_neuron[i]),
            "peak_activity": float(peak_per_neuron[i]),
        })


# ============================================================
# 4. EXPORT NETWORK SUMMARY
# ============================================================
with open(
    "controlled_stimulus_results.csv",
    "w",
    newline="",
    encoding="utf-8",
) as f:
    writer = csv.DictWriter(
        f,
        fieldnames=list(summary_rows[0].keys()),
    )
    writer.writeheader()
    writer.writerows(summary_rows)


# ============================================================
# 5. EXPORT ALL NEURON RANKINGS
# ============================================================
with open(
    "controlled_neuron_rankings.csv",
    "w",
    newline="",
    encoding="utf-8",
) as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "condition",
            "rank",
            "neuron",
            "integrated_activity",
            "peak_activity",
        ],
    )
    writer.writeheader()
    writer.writerows(ranking_rows)


# ============================================================
# 6. EXPORT FULL ACTIVITY MATRIX FOR EACH EXPERIMENT
# ============================================================
for label, activity in traces.items():
    filename = f"activity_{label}.csv"

    with open(
        filename,
        "w",
        newline="",
        encoding="utf-8",
    ) as f:
        writer = csv.writer(f)
        writer.writerow(["time_seconds"] + list(cells))

        for t in range(steps):
            writer.writerow(
                [time[t]] + activity[t].tolist()
            )


# ============================================================
# 7. PLOT NETWORK MEAN ACTIVITY
# ============================================================
plt.figure(figsize=(11, 6))

for label, activity in traces.items():
    plt.plot(
        time,
        activity.mean(axis=1),
        label=label,
    )

plt.xlabel("Time (seconds)")
plt.ylabel("Mean activity (arbitrary units)")
plt.title("Controlled Sensory Input Comparison")
plt.legend()
plt.tight_layout()
plt.savefig("12_controlled_stimuli.png", dpi=160)
plt.close()


# ============================================================
# 8. PLOT PEAK NETWORK ACTIVITY
# ============================================================
plt.figure(figsize=(11, 6))

for label, activity in traces.items():
    plt.plot(
        time,
        activity.max(axis=1),
        label=label,
    )

plt.xlabel("Time (seconds)")
plt.ylabel("Peak activity (arbitrary units)")
plt.title("Peak Network Response to Controlled Inputs")
plt.legend()
plt.tight_layout()
plt.savefig("12_peak_controlled_stimuli.png", dpi=160)
plt.close()


# ============================================================
# 9. PLOT TOP 10 NEURONS FOR EACH EXPERIMENT
# ============================================================
fig, axes = plt.subplots(
    2, 2,
    figsize=(15, 10),
    sharex=True,
)

for ax, (label, activity) in zip(
    axes.flat, traces.items()
):
    auc_per_neuron = np.trapezoid(
        activity, x=time, axis=0
    )
    top_indices = np.argsort(auc_per_neuron)[::-1][:10]

    for i in top_indices:
        ax.plot(
            time,
            activity[:, i],
            label=cells[i],
        )

    ax.set_title(label)
    ax.set_xlabel("Time (seconds)")
    ax.set_ylabel("Activity (arbitrary units)")
    ax.legend(fontsize=7, ncol=2)

fig.suptitle(
    "Top 10 Neurons by Integrated Activity",
    fontsize=14,
)
fig.tight_layout()
fig.savefig("12_top_neurons_by_condition.png", dpi=160)
plt.close(fig)


# ============================================================
# 10. FINISH
# ============================================================
print("\n" + "=" * 60)
print("ALL EXPERIMENTS COMPLETED")
print("=" * 60)

print("\nGenerated files:")
print("  controlled_stimulus_results.csv")
print("  controlled_neuron_rankings.csv")
print("  activity_ASHL_full.csv")
print("  activity_ASHR_full.csv")
print("  activity_ASH_pair_half.csv")
print("  activity_ASH_pair_full.csv")
print("  12_controlled_stimuli.png")
print("  12_peak_controlled_stimuli.png")
print("  12_top_neurons_by_condition.png")

print(
    "\nScientific note: chemical connections are treated as "
    "positive inputs in this exploratory model. The results "
    "are computational outputs, not validated biological "
    "predictions."
)
