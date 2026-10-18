
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from cect.readers import White_whole

# ============================================================
# C. elegans: temporal propagation analysis
# ============================================================

OUT = Path(".")
DT = 0.02
DURATION = 10.0
TAU = 0.5
CHEMICAL_STRENGTH = 0.8
GAP_STRENGTH = 0.25
DECAY = 0.15
STIMULUS_DURATION = 1.0
STIMULUS_AMPLITUDE = 1.0
ACTIVE_THRESHOLD = 1e-4

print("Loading White et al. 1986 connectome...")
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
    nonzero = sums > 0
    result[nonzero] /= sums[nonzero, None]
    return result

chemical_norm = normalize_rows(chemical)
gap_norm = normalize_rows(gap)

def simulate(stimulus, use_chemical=True, use_gap=True):
    steps = int(DURATION / DT)
    activity = np.zeros((steps, n), dtype=float)
    stimulus_index = index[stimulus]

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
            external[stimulus_index] = STIMULUS_AMPLITUDE

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

modes = {
    "full": (True, True),
    "chemical_only": (True, False),
    "gap_only": (False, True),
}

times = np.arange(int(DURATION / DT)) * DT
summary_rows = []
neuron_rows = []
timecourse_rows = []

# Use no-recurrent-coupling simulations as a matched baseline.
for stimulus in ["ASHL", "ASHR"]:
    baseline = simulate(stimulus, False, False)

    for mode, (use_chem, use_gap) in modes.items():
        print(f"Simulating {stimulus}: {mode}")
        activity = simulate(stimulus, use_chem, use_gap)
        excess = np.maximum(activity - baseline, 0.0)

        active_counts = np.count_nonzero(
            activity > ACTIVE_THRESHOLD, axis=1
        )
        excess_counts = np.count_nonzero(
            excess > ACTIVE_THRESHOLD, axis=1
        )

        summary_rows.append({
            "stimulus": stimulus,
            "mode": mode,
            "peak_network_activity": float(activity.max()),
            "peak_mean_activity": float(activity.mean(axis=1).max()),
            "integrated_network_activity": float(activity.sum() * DT),
            "maximum_active_neurons": int(active_counts.max()),
            "maximum_excess_active_neurons": int(excess_counts.max()),
            "neurons_above_threshold_at_end": int(
                np.count_nonzero(activity[-1] > ACTIVE_THRESHOLD)
            ),
        })

        # Per-neuron activity timing and magnitude.
        for i, name in enumerate(names):
            signal = activity[:, i]
            active = np.flatnonzero(signal > ACTIVE_THRESHOLD)

            neuron_rows.append({
                "stimulus": stimulus,
                "mode": mode,
                "neuron": name,
                "peak_activity": float(signal.max()),
                "integrated_activity": float(signal.sum() * DT),
                "first_active_time_s": (
                    float(times[active[0]]) if len(active) else np.nan
                ),
                "last_active_time_s": (
                    float(times[active[-1]]) if len(active) else np.nan
                ),
            })

        # Population activity by time.
        for t in range(len(times)):
            timecourse_rows.append({
                "stimulus": stimulus,
                "mode": mode,
                "time_s": float(times[t]),
                "mean_activity": float(activity[t].mean()),
                "active_neurons": int(active_counts[t]),
                "excess_active_neurons": int(excess_counts[t]),
            })

summary = pd.DataFrame(summary_rows)
neurons = pd.DataFrame(neuron_rows)
timecourses = pd.DataFrame(timecourse_rows)

summary.to_csv(OUT / "propagation_summary.csv", index=False)
neurons.to_csv(OUT / "propagation_neuron_details.csv", index=False)
timecourses.to_csv(OUT / "propagation_timecourses.csv", index=False)

# Plot 1: population mean activity over time.
fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)

for ax, stimulus in zip(axes, ["ASHL", "ASHR"]):
    subset = timecourses[timecourses["stimulus"] == stimulus]
    for mode in modes:
        part = subset[subset["mode"] == mode]
        ax.plot(part["time_s"], part["mean_activity"], label=mode)

    ax.axvline(STIMULUS_DURATION, linestyle="--", linewidth=1)
    ax.set_title(f"{stimulus} stimulation")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Mean activity across cells")
    ax.legend()
    ax.grid(alpha=0.25)

fig.tight_layout()
fig.savefig(OUT / "20_population_activity.png", dpi=180)
plt.close(fig)

# Plot 2: number of active neurons over time.
fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)

for ax, stimulus in zip(axes, ["ASHL", "ASHR"]):
    subset = timecourses[timecourses["stimulus"] == stimulus]
    for mode in modes:
        part = subset[subset["mode"] == mode]
        ax.plot(part["time_s"], part["active_neurons"], label=mode)

    ax.axvline(STIMULUS_DURATION, linestyle="--", linewidth=1)
    ax.set_title(f"{stimulus} stimulation")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Cells above activity threshold")
    ax.legend()
    ax.grid(alpha=0.25)

fig.tight_layout()
fig.savefig(OUT / "20_active_neurons.png", dpi=180)
plt.close(fig)

print("\n" + "=" * 65)
print("PROPAGATION ANALYSIS SUMMARY")
print("=" * 65)
print(summary.to_string(index=False, float_format=lambda x: f"{x:.6f}"))

print("\nTop 10 neurons by integrated activity per stimulus and mode:")
for stimulus in ["ASHL", "ASHR"]:
    for mode in modes:
        subset = neurons[
            (neurons["stimulus"] == stimulus)
            & (neurons["mode"] == mode)
        ]
        top = subset.nlargest(10, "integrated_activity")
        print(f"\n{stimulus} / {mode}")
        print(
            top[["neuron", "integrated_activity", "peak_activity",
                 "first_active_time_s"]].to_string(
                index=False,
                float_format=lambda x: f"{x:.6f}"
            )
        )

print("\nSaved files:")
for filename in [
    "propagation_summary.csv",
    "propagation_neuron_details.csv",
    "propagation_timecourses.csv",
    "20_population_activity.png",
    "20_active_neurons.png",
]:
    print(f"  {OUT / filename}")

print(
    "\nCaution: activity thresholds and model parameters are exploratory; "
    "these are simulated responses, not measured biological firing."
)
