
from collections import defaultdict
import csv

import matplotlib.pyplot as plt
import numpy as np
from cect.readers import White_whole

# --------------------------------------------------
# 1. Load the connectome
# --------------------------------------------------
ds = White_whole.get_instance()
cells, connections = ds.get_neuron_to_neuron_conns()

cell_index = {name: i for i, name in enumerate(cells)}
n = len(cells)

chemical = np.zeros((n, n), dtype=float)
gap = np.zeros((n, n), dtype=float)

# Chemical connections remain directed.
# Gap junctions are treated as undirected electrical coupling.
gap_weights = defaultdict(float)

for conn in connections:
    pre = conn.pre_cell
    post = conn.post_cell
    weight = float(conn.number)

    if pre not in cell_index or post not in cell_index:
        continue

    i = cell_index[pre]
    j = cell_index[post]

    if conn.syntype == "Chemical":
        chemical[i, j] += weight

    elif conn.syntype == "GapJunction":
        # Canonical pair prevents counting reciprocal records twice.
        pair = tuple(sorted((i, j)))
        gap_weights[pair] = max(gap_weights[pair], weight)

for (i, j), weight in gap_weights.items():
    gap[i, j] = weight
    gap[j, i] = weight

# Normalize each matrix so model parameters are easier to interpret.
def row_normalize(matrix):
    row_sums = matrix.sum(axis=1, keepdims=True)
    return np.divide(
        matrix,
        row_sums,
        out=np.zeros_like(matrix),
        where=row_sums != 0,
    )

chemical_norm = row_normalize(chemical)
gap_norm = row_normalize(gap)

print("=== NETWORK ===")
print("Cells:", n)
print("Chemical directed matrix entries:", np.count_nonzero(chemical))
print("Undirected gap-junction pairs:", len(gap_weights))
print("Chemical matrix shape:", chemical.shape)
print("Gap matrix shape:", gap.shape)

# --------------------------------------------------
# 2. Configure an exploratory rate-based model
# --------------------------------------------------
# This is a simplified computational model.
# Chemical effects are provisionally treated as positive;
# actual biological polarity is not inferred from this dataset.
dt = 0.02
duration = 10.0
steps = int(duration / dt)
time = np.arange(steps) * dt

tau = 0.5
chemical_strength = 0.8
gap_strength = 0.25
decay = 0.15

activity = np.zeros((steps, n), dtype=float)

# Stimulate ASHL and ASHR during the first second.
sensory_cells = ["ASHL", "ASHR"]
missing = [name for name in sensory_cells if name not in cell_index]

if missing:
    raise ValueError(f"Missing expected sensory cells: {missing}")

stimulus = np.zeros(n)
for name in sensory_cells:
    stimulus[cell_index[name]] = 1.0

for t in range(1, steps):
    previous = activity[t - 1]

    chemical_input = chemical_norm.T @ previous

    # Electrical coupling responds to activity differences.
    electrical_input = gap_norm @ previous - previous

    external_input = stimulus if time[t] < 1.0 else 0.0

    drive = (
        chemical_strength * chemical_input
        + gap_strength * electrical_input
        - decay * previous
        + external_input
    )

    # Leaky rate dynamics with a nonnegative activation function.
    target = np.maximum(drive, 0.0)
    activity[t] = previous + (dt / tau) * (target - previous)
    activity[t] = np.clip(activity[t], 0.0, 1.0)


# --------------------------------------------------
# 3. Summarize and export results
# --------------------------------------------------
import csv

mean_activity = activity.mean(axis=1)
peak_activity = activity.max(axis=1)

print("\n=== SIMULATION RESULTS ===")
print("Duration:", duration, "seconds")
print("Time step:", dt, "seconds")
print("Peak network activity:", float(peak_activity.max()))
print("Final mean activity:", float(mean_activity[-1]))

# Network-wide summary
with open("simulation_activity.csv", "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["time_seconds", "mean_activity", "peak_activity"])
    writer.writerows(zip(time, mean_activity, peak_activity))

# Full activity matrix: one row per time step, one column per cell.
with open("neuron_activity.csv", "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["time_seconds"] + list(cells))
    for t in range(steps):
        writer.writerow([time[t]] + activity[t].tolist())

# Rank neurons by their maximum activity and integrated activity.
peak_per_neuron = activity.max(axis=0)
auc_per_neuron = np.trapezoid(activity, x=time, axis=0)

ranking = sorted(
    range(n),
    key=lambda i: peak_per_neuron[i],
    reverse=True
)

print("\n=== TOP 15 NEURONS BY PEAK ACTIVITY ===")
for rank, i in enumerate(ranking[:15], start=1):
    print(
        f"{rank:2}. {cells[i]:8s} "
        f"peak={peak_per_neuron[i]:.6f} "
        f"integrated={auc_per_neuron[i]:.6f}"
    )

with open("neuron_activity_ranking.csv", "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["neuron", "peak_activity", "integrated_activity"])
    for i in ranking:
        writer.writerow([cells[i], peak_per_neuron[i], auc_per_neuron[i]])

# Plot network-wide activity.
plt.figure(figsize=(10, 5))
plt.plot(time, mean_activity, label="Mean activity")
plt.plot(time, peak_activity, label="Peak activity", alpha=0.8)
plt.xlabel("Time (seconds)")
plt.ylabel("Activity (arbitrary units)")
plt.title("Exploratory C. elegans Connectome Simulation")
plt.legend()
plt.tight_layout()
plt.savefig("09_network_activity.png", dpi=160)
plt.close()

# Plot sensory-neuron activity.
plt.figure(figsize=(10, 5))
for name in sensory_cells:
    plt.plot(time, activity[:, cell_index[name]], label=name)
plt.xlabel("Time (seconds)")
plt.ylabel("Activity (arbitrary units)")
plt.title("Stimulated Sensory Neurons")
plt.legend()
plt.tight_layout()
plt.savefig("09_sensory_activity.png", dpi=160)
plt.close()

# Plot the top 10 neurons by peak activity.
top_n = ranking[:10]
plt.figure(figsize=(11, 6))
for i in top_n:
    plt.plot(time, activity[:, i], label=cells[i])
plt.xlabel("Time (seconds)")
plt.ylabel("Activity (arbitrary units)")
plt.title("Top 10 Neurons by Peak Activity")
plt.legend(fontsize=8, ncol=2)
plt.tight_layout()
plt.savefig("09_top_neuron_activity.png", dpi=160)
plt.close()

print("\nSaved:")
print("  simulation_activity.csv")
print("  neuron_activity.csv")
print("  neuron_activity_ranking.csv")
print("  09_network_activity.png")
print("  09_sensory_activity.png")
print("  09_top_neuron_activity.png")
print("  09_chemical_connectivity.png (from earlier run)")
