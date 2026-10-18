
import csv
import numpy as np
from cect.readers import White_whole

# ---------------------------------------------------------
# 1. Load connectome
# ---------------------------------------------------------
ds = White_whole.get_instance()
cells, connections = ds.get_neuron_to_neuron_conns()
cell_index = {name: i for i, name in enumerate(cells)}
n = len(cells)

chemical = np.zeros((n, n), dtype=float)
gap_pairs = {}

for conn in connections:
    pre, post = conn.pre_cell, conn.post_cell

    if pre not in cell_index or post not in cell_index:
        continue

    i, j = cell_index[pre], cell_index[post]
    weight = float(conn.number)

    if conn.syntype == "Chemical":
        chemical[i, j] += weight
    elif conn.syntype == "GapJunction":
        pair = tuple(sorted((i, j)))
        gap_pairs[pair] = max(gap_pairs.get(pair, 0.0), weight)

gap = np.zeros((n, n), dtype=float)

for (i, j), weight in gap_pairs.items():
    gap[i, j] = weight
    gap[j, i] = weight


def normalize_rows(matrix):
    sums = matrix.sum(axis=1, keepdims=True)
    return np.divide(
        matrix, sums,
        out=np.zeros_like(matrix),
        where=sums != 0
    )


chemical = normalize_rows(chemical)
gap = normalize_rows(gap)

# ---------------------------------------------------------
# 2. Simulation function
# ---------------------------------------------------------
dt = 0.02
duration = 10.0
steps = int(duration / dt)
time = np.arange(steps) * dt

stimulus_index = {
    "ASHL": cell_index["ASHL"],
    "ASHR": cell_index["ASHR"],
}

def simulate(stimulated_cell, chemical_strength,
             gap_strength, decay):
    stimulus = np.zeros(n)
    stimulus[stimulus_index[stimulated_cell]] = 1.0

    activity = np.zeros((steps, n), dtype=float)
    tau = 0.5

    for t in range(1, steps):
        previous = activity[t - 1]

        chemical_input = chemical.T @ previous
        electrical_input = gap @ previous - previous

        external = stimulus if time[t] < 1.0 else 0.0

        drive = (
            chemical_strength * chemical_input
            + gap_strength * electrical_input
            - decay * previous
            + external
        )

        target = np.maximum(drive, 0.0)
        activity[t] = previous + (dt / tau) * (
            target - previous
        )
        activity[t] = np.clip(activity[t], 0.0, 1.0)

    # Integrated response per neuron.
    return np.trapezoid(activity, x=time, axis=0)


# ---------------------------------------------------------
# 3. Parameter sweep
# ---------------------------------------------------------
chemical_values = [0.40, 0.80, 1.20]
gap_values = [0.125, 0.25, 0.375]
decay_values = [0.075, 0.15, 0.225]

rows = []
downstream_indices = [
    i for i, name in enumerate(cells)
    if name not in ("ASHL", "ASHR")
]

for cs in chemical_values:
    for gs in gap_values:
        for decay in decay_values:
            left = simulate("ASHL", cs, gs, decay)
            right = simulate("ASHR", cs, gs, decay)

            left_norm = np.linalg.norm(left)
            right_norm = np.linalg.norm(right)

            similarity = (
                float(np.dot(left, right) /
                      (left_norm * right_norm))
                if left_norm > 0 and right_norm > 0
                else float("nan")
            )

            differences = right - left
            downstream_diff = differences[downstream_indices]

            right_higher = int(np.sum(downstream_diff > 1e-12))
            left_higher = int(np.sum(downstream_diff < -1e-12))

            # Identify the strongest downstream difference
            # in each direction.
            right_idx = max(
                downstream_indices,
                key=lambda i: differences[i]
            )
            left_idx = min(
                downstream_indices,
                key=lambda i: differences[i]
            )

            rows.append({
                "chemical_strength": cs,
                "gap_strength": gs,
                "decay": decay,
                "cosine_similarity": similarity,
                "right_higher_count": right_higher,
                "left_higher_count": left_higher,
                "strongest_right_higher_neuron": cells[right_idx],
                "strongest_right_difference": float(
                    differences[right_idx]
                ),
                "strongest_left_higher_neuron": cells[left_idx],
                "strongest_left_difference": float(
                    differences[left_idx]
                ),
            })

# ---------------------------------------------------------
# 4. Save and summarize
# ---------------------------------------------------------
output_file = "robustness_results.csv"

with open(output_file, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)

similarities = np.array([
    row["cosine_similarity"] for row in rows
])

print("=" * 65)
print("LEFT/RIGHT ROBUSTNESS TEST")
print("=" * 65)
print("Parameter combinations:", len(rows))
print(f"Mean cosine similarity: {np.nanmean(similarities):.6f}")
print(f"Minimum cosine similarity: {np.nanmin(similarities):.6f}")
print(f"Maximum cosine similarity: {np.nanmax(similarities):.6f}")

right_counts = [row["right_higher_count"] for row in rows]
left_counts = [row["left_higher_count"] for row in rows]

print("\nDownstream response differences across all parameter sets:")
print(
    "Mean neurons with ASHR > ASHL:",
    f"{np.mean(right_counts):.2f}"
)
print(
    "Mean neurons with ASHL > ASHR:",
    f"{np.mean(left_counts):.2f}"
)

from collections import Counter

right_top = Counter(
    row["strongest_right_higher_neuron"] for row in rows
)
left_top = Counter(
    row["strongest_left_higher_neuron"] for row in rows
)

print("\nMost frequent strongest right-higher neurons:")
for name, count in right_top.most_common(10):
    print(f"  {name}: {count}/27 parameter sets")

print("\nMost frequent strongest left-higher neurons:")
for name, count in left_top.most_common(10):
    print(f"  {name}: {count}/27 parameter sets")

print(f"\nSaved: {output_file}")
