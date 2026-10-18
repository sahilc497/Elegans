
import csv
import numpy as np
from cect.readers import White_whole

# ============================================================
# 1. LOAD CONNECTOME
# ============================================================
ds = White_whole.get_instance()
cells, connections = ds.get_neuron_to_neuron_conns()
cell_index = {name: i for i, name in enumerate(cells)}
n = len(cells)

chemical = np.zeros((n, n), dtype=float)
gap_pairs = {}

for conn in connections:
    if conn.pre_cell not in cell_index:
        continue
    if conn.post_cell not in cell_index:
        continue

    i = cell_index[conn.pre_cell]
    j = cell_index[conn.post_cell]
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
        where=sums != 0,
    )


chemical = normalize_rows(chemical)
gap = normalize_rows(gap)

# ============================================================
# 2. SIMULATION FUNCTION
# ============================================================
dt = 0.02
duration = 10.0
steps = int(duration / dt)
time = np.arange(steps) * dt
tau = 0.5

def simulate(stimulated_cell, chemical_strength,
             gap_strength, decay, bounded=True):

    stimulus = np.zeros(n)
    stimulus[cell_index[stimulated_cell]] = 1.0
    activity = np.zeros((steps, n), dtype=float)

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
        current = previous + (dt / tau) * (target - previous)

        # Compare the existing upper cap against no upper cap.
        if bounded:
            current = np.clip(current, 0.0, 1.0)
        else:
            current = np.maximum(current, 0.0)

        if not np.all(np.isfinite(current)):
            raise FloatingPointError(
                f"Non-finite activity: {stimulated_cell}, "
                f"chemical={chemical_strength}, "
                f"gap={gap_strength}, decay={decay}"
            )

        activity[t] = current

    integrated = np.trapezoid(activity, x=time, axis=0)

    return activity, integrated


def cosine_similarity(a, b):
    na = np.linalg.norm(a)
    nb = np.linalg.norm(b)

    if na == 0 or nb == 0:
        return float("nan")

    return float(np.dot(a, b) / (na * nb))


# ============================================================
# 3. RUN PARAMETER SWEEP
# ============================================================
chemical_values = [0.4, 0.8, 1.2]
gap_values = [0.125, 0.25, 0.375]
decay_values = [0.075, 0.15, 0.225]

results = []

for bounded in (True, False):
    for cs in chemical_values:
        for gs in gap_values:
            for decay in decay_values:

                left_activity, left_auc = simulate(
                    "ASHL", cs, gs, decay, bounded
                )
                right_activity, right_auc = simulate(
                    "ASHR", cs, gs, decay, bounded
                )

                similarity = cosine_similarity(
                    left_auc, right_auc
                )

                # Exclude directly stimulated neurons.
                downstream = [
                    i for i, name in enumerate(cells)
                    if name not in ("ASHL", "ASHR")
                ]

                left = left_activity[:, downstream]
                right = right_activity[:, downstream]

                results.append({
                    "bounded": bounded,
                    "chemical_strength": cs,
                    "gap_strength": gs,
                    "decay": decay,
                    "cosine_similarity": similarity,
                    "left_total_integrated_activity": float(
                        left_auc.sum()
                    ),
                    "right_total_integrated_activity": float(
                        right_auc.sum()
                    ),
                    "left_max_activity": float(
                        left_activity.max()
                    ),
                    "right_max_activity": float(
                        right_activity.max()
                    ),
                    "left_fraction_at_cap": float(
                        np.mean(left_activity >= 1.0 - 1e-12)
                    ) if bounded else float("nan"),
                    "right_fraction_at_cap": float(
                        np.mean(right_activity >= 1.0 - 1e-12)
                    ) if bounded else float("nan"),
                    "left_downstream_max": float(left.max()),
                    "right_downstream_max": float(right.max()),
                    "left_downstream_integrated": float(left.sum() * dt),
                    "right_downstream_integrated": float(right.sum() * dt),
                })

# ============================================================
# 4. SAVE AND SUMMARIZE
# ============================================================
output_file = "saturation_test_results.csv"

with open(output_file, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=results[0].keys())
    writer.writeheader()
    writer.writerows(results)

print("=" * 65)
print("BOUNDED VS UNBOUNDED DYNAMICS")
print("=" * 65)

for bounded in (True, False):
    subset = [
        row for row in results
        if row["bounded"] == bounded
    ]

    similarities = np.array([
        row["cosine_similarity"] for row in subset
    ])

    print("\nModel:", "BOUNDED [0, 1]" if bounded else "UNBOUNDED ABOVE")
    print("Parameter combinations:", len(subset))
    print("Mean cosine similarity:", f"{np.nanmean(similarities):.6f}")
    print("Min cosine similarity:", f"{np.nanmin(similarities):.6f}")
    print("Max cosine similarity:", f"{np.nanmax(similarities):.6f}")

    if bounded:
        fractions = [
            max(row["left_fraction_at_cap"],
                row["right_fraction_at_cap"])
            for row in subset
        ]
        print(
            "Largest fraction of all activity entries at cap:",
            f"{max(fractions):.6f}"
        )
    else:
        maxima = [
            max(row["left_max_activity"],
                row["right_max_activity"])
            for row in subset
        ]
        print(
            "Largest activity observed:",
            f"{max(maxima):.6f}"
        )

print("\n=== HIGH CHEMICAL STRENGTH: DETAILED COMPARISON ===")

for bounded in (True, False):
    subset = [
        row for row in results
        if row["bounded"] == bounded
        and row["chemical_strength"] == 1.2
    ]

    print("\nModel:", "BOUNDED" if bounded else "UNBOUNDED")
    for row in subset:
        print(
            f"gap={row['gap_strength']:.3f}, "
            f"decay={row['decay']:.3f}, "
            f"similarity={row['cosine_similarity']:.6f}, "
            f"left_max={row['left_max_activity']:.6f}, "
            f"right_max={row['right_max_activity']:.6f}"
        )

print("\nSaved:", output_file)
print(
    "Note: removing the upper cap is a numerical diagnostic, "
    "not a biologically validated alternative model."
)
