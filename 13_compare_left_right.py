
import numpy as np
import pandas as pd

# Load the individual-neuron activity data.
left = pd.read_csv("activity_ASHL_full.csv")
right = pd.read_csv("activity_ASHR_full.csv")

# Confirm both files use the same cell ordering.
if list(left.columns) != list(right.columns):
    raise ValueError("Neuron columns differ between experiment files.")

if not np.allclose(left["time_seconds"], right["time_seconds"]):
    raise ValueError("Time points differ between experiment files.")

neuron_names = list(left.columns[1:])

left_activity = left[neuron_names].to_numpy()
right_activity = right[neuron_names].to_numpy()
time = left["time_seconds"].to_numpy()

# Integrated activity for each neuron.
left_auc = np.trapezoid(left_activity, x=time, axis=0)
right_auc = np.trapezoid(right_activity, x=time, axis=0)

# Compare absolute and relative differences.
difference = right_auc - left_auc
mean_response = (right_auc + left_auc) / 2

relative_difference = np.divide(
    difference,
    mean_response,
    out=np.zeros_like(difference),
    where=mean_response > 1e-12,
)

comparison = pd.DataFrame({
    "neuron": neuron_names,
    "ASHL_integrated_activity": left_auc,
    "ASHR_integrated_activity": right_auc,
    "ASHR_minus_ASHL": difference,
    "relative_difference": relative_difference,
})

comparison["absolute_difference"] = np.abs(difference)

# Exclude the directly stimulated neurons when ranking
# downstream differences.
downstream = comparison[
    ~comparison["neuron"].isin(["ASHL", "ASHR"])
].copy()

downstream = downstream.sort_values(
    "absolute_difference",
    ascending=False,
)

comparison.to_csv(
    "left_right_neuron_comparison.csv",
    index=False,
)

print("=" * 65)
print("LEFT/RIGHT NEURAL RESPONSE COMPARISON")
print("=" * 65)

print("\nTOP 15 DOWNSTREAM NEURONS WITH THE LARGEST DIFFERENCES")

for rank, (_, row) in enumerate(
    downstream.head(15).iterrows(), start=1
):
    print(
        f"{rank:2}. {row['neuron']:8s} "
        f"ASHL={row['ASHL_integrated_activity']:.8f} "
        f"ASHR={row['ASHR_integrated_activity']:.8f} "
        f"difference={row['ASHR_minus_ASHL']:+.8f}"
    )

# Pattern similarity across all neurons.
left_norm = np.linalg.norm(left_auc)
right_norm = np.linalg.norm(right_auc)

if left_norm > 0 and right_norm > 0:
    cosine_similarity = np.dot(left_auc, right_auc) / (
        left_norm * right_norm
    )
else:
    cosine_similarity = float("nan")

print("\n=== OVERALL COMPARISON ===")
print(f"Cosine similarity: {cosine_similarity:.6f}")
print(
    "Mean absolute neuron-level difference:",
    f"{np.mean(np.abs(difference)):.8f}",
)

print("\nSaved: left_right_neuron_comparison.csv")
