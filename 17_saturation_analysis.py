
import pandas as pd

df = pd.read_csv("robustness_results.csv")

print("=" * 65)
print("MODEL SENSITIVITY SUMMARY")
print("=" * 65)

# Compare response similarity across chemical strengths.
print("\n=== CHEMICAL STRENGTH ===")
print(
    df.groupby("chemical_strength")["cosine_similarity"]
    .agg(["mean", "min", "max"])
    .round(6)
    .to_string()
)

# Report the parameter combinations for closer inspection.
for label, subset in [
    ("LOW CHEMICAL STRENGTH", df[df["chemical_strength"] == 0.4]),
    ("MEDIUM CHEMICAL STRENGTH", df[df["chemical_strength"] == 0.8]),
    ("HIGH CHEMICAL STRENGTH", df[df["chemical_strength"] == 1.2]),
]:
    print(f"\n=== {label} ===")
    print(
        subset[
            [
                "gap_strength",
                "decay",
                "cosine_similarity",
                "strongest_right_higher_neuron",
                "strongest_left_higher_neuron",
            ]
        ].sort_values("cosine_similarity").to_string(index=False)
    )

print(
    "\nNote: robustness_results.csv does not contain per-neuron "
    "activity maxima or total integrated activity. Saturation "
    "cannot be established from this file alone."
)
print(
    "\nNext, we should rerun the simulation with explicit "
    "saturation counters and activity-magnitude measurements."
)
