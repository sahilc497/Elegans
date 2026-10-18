
import pandas as pd

df = pd.read_csv("robustness_results.csv")

print("=" * 65)
print("ROBUSTNESS ANALYSIS")
print("=" * 65)

# 1. Overall statistics
print("\n=== OVERALL STATISTICS ===")
print(df["cosine_similarity"].describe().to_string())

# 2. Average similarity grouped by each parameter.
parameters = [
    "chemical_strength",
    "gap_strength",
    "decay",
]

for parameter in parameters:
    print(f"\n=== SIMILARITY BY {parameter.upper()} ===")

    summary = (
        df.groupby(parameter)["cosine_similarity"]
        .agg(["mean", "min", "max", "std"])
        .round(6)
    )

    print(summary.to_string())

# 3. Most and least similar profiles.
print("\n=== FIVE MOST SIMILAR RESPONSE PROFILES ===")

most_similar = df.nlargest(5, "cosine_similarity")

print(
    most_similar[
        [
            "chemical_strength",
            "gap_strength",
            "decay",
            "cosine_similarity",
            "strongest_right_higher_neuron",
            "strongest_left_higher_neuron",
        ]
    ].to_string(index=False)
)

print("\n=== FIVE LEAST SIMILAR RESPONSE PROFILES ===")

least_similar = df.nsmallest(5, "cosine_similarity")

print(
    least_similar[
        [
            "chemical_strength",
            "gap_strength",
            "decay",
            "cosine_similarity",
            "strongest_right_higher_neuron",
            "strongest_left_higher_neuron",
        ]
    ].to_string(index=False)
)

# 4. Save grouped summaries.
for parameter in parameters:
    summary = (
        df.groupby(parameter)["cosine_similarity"]
        .agg(["mean", "min", "max", "std"])
        .reset_index()
    )

    summary.to_csv(
        f"robustness_by_{parameter}.csv",
        index=False,
    )

print("\nSaved:")
print("  robustness_by_chemical_strength.csv")
print("  robustness_by_gap_strength.csv")
print("  robustness_by_decay.csv")
