
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(".")
OUT = ROOT / "research_report"
OUT.mkdir(exist_ok=True)

FILES = {
    "robustness": ROOT / "robustness_results.csv",
    "ablation": ROOT / "coupling_ablation_results.csv",
    "threshold": ROOT / "threshold_sensitivity_results.csv",
    "propagation": ROOT / "propagation_summary.csv",
}

for name, path in FILES.items():
    if not path.exists():
        raise FileNotFoundError(f"Missing {name} file: {path.resolve()}")

robustness = pd.read_csv(FILES["robustness"])
ablation = pd.read_csv(FILES["ablation"])
threshold = pd.read_csv(FILES["threshold"])
propagation = pd.read_csv(FILES["propagation"])

def save_table(df, filename):
    df.to_csv(OUT / filename, index=False)

def fmt(value):
    return f"{value:.6f}" if pd.notna(value) else "NA"

# ------------------------------------------------------------
# 1. Robustness summary
# ------------------------------------------------------------
robust_summary = robustness["cosine_similarity"].describe().to_frame(
    "cosine_similarity"
).reset_index(names="statistic")
save_table(robust_summary, "report_robustness_summary.csv")

robust_by_parameter = (
    robustness.groupby("chemical_strength")["cosine_similarity"]
    .agg(["count", "mean", "min", "max"])
    .reset_index()
)
save_table(robust_by_parameter, "report_robustness_by_chemical_strength.csv")

# ------------------------------------------------------------
# 2. Coupling ablation summary
# ------------------------------------------------------------
ablation_summary = (
    ablation.groupby("mode")
    .agg(
        combinations=("cosine_similarity", "count"),
        mean_similarity=("cosine_similarity", "mean"),
        min_similarity=("cosine_similarity", "min"),
        max_similarity=("cosine_similarity", "max"),
        mean_left_integrated=(
            "left_integrated_network_activity", "mean"
        ),
        mean_right_integrated=(
            "right_integrated_network_activity", "mean"
        ),
    )
    .reset_index()
)
save_table(ablation_summary, "report_coupling_ablation_summary.csv")

# ------------------------------------------------------------
# 3. Threshold sensitivity
# ------------------------------------------------------------
threshold_summary = threshold[
    [
        "stimulus", "mode", "threshold",
        "neurons_ever_above_threshold",
        "neurons_above_threshold_at_end",
        "neurons_with_integrated_activity_above_threshold",
    ]
].copy()
save_table(threshold_summary, "report_threshold_summary.csv")

# ------------------------------------------------------------
# 4. Propagation summary
# ------------------------------------------------------------
save_table(propagation, "report_propagation_summary.csv")

# ------------------------------------------------------------
# 5. Figure: robustness by chemical strength
# ------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(
    robust_by_parameter["chemical_strength"],
    robust_by_parameter["mean"],
    marker="o",
)
ax.set_xlabel("Chemical strength")
ax.set_ylabel("Mean left/right cosine similarity")
ax.set_title("Robustness across chemical-strength settings")
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(OUT / "figure_robustness.png", dpi=180)
plt.close(fig)

# ------------------------------------------------------------
# 6. Figure: coupling ablation
# ------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 4.5))
ax.bar(
    ablation_summary["mode"],
    ablation_summary["mean_similarity"],
)
ax.set_ylabel("Mean left/right cosine similarity")
ax.set_title("Coupling ablation comparison")
ax.tick_params(axis="x", rotation=20)
ax.grid(axis="y", alpha=0.3)
fig.tight_layout()
fig.savefig(OUT / "figure_coupling_ablation.png", dpi=180)
plt.close(fig)

# ------------------------------------------------------------
# 7. Figure: threshold sensitivity, full model
# ------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 4.5))
for stimulus in ["ASHL", "ASHR"]:
    subset = threshold[
        (threshold["stimulus"] == stimulus)
        & (threshold["mode"] == "full")
    ].sort_values("threshold")
    ax.plot(
        subset["threshold"],
        subset["neurons_ever_above_threshold"],
        marker="o",
        label=stimulus,
    )

ax.set_xscale("log")
ax.set_xlabel("Activity threshold (log scale)")
ax.set_ylabel("Neurons crossing threshold")
ax.set_title("Sensitivity to activity threshold: full model")
ax.legend()
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(OUT / "figure_threshold_sensitivity.png", dpi=180)
plt.close(fig)

# ------------------------------------------------------------
# 8. Research report
# ------------------------------------------------------------
robust_values = robustness["cosine_similarity"].dropna()
mode_order = ["full", "chemical_only", "gap_only", "no_recurrent_coupling"]
ablation_summary["mode"] = pd.Categorical(
    ablation_summary["mode"], mode_order, ordered=True
)
ablation_summary = ablation_summary.sort_values("mode")

lines = [
    "# C. elegans Connectome Simulation: Results Report",
    "",
    "## 1. Scope",
    "",
    "This report summarizes the existing simulation output files. "
    "No simulations are rerun by this script.",
    "",
    "The model uses the White et al. 1986 connectome loaded through CECT. "
    "It is an exploratory rate-based network model, not a validated "
    "physiological model of the worm.",
    "",
    "## 2. Robustness analysis",
    "",
    f"- Parameter combinations represented: {len(robustness)}",
    f"- Mean left/right cosine similarity: {fmt(robust_values.mean())}",
    f"- Minimum similarity: {fmt(robust_values.min())}",
    f"- Maximum similarity: {fmt(robust_values.max())}",
    "",
    "### By chemical strength",
    "",
    robust_by_parameter.to_markdown(index=False),
    "",
    "## 3. Coupling ablation",
    "",
    "The following table compares the full model with chemical-only, "
    "gap-only, and no-recurrent-coupling configurations.",
    "",
    ablation_summary.to_markdown(index=False),
    "",
    "## 4. Propagation summary",
    "",
    propagation.to_markdown(index=False),
    "",
    "## 5. Threshold sensitivity",
    "",
    "The number of neurons counted as active changes with the chosen "
    "numerical threshold. These thresholds are analytical cutoffs, "
    "not experimentally established firing thresholds.",
    "",
    threshold_summary.to_markdown(index=False),
    "",
    "## 6. Figures",
    "",
    "- `figure_robustness.png`: response similarity by chemical strength.",
    "- `figure_coupling_ablation.png`: response similarity across coupling modes.",
    "- `figure_threshold_sensitivity.png`: active-neuron count versus threshold.",
    "",
    "## 7. Interpretation and limitations",
    "",
    "1. Chemical-only coupling produced broader activity propagation than "
    "gap-only coupling under the tested equations and parameters.",
    "2. Left/right stimulation produced different neuron-level response "
    "profiles in this implementation.",
    "3. Similarity and active-neuron counts depend on model parameters and "
    "the activity threshold.",
    "4. The model treats chemical connections as positive; real chemical "
    "synapses cannot all be assumed excitatory.",
    "5. Gap-junction symmetrization, connection normalization, activity "
    "dynamics, and parameter values are modeling assumptions.",
    "6. These results establish behavior of this implementation, not "
    "experimental biological causality or validated worm behavior.",
    "",
    "## 8. Reproducibility",
    "",
    "The CSV tables in this folder are derived from the existing result files. "
    "The original simulation scripts and source CSVs should be retained "
    "alongside this report for full reproducibility.",
    "",
]

report_path = OUT / "research_report.md"
report_path.write_text("\n".join(lines), encoding="utf-8")

print("\nREPORT GENERATED")
print(f"Folder: {OUT.resolve()}")
print(f"Report: {report_path.resolve()}")
print("\nFiles:")
for path in sorted(OUT.iterdir()):
    print(f"  {path.name}")

print("\nRobustness:")
print(f"  Mean similarity: {fmt(robust_values.mean())}")
print(f"  Min similarity:  {fmt(robust_values.min())}")
print(f"  Max similarity:  {fmt(robust_values.max())}")

print("\nCoupling ablation:")
print(
    ablation_summary[
        ["mode", "mean_similarity", "mean_left_integrated",
         "mean_right_integrated"]
    ].to_string(index=False)
)
