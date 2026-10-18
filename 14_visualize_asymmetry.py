
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Load the left/right comparison produced by script 13.
df = pd.read_csv("left_right_neuron_comparison.csv")

# Exclude the directly stimulated sensory neurons to focus
# on downstream response differences.
downstream = df[
    ~df["neuron"].isin(["ASHL", "ASHR"])
].copy()

downstream = downstream.sort_values(
    "ASHR_minus_ASHL",
    ascending=True,
)

# Plot all downstream neurons.
values = downstream["ASHR_minus_ASHL"].to_numpy()
names = downstream["neuron"].to_numpy()

fig, ax = plt.subplots(figsize=(10, 12))

image = ax.imshow(
    values.reshape(-1, 1),
    aspect="auto",
    cmap="coolwarm",
    interpolation="nearest",
    vmin=-max(abs(values.min()), abs(values.max())),
    vmax=max(abs(values.min()), abs(values.max())),
)

ax.set_yticks(np.arange(len(names)))
ax.set_yticklabels(names, fontsize=6)
ax.set_xticks([0])
ax.set_xticklabels(["ASHR − ASHL"])

ax.set_title("Modeled Left/Right Neural Response Difference")
ax.set_ylabel("Neuron")

fig.colorbar(
    image,
    ax=ax,
    label="Difference in integrated activity",
)

fig.tight_layout()
fig.savefig("14_left_right_asymmetry.png", dpi=180)
plt.close(fig)

# Save the 20 strongest differences in each direction.
right_preferred = downstream.sort_values(
    "ASHR_minus_ASHL",
    ascending=False,
).head(20)

left_preferred = downstream.sort_values(
    "ASHR_minus_ASHL",
    ascending=True,
).head(20)

right_preferred.to_csv(
    "14_right_higher_response.csv",
    index=False,
)

left_preferred.to_csv(
    "14_left_higher_response.csv",
    index=False,
)

print("=== TOP 10 RIGHT-HIGHER RESPONSES ===")
print(
    right_preferred.head(10)[
        ["neuron", "ASHL_integrated_activity",
         "ASHR_integrated_activity", "ASHR_minus_ASHL"]
    ].to_string(index=False)
)

print("\n=== TOP 10 LEFT-HIGHER RESPONSES ===")
print(
    left_preferred.head(10)[
        ["neuron", "ASHL_integrated_activity",
         "ASHR_integrated_activity", "ASHR_minus_ASHL"]
    ].to_string(index=False)
)

print("\nSaved:")
print("  14_left_right_asymmetry.png")
print("  14_right_higher_response.csv")
print("  14_left_higher_response.csv")
