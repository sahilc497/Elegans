
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from cect.readers.White_whole import get_instance


# --------------------------------------------------
# 1. Load the real connectome
# --------------------------------------------------
ds = get_instance()
nodes = list(ds.nodes)

W_gj = np.asarray(ds.connections["Generic_GJ"], dtype=float)
W_cs = np.asarray(ds.connections["Generic_CS"], dtype=float)

n = len(nodes)

assert W_gj.shape == (n, n)
assert W_cs.shape == (n, n)

print("Neurons/cells:", n)
print("Chemical connections:", np.count_nonzero(W_cs))
print("Gap-junction connections:", np.count_nonzero(W_gj))


# --------------------------------------------------
# 2. Normalize the connectivity matrices
# Rows represent source cells; columns represent targets.
# Normalization prevents large weights dominating the model.
# --------------------------------------------------
def normalize_rows(W):
    W = W.copy()
    row_sums = W.sum(axis=1, keepdims=True)
    return np.divide(
        W,
        row_sums,
        out=np.zeros_like(W),
        where=row_sums != 0,
    )


Wcs = normalize_rows(W_cs)
Wgj = normalize_rows(W_gj)

# Symmetrize gap junctions for this simplified model.
Wgj = (Wgj + Wgj.T) / 2.0
Wgj = normalize_rows(Wgj)


# --------------------------------------------------
# 3. Choose sensory cells for an exploratory stimulus
# --------------------------------------------------
sensory_names = ["ASHL", "ASHR"]
sensory_indices = [
    nodes.index(name)
    for name in sensory_names
    if name in nodes
]

if not sensory_indices:
    raise ValueError("Requested sensory cells were not found.")

stimulus = np.zeros(n)
stimulus[sensory_indices] = 1.0

print("Stimulated cells:", [nodes[i] for i in sensory_indices])


# --------------------------------------------------
# 4. Run a simple leaky rate network
# --------------------------------------------------
dt = 0.01
duration = 5.0
steps = int(duration / dt)

tau = 0.15
chemical_gain = 1.0
gap_gain = 0.15

activity = np.zeros(n)
history = np.zeros((steps, n))

for t in range(steps):
    # Positive chemical coupling is a temporary modeling
    # assumption, not a claim about real synaptic polarity.
    chemical_input = activity @ Wcs

    # Approximate diffusive coupling through gap junctions.
    gap_input = activity @ Wgj - activity

    # Apply sensory input during the first second.
    external_input = stimulus if t * dt < 1.0 else np.zeros(n)

    target = (
        chemical_gain * chemical_input
        + gap_gain * gap_input
        + external_input
    )

    # Bounded continuous firing-rate proxy.
    target = np.maximum(target, 0.0)
    target = target / (1.0 + target)

    activity += (dt / tau) * (-activity + target)
    activity = np.clip(activity, 0.0, 1.0)

    history[t] = activity


# --------------------------------------------------
# 5. Plot activity over time
# --------------------------------------------------
times = np.arange(steps) * dt

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

axes[0].imshow(
    history.T,
    aspect="auto",
    origin="lower",
    extent=[times[0], times[-1], 0, n],
    interpolation="nearest",
    cmap="viridis",
)
axes[0].set_title("Activity across the connectome")
axes[0].set_xlabel("Time (s)")
axes[0].set_ylabel("Cell index")

axes[1].plot(times, history[:, sensory_indices].mean(axis=1))
axes[1].set_title("Mean activity of stimulated sensory cells")
axes[1].set_xlabel("Time (s)")
axes[1].set_ylabel("Activity (arbitrary units)")
axes[1].grid(True, alpha=0.3)

fig.tight_layout()

output = Path("neural_activity.png")
fig.savefig(output, dpi=160, bbox_inches="tight")
plt.show()

print("\nSimulation completed.")
print("Peak activity:", float(history.max()))
print("Final mean activity:", float(history[-1].mean()))
print("Saved:", output.resolve())
