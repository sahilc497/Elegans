
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from cect.readers.White_whole import get_instance


# Load the published White et al. connectome
ds = get_instance()

nodes = list(ds.nodes)
connections = ds.connections

print("\n=== C. ELEGANS CONNECTOME ===")
print("Total dataset nodes:", len(nodes))
print("First 20 nodes:", nodes[:20])
print("Connection types:", list(connections.keys()))

# Report matrix statistics
for connection_type, matrix in connections.items():
    matrix = np.asarray(matrix)

    print(f"\n--- {connection_type} ---")
    print("Shape:", matrix.shape)
    print("Non-zero entries:", np.count_nonzero(matrix))
    print("Sum of entries:", matrix.sum())
    print("Data type:", matrix.dtype)

    nonzero = matrix[matrix != 0]
    if nonzero.size:
        print("Minimum non-zero value:", nonzero.min())
        print("Maximum value:", nonzero.max())

# Visualize each connectivity matrix
fig, axes = plt.subplots(
    1, len(connections), figsize=(8 * len(connections), 7)
)

axes = np.atleast_1d(axes)

for ax, (connection_type, matrix) in zip(
    axes, connections.items()
):
    matrix = np.asarray(matrix)

    image = ax.imshow(
        matrix,
        cmap="viridis",
        interpolation="nearest",
        aspect="auto",
    )

    ax.set_title(connection_type)
    ax.set_xlabel("Target neuron index")
    ax.set_ylabel("Source neuron index")

    fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)

fig.suptitle("C. elegans Connectome — White et al. (1986)")
fig.tight_layout()

output = Path("connectome_matrices.png")
fig.savefig(output, dpi=160, bbox_inches="tight")
plt.show()

print("\nSaved visualization:", output.resolve())
