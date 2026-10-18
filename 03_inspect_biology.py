
import numpy as np
from cect.readers.White_whole import get_instance
from cect import ConnectomeDataset as CD

ds = get_instance()
nodes = list(ds.nodes)

print("=== CONNECTOME BIOLOGY INSPECTION ===")
print("Cells in dataset:", len(nodes))

print("\n=== CONNECTION TYPES ===")
for name, matrix in ds.connections.items():
    matrix = np.asarray(matrix)
    print(f"\n{name}")
    print("  Shape:", matrix.shape)
    print("  Non-zero entries:", np.count_nonzero(matrix))
    print("  Weight range:",
          (matrix[matrix > 0].min(), matrix.max())
          if np.any(matrix > 0) else "No positive entries")

print("\n=== CECT SYNAPSE CLASS DEFINITIONS ===")
print("Available synapse classes:")
for name in dir(CD):
    if "SYN_CLASS" in name or "SYN_TYPE" in name:
        value = getattr(CD, name)
        print(f"\n{name}:")
        print(value)

print("\n=== SAMPLE NEURON NAMES ===")
for name in nodes[:40]:
    print(name)
