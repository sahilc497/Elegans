
import numpy as np
import pandas as pd
from cect.readers import White_whole

# Reconstruct the same cell ordering used by the simulation.
ds = White_whole.get_instance()
cells, connections = ds.get_neuron_to_neuron_conns()

# Load the simulation's time-series summary.
df = pd.read_csv("simulation_activity.csv")

print("=== SAVED ACTIVITY SUMMARY ===")
print(df.describe().to_string())

print("\nThe CSV contains network-wide mean and peak activity.")
print("To rank individual neurons, the simulation must also export")
print("the full activity matrix; the current CSV does not contain it.")

print("\n=== SENSORY CELL INDICES ===")
for name in ("ASHL", "ASHR"):
    if name in cells:
        print(f"{name}: index {cells.index(name)}")
    else:
        print(f"{name}: not found")
