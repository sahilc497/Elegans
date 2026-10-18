
from collections import Counter
from cect.readers import White_whole

ds = White_whole.get_instance()

# The method returns two collections:
# (cell_names, connection_records)
cells, connections = ds.get_neuron_to_neuron_conns()

print("=== CONNECTOME SUMMARY ===")
print("Cell names:", len(cells))
print("Connection records:", len(connections))

print("\n=== CONNECTION TYPE COUNTS ===")
type_counts = Counter(
    getattr(conn, "syntype", "Unknown")
    for conn in connections
)
for name, count in type_counts.items():
    print(f"{name}: {count}")

print("\n=== SYNAPSE CLASS COUNTS ===")
class_counts = Counter(
    getattr(conn, "synclass", "Unknown")
    for conn in connections
)
for name, count in class_counts.items():
    print(f"{name}: {count}")

print("\n=== FIRST 10 CHEMICAL SYNAPSES ===")
chemical = [
    conn for conn in connections
    if getattr(conn, "syntype", "") == "Send"
    or getattr(conn, "synclass", "") == "Generic_CS"
]

for conn in chemical[:10]:
    print(
        f"{conn.pre_cell} -> {conn.post_cell} | "
        f"weight/count={conn.number} | "
        f"class={conn.synclass} | "
        f"type={conn.syntype}"
    )

print("\n=== FIRST 10 GAP JUNCTIONS ===")
gap_junctions = [
    conn for conn in connections
    if getattr(conn, "syntype", "") == "GapJunction"
    or getattr(conn, "synclass", "") == "Generic_GJ"
]

for conn in gap_junctions[:10]:
    print(
        f"{conn.pre_cell} <-> {conn.post_cell} | "
        f"weight/count={conn.number} | "
        f"class={conn.synclass} | "
        f"type={conn.syntype}"
    )

print("\n=== AVAILABLE RECORD ATTRIBUTES ===")
if connections:
    first = connections[0]
    for field in (
        "pre_cell", "post_cell", "number",
        "synclass", "syntype", "neurotransmitter"
    ):
        if hasattr(first, field):
            print(f"{field}: {getattr(first, field)}")
