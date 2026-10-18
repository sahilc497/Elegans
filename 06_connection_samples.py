
from cect.readers import White_whole
from cect.ConnectomeReader import ConnectionInfo

ds = White_whole.get_instance()

print("=== CONNECTION RECORD TYPE ===")
print(ConnectionInfo)

connections = ds.get_neuron_to_neuron_conns()

print("\n=== CONNECTION COLLECTION ===")
print("Type:", type(connections))
print("Count:", len(connections))

print("\n=== FIRST 15 CONNECTIONS ===")
for conn in connections[:15]:
    print(conn)

if connections:
    first = connections[0]
    print("\n=== FIRST RECORD FIELDS ===")
    print(vars(first))

print("\n=== AVAILABLE SYNAPSE CLASSES ===")
print(ds.get_synclasses_with_connections())
