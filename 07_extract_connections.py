
from cect.readers import White_whole

ds = White_whole.get_instance()

print("=== CONNECTION COLLECTIONS ===")

result = ds.get_neuron_to_neuron_conns()

print("Returned type:", type(result))
print("Number of returned collections:", len(result))

for i, collection in enumerate(result):
    print(f"\nCollection {i}:")
    print("Type:", type(collection))
    print("Length:", len(collection))

    if isinstance(collection, (list, tuple)) and collection:
        first = collection[0]
        print("First item type:", type(first))
        print("First item:", first)

        print("\nRecord fields:")
        for field in (
            "pre_cell",
            "post_cell",
            "number",
            "synclass",
            "syntype",
            "neurotransmitter",
        ):
            if hasattr(first, field):
                print(f"  {field} = {getattr(first, field)}")

print("\n=== AVAILABLE SYNAPSE CLASSES ===")
print(ds.get_synclasses_with_connections())
