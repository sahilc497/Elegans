
import inspect
from cect.readers import White_whole

print("=== WHITE ET AL. READER ===")
print("Source file:", inspect.getfile(White_whole))

print("\n=== get_instance() IMPLEMENTATION ===")
print(inspect.getsource(White_whole.get_instance))

print("\n=== CONNECTION-RELATED FUNCTIONS ===")
for name, obj in inspect.getmembers(White_whole):
    if inspect.isfunction(obj) and any(
        term in name.lower()
        for term in ("connection", "synapse", "neuron", "cell")
    ):
        print(f"\n--- {name} ---")
        try:
            print(inspect.getsource(obj))
        except (OSError, TypeError):
            print("Source unavailable")
