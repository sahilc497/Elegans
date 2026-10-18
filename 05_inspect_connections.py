
import inspect
from cect.readers import White_whole

# Inspect the reader class
Reader = White_whole.WhiteDataReader

print("=== READER METHODS ===")
for name, obj in inspect.getmembers(Reader):
    if not name.startswith("_") and callable(obj):
        try:
            print(f"{name}{inspect.signature(obj)}")
        except (ValueError, TypeError):
            print(name)

# Load the dataset
ds = White_whole.get_instance()

print("\n=== DATASET OBJECT TYPE ===")
print(type(ds))

print("\n=== DATASET ATTRIBUTES ===")
for name in dir(ds):
    if any(term in name.lower() for term in (
        "connection", "neuron", "cell", "synapse", "muscle"
    )):
        if not name.startswith("__"):
            print(name)

# Inspect the connection record class
ConnectionInfo = White_whole.ConnectionInfo

print("\n=== CONNECTION RECORD FIELDS ===")
print("Class:", ConnectionInfo)
try:
    print(inspect.signature(ConnectionInfo))
except (ValueError, TypeError):
    pass

print("\n=== CONNECTION RECORD SOURCE ===")
try:
    print(inspect.getsource(ConnectionInfo))
except (OSError, TypeError):
    print("Source unavailable")
