import numpy as np


def print_nested(d, indent=0):
    for key, value in d.items():
        prefix = "  " * indent
        if isinstance(value, dict):
            print(f"{prefix}{key}:")
            print_nested(value, indent + 1)
        elif (
            isinstance(value, list)
            and value
            and all(isinstance(i, dict) for i in value)
        ):
            print(f"{prefix}{key}: list ({len(value)})")
            print_nested(value[0], indent + 1)
        elif isinstance(value, np.ndarray):
            print(f"{prefix}{key}: ndarray {value.shape}")
        else:
            print(f"{prefix}{key}: {type(value).__name__}")
