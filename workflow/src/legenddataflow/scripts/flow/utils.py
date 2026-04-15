from pathlib import Path

def replace_path(d, old_path, new_path):
    if isinstance(d, dict):
        for k, v in d.items():
            d[k] = replace_path(v, old_path, new_path)
    elif isinstance(d, list):
        for i in range(len(d)):
            d[i] = replace_path(d[i], old_path, new_path)
    elif isinstance(d, str) and old_path in d:
        d = d.replace(old_path, new_path)
        d = d.replace(new_path, f"$_/{Path(new_path).name}")
    return d