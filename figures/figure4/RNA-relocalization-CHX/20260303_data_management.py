from pathlib import Path

root = Path("/Users/julian/local_files/chx_kim")
out_dir = root / "data"
out_dir.mkdir(exist_ok=True)

list_of_days = [path for path in root.iterdir() if path.is_dir() and not path.name.startswith(".")]
list_of_repeats = []


for day in list_of_days:
    list_of_repeats += [path for path in day.iterdir() if path.is_dir() and not path.name.startswith(".")]
list_of_image_files = []


for repeat in list_of_repeats:
    list_of_image_files += [path for path in repeat.iterdir() if path.is_file() and path.suffix == ".lsm"]

print(f"Found {len(list_of_image_files)} image files to move.")

for file_path in list_of_image_files:

    if file_path.name.startswith("."):
        continue
    if file_path.exists():
        condition = "chx" if "chx" in file_path.parent.name.lower() else "ctrl"
        out_file_path = out_dir / f"{condition}__{file_path.parent.parent.name}_{file_path.parent.name}_{file_path.name}"
        print(f"Moving {file_path} to {out_file_path}")

        # Move the file to data output directory
        out_file_path.write_bytes(file_path.read_bytes())


