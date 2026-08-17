from pathlib import Path
import numpy as np

data_root = Path('/Users/julian/local_files/pata_kim_data/RNAscope_tdrd7_nanos/collected_data')
all_lsm_files = [f for f in list(data_root.glob('*.lsm')) if not f.name.startswith('.') and f.name.endswith('.lsm')]

# Rearrange files in folders by sample
for lsm_file in all_lsm_files:
    condition = 'ctrl' if "F125" in lsm_file.stem.lower() or "dmso" in lsm_file.stem.lower() else 'exp'
    new_name = f"{hex(np.random.randint(0, 2 ** 32 - 1))}__{condition}__{lsm_file.stem}"
    sample_folder = data_root / new_name
    sample_folder.mkdir(exist_ok=True, parents=True)
    new_path = sample_folder / ''.join([new_name, lsm_file.suffix])
    #copy the file to the new folder
    lsm_file.rename(new_path)
