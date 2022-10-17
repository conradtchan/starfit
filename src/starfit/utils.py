import os
from pathlib import Path
from . import DATA_DIRS


def find_data(filename):
    fullpath = Path(filename).expanduser().resolve()
    if os.path.isfile(fullpath):
        return fullpath
    for data_dir in DATA_DIRS:
        files = list(Path(data_dir).glob(f"**/{filename}"))
        if len(files) > 1:
            raise RuntimeError(f"Multiple files with matching names found: {files}")
        elif len(files) == 1:
            return files[0]
    raise IOError(f"file {filename} not found")
