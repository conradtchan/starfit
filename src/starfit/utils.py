from pathlib import Path

from . import DATA_DIRS


def find_data(filename, dtype):
    """
    Search through all the data directories for the filename.
    Return the abosulte path if found. Precedence is:
    1) absolute path provided
    2) STARFIT_DATA env var
    3) Installation DATA_DIR
    """
    if dtype in ("stars", "star"):
        dtype = "stars"
    elif dtype in ("data_base", "database", "db", "dbs"):
        dtype = "db"
    elif dtype in ("ref", "refs", "reference", "references"):
        dtype = "ref"
    else:
        raise ValueError

    fullpath = Path(filename).expanduser().resolve()
    if fullpath.is_file():
        return fullpath
    for data_dir in DATA_DIRS:
        files = list(Path(data_dir).glob(f"{dtype}/{filename}"))
        if len(files) > 1:
            raise RuntimeError(f"Multiple files with matching names found: {files}")
        elif len(files) == 1:
            return files[0]
    raise IOError(f"file {filename} not found")
