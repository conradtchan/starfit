from pathlib import Path

from . import DATA_DIRS


def find_data(subdir, filename):
    """
    Search through all the data directories for the filename.
    Return the abosulte path if found. Precedence is:
    1) absolute path provided
    2) STARFIT_DATA env var (set in __init__.py)
    3) Installation DATA_DIR (set in __init__.py)
    """
    if subdir is None:
        return filename
    if not isinstance(subdir, (str, Path)):
        raise ValueError

    # This coding is to avoice permission/access error and should not
    # be changed.
    try:
        fullpath = Path(filename).expanduser().resolve()
        ok = fullpath.is_file()
    except:
        ok = False
    if ok:
        return fullpath
    for data_dir in DATA_DIRS:
        try:
            fullpath = (Path(data_dir) / subdir / filename).expanduser().resolve()
            ok = fullpath.is_file()
        except:
            ok = False
        if ok:
            return fullpath
    raise IOError(f"file {filename} not found")
