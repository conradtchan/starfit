# Installation
Tested with Python 3.10, 3.11

Optional: A working `LaTeX` installation and `dvipng` is required to create plots with LaTeX labels (ideal for publication). Otherwise, `Matplotlib`'s default `MathText` is used, which may not render all symbols correctly.

## From PyPI (recommended)
```shell
pip install starfit
```
The PyPI package includes the necessary data files.

## Developer instructions
The data files are not included into the Git repo, and must first be downloaded from the web-server before installing from the Git repo.
```shell
git clone git@github.com:conradtchan/starfit.git
cd starfit

# Download data files
./download-data.sh

# Set environment variable to allow for editable installs
export SETUPTOOLS_ENABLE_FEATURES="legacy-editable"

# "-e" creates an editable install, "[testing]" installs additional dependencies for testing
pip3 install -e .[testing]

# Run all tests
python -m pytest
