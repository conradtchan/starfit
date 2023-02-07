# Installation
Tested with Python [3.10](https://www.python.org/downloads/release/python-3100/), [3.11](https://www.python.org/downloads/release/python-3110/)

Optional: A working [LaTeX](https://www.latex-project.org>) installation and [dvipng](https://ctan.org/pkg/dvipng) is required to create plots with LaTeX labels (ideal for publication). Otherwise, [Matplotlib](https://matplotlib.org>)'s default [MathText](https://matplotlib.org/stable/tutorials/text/mathtext.html) is used, which may not render all symbols correctly.

## From PyPI (recommended)
```shell
pip install starfit
```
The [PyPI package](https://pypi.org/project/starfit/) includes the necessary data files.

## Developer instructions
The data files are not included into the Git repo, and must first be downloaded from the web-server ([starfit.org/data](https://starfit.org/data)) before installing from the Git repo.
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
