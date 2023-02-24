# Contributing to StarFit
Contributions to the StarFit code are welcome. The `main` branch is protected and cannot be committed to directly. Instead, please create a Pull Request with your proposed contributions.  To make a new branch and set to track `origin`
```shell
git checkout -b <new_branch>
git push --set-upstream origin <new_branch>
```

1. If you changed the Fortran code and want to test locally, remember to re-compile / re-install the package.  First, set the legacy environment variable, and then install as editable packages (see instructions above):
```shell
# Set environment variable to allow for editable installs
export SETUPTOOLS_ENABLE_FEATURES="legacy-editable"

# remove artefacts from previous build
rm -rf ./build
make -C ./src/starfit/fitness clean

# "-e" creates an editable install, "[testing]" installs additional dependencies for testing
pip3 install -e .[testing]
```

To make this step more convenient we provide a `Makefile` in the root directory that does all three steps:
```shell
make
```

If there are issues with the `fitness` sub-module, there is a `Makefile` in its source directory that can be used to compile a test program outside of the python package build process.

Two automated checks (on Github Actions) must be passed (Items 2 and 3):

2. Code formatting using pre-commit. To ensure your changes are compliant with this project's linters, we recommend installing pre-commit prior to making any commits locally.
```shell
pip install pre-commit
pre-commit install
```
If you have already made non-compliant commits prior to installing pre-commit, then the pre-commit check on `GitHub` will fail. To make the code compliant again, run
```shell
pre-commit run --all
```

3. Code tests using `pytest`.  New tests can be added to the `tests/` directory.

Run these tests as a check
```shell
python -m pytest
```
and include any necessary changes in the commit.


## Development branch

Development branches are generated and uploaded to Test PyPI if the version number ends in `.dev*` where `*` can be blank or a optional number.  For example, '`0.3.11.dev22`.
They may also be flagged as pre-releases.

To install packages from Test PyPI use
```shell
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ starfit
```
You may include the `-pre` flag or specify a specific version.

## Adding new database files
Database files specified in the `.hashlist` files in `src/starfit/data/db,ref,stars` are downloaded from the web server. To add new data files:
1. Add the new files to the web server hosting the data files at `/var/www/html/data`
2. Generate the hash using `shasum -a 256` (or `sha256sum`)
3. Add an entry into the hash list

When adding new databases into `data/db`, add corresponding labels into the file `data/db/labels` and a description of the data base into the file `data/db/databases` on the web server.

## Creating database files
New database files can be made using the `StarDB` class in `autils/stardb.py`.  A demonstration may be found at  `src/starfit/example/lc12_stardb.py`.  This file serves as a demonstration only and will not work as is.

# Publishing to PyPI
Github releases will automatically be published to [pypi.org/project/starfit/](https://pypi.org/project/starfit/)
