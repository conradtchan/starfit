![x](https://github.com/conradtchan/starfit/actions/workflows/test.yml/badge.svg)
![x](https://github.com/conradtchan/starfit/actions/workflows/pre-commit.yml/badge.svg)
![x](https://github.com/conradtchan/starfit/actions/workflows/publish.yml/badge.svg)

![GitHub commits since latest release (by SemVer including pre-releases)](https://img.shields.io/github/commits-since/conradtchan/starfit/latest/main?include_prereleases)
![GitHub Release Date](https://img.shields.io/github/release-date/conradtchan/starfit)

Python package for matching stellar abundance measurements against a database of model stellar explosions. Based on the [old IDL code](https://2sn.org/starfit/) by [Alexander Heger](https://2sn.org).

StarFit can match combined abundances of multiple models. For single stars and combinations of multiple stars, a complete search can be found.  For three or more stars, the problem is extremely expensive, so a [Genetic Algorithm](https://en.wikipedia.org/wiki/Genetic_algorithm) has been implemented by Conrad Chan to efficiently find an approximate solution.

An online interface (with a subset of functionality) is available at [starfit.org](https://starfit.org).

# Installation
Tested with Python 3.10, 3.11

Optional: A working LaTeX installation and dvipng is required to create plots with LaTeX labels (ideal for publication). Otherwise, Matplotlib's default MathText is used, which may not render all symbols correctly.

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
```

# Usage

## Single star matches

`starfit.Single` fits an abundance pattern to a single model from the database.

Required arguments:
- `filename`: filename of star.  Can be absolute or relative path.  The files will also be searched for in the distribution files and in the search path specified by environment variable `STARFIT_DATA` in subdirectory `stars`.
- `db`: database file or tuple of data base files.  String or `Path` object.  Can be absolute or relative path.  Files will also be searched in the distribution files and in the search path specified by environment variable `STARFIT_DATA` in subdirectory `db`.  You may also use the "jokers" (`"*"`) in the data base name.  The code will then try to resolve all matching data bases in the first source directory that contains any matching file.  The plain `*` argument will include all data bases in the first source that contains any data base; the matching is done against the pattern `*.stardb.*`.  The `Ellipis` (`...` Python object, not in quotation marks) will do the same as the plain `*` argument, but will continue searching through all data souces.  This allows for an easy way to search across all model data bases available.

Optional arguments:
- `combine`: a list of lists of element charge numbers to treat as combined abundances (e.g. combine the CNO elements)
- `z_max`: highest element charge number to fit
- `z_exclude`: element charge numbers to exclude from fit
- `z_lolim`: elements that are *model* lower limits (effectively the same as *observational* upper limits)
- `upper_lim`: include observational upper limits in data fitting
- `cdf`: use the uncertainty of upper limits to calculate a cumulative distribution function when calculating error contribution (otherwise treat the upper limit as a simple one-sided &#x1D6D8;&sup2; error)
- `det`: use the detection limits when calculating error contribution (experimental).
- `cov`: use the error covariances when calculating error contribution (experimental).
- `limit_solver`: solver/search will only allow solutions for each star that contribute no more than 100%
- `limit_solution`: solver/search will only allow solutions where the total adds up to no more than 100% contributions from all stars.  Results from the search are renormalised accordingly.
- `y_floor`: floor value for abundances to assume in models (default: `1e.0e-99`).  This is useful for elements not produced in a model, otherwise &#x1D6D8;&sup2; of -&infin; may result.
- `db_label`: a list of labels for the data bases to be used in plots and tables.  Will only be shown if there are more than one databases specified.  If present, needs to match the number of databases specified.  If not present, databases will be numbbeded starting with `0`, unless the `StarDB` has a `label` field that will be used instead.  The maximum label length currently allowed is `8`.
- `show`: show list of loaded databases with label/number and name, then quit.
```python
import starfit

s = starfit.Single(
    filename = 'HE1327-2326.dat',
    db = 'znuc2012.S4.star.el.y.stardb.gz',
    combine = [[6, 7, 8]],
    z_max = 30,
    z_exclude = [3, 24, 30],
    z_lolim = [21, 29],
    upper_lim = True,
    cdf = True,
    )

s.print()
```
the `print` method allows to specify the number of lines to be printed (`n`), the offset for the first entry to print (`n0`, default is `0`) and the maximum number of columns to use as a "wide" table (`wide`, default `12`).  A `format` and be specified as `"html"` or `"unicode"`, plain text otherwise.  Default is `unicode`.
```python
s.print(n0=3, n=1, wide=8, format=None)
```

The `info` method allows to print information about individual table entries, starting with index `0` (default).  For example, to print the third model info uee
```python
s.info(2)
```

The **database indices** of the best fitting models (sorted from best to worst) are given by:
```python
s.sorted_stars['index']
```

The corresponding **reduced &#x1D6D8;&sup2;** values are:
```python
s.sorted_fitness
```

The physical properties of the models corresponding to these indices can be accessed using the database:
```python
i_bestfit = s.sorted_stars['index'][0]
s.db.fielddata[i_bestfit]
```

The chemical yield of the models (and the respective element names) are:
```python
list(zip(s.list_db, s.full_abudata[:, i_bestfit]))
```

### Plots

To make the same plots as the web version:
```python
s.plot()
```
If you want to plot a solution other than the best one, use the parameter `num` (default: `0`)  To plot the 5th best solution, skipping the first `4`, use
```python
s.plot(num=4)
```
The legend as well as the star name and copyright string can be moved (dragged).
Plot parameters iclude
- `num`: Number of solution, from the top (default: `0`).
- `yscale`: select the y-scale of the plot.  Numerical value identical to those used for the star data formats.
- `ynorm`: elemets to use a norm for `[X/Y]' plots (`yscale=3`).
- `multi`: plot this many best solutions as gray lines (default: `0`).  For `multi=-1` lines will be shaded according to relative data point probability based on &#x1D6D8;&sup2; and assuming multi-dimensional Gaussian error.
- `save`: filename to save plot.
- `range_det`: adjust range to include detection thresholds (default: `False`)
- `range_lim`: adjust range to incoude detection limits (default: `True`)
- `range_zmin`: in Z ro conider for determining y range (default: `3`)
- `pad_abu`: fraction of plot range to use at bottom/top boundary (default: `0.1`).
- `pad_det`: fraction of plot range to pad detection tresholds (default: `0.05`).
- `figsize`: dimensions of figure in inches (dedault: `(10, 6)`).
- `dpi`: resolution of image (default: `102`).
- `xlim`: overwrite x range (low, high).
- `ylim`: overwrite y range (low, high).
- `data_size`: Size of data lines and symols (default: `3`).
- `fontsize`: size used for axis labels (default: `12`).
- `annosize`: size used for element symbols (default: `small`).
- `dist`:  distance of lables from data points.
- `fig`: figure object to use as canvas, otherwise as new figure is created.
- `ax`: axis objects to use for drawing, otherwise axis and parent figure are created as needed.
- `xlabel`: overwrite label for x-axis.
- `ylabel`: overwrite label for y-axis.

## Full multi-star search

`starfit.Multi` fits an abundance pattern to a combination of models from the database(s).  This can take a long time as there can be many combinations.

Additional arguments:
- `fixed_offsets`: Use dilution factors based on the ejecta mass, rather than solving for the optimal dilution ratio of each explosion independently (decreases solve time)
- `threads`: Number of threads to use.  Default is to use the CPU count (including hyperthreading)
- `nice`: Nice level of background threads.  Default is 19 (lowest priority on unix systems).
- `group`: by default, all data are merged in one big list and all possible combinations (excluding duplicates) are explored.  If `group` is specified, only combinations form different databases are considered.  This can significantly reduce the cost and often may be more what is intended.  In this case, the number of data base partitons needs to match the number of stars (`sol_size`).  `group` can be a vector with number of data bases to group into each group.  The number of groups needs match the `sol_size` vector.

Changed arguments:
- sol_size can now be a vector with one entry for each partiton.  The number of entries need to match the number of groups.  A scalar value is equivalent to vector with that many `1`s.  All combinations in each group (without repetitions) are tested.
```python
s = starfit.Multi(
    filename = 'SMSS2003-1142.dat',
    db = (
        'he2sn.HW02.star.el.y.stardb.gz',
        'rproc.just15.star.el.y.stardb.xz',
	'rproc.wu.star.el.y.stardb.xz',
        ),
    z_max = 999,
    z_exclude = [3, 24, 30],
    z_lolim = [21, 29],
    upper_lim = True,
    cdf = True,
    fixed_offsets = False,
    sol_size = [2,1],
    group = [1,2],
    )
```

## Genetic algorithm

`starfit.Ga` fits an abundance pattern to a combination of two or more models from the database. The solution is approximate, but approaches the best solution with increased run time.

Additional arguments:
- `gen`: maximum number of generations (iterations) to search; no limit if `0` or `None` (default: `1000`).
- `time_limit`: maximum amount of time (in seconds) to search for solution.  Infinite if `None` (default: `20 s`).
- `sol_size`: number of nucleosynthesis models to combinatine for solution (default: `2`).
- `pop_size`: GA parameter - number of solutions in the population (default: `200`).
- `tour_size`: GA parameter - number of solutions per tournament selection (default: `2`).
- `frac_mating_pool`: GA parameter - fraction of solutions in the mating pool (default: `1`).
- `frac_elite`: GA parameter - top fraction of elite solutions (default: `0.5`).
- `mut_rate_index`: GA parameter - mutation rate of the star index in the databases (default: `0.2`).
- `mut_rate_offset`: GA parameter - mutation rate of the dilution factor (default: `0.1`).
- `mut_offset_magnitude`: GA parameter - size of the mutation of the dilution factor (default `1`).
- `local_search`: GA parameter - solve for the best dilution factors rather than relying on the GA (default: `True`).
- `spread`: GA parameter - ensure no database sources are skipped unless there are fewer stars than data bases.  This can be useful if there is a large disparity in the number of models between the different data bases and if you have a prior that all data bases should be used.  Eventually, the genetic algorithm should find all combinations that match best anyway, however.
- `group`: grouping of data bases, for use with `spread`: try to cover each group but not each database within it separately.  Provide a vector of group length or of tuples with database indices (`0`-based), no duplications allowed.  Same rules as above apply: if group is specified, you need to a provide grouping that covers each database listed by index.
- `pin`: number or list of groups to require to be included.  Repetirions are allowed to enforce multiple selections from that group.

The default GA parameters should be used unless you really know what you are doing.
```python
s = starfit.Ga(
    filename = 'HE1327-2326.dat',
    db = (
        'rproc.just15.star.el.y.stardb.xz',
        'znuc2012.S4.star.el.y.stardb.gz',
	'rproc.wu.star.el.y.stardb.xz',
        ),
    combine = [[6, 7, 8]],
    z_max = 30,
    z_exclude = [3, 24, 30],
    z_lolim = [21, 29],
    upper_lim = True,
    cdf = True,
    time_limit = 20,
    sol_size = 2,
    spread = True,
    group=[[0,2],[1]]
    )
```
The execution can be terminated pressing the `<Enter>` key.

### Evolutionary History
The history of fitness evolution can be plotted using the `plot_fitness` method.

Additional arguments:
- `gen`: when set to `True`, plot as a function of generation number.  Otheriwse, plot as a function of computational time (default).
```python
s.plot_fitness(gen=True)
```

## Matching specific star combinations

`starfit.Direct` allows to find the best fit to pre-selected group, or groups of stars.

Additional arguments:
- `stars`: List of lists of models.  For each model, specify a list of database and index.  Both, the database index and the star index, are `0`-based.

The following selects two groups of models: the first selects model index `0` from the first database with index `0`, and the second model (index `1`) from the second database (index `1`); the second selects the third model (index `2`) from the first database (index `0`) and the fourth model (index `3`) from the second database (index `1`):
```python
s = starfit.Direct(
    filename = 'HE1327-2326.dat',
    db = (
        'he2sn.HW02.star.el.y.stardb.gz',
        'rproc.just15.star.el.y.stardb.xz',
        ),
    stars = [
        [[0,0], [1,1]],
        [[0,2], [1,3]]],
    )
```
The results are sorted by fitness and stored in the returned object as usual, allowing to print and plot the results.

## Multiple databases info

By default, data bases are numbered in the order provided.  The database numbers are only listed when there is more than one database provided.  Full database information can be printed using the `print_comments` method of the solution object:
```python
s.print_comments()
```
or if the `full` parameter is specified to the `print` method
```python
s.print(full=True)
```

## Error matrix plots

The `StarFit` object provide three functions analyse error and plot error contributions.

### Error matrix of data

`plot_star_matrix` plots the error as computed/used by the fitter.

Arguments are
- `zoom`: How much to zoom in around zero.  Use `False` to disable. (default: `1000`)
- `nlab`: How many labels to draw on the colorbar (default: `9`).
- `compress`: Whether to skip elements for which there are no measurements (default: True).

### Inverse of error matrix of data

`plot_star_inverse` plots the inverted error matrix as computed/used by the fitter.

Arguments are
- `zoom`: How much to zoom in around zero.  Use `False` to disable. (default: `0.1`)
- `nlab`: How many labels to draw on the colorbar (default: `9`).
- `compress`: Whether to skip elements for which there are no measurements (default: True).

### Error contributions as computed

`plot_error_matrix` plots the error contributions as computed by the fitter for a given star.

Arguments are
- `num`: Number of solution to plot, counted from the best (default: `0`).
- `zoom`: How much to zoom in around zero.  Use `False` to disable. (default: `1000`)
- `nlab`: How many labels to draw on the colorbar (default: `9`).


# Custom data directory
Custom stellar data and model database files can always be used by providing a full path in the argument. However, users may optionally specify their own data directory using the environment variable `STARFIT_DATA` for convenience:
```shell
export STARFIT_DATA='/your/custom/data'
```
Files found in the custom data directory will take precedence over the default data directory.
Your custom data directory must have the same structure as `src/starfit/data`, i.e. it should contain the `db`, `ref`, and `stars` directories:
```shell
❯ ls
db
ref
stars
```

# Contributing to StarFit
Contributions to the StarFit code are welcome. The `main` branch is protected and cannot be committed to directly. Instead, please create a Pull Request with your proposed contributions.  To make a new branch and set to track `origin`
```shell
git checkout -b <new_branch>
git push --set-upstream origin <new_branch>
```

1. If you changed the Fortran code and want to test locally, remeber to re-compile / re-install the package.  First, set the legacy environment variable, and then install as editable packages (see instructions above):
```shell
# Set environment variable to allow for editable installs
export SETUPTOOLS_ENABLE_FEATURES="legacy-editable"

# remove artefacts from previous build
rm -rf ./build
make -C ./src/starfit/fitness clean

# "-e" creates an editable install, "[testing]" installs additional dependencies for testing
pip3 install -e .[testing]
```

To make this step more convenient we povide a `Makefile` in the root directory that does all three steps:
```shell
make
```

If there are issues with the `fitness` submodule, there is a `Makefile` in its source directory that can be used to compile a test program outside of the python package build process.

Two automated checks (on Github Actions) must be passed (Items 2 and 3):

2. Code formatting using pre-commit. To ensure your changes are compliant with this project's linters, we recommend installing pre-commit prior to making any commits locally.
```shell
pip install pre-commit
pre-commit install
```
If you have already made non-compliant commits prior to installing pre-commit, then the pre-commit check on Github will fail. To make the code compliant again, run
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
Github releases will automatically be published to https://pypi.org/project/starfit/
