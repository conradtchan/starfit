Python package for matching stellar abundance measurements against a database of model stellar explosions. Based on the [old IDL code](https://2sn.org/starfit/).

StarFit can match combined abundances of multiple models. For single stars and combinations of two stars, a complete search can be found. For three or more stars, the problem is extremely expensive, so a [Genetic Algorithm](https://en.wikipedia.org/wiki/Genetic_algorithm) is used to find an approximate solution.

# Installation
Tested with Python 3.10

Optional: A working LaTeX installation and dvipng is required to create plots with LaTeX labels (ideal for publication). Otherwise, Matplotlib's default MathText is used, which may not render all symbols correctly.

## From PyPI (recommended)
```
pip install starfit
```
The PyPI package includes the necessary data files.

## From git repo
The data files are not included into the Git repo, and must first be downloaded from the web-server before installing from the Git repo.
```
./download-data.sh
pip install .
```

# Usage
`starfit.Single` fits an abundance pattern to a single model from the database.

Optional arguments:
- `combine`: a list of lists of element charge numbers to treat as combined abundances (e.g. combine the CNO elements)
- `z_max`: highest element charge number to fit
- `z_exclude`: element charge numbers to exclude from fit
- `z_lolim`: elements that are *model* lower limits (effectively the same as *observational* upper limits)
- `upper_lim`: include observational upper limits in data fitting
- `cdf`: use the uncertainty of upper limits to calculate a cumulative distribution function when calculating error contribution (otherwise treat the upper limit as a simple one-sided chi-squared error)
```
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

print(s)
```

The **database indices** of the best fitting models (sorted from best to worst) are given by:
```
s.sorted_stars['index']
```

The corresponding **reduced chi^2** values are:
```
s.sorted_fitness
```

The physical properties of the models corresponding to these indices can be accessed using the database:
```
i_bestfit = s.sorted_stars['index'][0]
s.db.fielddata[i_bestfit]
```

The chemical yield of the models (and the respective element names) are:
```
list( zip(s.list_db, s.full_abudata[:,i_bestfit]) )
```

To make the same plots as the web version:
```
s.plot()
```

`starfit.Double` fits an abundance pattern to a combination of two models from the database. This takes approximately 1 hour, depending on your machine.

Additional arguments:
- `fixed`: Use dilution factors based on the ejecta mass, rather than solving for the optimal dilution ratio of each explosion independently (decreases solve time)
```
s = starfit.Single(
    filename = 'HE1327-2326.dat',
    db = 'znuc2012.S4.star.el.y.stardb.gz',
    combine = [[6, 7, 8]],
    z_max = 30,
    z_exclude = [3, 24, 30],
    z_lolim = [21, 29],
    upper_lim = True,
    cdf = True,
    fixed = False,
)
```

`starfit.Ga` fits an abundance pattern to a combination of two or more models from the database. The solution is approximate, but approaches the best solution with increased run time.

Additional arguments:
- `time_limit`: amount of time (in seconds) to search for solution
- `sol_size`: number of explosion models to find combinations of
- `pop_size`: GA parameter - number of solutions in the population
- `tour_size`: GA parameter - number of solutions per tournament selection
- `frac_mating_pool`: GA parameter - fraction of solutions in the mating pool
- `frac_elite`: GA parameter - top fraction of elite solutions
- `mut_rate_index`: GA parameter - mutation rate of the database index
- `mut_rate_offset`: GA parameter - mutation rate of the dilution factor
- `mut_offset_magnitude`: GA parameter - size of the mutation of the dilution factor
- `local_search`: GA parameter - solve for the best dilution factors rather than relying on the GA

The default GA parameters should be used unless you really know what you are doing.

```
s = starfit.Ga(
    filename = 'HE1327-2326.dat',
    db = 'znuc2012.S4.star.el.y.stardb.gz',
    combine = [[6, 7, 8]],
    z_max = 30,
    z_exclude = [3, 24, 30],
    z_lolim = [21, 29],
    upper_lim = True,
    cdf = True,
    time_limit=20,
    sol_size=3,
)
```

# Contributing to StarFit
Contributions to the StarFit code are welcome. The `master` branch is protected and cannot be committed to directly. Instead, please create a Pull Request with your proposed contributions. Two automated checks (on Github Actions) must be passed:
1. Code formatting using pre-commit. To ensure your changes are compliant with this project's linters, we recommend installing pre-commit prior to making any commits locally.
```
pip install pre-commit
pre-commit install
```
If you have already made non-compliant commits prior to installing pre-commit, then the pre-commit check on Github will fail. To make the code compliant again, run
```
pre-commit run --all
```
and include these changes in a follow-up commit.

2. Code tests using `pytest`. New tests can be added to the `tests/` directory.

## Adding new data files
Data files specified in the `.hashlist` files in `src/starfit/data/db,ref,stars` are downloaded from the web server. To add new data files:
1. Add the new files to the web server hosting the data files at `/var/www/html/data`
2. Generate the hash using `shasum -a 256` (or `sha256sum`)
3. Add an entry into the hash list

When adding new databases into `data/db`, add corresponding labels into the file `data/db/labels`. This label is used in the web application.

# Publishing to PyPI
Github releases will automatically be published to https://pypi.org/project/starfit/
