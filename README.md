Python package for matching stellar abundance measurements against a database of model stellar explosions. Can match combined abundances of multiple models. For single stars and combinations of two stars, a complete search can be found. For three or more stars, the problem is extremely expensive, so a [Genetic Algorithm](https://en.wikipedia.org/wiki/Genetic_algorithm) is used to find an approximate solution.

# Installation
## From PyPI (recommended)
```
pip install starfit
```
The PyPI includes the necessary data files.

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
Data files specified in the `.hashlist` files in `src/starfit/data/db,ref,stars` are downloaded from the web server. To add new data files, add them to the web server hosting the data files, generate the hash using `shasum`, and add an entry into the hash list.

# Publishing to PyPI
Github releases will automatically be published to https://pypi.org/project/starfit/
