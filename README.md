# StarFit
Match stellar abundance measurements against a database of stellar models.

# Installation
## From PyPI
```
pip install starfit
```

## From git repo
Download data files before installing from the repo
```
./download-data.sh
pip install .
```

# Usage
`starfit.Single` fits an abundance pattern to a single model from the database.

Optional arguments:
- `z_max`: highest element charge number to fit
- `z_exclude`: element charge numbers to exclude from fit
- `z_lolim`: elements that are *model* lower limits (effectively the same as *observational* upper limits)

```
import starfit

s = starfit.Single(
    filename = 'HE1327-2326.dat',
    db = 'znuc2012.S4.star.el.y.stardb.gz',
    z_max = 30,
    z_exclude = [3, 24, 30],
    z_lolim = [21, 29],

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

# Publishing to PyPI
Github releases will automatically be published to https://pypi.org/project/starfit/
