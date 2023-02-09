
# Usage

## Single star matches

{py:class}`starfit.Single` fits an abundance pattern to a single model from the database. Example usage:

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
    constraints = 'energy <= 5',
    )

s.print()
```
The {py:meth}`.StarFit.print` method allows to specify the number of lines to be printed (`n`), the offset for the first entry to print (`n0`, default is `0`) and the maximum number of columns to use as a "wide" table (`wide`, default `12`).  A `format` and be specified as `"html"` or `"unicode"`, plain text otherwise.  Default is `unicode`.
```python
s.print(n0=3, n=1, wide=8, format=None)
```

The {py:meth}`.StarFit.info` method allows to print information about individual table entries, starting with index `0` (default).  For example, to print the third model info use
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

{py:meth}`.StarFit.plot` can be used to make the same plots as the web version:
```python
s.plot()
```
If you want to plot a solution other than the best one, use the parameter `num` (default: `0`)  To plot the 5th best solution, skipping the first `4`, use
```python
s.plot(num=4)
```
The legend as well as the star name and copyright string can be moved (dragged).

## Full multi-star search

`starfit.Multi` fits an abundance pattern to a combination of models from the database(s).  This can take a long time as there can be many combinations.

Additional arguments:
- `fixed_offsets`: Use dilution factors based on the ejecta mass, rather than solving for the optimal dilution ratio of each explosion independently (decreases solve time)
- `threads`: Number of threads to use.  Default is to use the CPU count (including hyper-threading)
- `nice`: Nice level of background threads.  Default is 19 (lowest priority on Unix systems).
- `group`: by default, all data are merged in one big list and all possible combinations (excluding duplicates) are explored.  If `group` is specified, only combinations form different databases are considered.  This can significantly reduce the cost and often may be more what is intended.  In this case, the number of data base partitions needs to match the number of stars (`sol_size`).  `group` can be a vector with number of data bases to group into each group.  The number of groups needs match the `sol_size` vector.

Changed arguments:
- sol_size can now be a vector with one entry for each partition.  The number of entries need to match the number of groups.  A scalar value is equivalent to vector with that many `1`s.  All combinations in each group (without repetitions) are tested.
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
- `sol_size`: number of nucleosynthesis models to combine for the solution (default: `2`).
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
- `pin`: number or list of groups to require to be included.  Repetitions are allowed to enforce multiple selections from that group.

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
- `gen`: when set to `True`, plot as a function of generation number.  Otherwise, plot as a function of computational time (default).
```python
s.plot_fitness(gen=True)
```

## Matching specific star combinations

`starfit.Direct` allows to find the best fit to pre-selected group, or groups of stars.

Additional arguments:
- `stars`: Nested list of lists of models.  For each model, specify a list of database and index.  Both, the database index and the star index, are `0`-based.
- `offsets`: Nested list of offsets .  For each model, specify a list of offsets.  If not provided, a default (starting) value (`1e-4` total) will be assumed.
- `optimize`: Whether to find best matching offset or use offsets as is (default: `True`).

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
