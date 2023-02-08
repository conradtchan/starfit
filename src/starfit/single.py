import time

import numpy as np

from .autils.human import time2human
from .starfit import StarFit, get_fitness


class Single(StarFit):
    """
    Find the best fit for single stars (complete search)

    Example::

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

    """

    #: ``Single`` has ``sol_size`` hardcoded to 1.
    #: All other arguments are inherited from :clas:`starfit.StarFit`
    sol_size = 1

    def __init__(self, *args, **kwargs):
        """
        Create an instance of the single-star fitting object and calculate the
        best fit.
        """
        super().__init__(*args, **kwargs)
        if self.show:
            return

        # Initialize array
        stars = np.recarray(
            (self.db_size, 1), dtype=[("index", np.int64), ("offset", np.float64)]
        )
        stars.offset = 1.0e-4
        stars.index[:, 0] = np.arange(self.db_size)

        t_start = time.perf_counter()
        fitness = get_fitness(
            trimmed_db=self.trimmed_db,
            eval_data=self.eval_data,
            z_exclude_index=self.exclude_index,
            sol=stars,
            cdf=self.cdf,
            dst=self.dst,
            limit_solver=self.limit_solver,
            limit_solution=self.limit_solution,
            local_search=False,
        )
        self.logger.info("Calculation time:")
        self.logger.info(f" {time2human(time.perf_counter() - t_start)}")

        # Identify the best fit by finding the entry with the smallest error
        self.unsorted_fitness = fitness
        self.unsorted_stars = stars
        sort_index = np.argsort(fitness)
        self.sorted_fitness = fitness[sort_index]
        self.sorted_stars = stars[sort_index]
        self.rank = sort_index

        self.logger.info("Best solution:")
        self.logger.info(" " + str(self.sorted_stars[0, 0]))
        self.logger.info("Fitness:")
        self.logger.info(" " + str(self.sorted_fitness[0]))

        self.close_logger(timing="Finished in")
