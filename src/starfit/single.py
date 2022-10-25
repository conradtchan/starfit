import time

import numpy as np

from .autils.human import time2human
from .starfit import StarFit, get_fitness


class Single(StarFit):
    """Find the best fit for single stars (complete search)"""

    sol_size = 1

    def __init__(
        self,
        filename,
        db="znuc2012.S4.star.el.y.stardb.gz",
        *,
        combine=[[]],
        z_exclude=[],
        z_min=1,
        z_max=999,
        z_lolim=[],
        upper_lim=True,
        silent=False,
        cdf=True,
        y_floor=1.0e-99,
    ):
        self.silent = silent
        super().__init__()

        self._setup(
            filename=filename,
            database=db,
            combine=combine,
            z_exclude=z_exclude,
            z_min=z_min,
            z_max=z_max,
            upper_lim=upper_lim,
            z_lolim=z_lolim,
            y_floor=y_floor,
        )

        self.gene_size = 1
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
            cdf=cdf,
            ls=True,
        )
        self.logger.info("Calculation time:")
        self.logger.info(f" {time2human(time.perf_counter() - t_start)}")

        # Identify the best fit by finding the entry with the smallest error
        self.unsorted_fitness = fitness
        sort_index = np.argsort(fitness)
        self.sorted_fitness = fitness[sort_index]
        self.sorted_stars = stars[sort_index]

        self.logger.info("Best solution:")
        self.logger.info(" " + str(self.sorted_stars[0, 0]))
        self.logger.info("Fitness:")
        self.logger.info(" " + str(self.sorted_fitness[0]))

        self.close_logger(timing="Finished in")