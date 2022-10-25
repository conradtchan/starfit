"""find best fot for a selection of stars"""

import numpy as np

from .starfit import StarFit


class Direct(StarFit):
    """Find the best fit for stars you specify

    when multiple DBs are specified, provide tuples of (DB, index) where DB is 1-based
    """

    def __init__(
        self,
        filename,
        db="znuc2012.S4.star.el.y.stardb.gz",
        stars=[[114, 3604]],  # A list of the indices of the stars
        offsets=None,
        z_min=1,
        z_max=999,
        y_floor=1.0e-99,
        fig=None,
        silent=False,
        combine=[[]],
        z_exclude=[],
        z_lolim=[],
        upper_lim=True,
        fixed=False,
        cdf=True,
    ):
        self.silent = silent
        super().__init__()

        assert stars is not None, "require stars"

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

        stars = np.array(stars)
        if stars.ndim == 1:
            stars = stars.reshape((1, -1))
        if stars.ndim == 2 and self.db_n == 1:
            stars = np.stack((np.full_like(stars, 1), stars), axis=-1)
        assert stars.ndim == 3, "require database info"
        assert (
            stars.shape[-1] == 2
        ), "expecting sequence of 2 for DB and index as last dimension"
        assert np.all(stars[:, :, 0] > 0), "DB index starts at 1"
        assert np.all(
            stars[:, :, 0] <= self.db_n
        ), f"DB index needs to be <= {self.db_n}"

        # convert to actual indices
        stars = self.db_off[stars[:, :, 0] - 1] + stars[:, :, 1]

        self.cdf = cdf
        sol, fitness = self.run(stars, offsets, fixed)

        sort = np.argsort(fitness)
        self.sorted_fitness = fitness[sort]
        self.sorted_stars = sol[sort]
        self.sol_size = stars.shape[1]

        text = self.format().splitlines()
        n = len(text[2])
        self.logger.info("")
        self.logger.info("-" * n)
        for t in text:
            self.logger.info(t)
        self.logger.info("-" * n)

        self.close_logger(timing="Finished in")
