"""find best fot for a selection of stars"""

import numpy as np

from .starfit import StarFit


class Direct(StarFit):
    """Find the best fit for stars you specify

    when multiple DBs are specified, provide tuples of (DB, index) where DB is 1-based
    """

    def __init__(
        self,
        *args,
        stars=None,  # A list of the indices of the stars, such as [[114, 3604]]
        offsets=None,
        fig=None,
        fixed_offsets=False,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        if self.show:
            return

        assert stars is not None, "require stars"

        stars = np.array(stars)
        if stars.ndim == 1:
            stars = stars.reshape((1, -1))
        if stars.ndim == 2 and self.db_n == 1:
            stars = np.stack((np.full_like(stars, 0), stars), axis=-1)
        assert stars.ndim == 3, "require database info"
        assert (
            stars.shape[-1] == 2
        ), "expecting sequence of 2 for DB and index as last dimension"
        assert np.all(stars[:, :, 0] >= 0), "DB index starts at 0"
        assert np.all(
            stars[:, :, 0] < self.db_n
        ), f"DB index needs to be < number of data bases ({self.db_n})"

        # convert to actual indices
        stars = self.db_off[stars[:, :, 0]] + stars[:, :, 1]

        sol, fitness = self.run(
            stars=stars, offsets=offsets, fixed_offsets=fixed_offsets
        )

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
