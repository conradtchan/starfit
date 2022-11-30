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
        fixed_offsets=False,
        optimize=True,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        if self.show:
            return

        assert stars is not None, "require stars"

        stars = np.atleast_2d(np.array(stars))
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

        if offsets is not None:
            if np.ndim(offsets) == 0 and stars.shape[:] == (1, 1):
                offsets = np.asarray([[offsets]])
            elif np.ndim(offsets) == 1:
                n_offsets = len(offsets)
                if stars.shape[:2] == (n_offsets, 1):
                    offsets = np.asarray(offsets).reshape((n_offsets, 1))
                elif stars.shape[:2] == (1, n_offsets):
                    offsets = np.asarray(offsets).reshape((1, n_offsets))
            offsets = np.asarray(offsets)
            assert (
                offsets.shape == stars.shape[:2]
            ), "offsets dimensions need to match stars dimensions"

        sol, fitness = self.run(
            stars=stars,
            offsets=offsets,
            fixed_offsets=fixed_offsets,
            optimize=optimize,
        )

        self.unsorted_fitness = fitness
        self.unsorted_stars = sol

        sort = np.argsort(fitness)
        self.sorted_fitness = fitness[sort]
        self.sorted_stars = sol[sort]
        self.sol_size = stars.shape[1]

        text = self.format(best=False, n=2**12).splitlines()
        n = len(text[2])
        self.logger.info("")
        self.logger.info("-" * n)
        for t in text:
            self.logger.info(t)
        self.logger.info("-" * n)

        self.close_logger(timing="Finished in")
