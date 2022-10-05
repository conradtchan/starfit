import numpy as np
from scipy.interpolate import interp2d

from .autils.stardb import StarDB


class TrimDB(StarDB):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        try:
            parameter = args[1]
            self.trim(parameter)
        except:
            pass

    def trim(self, parameter):
        """
        get_star_slice is probably better than this
        """
        # Format:
        # parameter = [['mixing', 0, 0.1], ['mass', 10, 20]]

        mask = np.ones_like(self.data.transpose()[1]).astype(int)
        for p in parameter:
            mask *= p[1] <= self.field_data[p[0]]
            mask *= p[2] >= self.field_data[p[0]]

        ind = np.where(mask)[0]
        abu_data = self.data.transpose()[:, ind]
        field_data = self.field_data[ind]

        self.data = np.copy(abu_data.transpose())
        self.field_data = np.copy(field_data)

    def grid(
        self,
    ):
        """
        return a grid arranged by parameters
        """

        param_shape = self.nvalues[2::-1]

        param_grid = self.field_data.reshape(param_shape)
        index_grid = np.arange(self.nstar).reshape(param_shape)
        abu_grid = self.data.transpose().reshape((83,) + tuple(param_shape))

        return abu_grid, index_grid, param_grid, param_shape

    def inter_db(self, new_mixing=100, new_energy=100):
        # [mixing, energy, mass]
        abu_grid, index_grid, param_grid, shape = self.grid()

        n_mass = shape[2]
        n_el = self.data.transpose().shape[0]

        interpolators = np.ndarray((n_mass, n_el), dtype=object)
        rem_interpolator = np.ndarray((n_mass,), dtype=object)

        # +1 to prevent infinities, -1 later on
        x = param_grid["energy"][0, :, 0]
        y = param_grid["mixing"][:, 0, 0]
        for m in range(n_mass):
            for el in range(n_el):
                z = abu_grid[el, :, :, m]
                rem = param_grid["remnant"][:, :, m]
                interpolators[m, el] = interp2d(x, y, z, "linear")
                rem_interpolator[m] = interp2d(x, y, rem, "linear")

        new_size = {"mixing": new_mixing, "energy": new_energy}
        names = ["mixing", "energy"]
        new_params = {}
        for n in names:
            p = np.sort(np.unique(param_grid[n]))
            if p[0] == 0:
                logmin = np.log10(p[1])
                inf_start = True
            else:
                logmin = np.log10(p[0])
                inf_start = False
            logmax = np.log10(p[-1])

            new_params[n] = np.logspace(
                logmin, logmax, new_size[n] - inf_start, base=10
            )
            if inf_start:
                new_params[n] = np.concatenate(([0], new_params[n]))

        inter_grid = np.ndarray((n_el, new_size["mixing"], new_size["energy"], n_mass))
        field_grid = np.ndarray(
            (new_size["mixing"], new_size["energy"], n_mass),
            dtype=[
                ("mass", "f8"),
                ("energy", "f8"),
                ("mixing", "f8"),
                ("remnant", "f8"),
            ],
        )

        field_grid["mass"] = param_grid["mass"][0, 0, :]
        field_grid["energy"] = new_params["energy"][np.newaxis, :, np.newaxis]
        field_grid["mixing"] = new_params["mixing"][:, np.newaxis, np.newaxis]

        for m in range(n_mass):
            for el in range(n_el):
                x = new_params["energy"]
                y = new_params["mixing"]
                inter_grid[el, :, :, m] = interpolators[m, el](x, y)
                field_grid["remnant"][:, :, m] = rem_interpolator[m](x, y)

        self.data = inter_grid.reshape(n_el, -1).transpose()
        self.field_data = field_grid.flatten()
        self.nstar = self.data.transpose().shape[1]
        self.nvalues = np.array(
            [n_mass, new_size["energy"], new_size["mixing"], -1]
        ).astype("int")
