"""
define configurations

TODO: read from config file?
"""

from collections import OrderedDict
from collections.abc import Iterable, Mapping

import numpy as np
import yaml

from .prefix import Prefix, _BinPrefixes, _binPrefixes
from .prefix import _Prefixes2 as _Prefixes
from .prefix import _prefixes2 as _prefixes
from .unit import Unit, Units


class Config(object):
    """ """

    def __init__(self, config=None, /, **kwargs):
        if config is not None:
            if isinstance(config, str):
                try:
                    config = globals()[config]
                except KeyError:
                    pass
            if isinstance(config, str) and config.endswith(".yaml"):
                try:
                    with open(config, "rt") as f:
                        data = f.read()
                    config = yaml.safe_load(data)
                except:
                    raise AttributeError(f'Could not load "{config}".')
                assert isinstance(config, dict)
            if isinstance(config, str) and config.endswith(".py"):
                config = eval(config)
            self._config = OrderedDict(config)
        else:
            self._config = OrderedDict()
        self._config.update(kwargs)

        # some defaults
        self._config.setdefault("name", "value")
        self._config.setdefault("latex", self.name)
        self._config.setdefault("unicode", self.name)

        self._config.setdefault("convert", None)
        self._config.setdefault("integer", False)

        # maybe do some processing in units is no type Units
        units = self._config.get("units", None)
        if units is None:
            units = Units()
        elif isinstance(units, (list, tuple, set, frozenset)):
            units = Units(*units)
        elif isinstance(units, Unit):
            units = Units(units)
        elif isinstance(units, str):
            prefix = self._config.pop("prefix", None)
            units = Units(Unit(units, 1, prefix=prefix))
        self._config["units"] = units

        hi_prefixes = self._config.get("hi_prefixes", _Prefixes)
        if hi_prefixes in (
            None,
            False,
        ):
            hi_prefixes = tuple()
        elif hi_prefixes in (
            Ellipsis,
            True,
        ):
            hi_prefixes = _Prefixes
        hi_prefixes = (
            Prefix(**p)
            if isinstance(p, Mapping)
            else Prefix(*p)
            if isinstance(p, Iterable)
            else p
            for p in hi_prefixes
        )
        hi_prefixes = sorted(hi_prefixes, key=lambda x: x.value, reverse=False)
        self._config["hi_prefixes"] = tuple(hi_prefixes)

        lo_prefixes = self._config.get("lo_prefixes", _prefixes)
        if lo_prefixes in (
            None,
            False,
        ):
            lo_prefixes = tuple()
        elif lo_prefixes in (
            Ellipsis,
            True,
        ):
            lo_prefixes = _prefixes
        lo_prefixes = (
            Prefix(**p)
            if isinstance(p, Mapping)
            else Prefix(*p)
            if isinstance(p, Iterable)
            else p
            for p in lo_prefixes
        )
        lo_prefixes = sorted(lo_prefixes, key=lambda x: x.value, reverse=True)
        self._config["lo_prefixes"] = tuple(lo_prefixes)

        self._config.setdefault("zero_unit", self.units.base)

        prefix_base = self._config.get("prefix_base", None)
        self._config.setdefault("hi_prefix_base", prefix_base)
        self._config.setdefault("lo_prefix_base", prefix_base)

        power_base = self._config.get("power_base", None)
        self._config.setdefault("hi_power_base", power_base)
        self._config.setdefault("lo_power_base", power_base)

        # self._config.setdefault

    def __getitem__(self, index):
        return self._config[index]

    def __len__(self):
        return len(self._config)

    def __getattr__(self, key):
        try:
            return self._config[key]
        except:
            raise AttributeError()


import datetime

SEC = 31556952


def time_convert(val):
    if isinstance(val, datetime.timedelta):
        val = val.total_seconds()
    return val


time_config = {
    "name": "time",
    "convert": time_convert,
    "units": Units(
        Unit("s", 1, div_lim=100),
        Unit("min", 60, div_lim=100),
        Unit("h", 3600),
        Unit("d", 86400),
        Unit("wk", 86400 * 7, prefix="+", output=False),
        Unit("yr", 31556952),
    ),
    # 'hi_prefix_base' : 'yr',  # default
    "hi_power_base": "yr",
}

byte_config = {
    "name": "capacity",
    "integer": True,
    "units": Units(
        Unit("B", 1, prefix="+"),
    ),
    "hi_prefixes": _BinPrefixes,
    "lo_prefixes": _binPrefixes,
}

bytefrac_config = {
    "name": "capacity",
    "units": Units(
        Unit("B", 1),
    ),
    "hi_prefixes": _BinPrefixes,
    "lo_prefixes": _binPrefixes,
}

length_config = {
    "name": "length",
    "units": Units(
        # Unit('mm', 0.1, prefix = None),
        Unit("cm", 1, prefix=None),
        Unit(
            "m",
            100,
            prefix="-c",
            out_prefix=(
                "-",
                "k",
            ),
            hi_power_limit=1e6,
        ),
        # Unit('km', 100000, prefix = None, hi_power_limit = 6.957e7),
        # Unit('Mm', 100000000, prefix = None),
        Unit("Rsun", 6.957e10, latex=r"\mathrm{R}_\odot", out_prefix=None),
        Unit("AU", 14959787070000, out_prefix=None),
        # Unit('kAU', 14959787070000000, prefix = None),
        Unit("ly", 29979245800 * 31556952, output=False),
        Unit("pc", 14959787070000 / np.tan(np.pi / 648000), out_prefix="+"),
    ),
    "lo_prefix_base": "m",
    # 'hi_prefix_base' : 'm',
    "hi_prefix_base": "pc",
    "power_base": "cm",
    "zero_unit": "cm",
}

mass_config = {
    "name": "mass",
    "units": Units(
        Unit("g", 1, hi_power_limit=1e3),
        # Unit('kg', 1000, prefix = None),
        # Unit('t', 1e6),
        Unit("Msun", 1.98847e33, latex=r"\mathrm{M}_\odot"),
    ),
    "hi_prefix_base": "Msun",
}

erg_config = {
    "name": "energy",
    "units": Units(
        Unit("erg", 1),
    ),
    "lo_prefix_base": None,
    "hi_prefix_base": None,
}

# TODO - should use cgs internally?
eV_config = {
    "name": "energy",
    "units": Units(
        Unit("eV", 1),
    ),
}

energy_config = {
    "name": "energy",
    "units": Units(
        Unit(
            "eV",
            1.602176634e-12,
            out_prefix=("-", "k", "M", "G"),
            hi_power_limit=-1e-3,
        ),
        Unit("erg", 1, out_prefix=None),
        Unit("B", 1e51, out_prefix=None),
    ),
    "lo_prefix_base": "eV",
    "hi_prefix_base": None,
    "power_base": "erg",
}

frequency_config = {
    "name": "frequency",
    "units": Units(
        Unit("Hz", 1),
    ),
}

temperature_config = {
    "name": "temperature",
    "units": Units(
        Unit("K", 1),
    ),
}

density_config = {
    "name": "dnsity",
    "units": Units(
        dict(
            name="g/cm**3",
            value=1,
            prefix=False,
            latex=r"g\,cm{^{-3}}",
            unicode="g cm\u00B3",
            alt="gcc",
        ),
    ),
}

column_config = {
    "name": "column density",
    "units": Unit(
        "g/cm**2",
        prefix=False,
        latex=r"g\,cm{^{-2}}",
        unicode="g cm\u207b\u00B2",
    ),
}

volume_config = dict(
    name="volume",
    units=(
        dict(
            name="cm**3",
            value=1,
            prefix=False,
            latex="cm{^3}",
            unicode="cm\u00B3",
            alt="cc",
        ),
        dict(
            name="l",
            value=1_000,
            prefix=False,
            latex=r"{\ell}",
            unicode="\u2113",
        ),
        dict(
            name="m**3",
            value=1_000_000,
            prefix="-c",
            latex="m{^3}",
            unicode="m\u00B3",
            prefix_dimension=3,
        ),
        # dict(
        #     'km**3',
        #     value = 1_000_000_000_000_000,
        #     prefix = False,
        #     latex = 'km{^3}',
        #     unicode = "km\u00B3",
        #     ),
        dict(
            name="Vsun",
            value=6.957e10**3 * np.pi * 4 / 3,
            prefix=False,
            latex=r"V{_\odot}",
        ),
        dict(
            name="pc**3",
            value=(14959787070000 / np.tan(np.pi / 648000)) ** 3,
            prefix="+",
            latex="pc{^3}",
            unicode="pc\u00B3",
            prefix_dimension=3,
        ),
    ),
)
