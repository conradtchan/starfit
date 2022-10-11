# OLD:
# from .time import time2human
# from .byte import byte2human

from .version import version2human

__all__ = [
    "byte2humen",
    "time2human",
    "version2human",
]

# Lazy creation using module __getattr__
from . import config
from .template import Human

_map = {
    "time2human": "time_config",
    "byte2human": "byte_config",
    "length2human": "length_config",
    "volume2human": "volume_config",
    "temperature2human": "temperature_config",
    "density2human": "density_config",
    "mass2human": "mass_config",
}

_dynamic = True

if _dynamic:

    def __getattr__(name):
        if name in _map:
            if isinstance(_map[name], str):
                h = Human(getattr(config, _map[name]))
                _map[name] = h
            else:
                h = _map[name]
            return h
        raise AttributeError(f'[{__name__}] "{name}" not found.')

else:
    from .config import byte_config, length_config, time_config

    time2human = Human(time_config)
    byte2human = Human(byte_config)
    length2human = Human(length_config)
