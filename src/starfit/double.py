"""Fit double stars - compatibility module"""

from .multi import Multi


class Double(Multi):
    sol_size = 2

    def __init__(self, *args, **kwargs):
        sol_size = kwargs.setdefault("sol_size", self.sol_size)
        assert sol_size == 2, f"require sol_size == 2, provided: {sol_size=}"
        super().__init__(*args, **kwargs)
