"""
Module for reaction based on Ion class

TDOD - put on PyPi
"""
import copy
import io
import os
import re
import sys
import time
from collections.abc import Iterable, MappingView
from pathlib import Path
from textwrap import wrap

import numpy as np
from numpy.linalg import matrix_power
from scipy.linalg import solve_banded

from .human import version2human
from .isotope import VOID, VOID_IDX, Elements, elements, get_ufunc, get_ufunc_idx
from .isotope import ion as I
from .isotope import (
    ioncacheidx,
    ision,
    isionclass,
    ufunc_A,
    ufunc_E,
    ufunc_element,
    ufunc_element_idx,
    ufunc_idx,
    ufunc_ion_from_idx,
    ufunc_is_element,
    ufunc_is_ion,
    ufunc_is_isobar,
    ufunc_is_isomer,
    ufunc_is_isotone,
    ufunc_is_isotope,
    ufunc_isobar,
    ufunc_isobar_idx,
    ufunc_isomer,
    ufunc_isomer_idx,
    ufunc_isotone,
    ufunc_isotone_idx,
    ufunc_isotope,
    ufunc_isotope_idx,
    ufunc_N,
    ufunc_Z,
)
from .logged import Logged
from .physconst import MSUN
from .utils import CachedAttribute, cachedmethod, index1d, magic, project, stuple


class DuplicateIons(Exception):
    """
    Exception raise when ions list or sets illegally contain duplicate ions
    """

    def __init__(self, ions=None):
        if ions is None:
            super().__init__("Duplicate ions.")
        else:
            super().__init__(f"Duplicate ions: {ions}")


class IonList(object):
    """
    Unsorted list of ions - fixed positions.

    TODO - Maybe not allow KepIon and Ion mixed - all things with extra flags.

    TODO - add some of AbuSet functionality, e.g., return all
    Elements, maybe mapping matrices (move functionality here?)
    """

    def __init__(self, ions=None, duplicates=None):
        """
        duplicates:
          True
          False - drop
          None - raise exception
        """
        self._ions = np.array((), dtype=object)
        self.duplicates = duplicates
        self.add(ions)

    def add(self, ions=None):
        """
        Add ion(s) to list.
        """
        if ions is None:
            return
        if isinstance(ions, MappingView):
            ions = tuple(ions)
        if isinstance(ions, self.__class__):
            self._ions = np.append(self._ions, ions)
            return
        ions = np.atleast_1d(ions)
        ions = np.array([i if ision(i) else I(i) for i in ions])
        if self.duplicates is True:
            self._ions = np.append(self._ions, ions)
        elif self.duplicates is False:
            _, ii = np.unique(ions, return_index=True)
            ions = ions[sorted(ii)]
            ii = np.in1d(ions, self._ions)
            self._ions = np.append(self._ions, ions[~ii])
        elif self.duplicates is None:
            ions = np.append(self._ions, ions)
            ix, nx = np.unique(ions, return_counts=True)
            ix = ix[nx > 1]
            if len(ix) > 0:
                raise DuplicateIons(ix)
            self._ions = ions
        else:
            raise Exception(f'Invalid mode for "duplicates": {self.duplicates}')

    def index(self, ion):
        if not ision(ion):
            ion = I(ion)
        i = np.argwhere(self._ions == ion)
        if len(i) == 0:
            return -1
        else:
            return i.flatten()[0]

    def __eq__(self, other):
        # assert isinstance(other, IonList)
        if len(self) != len(other):
            return False
        for i, j in zip(self, other):
            if not i == j:
                return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def ions(self):
        return copy.deepcopy(self._ions)

    def copy(self):
        new = self.__class__()
        new._ions = self._ions.copy()
        return new

    # copy from IonList:
    def __str__(self):
        return f"[{','.join(str(i) for i in self._ions)}]"

    def __repr__(self):
        return f"{self.__class__.__name__}({str(self)})"

    def __getitem__(self, index):
        return self._ions[index]

    def __len__(self):
        """
        Return number of isotopes.
        """
        return len(self._ions)

    def __iter__(self):
        for ion in self._ions:
            yield ion

    def __add__(self, other):
        new = self.copy()
        new.add(other)
        return new

    def __call__(self):
        return self._ions

    @property
    def A(self):
        return ufunc_A(self._ions)

    @property
    def Z(self):
        return ufunc_Z(self._ions)

    @property
    def N(self):
        return ufunc_N(self._ions)

    @property
    def E(self):
        return ufunc_E(self._ions)

    @property
    def idx(self):
        return ufunc_idx(self._ions)


class IonSet(object):
    """
    Provide sorted list of isotopes.

    TODO: efficient additions of arrays of ions.

    I suppose in principle we do want :set: properties but it also is
    supposed to be sorted by index at all times if _sort is set.

    Maybe there need be 2 kinds to allow interface with ionmap:
    ionmap also allows to have elements, etc., mixed.

    TODO - Maybe not KepIon and Ion mixed - things with extra flags.

    TODO - Do not allow duplicates

    NOTE - I think there should be only one, IonSet or IonList
    """

    def __init__(self, *args, **kwargs):
        self._ions = np.array((), dtype=object)
        self._sort = kwargs.get("sort", True)
        self._type = kwargs.get("type", None)
        assert self._type is None or not isionclass(
            self._type
        ), "type needs to be subclass of Ion"
        self.add(args)

    def copy(self):
        """
        return a copy
        """
        return copy.copy(self)

    def add(self, ions=None):
        """
        Add ion(s) to list.

        Ions need to be same "type"
        """
        if isinstance(ions, str):
            ions = tuple(I(i) for i in re.split(r"['\",;\s]+", ions))
        if len(ions) == 1:
            try:
                ions = ions[0]
            except:
                pass
        if ions is not None:
            if np.isscalar(ions):
                self._add_one(ions)
            elif len(ions) == 1:
                while isinstance(ions, Iterable) and not isinstance(ions, str):
                    ions = ions[0]
                self._add_one(ions)
            else:
                for ion in ions:
                    self._add_one(ion, update=False)
                self._update()

    def _add_one(self, ix, update=True):
        """
        add one ion and check compatibility
        """
        if self._type is None:
            self._type = type(ix)
            if not isionclass(self._type):
                self._type = type(I(ix))
        if ision(ix) and not isinstance(ix, self._type):
            raise TypeError(
                f"All ions need compatible type: {str(self._type):s}, {str(type(ix)):s}"
            )
        if not isinstance(ix, self._type):
            ix = self._type(ix)
        self._ions = np.append(self._ions, ix)
        if update:
            self._update()

    def _update(self):
        if self._sort:
            self._ions.sort()
        self._ions = np.unique(self._ions)

    def __str__(self):
        return f"[{','.join(str(i) for i in self._ions)}]"

    def __repr__(self):
        return f"{self.__class__.__name__}({str(self)})"

    def __getitem__(self, index):
        return self._ions[index]

    def __len__(self):
        """
        Return number of isotopes.
        """
        return len(self._ions)

    def __iter__(self):
        i = 0
        while i < self.__len__():
            yield self._ions[i]
            i += 1

    def __add__(self, other):
        new = self.copy()
        new.add(other)
        return new


class Abu(object):
    """
    Ion plus single abundance.  Provide some general functionality.

    TODO - add most of AbuSet functionality

    TODO  Contain mass information?
    """

    def __init__(self, ion, abu=0):
        self.ion = ion
        self.abu = np.float64(abu)

    def __str__(self):
        return f"{self.ion!s}:{self.abu:>12.5f}"

    def X(self):
        return self.abu

    def Y(self):
        return self.abu / self.ion.A


class AbuDist(object):
    """
    Ion and distribution information, e.g., structure.

    TDOD - implement
    TDOD - add Abu(Set) and AbuData functionalities

    TODO - internam structures of AbuData should be (ions, dimensions) as we want to logically write
    AbuData['he4'] returning an ndarray
    """

    # should be initialized with Ion and abundace vector,
    # mass and radius would be good.
    pass


class AbuData(object):
    """
    n-D abu data plus extra info.
    ions are in shape[-1]

    Probably base class (if not identical) to starfit data.

    Later, this should be combined with AbuSet toward a multi-D version of that.
    """

    def __init__(
        self,
        data=None,
        ions=None,
        *,
        molfrac=None,
        # fields = None,
        # formats = None,
        # values = None,
        # info = None,
        copy=True,
        mass_table=None,
        silent=False,
    ):

        self.silent = silent
        self.ions = IonList(ions)
        assert isinstance(data, np.ndarray), "require numpy data array"
        if copy:
            self.data = data.copy()
        else:
            self.data = data
        nions = len(self.ions)
        if self.data.shape[-1] != nions:
            try:
                ii = self.data.shape.index(nions)
                self.data = np.swapaxes(self.data, ii, -1)
            except ValueError:
                pass

        assert len(self.ions) == self.data.shape[-1]
        if molfrac is None:
            x = np.sum(data.reshape(-1, (len(self.ions)))[0, :])
            molfrac = np.abs(x - 1) > 1e-4
            if not self.silent:
                print(f" [AbuData] Setting molfrac = {molfrac}")
        self.molfrac = molfrac

        # self.fields = fields
        # self.formats = formats
        # self.values = values
        self._mass_table = mass_table
        self.i0 = 0
        self.i1 = self.data.shape[0]

    def as_molfrac(self):
        if self.molfrac:
            return self
        out = self.copy()
        out.molfrac = True
        out.data /= self.ions.A
        return out

    def as_massfrac(self):
        if not self.molfrac:
            return self
        out = self.copy()
        out.molfrac = False
        out.data *= self.ions.A
        return out

    def __len__(self):
        return np.product(self.data.shape[:-1])

    def __iadd__(self, other):
        new = self.join(self, other)
        self.__dict__.update(new.__dict__)
        return self

    def __radd__(self, other):
        return self.join(other, self)

    def __add__(self, other):
        return self.join(self, other)

    @classmethod
    def join(cls, objects, axis=0, along=True, silent=False):
        if len(objects) == 0:
            return None
        for o in objects:
            assert isinstance(o, cls)
        assert np.all([objects[0].molfrac == o.molfrac for o in objects])

        # merge (including sort) if needed
        if not np.all([objects[0].ions == o.ions for o in objects]):
            if not silent:
                print(f" [{cls.__name__}] merging ion lists.")
            ions = set()
            ionsa = list()
            for o in objects:
                idx = ufunc_idx(o.ions)
                ions |= set(idx)
                ionsa.append(idx)
            ions = np.array(sorted(ions))
            nions = len(ions)
            ionlist = IonList(ioncacheidx(ions))
            for i, o in enumerate(objects):
                o = o.copy()
                data = np.zeros(o.data.shape[:-1] + (nions,))
                ii = index1d(ionsa[i], ions)
                data[..., ii] = o.data
                o.data = data
                o.ions = ionlist
                objects[i] = o

        assert np.all([objects[0].data.ndim == o.data.ndim for o in objects])
        if along:
            assert axis < objects[0].data.ndim - 1
        else:
            assert axis < objects[0].data.ndim
        dims = list(range(objects[0].data.ndim))
        if along:
            del dims[axis]
        assert np.all(
            [
                np.all(
                    np.array(objects[0].data.shape)[dims]
                    == np.array(o.data.shape)[dims]
                )
                for o in objects
            ]
        )

        mass_table = None
        for o in objects:
            if mass_table is None:
                mass_table = o._mass_table
            else:
                if o._mass_table is not None:
                    assert mass_table == o._mass_table

        if axis == 0 and along:
            assert np.all([o.i0 == 0 for o in objects[1:]])
            assert np.all([o.i1 == o.data.shape[0] for o in objects[:-1]])
        else:
            assert np.all([objects[0].i0 == o.i0 for o in objects])
            assert np.all([objects[0].i1 == o.i1 for o in objects])

        new = objects[0].copy()
        if along:
            new.data = np.concatenate(tuple(o.data for o in objects), axis=axis)
        else:
            new.data = np.array([o.data for o in objects])
            if axis != 0:
                new.data = np.moveaxis(new.data, 0, axis)
        if axis == 0 and along:
            new.i1 = objects[-1].i1
        return new

    @classmethod
    def from_abusets(cls, abu, molfrac=None):
        if isinstance(abu, (list, tuple)):
            s = (len(abu),)
            x = np.ndarray(s, dtype=object)
            x[:] = abu
            abu = x
        elif isinstance(abu, np.ndarray):
            s = abu.shape
            abu = abu.flatten()
        else:
            raise AttributeError(f"Unsupported type {type(abu)=}")
        ions = set()
        ionsa = list()
        for a in abu:
            idx = ufunc_idx(a.iso)
            ions |= set(idx)
            ionsa.append(idx)
        ions = np.array(sorted(ions))
        nions = len(ions)
        abub = np.zeros((len(abu), nions))
        for i, a in enumerate(abu):
            ii = index1d(ionsa[i], ions)
            abub[i, ii] = a.abu
        ions = ioncacheidx(ions)
        abub = np.reshape(abub, s + (nions,))
        return cls(abub, ions, molfrac=molfrac)

    def copy(self, molfrac=None):
        new = copy.deepcopy(self)
        if molfrac is None or molfrac == self.molfrac:
            return new
        if new.molfrac:
            new.data *= new.ions.A
        else:
            new.data /= new.ions.A
        new.molfrac = not new.molfrac
        return new

    def abu(self, molfrac=None):
        """
        abundace
        """
        if molfrac is None:
            molfrac = self.molfrac

        if molfrac == self.molfrac:
            yfac = np.ones(self.data.shape[1])
        else:
            yfac = ufunc_A(self.ions)
        if molfrac and not self.molfrac:
            yfac = 1 / yfac

        value = np.ndarray(self.data.shape, dtype=np.float64)
        value[..., :] = self.data[..., :] * yfac
        return value

    @property
    def X(self):
        return self.abu(molfrac=False)

    @property
    def Y(self):
        return self.abu(molfrac=True)

    @property
    def mass_nc(self):
        return np.sum(self.X, axis=-1) - 1

    def updated(
        self,
        data=None,
        ions=None,
        molfrac=None,
    ):
        """
        return updated copy of self
        """
        new = copy.deepcopy(self)
        new.data = data.copy()
        new.ions = IonList(ions.copy())
        if molfrac is not None:
            new.molfrac = molfrac
        return new

    def decayed(self, decay=None, **decpar):
        """
        Return copy of self with yields replaced by decayed yields.
        """
        from ionmap import Decay

        return Decay.Map(self, **decpar)

    def cleaned(self, threshold=0.0, debug=False):
        """
        return new data set with exess ions removed
        """
        if threshold is None:
            threshold = -1.0
        maxima = np.max(
            self.data,
            axis=tuple(np.arange(self.data.ndim - 1, dtype=np.int64).tolist()),
        )
        (ii,) = np.where(maxima > threshold)
        save_data = self.data
        save_ions = self.ions
        self.data = self.data[:, ii]
        self.ions = self.ions.__class__(self.ions[ii])
        if debug:
            delions = set(save_ions) - set(self.ions)
            n = len(delions)
            if n > 0:
                s = ", ".join(str(i) for i in delions)
                print(f"[{self.__class__.__name__}] removing {n} ions: {s}.")
        new = copy.deepcopy(self)
        self.data = save_data
        self.ions = save_ions
        return new

    def sliced(self, index):
        """
        return data-sliced copy of self
        """
        new = copy.deepcopy(self)
        new.data = new.data[index, :].copy()
        return new

    @CachedAttribute
    def idx(self):
        """
        return ion indices
        """
        return ufunc_idx(self.ions)

    @property
    def A(self):
        return self.ions.A

    @property
    def Z(self):
        return self.ions.Z

    @property
    def N(self):
        return self.ions.N

    @property
    def E(self):
        return self.ions.E

    def __getitem__(self, key):
        i = self.ions.index(key)
        return self.data[..., i]

    def zone_abu(self, i):
        return AbuSet(
            self.ions,
            self.data[i],
            molfrac=self.molfrac,
        )

    @CachedAttribute
    def abar(self):
        """
        Compute Abar as a function of coordinate
        """
        if not self.molfrac:
            yf = 1 / ufunc_A(self.ions)
        else:
            yf = np.ones(self.data.shape[1])
        xY = np.tensordot(self.data, yf, axes=(-1, 0))
        return 1 / xY

    muI = abar

    @CachedAttribute
    def zbar(self):
        """
        mean charge number
        """
        Z = ufunc_Z(self.ions)
        if not self.molfrac:
            yf = 1 / ufunc_A(self.ions)
            Z = Z * yf
        else:
            yf = np.ones(self.data.shape[1])
        xZ = np.tensordot(self.data, Z, axes=(-1, 0))
        xY = np.tensordot(self.data, yf, axes=(-1, 0))
        return xZ / xY

    @CachedAttribute
    def Ye(self):
        """
        Compute Ye as a function of coordinate
        """
        Z = ufunc_Z(self.ions)
        if not self.molfrac:
            Z = Z / ufunc_A(self.ions)
        Ye = np.tensordot(self.data, Z, axes=(-1, 0))
        return Ye

    @CachedAttribute
    def eta(self):
        """
        neutron excess eta = 1-2*Ye

        this is approximate

        really should be

        sum_i(Y_i * (N_i - Z_i)) / sum_i(Y_i * (N_i + Z_i))

        or

        ~ sum_i(Y_i * (N_i - Z_i) / (N_i + Z_i))

        to account for true atiomic masses
        """
        return 1 - 2 * self.Ye

    @CachedAttribute
    def mu(self):
        """
        mean molecular weight
        """
        Z1 = 1 + ufunc_Z(self.ions)
        if not self.molfrac:
            Z1 = Z1 / ufunc_A(self.ions)
        xmu = np.tensordot(self.data, Z1, axes=(-1, 0))
        return 1 / xmu

    @CachedAttribute
    def mue(self):
        """
        mean molecular weight per electron
        """
        return 1 / self.Ye

    @property
    def mass_table(self):
        import bdat
        import mass_table

        if not hasattr(self._mass_table, "mass_excess"):
            if self._mass_table is None:
                self._mass_table = mass_table.MassTable()
            elif isinstance(self._mass_table, str) and self._mass_table.endswith(
                "bdat"
            ):
                self._mass_table = bdat.BDat(self._mass_table)
            elif isinstance(self._mass_table, str):
                self._mass_table = mass_table.MassTable(self._mass_table)
            else:
                print(" [AbuData] invalid mass table")
                raise AttributeError("invalid mass table")
        return self._mass_table

    @CachedAttribute
    def ME(self):
        """
        mass excess in MeV
        """
        m = self.mass_table
        v = m.mass_excess(self.ions)
        if not self.molfrac:
            v /= ufunc_A(self.ions)
        return np.tensordot(self.data, v, axes=(-1, 0))

    @CachedAttribute
    def BE(self):
        """
        Binding energy in MeV
        """
        m = self.mass_table
        v = m.binding_energy()
        if not self.molfrac:
            v /= ufunc_A(self.ions)
        return np.tensordot(self.data, v, axes=(-1, 0))

    def write_burngen(
        self,
        filename,
        comment=list(),
        i0=None,
        i1=None,
        zone0=1,
    ):

        # version 1.1 adds 'abu' card
        version = 10100

        card_cmt = "c"
        card_abu = "abu"
        card_net = "net"
        net = 1

        iions, imap, oions, omap = get_burn_map(self.ions)
        abu = np.zeros(self.data.shape[:-1] + (len(oions),), dtype=np.float)
        abu[..., omap] = self.abu(molfrac=True)[..., imap]

        if i0 is None:
            i0 = 0
        if i1 is None:
            i1 = self.data.shape[0] - 1
        filename = Path(filename).expanduser()
        with open(filename, "wt") as f:
            f.write(f"{card_cmt} COMPUTER-GENERATED BURN GENERATOR FILE\n")
            f.write(f"{card_cmt} VERSION {version2human(version)}\n")
            f.write(f"{card_cmt} {time.asctime(time.gmtime())} UTC\n")
            f.write(f"{card_cmt}\n")

            for c in comment:
                f.write(f"{card_cmt} {c}\n")
            f.write(f"{card_cmt}\n")

            indent = f"{card_net} {net:d} "
            for line in wrap(
                " ".join([iso.Kepler for iso in oions]),
                initial_indent=indent,
                subsequent_indent=indent,
            ):
                f.write(f"{line}\n")
            f.write(f"{card_cmt}\n")
            f.write(f"{card_cmt} raw ppnb data (as needed in ppnb)\n")
            s = "".join([f"{i.name():>25s}" for i in oions])
            f.write(f"{card_cmt} {s[2:]}\n")
            f.write(f"{card_abu} {net:d} {zone0} {zone0 + i1 - i0}\n")
            for i in range(i0, i1 + 1):
                x = abu[i, :]
                x = np.maximum(x, 0.0)
                x[x < 1.00001e-99] = 0.0
                f.write("".join([f"{float(a):25.17e}" for a in x]) + "\n")

    # really slow
    # def __iter__(self, **kw):
    #     for i,ion in enumerate(self.ions):
    #         yield (ion, self.ion_abu(ion))

    # fastest, but not as flexible or consistent with returned arrays
    # def __iter__(self):
    #     for i,ion in enumerate(self.ions):
    #         yield (ion, self.data[self.i0:self.i1,i])

    # not quite as fast but consistent and flexible
    def __iter__(self, **kw):
        for i, ion in enumerate(self.ions):
            yield (ion, self.ion_abu(index=i))

    def ion_abu(
        self,
        ion=None,
        molfrac=False,
        missing=np.nan,
        index=None,
        ionmissing=None,
    ):
        """
        Return isotope abundace
        """
        if index is not None:
            ion = self.ions[index]
        elif not ision(ion):
            ion = I(ion)
            if ion == VOID:
                value = np.ndarray(self.data.shape[0], dtype=np.float64)
                value.fill(missing)
                return value
            # return self.ion_abu(ion, molfrac = molfrac, missing = missing)
        yfac = max(1, ion.A)
        if molfrac == self.molfrac:
            yfac = 1
        if molfrac and not self.molfrac:
            yfac = 1 / yfac
        value = np.ndarray(self.data.shape[0], dtype=np.float64)
        if index is None:
            # this is really bad and should not be used - rewrite to use ion.idx
            index = self.ions.index(ion)
        if index == -1:
            if ionmissing is not None:
                value[self.i0 : self.i1] = ionmissing
            else:
                raise AttributeError(f"ion {ion} not in data set")
        else:
            value[self.i0 : self.i1] = self.data[self.i0 : self.i1, index] * yfac

        value[: self.i0] = missing
        value[self.i1 :] = missing
        return value


def get_burn_map(ions):
    """
    make minimum net and set alpha particle as 3rd species

    returns:
        iions:
            the subset of ions that can be used.
        imap:
            llist of indices into original array
        oions:
            list of ions to be mapped to, including minnet
        omap:
            indices of iions in oions list
    """

    minnet = [
        [[1, 1]],
        [[1, 3]],
        [[3, 4]],
        [[6, 7]],
        [[7, 7], [9, 9]],
        [[8, 8], [10, 11]],
        [[11, 13]],
        [[13, 15]],
        [[14, 18]],
        [[17, 19]],
        [[19, 22]],
        [[21, 23]],
        [[23, 26]],
        [[25, 27]],
        [[27, 30]],
        [[30, 31]],
        [[31, 36]],
    ]

    zmax = max(np.max(ions.Z), len(minnet))
    netw = np.ndarray((zmax + 1, 2, 2), dtype=np.int64)
    netw_count = np.zeros(zmax + 1, dtype=np.int64)

    for iz, aa in enumerate(minnet):
        for iar, ar in enumerate(aa):
            netw[iz, iar, :] = ar
        netw_count[iz] = len(aa)

    # assume ions are sorted
    imap = list()
    iions = list()
    for i, iso in enumerate(ions):
        Z = iso.Z
        if Z > zmax:
            break
        A = iso.A
        n = netw_count[Z]
        if n == 0:
            netw_count[Z] = 1
            netw[Z, 0, :] = A
        else:
            netw[Z, 0, 0] = min(A, netw[Z, 0, 0])
            netw[Z, n - 1, 1] = max(A, netw[Z, n - 1, 1])
        add = True
        if n == 2:
            if (A > netw[Z, 0, 1]) and (A < netw[Z, 1, 0]):
                print(f"{iso.name()} inside required gap.")
                add = False
        if add:
            imap.append(i)
            iions.append(iso)

    oions = list()
    omap = list()
    noions = 0
    for iz, aa in enumerate(netw):
        for iar, ar in enumerate(aa[: netw_count[iz]]):
            for a in range(ar[0], ar[1] + 1):
                iso = I(A=a, Z=iz)
                oions.append(iso)
                if iso in ions:
                    omap.append(noions)
                noions += 1

    # fix alpha particle as 3rd particle
    alpha = I.He4
    for ia, ion in enumerate(oions):
        if ion == alpha:
            break
    oions[3 : ia + 1] = oions[2:ia]
    oions[2] = alpha

    for im, iom in enumerate(omap):
        if iom == ia:
            break
    for i in range(2, im):
        omap[i] += 1
    omap[im] = 2

    iions = IonList(iions)
    oions = IonList(oions)

    return iions, imap, oions, omap


class AbuDump(AbuData):
    """
    KEPLER BURN abundance set.  To hold KEPLER ppnb data (plus wind)

    Should be combined with AbuSet (subclass of)
    """

    def __init__(
        self,
        data=None,
        ions=None,
        windb=None,
        molfrac=None,
        bottom=None,
        top=None,
        xm=None,
        zm=None,
        mass_table="bdat",
    ):
        """
        xm is the mass of zone in g
        zm is the mass coordinate of outer zone boundary

        if there is n star data, it should be in zones 1..n zone 0 is
        the inner boundady condition (could become having non-zero data
        in accretion problems); zone n+1 holds the wind data

        abundace data should be mol fractions or mass fractions

        data on import has format [isotopes, zones+2] as in Kepler
        ppnb, but is stored internally as [zones+2, isotopes]

        molfrac determines whether the input is in mol fractions;
        wise to set if known; otherwise the code tries automatic determination,
        which may fail.
        """
        if isinstance(data, AbuData):
            self.data = data.data.copy()
            self.ions = data.ions.copy()
            self.molfrac = data.molfrac
            self._mass_table = data._mass_table
            assert ions is None
            assert molfrac is None
            if bottom is None:
                bottom = 1
            self.i0 = bottom
            if top is None:
                top = data.data.shape[1]
            self.i1 = top
            # check for wind - if abu in wind is 0, reduce i1 by 1
            if np.sum(self.data[-1, :]) == 0:
                self.has_wind = False
                self.i1 -= 1
            else:
                self.has_wind = True
        else:
            self.ions = IonList(ions)
            self.data = data.copy()
            if data.shape[0] == len(self.ions):
                self.data = self.data.transpose()
            self.ions = IonList(ions)
            if molfrac is None:
                molfrac = (1 - np.sum(self.data[1, :])) > 0.1
            self.molfrac = molfrac
            if bottom is None:
                bottom = 1
            self.i0 = bottom
            if top is None:
                top = self.data.shape[0] - 1
            self.i1 = top
            if windb is not None:
                self.data[-1, :] = windb
                self.i1 += 1
                self.has_wind = True
            else:
                self.has_wind = False
            self._mass_table = mass_table
        self.xm = xm
        self.zm = zm

    @classmethod
    def from_abusets(cls, abu, molfrac=None, **kwargs):
        return cls(AbuData.from_abusets(abu, molfrac=molfrac), molfrac=None, **kwargs)

    @cachedmethod
    def __call__(self, *args, **kw):
        return self.ion_abu(*args, **kw)

    @cachedmethod
    def __getitem__(self, ion):
        return self.ion_abu(ion)

    def __setitem__(self, ion, value):
        """we cannot check molfrac"""
        ii = np.where(self.ions == ion)[0][0]
        if len(ii) == 1:
            self.data[ii[0]] = value
            return
        raise AttributeError(ion)

    def __setattr__(self, attr, value):
        """we cannot check molfrac"""
        if "ions" in self.__dict__ and "data" in self.__dict__:
            ii = np.where(self.ions == attr)[0]
            if len(ii) == 1:
                self.data[ii[0]] = value
                return
        super().__setattr__(attr, value)

    def __getattr__(self, attr):
        if "ions" in self.__dict__:
            if attr in self.ions:
                return self.ion_abu(attr)
        raise AttributeError(attr)

    def abu(self, molfrac=None, missing=np.nan, wind=True, bottom=True):
        """
        abundace
        """
        i1 = self.data.shape[0] - (0 if (wind and self.has_wind) else 1)
        i0 = 0 if bottom else 1

        if molfrac is None:
            molfrac = self.molfrac

        if molfrac == self.molfrac:
            yfac = np.ones(self.data.shape[1])
        else:
            yfac = ufunc_A(self.ions)
        if molfrac and not self.molfrac:
            yfac = 1 / yfac
        value = np.ndarray(self.data.shape, dtype=np.float64)
        value[self.i0 : self.i1, :] = (
            self.data[self.i0 : self.i1, :] * yfac[np.newaxis, :]
        )
        value[: self.i0, :] = missing
        value[self.i1 :, :] = missing
        return value[i0:i1]

    def wind_abu(self):
        """
        return wind abu if wind is present, None otherwise
        """
        if self.has_wind:
            return self.zone_abu(-1)
        return None

    def project(
        self,
        output="massfrac",
        xm=None,
        zones=None,
    ):
        """
        Return projected yields as AbusSet object.

        output = massfrac | molfrac | g | mol | Msun
        """
        if xm is None:
            try:
                xm = self.xm
            except AttributeError:
                raise Exception("Need to have zone mass defined.")
        assert xm.shape[0] == self.data.shape[0]

        if zones is None:
            # exclude "phony" zones but leave in wind
            zones = slice(self.i0, self.i1)
        elif isinstance(zones, int):
            zones = slice(zones, zones + 1)

        norm = 1
        if output == "massfrac":
            molfrac = False
            norm = 1 / np.sum(xm[zones])
        elif output == "molfrac":
            molfrac = True
            norm = 1 / np.sum(xm[zones])
        elif output == "mol":
            molfrac = True
        elif output == "g":
            molfrac = False
        elif output == "Msun":
            molfrac = False
            norm = 1 / MSUN
        else:
            molfrac = None

        if molfrac is None:
            molfrac = self.molfrac

        if molfrac == self.molfrac:
            yfac = np.ones(self.data.shape[1])
        else:
            yfac = ufunc_A(self.ions)
        if molfrac and not self.molfrac:
            yfac = 1 / yfac

        if xm[zones].shape == ():
            value = self.data[zones, :] * xm[zones] * yfac * norm
        else:
            value = (
                np.tensordot(self.data[zones, :], xm[zones], axes=(0, 0)) * yfac * norm
            )
        # add proper norm
        return AbuSet(iso=self.ions, abu=value, normalize=False, unit="?")

    def mix(
        self,
        mbox,
        iterations=4,
        xm=None,
        zones=None,
    ):
        """
        Mix abundance dump and return new mixed object.
        """
        if xm is None:
            try:
                xm = self.xm
            except AttributeError:
                raise Exception("Need to have zone mass defined.")
        assert (
            xm.shape[0] == self.data.shape[0]
        ), "require same number of zones in xm and data set"

        if zones is None:
            # exclude "phony" and surface zones
            if self.has_wind:
                zones = slice(self.i0, self.i1 - 1)
            else:
                zones = slice(self.i0, self.i1)

        # convert to list of zones
        if isinstance(zones, slice):
            zones = np.array(range(*zones.indices(len(xm))))
        if not isinstance(zones, np.ndarray):
            zones = np.array(zones)
        nzones = len(zones)

        mat = np.zeros((nzones, nzones), dtype=np.float64)
        xm = xm[zones]

        # compute contributions
        # mat of [i,j] contains contributions of zone j to zone i
        xmt = 0.0
        j = -1
        i = 0
        for xmi in xm:
            j += 1
            xmt += xmi
            if xmt >= mbox:
                break
        if xmt < mbox:
            # complete mixing
            mat[:, :] = xm[np.newaxis, :] / xmt
        else:
            mat[i, i : j + 1] = xm[i : j + 1] / xmt

            for k in range(1, nzones):
                j0 = j
                xmt -= xm[k - 1]
                xmt0 = xmt
                while xmt <= mbox:
                    j += 1
                    xmt += xm[j]
                    if j == nzones - 1:
                        break
                if k <= j0:
                    mat[k, : j0 + 1] = mat[k - 1, : j0 + 1] * (xmt0 / xmt)
                if j > j0:
                    mat[k, j0 + 1 : j + 1] = xm[j0 + 1 : j + 1] / xmt
                i = k
                if j == nzones - 1:
                    break
            # add flat tail
            mat[i + 1 :, :] = mat[i, :]
        mat = matrix_power(np.transpose(mat), iterations)

        result = copy.deepcopy(self)
        result.data[zones, ...] = np.tensordot(
            mat, result.data[zones, ...], axes=(0, 0)
        )
        return result

    def mix_diff(
        self,
        diff,
        xm=None,
        zones=None,
    ):
        """
        Mix abundance dump and return new mixed object.

        'diff' is the mix width using diffusion-like approximation
        (Gaussian-like mixing with widths mbox).  In a diffusion
        apprximation diff would correstpond to sqrt(D * dt) where D is
        the diffusion coefficient and dt is the time [step].  Hence
        the unit of diff here is 'mass' (g) for conveniece.

        'diff' may also be an array of same dimension as xm.  The
        mixing magnitudes are defined on the zone interfaces similar
        to the mass coordinates zm.

        It appears similar results in the near field are obtained with
        a 'diff' value of about 2x the 'mix' 'mbox' value; in the far
        field, a value of sqrt(2) gives more similar results (on log
        scale).
        """
        if xm is None:
            try:
                xm = self.xm
            except AttributeError:
                raise Exception("Need to have zone mass defined.")
        assert (
            xm.shape[0] == self.data.shape[0]
        ), "require same number of zones in xm and data set"

        if zones is None:
            # exclude "phony" and surface zones
            if self.has_wind:
                zones = slice(self.i0, self.i1 - 1)
            else:
                zones = slice(self.i0, self.i1)

        # convert to list of zones
        if isinstance(zones, slice):
            zones = np.array(range(*zones.indices(len(xm))))
        if not isinstance(zones, np.ndarray):
            zones = np.array(zones)
        nzones = len(zones)

        xm = xm[zones]

        # compute forward matrix
        # mat of [i,j] contains contributions of zone j to zone i
        if np.isscalar(diff):
            diff = np.tile(diff, len(xm))
        elif len(diff) == len(self.xm):
            diff = diff[zones]

        assert len(diff) == len(xm)
        xmim = np.zeros_like(xm)
        xmip = np.zeros_like(xm)
        ip = np.arange(1, nzones)
        im = np.arange(0, nzones - 1)
        ia = np.arange(0, nzones)

        xmi = -2 / xm
        xmip[im] = diff[im] ** 2 / (xm[im] + xm[im + 1])
        xmim[ip] = xmip[im]

        mat = np.ndarray((3, nzones), dtype=np.float64)
        mat[0, ip] = xmi[im] * xmip[im]
        mat[1, ia] = -xmi[ia] * (xmip[ia] + xmim[ia]) + 1
        mat[2, im] = xmi[ip] * xmim[ip]

        result = copy.deepcopy(self)
        result.data[zones, ...] = solve_banded((1, 1), mat, result.data[zones, ...])
        return result

    @CachedAttribute
    def zm_sun(self):
        """
        mass shell interfaces in solar masses
        """
        return self.zm / MSUN


########################################################################


class AbuSet(Logged):
    """
    Prototype Class for single abundace set like solar abundance

    Should there be a subclass or a field for normalized/partial abundances?
    - Normalization is difficult if addition occurs piecewise ...

    How about support for different formats, like [] ... ?

    *** TODO ***
    - derive using common functions with IonSet
    - replace internal variables with "_"-names,
      most prominently: iso --> _ions
                        abu --> _abu

    *** Add a way to hold non-normalized values?
    """

    _init = True

    def __init__(
        self,
        iso=None,
        abu=None,
        comment=None,
        mixture=None,
        dat_file=None,
        bg_file=None,
        silent=False,
        # todo
        normalize=False,  # only useful for fractions
        normalized=False,  # only useful for fractions
        molfrac=False,  # ?
        unit="massfrac",  # set up proper unit system
        # - g, mol, massfrac, molfrac, M_sun
        sorted=False,
        sort=False,
        sentinel=None,
        comp=None,
        attr=None,  # dictionary of extra attributes
        allow_negative=False,
        **kwargs,
    ):
        """
        Initialize abundance set.

        Currently allow call with
         * list of isotopes only (abu set to 0)
         * list of isotopes and list of abundances
         * keyword iso=abu
         * dictionary {'iso':abu, ...}
           negative values collet up missing mass fraction, by relative ratio
           dict(c12=.1, h1=1e-3, fe56=1e-5, he4=-1,o16=-2)
           -->
           C12: 0.1, H1: 0.001, Fe56: 1E-05, He4: 0.299663, O16: 0.599327

        Initialize from data file:
         * dat_file (KEPLER abundance data file, like solabu.dat)
         * bg_file (KEPLER BURN generator data file)

        TODO:
          * sorting implementation
          * implement molfrac
          * proper deal with normalization on/off
          * make A, Z, etc., properties rather than functions
          * can initialise from another composition for derived classes?
          * axis parameter for multi-D data, unite with AbuData
          * use IonList or implement IonList functions
          * isotope caching
        """
        self.comment = stuple(comment)
        self.mixture = mixture
        self.sentinel = sentinel
        self.comp = comp
        self.unit = unit
        self.is_sorted = sorted
        if attr is not None:
            self.__dict__.update(**attr)
        if dat_file is not None:
            self._from_dat(dat_file, silent=silent)
            self._init = False
            self.molfrac = False
            return
        if bg_file is not None:
            self._from_bg(bg_file, silent=silent)
            self.molfrac = False
            self._init = False
            return
        if isinstance(iso, IonList):
            iso = iso.ions()
        if isinstance(iso, AbuSet):
            assert abu is None
            self.iso = iso.iso.copy()
            self.abu = iso.abu.copy()
        elif isinstance(iso, dict):
            assert abu is None
            self.iso, self.abu = self._ion_abu_from_dict(
                iso, allow_negative=allow_negative
            )
        elif iso is None:
            assert abu is None, "Need isotope name"
            self.iso = np.array([], dtype=object)
            self.abu = np.array([], dtype=np.float64)
        else:
            if not (
                isinstance(iso, np.ndarray) and (iso.ndim > 0) and (iso.dtype == object)
            ):
                iso = np.atleast_1d(I(iso))
            self.iso = iso
            if abu is not None:
                self.abu = np.array(np.atleast_1d(abu), dtype=np.float64)
                assert len(self.abu) == len(self.iso), "Need equal number of elements."
            else:
                self.abu = np.zeros_like(self.iso, dtype=np.float64)
        if len(kwargs) > 0:
            self._append(*self._ion_abu_from_dict(kwargs))
        if sort and not self.is_sorted:
            self.sort()
            self.is_sorted = True
        # need to update for general use, e.g., mapped elemental abundances
        # TODO - implement
        # assert molfrac == False
        # self.molfrac = molfrac
        if molfrac:
            self.abu *= self.A()
            molfrac = False
        self.molfrac = molfrac
        self._init = False
        self.silent = silent

    def setattr(self, attr, value):
        self.__dict__[attr] = value

    @staticmethod
    def _ion_abu_from_dict(abudict, allow_negative=False):
        # todo: add sorting?
        # todo: check uniqueness

        # todo - add isotope caching
        ions = np.array([I(i) for i in abudict.keys()], dtype=object)
        abu = np.array(list(abudict.values()), dtype=np.float64)
        if not allow_negative:
            jj = abu < 0.0
            if np.count_nonzero(jj) > 0:
                abu[jj] *= (1 - np.sum(abu[~jj])) / np.sum(abu[jj])
        return ions, abu

    @classmethod
    def from_kepabu_string(cls, line):
        tokens = line.split()
        assert len(tokens) % 2 == 0
        abu = np.array([float(t) for t in tokens[0::2]])
        iso = I(tokens[1::2])
        return cls(iso, abu)

    def cleaned(self, threshold=0.0, debug=False):
        """
        return new data set with exess ions removed
        """
        if threshold is None:
            threshold = -1.0
        (ii,) = np.where(self.abu > threshold)
        abu = self.abu[ii]
        iso = self.iso[ii]
        if debug:
            s = ", ".join(str(i) for i in sorted(set(self.iso) - set(iso)))
            print(f"[{self.__class__.__name__}] removing {s}.")
        new = self.copy()
        new.abu = abu
        new.iso = iso
        return new

    def _apply_limits(
        self,
        low_X_truncate=None,
        high_Z_truncate=None,
        high_A_truncate=None,
        high_N_truncate=None,
    ):
        if low_X_truncate is not None:
            ii = self.abu >= low_X_truncate
            self.iso = self.iso[ii]
            self.abu = self.abu[ii]
        if high_Z_truncate is not None:
            ii = ufunc_Z(self.iso) < high_Z_truncate
            self.iso = self.iso[ii]
            self.abu = self.abu[ii]
        if high_A_truncate is not None:
            ii = ufunc_A(self.iso) < high_A_truncate
            self.iso = self.iso[ii]
            self.abu = self.abu[ii]
        if high_N_truncate is not None:
            ii = ufunc_N(self.iso) < high_N_truncate
            self.iso = self.iso[ii]
            self.abu = self.abu[ii]
        self.normalize()
        if not self.silent:
            print(
                f" [{self.__class__.__name__}] Returning set with {len(self.abu)} isotopes."
            )

    def _append(self, iso, abu):
        # todo: add sorting? - YES
        # todo: check uniqueness
        self.iso = np.append(self.iso, iso)
        self.abu = np.append(self.abu, abu)

    def __str__(self):
        return (
            f"{self.__class__.__name__}(["
            + self.unit
            + "] "
            + ", ".join(
                [
                    f"{iso.Name():s}: " + f"{abu:8G}".strip()
                    for iso, abu in zip(self.iso, self.abu)
                ]
            )
            + ")"
        )

    def __repr__(self):
        if len(self.iso) < 5:
            return (
                f"{self.__class__.__name__}("
                + ",".join(
                    [
                        f"{iso.Name():s}=" + f"{abu:8G}".strip()
                        for iso, abu in zip(self.iso, self.abu)
                    ]
                )
                + ")"
            ).replace(" ", "")
        xyz = self.XYZ()
        s = ",".join(f"{x}={y:<g}".strip() for x, y in zip("XYZ", xyz))
        return f"{self.__class__.__name__}({s})"

    def _delete(self, iso):
        """
        remove isotope
        """
        ii = np.where(self.iso == iso)[0]
        assert len(ii) == 1
        self.iso = np.delete(self.iso, ii)
        self.abu = np.delete(self.abu, ii)

    def __delitem__(self, iso):
        """
        remove isotope
        """
        self._delete(iso)

    def massnc(self):
        """
        Print mass non-conservation.
        """
        return 1 - self.abu.sum()

    def normalize(self, total=None):
        """
        Normalize abundances to one.

        If sum == 0 just return.
        """
        abusum = self.abu.sum()
        if abusum == 0.0:
            return
        self.abu /= abusum
        if total is not None:
            self.abu *= total

    def normalized(self, total=None):
        """
        Return normalized copy of abu.
        """
        new = self.copy()
        new.normalize(total=total)
        return new

    def scaled_by(self, X=None, Y=None, Z=None):
        """
        return scaled abundance set
        """
        assert X is None or Y is None or Z is None
        if X is None and Y is None and Z is None:
            return self.copy()
        comment = list()
        if X is None:
            X = 1
        else:
            comment.append(f"X={X:12.5g}")
        if Y is None:
            Y = 1
        else:
            comment.append(f"Y={Y:12.5g}")
        if Z is None:
            Z = 1
        else:
            comment.append(f"Z={Z:12.5g}")
        comment = ", ".join(comment)

        ix = self.Z() == 1
        iy = self.Z() == 2
        iz = self.Z() >= 3

        f = 1 / np.sum(np.array([X, Y, Z] * self.XYZ()))

        new = self.copy()
        new.abu[ix] *= X * f
        new.abu[iy] *= Y * f
        new.abu[iz] *= Z * f

        new.comment += (f"Scaled by {comment} and normalized.",)

        return new

    def log_scaled_by(self, X=None, Y=None, Z=None):
        """
        return scaled abundance set for log offsets

        for example, Lodders 2020 uses log_scaled_by(Y=-0.07,Z=-0.088)
        """
        if X is not None:
            X = 10**X
        if Y is not None:
            Y = 10**Y
        if Z is not None:
            Z = 10**Z
        return self.scaled_by(X, Y, Z)

    def copy(self):
        """
        make a copy
        """
        iso = self.iso
        self.iso = None
        new = copy.deepcopy(self)
        self.iso = iso
        new.iso = iso.copy()
        return new

    def _from_bg(self, filename, silent=False):
        """
        Generate abundace set from BURN gen file.

        TODO - return tuple with multiple mixtures if file has several

               Actually, this should be a module function, not part of
               the constructor.

        TODO - add option to show comment
        """
        self.setup_logger(silent=silent)
        self.comment = stuple(self.comment)
        xre = re.compile("[-+a-zA-Z0-9.]+")
        self.iso = np.array([], dtype=object)
        self.abu = np.array([], dtype=np.float64)
        filename = Path(filename).expanduser()
        with open(filename, "r") as f:
            self.logger_file_info(f)
            for line in f:
                if line.startswith("c "):
                    self.comment += (line[2:].rstrip(),)
                elif line.startswith("m "):
                    xdata = xre.findall(line)
                    xnum = len(xdata)
                    if xnum < 2:
                        continue
                    mixture = xdata[1]
                    if self.mixture is None:
                        self.mixture = mixture
                        self.logger.info(f'Loading mixture "{self.mixture:s}".')
                    if self.mixture == mixture:
                        if xnum < 3:
                            continue
                        assert xnum % 2 == 0
                        xion = xdata[3::2]
                        xabu = xdata[2::2]
                        for i, a in zip(xion, xabu):
                            self._append(
                                I(i), np.double(a.replace("D", "E").replace("d", "e"))
                            )
        self.close_logger(timing="Data loaded in")

    def _from_dat(self, filename, silent=False):
        """
        Load abundace set from "dat" file.

        TODO - add option to show comment
        """
        self.setup_logger(silent=silent)
        self.comment = stuple(self.comment)
        xre = re.compile("[-+a-zA-Z0-9.]+")
        self.iso = np.array([], dtype=object)
        self.abu = np.array([], dtype=np.float64)
        filename = Path(filename).expanduser()
        with open(filename, "r") as f:
            self.logger_file_info(f)
            self.comment += (
                "",
                f'Generated from file "{filename}"',
                "Original file comments follow:",
                "",
            )
            for lineno, line in enumerate(f):
                if not line.startswith((";", "#")):
                    xdata = xre.findall(line)
                    xnum = len(xdata)
                    if xnum == 0:
                        continue
                    if xnum == 2:
                        xion, xabu = tuple(xdata)
                    else:
                        raise IOError(f'Bad format in line {lineno}: "{line}"')
                    self._append(I(xion), np.double(xabu))
                else:
                    self.comment += (line[2:].rstrip(),)
        self.filename = filename
        message = f"{len(self.iso):3d} isotopes loaded in"
        self.close_logger(timing=message)

    def _new_order(self):
        # reset things that depend on order
        # this should be called every time things are added or removed.
        # try:
        #     del self.approx_map
        # except:
        #     pass
        pass

    def ionset(self):
        return IonSet(self.ions)

    def sort(self):
        """
        Sort ions.
        """
        if self.is_sorted:
            return
        if self.iso.size == 0:
            return
        sort = self.iso.argsort()
        self.iso = self.iso[sort]
        self.abu = self.abu[sort]
        self._new_order()

    def write_bg(
        self,
        outfile=None,
        net=1,
        mixture=None,
        zmax=83,
        overwrite=False,
        silent=False,
        write_net=True,
    ):
        """
        Write out BURN generator.

        If outfile is file use this.
        If outfile is a filename open outfile for writing.

        We need to assert gaps around be8, b9
        (though it is not clear why KEPLER needs these gaps)

        We do need to assert, however, because ADAPNET preserves gaps
        for this reason, not to introduce gaps in new network.  Which
        also makes the work of ADAPNET easier.

        Eventually KEPLER should check for these gaps and issue errors
        if present.  Maybe it already does?
        """

        self.setup_logger(silent=silent)

        # minnet from adapnet general purpose network
        # we could load this, but for now, we just hard-code.
        # note the gaps for Be8 and B9 that appear to be requires
        # by KEPLER - though it is not clear why
        # maybe because these istopes are not in bdat
        # TODO - ask RDH to add to future versions of bdat.
        minnet = [
            [[1, 1]],
            [[1, 3]],
            [[3, 4]],
            [[6, 7]],
            [[7, 7], [9, 9]],
            [[8, 8], [10, 11]],
            [[11, 13]],
            [[13, 15]],
            [[14, 18]],
            [[17, 19]],
            [[19, 22]],
            [[21, 23]],
            [[23, 26]],
            [[25, 27]],
            [[27, 30]],
            [[30, 31]],
            [[31, 36]],
        ]

        version = 10000

        default_mixture = "x"

        card_cmt = "c"
        card_mix = "m"
        card_grd = "gg"
        card_ntw = "netw"

        zmax = min(zmax + 1, len(elements))
        netw = np.ndarray((zmax, 2, 2), dtype=np.int64)
        netw_count = np.zeros(zmax, dtype=np.int64)
        for iz, z in enumerate(minnet):
            for ia, a in enumerate(z):
                netw[iz, ia, :] = a
            netw_count[iz] = len(z)

        niso = 0
        iso = np.zeros_like(self.iso)
        abu = np.zeros_like(self.abu)
        for i, a in zip(self.iso, self.abu):
            Z = i.Z
            if Z < zmax:
                A = i.A
                n = netw_count[Z]
                if n == 0:
                    netw_count[Z] = 1
                    netw[Z, 0, :] = A
                else:
                    netw[Z, 0, 0] = min(A, netw[Z, 0, 0])
                    netw[Z, n - 1, 1] = max(A, netw[Z, n - 1, 1])
                add = True
                if n == 2:
                    if (A > netw[Z, 0, 1]) and (A < netw[Z, 1, 0]):
                        self.logger.error(f"{i.name():s} inside required gap.")
                        add = False
                if add:
                    iso[niso] = i
                    abu[niso] = a
                    niso += 1
        iso = iso[:niso]
        abu = abu[:niso]

        if mixture is None:
            try:
                mixture = self.mixture
            except:
                pass
        if mixture is None:
            mixture = default_mixture

        if outfile is None:
            f = sys.stdout
        if not isinstance(outfile, io.IOBase):
            filename = Path(outfile).expanduser()
            assert overwrite or not filename.exists(), f"file exists: {filename}"
            f = open(filename, "w")
        else:
            f = outfile

        if write_net:
            f.write(card_cmt + " COMPUTER-GENERATED BURN GENERATOR FILE\n")
            f.write(card_cmt + f" VERSION {version2human(version):s}" + "\n")
            f.write(card_cmt + " " + time.asctime(time.gmtime()) + " UTC\n")
            f.write(card_cmt + "\n")

            for c in self.comment:
                f.write(f"{card_cmt:s} {c:s}\n")
            f.write(card_cmt + "\n")
            f.write(card_cmt + f" define network (Z_max = {zmax - 1:d})\n")
            for i in range(zmax):
                nc = netw_count[i]
                if nc > 0:
                    c = " ".join(
                        ["{:3d} {:3d}".format(*netw[i, j, :]) for j in range(nc)]
                    )
                    f.write(f"{card_ntw:s} {net:d} {elements[i]:2s} {c:s}\n")
        f.write(card_cmt + "\n")
        f.write(card_cmt + " define composition " + mixture + "\n")
        for i, a in zip(iso, abu):
            f.write(f"{card_mix:s} {mixture:s} {a:14.8E} {i.Kepler:s}\n")
        if write_net:
            f.write(card_cmt + "\n")
            f.write(card_cmt + " specify grid composition (homogeneous star)\n")
            f.write(
                card_cmt + ' NOTE: call only after all "g" cards in main generator\n'
            )
            f.write(f"{card_grd:s} {net:d} {mixture:s}\n")

        if not isinstance(outfile, io.IOBase) and outfile is not None:
            f.close()

        if write_net:
            self.close_logger(timing=f'BURN generator written to "{f.name:s}" in')
        else:
            self.close_logger(timing=f"BURN mixture {mixture:s} written in")

    def write_dat(self, outfile=None, overwrite=False, silent=False):
        """
        Write out dat file.

        If file is file use this.
        If file is file name open file for write.
        """
        self.setup_logger(silent=silent)

        version = 10000

        card_cmt = ";"

        if outfile is None:
            f = sys.stdout
        elif not isinstance(outfile, io.IOBase):
            filename = Path(outfile).expanduser()
            assert overwrite or not filename.exists(), " file exists: " + filename
            f = open(filename, "wt")
        else:
            f = outfile

        f.write(card_cmt + " COMPUTER-GENERATED ABUNDANCE DATA FILE\n")
        f.write(card_cmt + f" VERSION {version2human(version):s}" + "\n")
        f.write(card_cmt + " " + time.asctime(time.gmtime()) + " UTC\n")
        f.write(card_cmt + "\n")

        for c in self.comment:
            f.write(f"{card_cmt:s} {c:s}\n")
        f.write(card_cmt + "\n")
        f.write(
            card_cmt
            + "----------------------------------------------------------------------\n"
        )
        for i, a in zip(self.iso, self._X()):
            f.write(f"{i.name():6s} {a:25.17e}\n")
        if not isinstance(outfile, io.IOBase) and outfile is not None:
            f.close()
        self.close_logger(timing=f'"{f.name:s}" written in')

    def write_compsurb(self, outfile, overwrite=False, silent=False):
        """
        Write out data file for setting compsurf BURN values.

        If file is file use this.
        If file is file name open file for write.
        """
        self.setup_logger(silent=silent)

        version = 10000

        card_cmt = "c"

        if not isinstance(outfile, io.IOBase):
            filename = Path(outfile).expanduser()
            assert overwrite or not filename.exists(), " file exists: " + filename
            f = open(os.path.expanduser(os.path.expandvars(filename)), "w")
        else:
            f = outfile

        f.write(card_cmt + " COMPUTER-GENERATED BURN COMPSURF LINK FILE\n")
        f.write(card_cmt + f" VERSION {version2human(version)}\n")
        f.write(card_cmt + f" {time.asctime(time.gmtime())} UTC\n")
        f.write(card_cmt + "\n")

        for c in self.comment:
            f.write("{card_cmt} {c}\n")
        f.write(card_cmt + "\n")
        f.write(
            card_cmt
            + "----------------------------------------------------------------------\n"
        )
        f.write("compsurb clear\n")
        for i, a in zip(self.iso, self._X()):
            f.write(f"compsurb {a:13.7E} {i.name():s}\n")
        f.write("mapsurfb\n")
        if not isinstance(outfile, io.IOBase):
            f.close()
        self.close_logger(timing=f'"{f.name:s}" written in')

    # some property routines

    @magic
    def A(self):
        return ufunc_A(self.iso)

    def Z(self):
        return ufunc_Z(self.iso)

    def N(self):
        return ufunc_N(self.iso)

    def E(self):
        return ufunc_E(self.iso)

    def isomer(self):
        return ufunc_isomer(self.iso)

    def isotope(self):
        return ufunc_isotope(self.iso)

    def element(self):
        return ufunc_element(self.iso)

    def isobar(self):
        return ufunc_isobar(self.iso)

    def isotone(self):
        return ufunc_isotone(self.iso)

    def element_name(self):
        return np.array([x.element_name() for x in self.iso])

    def element_symbol(self, upcase=True):
        return np.array([x.element_symbol(upcase=upcase) for x in self.iso])

    def idx(self):
        return ufunc_idx(self.iso)

    def isotone_idx(self):
        return ufunc_isotone_idx(self.iso)

    def isobar_idx(self):
        return ufunc_isobar_idx(self.iso)

    def element_idx(self):
        return ufunc_element_idx(self.iso)

    def isotope_idx(self):
        return ufunc_isotope_idx(self.iso)

    def isomer_idx(self):
        return ufunc_isomer_idx(self.iso)

    def isotones_idx(self):
        return np.unique(self.isotone_idx())

    def isobars_idx(self):
        return np.unique(self.isobar_idx())

    def elements_idx(self):
        return np.unique(self.element_idx())

    def isotopes_idx(self):
        return np.unique(self.isotope_idx())

    def isomers_idx(self):
        return np.unique(self.isomer_idx())

    def XYZ(self):
        """
        Return 'astronomical' X, Y, Z of composition by mass fraction
        """
        x = sum(self.X()[self.Z() <= 1])
        y = sum(self.X()[self.Z() == 2])
        z = sum(self.X()[self.Z() >= 3])
        return np.array([x, y, z])

    def xyz(self):
        """
        Return 'astronomical' X, Y, Z of composition by mol fraction
        """
        x = sum(self.Y()[self.Z() <= 1])
        y = sum(self.Y()[self.Z() == 2])
        z = sum(self.Y()[self.Z() >= 3])
        return np.array([x, y, z])

    def CNO(self):
        """
        Return 'CNO' mass fraction
        """
        Z = self.Z()
        return sum(self.X()[(6 <= Z) & (Z <= 8)])

    def cno(self):
        """
        Return 'CNO' mol fraction
        """
        Z = self.Z()
        return sum(self.Y()[(6 <= Z) & (Z <= 8)])

    def metallicity(self):
        """
        Return 'metallicity' Z of composition.
        """
        return sum(self.X()[self.Z() >= 3])

    @property
    def metals(self):
        """
        return metal mass fraction
        """
        return self.metallicity()

    def Ye(self):
        """
        Return electron to baryon ratio.
        """
        return np.sum(self.Z() * self.Y())

    def mue(self):
        """
        Return mean molecular weight per electron of composition.
        """
        return 1 / self.Ye()

    def eta(self):
        """
        Return neutron excess of mixture.
        """
        return 1 - 2 * self.Ye()

    def mu(self):
        """
        Return mean molecular weight of composition.
        """
        return 1 / np.sum((self.Z() + 1) * self.Y())

    def muI(self):
        """
        Return mean molecular weight of composition.
        """
        return 1 / np.sum(self.Y())

    abar = muI

    def zbar(self):
        """
        Return mean charge number of composition.
        """
        return self.Ye() * self.abar()

    def zmom(self, mom=1):
        """
        Return mean moment of charge number of composition.
        """
        return np.sum(self.Z() ** mom * self.Y()) * self.abar()

    def amom(self, mom=1):
        """
        Return mean moment of mass number of composition.

        Trivially, mom==1 is just abar.
        """
        return np.sum(self.A() ** mom * self.Y()) * self.abar()

    # selection functions

    def get_isotopes(self, ion):
        ion = I(ion)
        ii = np.where(ufunc_Z(self.iso) == ion.Z)[0]
        return self.iso[ii]

    def get_isotones(self, ion):
        ion = I(ion)
        ii = np.where(ufunc_N(self.iso) == ion.N)[0]
        return self.iso[ii]

    def get_isobars(self, ion):
        ion = I(ion)
        ii = np.where(ufunc_A(self.iso) == ion.A)[0]
        return self.iso[ii]

    def get_isomers(self, ion):
        ion = I(ion)
        ii = np.where(
            np.logical_and(ufunc_Z(self.iso) == ion.Z, ufunc_N(self.iso) == ion.N)
        )[0]
        return self.iso[ii]

    # here to add some general routines that return data for any Ion type
    # (isotope, isotone, isobar, element)
    def _get(self, selection=None, data="X", invalid=np.nan, exception=False):
        """
        Return data for selected (isotope, isotone, isobar, element).

        If nothing is specified, all isotopes will be returned.
        Default is mass fraction 'X'.
        """

        def output(x):
            nonlocal shape, exception, selection
            if exception:
                (index,) = np.where(np.isnan(x))
                print(index)
                if len(index) > 0:
                    raise AttributeError(
                        "Isotope(s) not found: "
                        + ", ".join([str(ion) for ion in selection[index]])
                    )
            if len(shape) == 0:
                x = x[0]
            elif len(shape) != 1:
                x = x.reshape(shape)
            return x

        # the following should go into class definition or init routine

        # setup
        uis_func = [
            ufunc_is_isomer,
            ufunc_is_isotope,
            ufunc_is_element,
            ufunc_is_isobar,
            ufunc_is_isotone,
        ]

        # we would only have to switch out this one definition to get them all...
        if data == "X":
            val_func = [
                self.isomers_X,
                self.isotopes_X,
                self.elements_X,
                self.isobars_X,
                self.isotones_X,
            ]
        elif data == "Y":
            val_func = [
                self.isomers_Y,
                self.isotopes_Y,
                self.elements_Y,
                self.isobars_Y,
                self.isotones_Y,
            ]
        elif data == "log_eps":
            val_func = [
                self.isomers_log_eps,
                self.isotopes_log_eps,
                self.elements_log_eps,
                self.isobars_log_eps,
                self.isotones_log_eps,
            ]
        elif data == "Abu":
            val_func = [
                self.isomers_Abu,
                self.isotopes_Abu,
                self.elements_Abu,
                self.isobars_Abu,
                self.isotones_Abu,
            ]
        else:
            raise ValueError(f"Invalid Data request: '{data}'.")

        # selection setup
        if selection is None:
            selection = self.iso
        shape = np.shape(selection)
        if shape == ():
            selection = np.array([selection])
        if not isinstance(selection, np.ndarray):
            selection = np.array(selection)
        if len(shape) > 1:
            selection = selection.reshape(-1)
        if selection.dtype == np.dtype(bool):
            assert len(selection) == len(self.iso)
            selection = self.iso[selection]
        else:
            # use factory function instead of Ion
            # TODO - add/use ion cache
            selection = np.array([I(ix) for ix in selection])

        selections = np.ndarray([len(uis_func), len(selection)], dtype=bool)

        # pure cases
        # here we test the most likely cases first ...
        for i, (uis, val) in enumerate(zip(uis_func, val_func)):
            selections[i, :] = uis(selection)
            if np.all(selections[i, :]):
                return output(val(selection))

        # now mixed cases
        if exception:
            invalid = np.nan
        x = np.tile(np.array(invalid, dtype=np.float64), len(selection))
        for i, val in enumerate(val_func):
            if np.any(selections[i, :]):
                (index,) = np.nonzero(selections[i, :])
                x[index] = val(selection[index])

        return output(x)

    def _X(self):
        return self.abu

    def _Y(self):
        return self.abu / self.A()

    def ppm(self, selection=None):
        """
        Return ppm (10,000 ppm = 1 mass%) for selected(isomer, isotope, isotone, isobar, element).

        If nothing is specified, all isotopes will be returned.
        """
        return self.X(selection) * 1e6

    @magic
    def X(self, selection=None):
        """
        Return mass fraction for selected (isomer, isotope, isotone, isobar, element).

        If nothing is specified, all isotopes will be returned.
        """
        return self._get(selection, data="X")

    def Y(self, selection=None):
        """
        Return mole/g for selected (isomer, isotope, isotone, isobar, element).

        If nothing is specified, all isotopes will be returned.
        """
        return self._get(selection, data="Y")

    def log_eps(self, selection=None):
        """
        Return log(eps) for selected (isomer, isotope, isotone, isobar, element).

        If nothing is specified, all isotopes will be returned.
        """
        return self._get(selection, data="log_eps")

    def Abu(self, selection=None):
        """
        Return A with A(Si) = 10**6 for selected (isomer, isotope, isotone, isobar, element).

        If nothing is specified, all isotopes will be returned.
        """
        return self._get(selection, data="Abu")

    # add functions to compute [] , \delta - those need to be in SolAbu

    # generic selection routine
    def _selection(
        self,
        x,
        idx,
        selection=None,
        missing=np.nan,
        check=None,
        convert=None,
    ):
        if selection is None:
            return x
        # we should really use IonList objects
        if np.shape(selection) == ():
            selection = np.array([selection])
        selection = I(selection)
        if not np.all(check(selection)):
            raise KeyError("Selection contained invalid entries.")
        selection = convert(selection)
        y = np.ndarray(selection.shape, dtype=x.dtype)
        ii = np.argsort(idx)
        sel = np.minimum(np.searchsorted(idx[ii], selection), len(ii) - 1)
        jj = idx[ii][sel] == selection
        y[jj] = x[ii][sel][jj]
        y[~jj] = missing
        return y

    def _generic_selection(
        self,
        data=None,
        selection=None,
        mode=None,
        missing=np.nan,
        return_selection=False,
        **kwargs,
    ):
        if data is None:
            data = self.abu
        elif data == "X":
            data = self._X()
        elif data == "Y":
            data = self._Y()
        if mode is None:
            convert = ufunc_idx
            check = ufunc_is_ion
            keys = self.idx()
        else:
            keys = (self.__getattribute__(mode + "_idx")(),)
            check = get_ufunc(mode)
            convert = get_ufunc_idx(mode)
        x, idx = project(
            data,
            keys,
            return_values=True,
        )
        if selection is None:
            result = x
        else:
            result = self._selection(
                x,
                idx,
                selection=selection,
                check=check,
                convert=convert,
                missing=missing,
                **kwargs,
            )
        if return_selection:
            if selection is None:
                selection = ufunc_ion_from_idx(idx)
            return result, selection
        return result

    # ------------------------------
    # isomer routines
    # ------------------------------
    def _isomers_selection(self, **kwargs):
        return self._generic_selection(mode="isomer", **kwargs)

    def isomers_X(self, selection=None):
        """
        Return isomer mass fractions.
        """
        return self._isomers_selection(data="X", selection=selection)

    def isomers_Y(self, selection=None):
        """
        Return isomer number fractions (mol per gram).
        """
        return self._isomers_selection(data="Y", selection=selection)

    def isomers_set(self, selection=None):
        """
        Return set of projected isomers.
        """
        abu, iso = self._isomers_selection(
            data="X", selection=selection, return_selection=True
        )
        return AbuSet(
            abu=abu,
            iso=iso,
            comment=self.comment,
        )
        # we may want to add something signifying isomers ...

    def isomers_log_eps(self, selection=None):
        """
        Return isomere log(eps).

        log10(# of atoms relative to H) + 12
        """
        h = self.Z() == 1
        A = self.A()
        x = sum(self._X()[h] / A[h])
        y = self.isomers_Y(selection)
        return np.log10(y / x) + 12

    def isomers_Abu(self, selection=None):
        """
        Return isomere Abundance A(Si) = 10**6.

        number of atoms relative to Si
        """
        si = self.Z() == 14
        A = self.A()
        x = sum(self._X()[si] / A[si])
        y = self.isomers_Y(selection)
        return y / x * 1e6

    def isomers(self):
        """
        Return list of all isomers in abundance set.
        """
        idx = np.unique(self.isomer_idx())
        ii = np.where(idx != VOID_IDX)
        return ufunc_ion_from_idx(idx[ii])

    # ------------------------------
    # isotope routines
    # ------------------------------
    def _isotopes_selection(self, **kwargs):
        return self._generic_selection(mode="isotope", **kwargs)

    def isotopes_X(self, selection=None):
        """
        Return isotope mass fractions.
        """
        return self._isotopes_selection(data="X", selection=selection)

    def isotopes_Y(self, selection=None):
        """
        Return isotope number fractions (mol per gram).
        """
        return self._isotopes_selection(data="Y", selection=selection)

    def isotopes_log_eps(self, selection=None):
        """
        Return isotope log(eps).

        log10(# of atoms relative to H) + 12
        """
        h = self.Z() == 1
        A = self.A()
        x = sum(self._X()[h] / A[h])
        y = self.isotopes_Y(selection)
        return np.log10(y / x) + 12

    def isotopes_Abu(self, selection=None):
        """
        Return isotope Abundance A(Si) = 10**6.

        number of atoms relative to Si
        """
        si = self.Z() == 14
        A = self.A()
        x = sum(self._X()[si] / A[si])
        y = self.isotopes_Y(selection)
        return y / x * 1e6

    def isotopes(self):
        """
        Return list of all isotopes in abundance set.
        """
        idx = np.unique(self.isotope_idx())
        ii = np.where(idx != VOID_IDX)
        return ufunc_ion_from_idx(idx[ii])

    def isotopes_set(self, selection=None):
        """
        Return set of projected isotopes.
        """
        abu, iso = self._isotopes_selection(
            data="X", selection=selection, return_selection=True
        )
        return AbuSet(
            abu=abu,
            iso=iso,
            comment=self.comment,
        )
        # we may want to add something signifying isotopes ...

    # ------------------------------
    # element routines
    # ------------------------------
    def _elements_selection(self, **kwargs):
        return self._generic_selection(mode="element", **kwargs)

    def elements(self):
        """
        Return list of all elements in abundance set.
        """
        idx = np.unique(self.element_idx())
        ii = np.where(idx != VOID_IDX)
        return ufunc_ion_from_idx(idx[ii])

    def elements_Z(self):
        """
        Return charge number of all elements.
        """
        return np.unique(self.Z())

    def elements_X(self, selection=None):
        """
        Return elemental mass fractions.
        """
        return self._elements_selection(data="X", selection=selection)

    def elements_Y(self, selection=None):
        """
        Return elemental number fractions (mol per gram).
        """
        return self._elements_selection(data="Y", selection=selection)

    def elements_log_eps(self, selection=None):
        """
        Return elemental log(eps).

        log10(# of atoms relative to H) + 12
        """
        h = self.Z() == 1
        A = self.A()
        x = sum(self._X()[h] / A[h])
        y = self.elements_Y(selection)
        return np.log10(y / x) + 12

    def elements_Abu(self, selection=None):
        """
        Return element abundance A(Si) = 10**6.

        number of atoms relative to Si
        """
        si = self.Z() == 14
        A = self.A()
        x = sum(self._X()[si] / A[si])
        y = self.elements_Y(selection)
        return y / x * 1e6

    def elements_name(self):
        """
        Return element name of each isotope.
        """
        return np.array([Elements[x] for x in self.elements_Z()])

    def elements_set(self, selection=None):
        """
        Return set of projected elements.
        """
        abu, iso = self._elements_selection(
            data="X", selection=selection, return_selection=True
        )
        return AbuSet(
            abu=abu,
            iso=iso,
            comment=self.comment,
        )
        # we may want to add something signifying elements ...

    # ------------------------------
    # isobar routines
    # ------------------------------
    def _isobars_selection(self, **kwargs):
        return self._generic_selection(mode="isobar", **kwargs)

    def isobars(self):
        """
        Return list of all elements in abundance set.
        """
        idx = np.unique(self.isobar_idx())
        ii = np.where(idx != VOID_IDX)
        return ufunc_ion_from_idx(idx[ii])

    def isobars_A(self):
        """
        Return mass number of all isobars.
        """
        return np.unique(self.A())

    def isobars_X(self, selection=None):
        """
        Return isobar mass fractions.
        """
        return self._isobars_selection(data="X", selection=selection)

    def isobars_Y(self, selection=None):
        """
        Return isobar number fractions (mol per gram).
        """
        return self._isobars_selection(data="Y", selection=selection)

    def isobars_log_eps(self, selection=None):
        """
        Return isobar log(eps).

        log10(# of atoms relative to H) + 12
        """
        h = self.Z() == 1
        A = self.A()
        x = sum(self._X()[h] / A[h])
        y = self.isobars_Y(selection)
        return np.log10(y / x) + 12

    def isobars_Abu(self, selection=None):
        """
        Return isobar abundance A(Si) = 10**6.

        number of atoms relative to Si
        """
        si = self.Z() == 14
        A = self.A()
        x = sum(self._X()[si] / A[si])
        y = self.isobars_Y(selection)
        return y / x * 1e6

    def isobars_set(self, selection=None):
        """
        Return set of projected isobars.
        """
        abu, iso = self._isobars_selection(
            data="X", selection=selection, return_selection=True
        )
        return AbuSet(abu=abu, iso=iso, comment=self.comment)
        # we may want to add something signifying isobars ...

    # ------------------------------
    # isotone routines
    # ------------------------------
    def _isotones_selection(self, **kwargs):
        return self._generic_selection(mode="isotone", **kwargs)

    def isotones(self):
        """
        Return list of all isotones in abundance set.
        """
        idx = np.unique(self.isotone_idx())
        ii = np.where(idx != VOID_IDX)
        return ufunc_ion_from_idx(idx[ii])

    def isotones_N(self):
        """
        Return neutron number of all isotones.
        """
        return np.unique(self.N())

    def isotones_X(self, selection=None):
        """
        Return isotone mass fractions.
        """
        return self._isotones_selection(data="X", selection=selection)

    def isotones_Y(self, selection=None):
        """
        Return isotone number fractions (mol per gram).
        """
        return self._isotones_selection(data="Y", selection=selection)

    def isotones_log_eps(self, selection=None):
        """
        Return isotone log(eps).

        log10(# of atoms relative to H) + 12
        """
        h = self.Z() == 1
        A = self.A()
        x = sum(self._X()[h] / A[h])
        y = self.isotones_Y(selection)
        return np.log10(y / x) + 12

    def isotones_Abu(self, selection=None):
        """
        Return isotone abundance A(Si) = 10**6.

        number of atoms relative to Si
        """
        si = self.Z() == 14
        A = self.A()
        x = sum(self._X()[si] / A[si])
        y = self.isotones_Y(selection)
        return y / x * 1e6

    def isotones_set(self, selection=None):
        """
        Return set of projected isotones.
        """
        abu, iso = self._isotones_selection(
            data="X", selection=selection, return_selection=True
        )
        return AbuSet(abu=abu, iso=iso, comment=self.comment)
        # we may want to add something signifying isotones ...

    def as_dict(self, selection=None, massfrac=True):
        if massfrac is True:
            data = "X"
        else:
            data = "Y"
        abu, iso = self._generic_selection(
            data=data,
            selection=selection,
            return_selection=True,
        )
        return {str(i): a for i, a in zip(iso, abu)}

    # =========================
    # general access interfaces
    # =========================
    def __getitem__(self, index):
        try:
            return self._get(index)
        except:
            raise AttributeError(f"Isotope not found: {index}")

    def __setitem__(self, index, item):
        # TODO add isotope if not in list?
        # maybe add parameter allow_new
        try:
            # this does not work for numpy
            if isinstance(index, str):
                index = I(index)
            if ision(index):
                index = np.where(self.iso == index)[0][0]
            if not np.isreal(item):
                raise AttributeError("Abundance needs to be a real number.")
            self.abu[index] = item
        except:
            raise AttributeError(f"Isotope to set not found: {index}")

    def __len__(self):
        """
        Return number of isotopes.
        """
        return len(self.iso)

    def __getattr__(self, attr):
        try:
            x = self._get(attr, invalid=np.nan)
        except AttributeError:
            pass
        else:
            if not np.isnan(x):
                return x
        return super().__getattribute__(attr)

    # maybe this should not be done as it makes adding new attributes rather slow...
    def __setattr__(self, attr, value):
        if not self._init and attr not in self.__dict__ and not attr == "__class__":
            if "iso" in self.__dict__:
                x = np.where(self.iso == attr)[0]
                if len(x) == 1:
                    i = x[0]
                    fac = 1 - value
                    self.abu[i] = 0.0
                    self._normalize(fac)  # not implemented
                    self.abu[i] = value
                    return
        super().__setattr__(attr, value)

    def _normalize(self, fac=None):
        if fac is None:
            fac = 1 - np.sum(self.X())
        self.abu *= fac

    def __iter__(self):
        for i, a in zip(self.iso, self.abu):
            yield i, a

    def __call__(self, *args, **kwargs):
        """
        This routine needs a lot more work ...

        TODO:
        1) add new isotopes to abundance pattern
        2) renormailzie
        3) check value is in Range 0<v<1
        """
        for k in args:
            x = np.where(self.iso == k)[0]
            if len(x) == 1:
                return self.abu[x[0]]
            else:
                raise AttributeError("Isotope not in list")
        for k, v in kwargs.items():
            if not 0 <= v <= 1:
                raise ValueError("Abundance for " + k + " is out of range.")
            x = np.where(self.iso == k)[0]
            if len(x) == 1:
                self.abu[x[0]] = v
            else:
                # TODO - add isotope instead (if valid)
                raise AttributeError("Isotope not in list")

    def __contains__(self, other):
        """
        Determine whether AbuSet contains other (iso) or AbuSet (vector).
        """
        if isinstance(other, AbuSet):
            return self.contains(other).all()
        ll = np.in1d(other, self.iso)
        if np.all(ll):
            return True
        return False

    def contains(self, iso):
        """
        Determine whether AbuSet contains iso (scalar) or which isotopes in iso are present (vector).

        If scalar argument: return True/False.
        If vector argument return array with result for each element.
        """
        if np.isscalar(iso):
            return iso in self
        else:
            if isinstance(iso, AbuSet):
                iso = iso.iso
            else:
                iso = I(iso)
            return self._in1d(iso, self.iso)

    @staticmethod
    def _in1d(ar1, ar2):

        """
        replace non-working in1d that uses mergesort

        we convert to index and then call the original in1d
        ... which is sort of what would have happened anyway ...
        """
        ar1 = ufunc_idx(ar1)
        ar2 = ufunc_idx(ar2)
        return np.in1d(ar1, ar2)
        # ar = np.concatenate( (ar1, ar2) )
        # order = ar.argsort()
        # sar = ar[order]
        # equal_adj = (sar[1:] == sar[:-1])
        # flag = np.concatenate( (equal_adj, [False] ) )
        # indx = order.argsort()[:len( ar1 )]
        # return flag[indx]

    def ions(self):
        """
        return IonList
        """
        return IonList(self.iso)

    def finite(self):
        """
        return finite subset (useful for ratios)
        """
        ii = np.isfinite(self.abu)
        return AbuSet(self.iso[ii], self.abu[ii])

    def ztruncated(self, Zmax):
        """
        return subset with Z <= Zmax
        """
        ii = self.ions().Z <= Zmax
        return AbuSet(self.iso[ii], self.abu[ii])

    def norm(self):
        """
        return total of abu
        """
        return np.sum(self.abu)

    def subset(self, mask):
        """
        return abuset subset based on selection
        """
        mask = np.atleast_1d(np.asarray(mask)).reshape(-1)
        assert len(mask) == len(self.iso)
        return AbuSet(
            iso=self.iso[mask],
            abu=self.abu[mask],
        )

    def __and__(self, other):
        return self.subset(other)

    # ------------------------------
    # arithmetic operations
    # ------------------------------

    # TODO - add hash for speed-up

    _zero_missing = None

    @staticmethod
    def _return_matching(
        iso1, abu1, iso2, abu2, missing=None, missing1=None, missing2=None
    ):
        """
        if missing is None: return intersection of isotopes, otherwise union

        if only one missing value is provided, use other set as basis
        """
        if missing1 is None:
            missing1 = missing
        if missing2 is None:
            missing2 = missing
        history = ionmatch[(iso1, missing1, iso2, missing2)]
        if history is not None:
            iso, d1, d2 = history
            if missing1 is None and missing2 is None:
                ii1 = d1
                ii2 = d2
                abu_1 = abu1[ii1]
                abu_2 = abu2[ii2]
            elif missing1 is None and missing2 is not None:
                ii1 = d1
                kk, jj = d2
                abu_1 = abu1[ii1]
                abu_2 = np.full(len(iso), missing2, dtype=np.float64)
                abu_2[kk] = abu2[jj]
            elif missing1 is not None and missing2 is None:
                ii2 = d2
                kk, jj = d1
                abu_2 = abu1[ii2]
                abu_1 = np.full(len(iso), missing1, dtype=np.float64)
                abu_1[kk] = abu2[jj]
            else:
                ii1 = d1
                ii2 = d2
                abu_1 = np.full(len(iso), missing1, dtype=np.float64)
                abu_2 = np.full(len(iso), missing2, dtype=np.float64)
                abu_1[ii1] = abu1
                abu_2[ii2] = abu2
            return iso, abu_1, abu_2
        idx1 = ufunc_idx(iso1)
        idx2 = ufunc_idx(iso2)
        if missing1 is None and missing2 is None:
            idx = np.sort(np.intersect1d(idx1, idx2))
            ii1 = np.argsort(idx1)
            ii2 = np.argsort(idx2)
            idx1 = idx1[ii1]
            idx2 = idx2[ii2]
            ii1 = ii1[np.searchsorted(idx1, idx)]
            ii2 = ii2[np.searchsorted(idx2, idx)]
            abu_1 = abu1[ii1]
            abu_2 = abu2[ii2]
            d1, d2 = ii1, ii2
        elif missing1 is None and missing2 is not None:
            ii1 = np.argsort(idx1)
            abu_1 = abu1[ii1]
            idx = idx1[ii1]
            abu_2 = np.full(len(idx), missing2, dtype=np.float64)
            sel = np.minimum(np.searchsorted(idx, idx2), len(idx) - 1)
            jj = idx[sel] == idx2
            kk = sel[jj]
            abu_2[kk] = abu2[jj]
            d1, d2 = ii1, (kk, jj)
        elif missing1 is not None and missing2 is None:
            ii2 = np.argsort(idx2)
            abu_2 = abu2[ii2]
            idx = idx2[ii2]
            abu_1 = np.full(len(idx), missing1, dtype=np.float64)
            sel = np.minimum(np.searchsorted(idx, idx1), len(idx) - 1)
            jj = idx[sel] == idx1
            kk = sel[jj]
            abu_1[kk] = abu1[jj]
            d1, d2 = (kk, jj), ii2
        else:
            idx = np.sort(np.union1d(idx1, idx2))
            abu_1 = np.full(len(idx), missing1, dtype=np.float64)
            abu_2 = np.full(len(idx), missing2, dtype=np.float64)
            ii1 = np.searchsorted(idx, idx1)
            ii2 = np.searchsorted(idx, idx2)
            abu_1[ii1] = abu1
            abu_2[ii2] = abu2
            d1, d2 = ii1, ii2
        iso = ufunc_ion_from_idx(idx)
        ionmatch[(iso1, missing1, iso2, missing2)] = (iso, d1, d2)
        return iso, abu_1, abu_2

    def __truediv__(self, other):
        """
        return abundance ratio

        TODO - check for X or Y, unit
        """
        if isinstance(other, AbuSet):
            iso, abu1, abu2 = self._return_matching(
                self.iso, self.X(), other.iso, other.X(), missing1=self._zero_missing
            )
            new = abu1 / abu2
            comment = other.comment
        else:
            try:
                new = self.abu / other
            except:
                return NotImplemented
            iso = self.iso
            comment = str(other)
        return AbuSet(
            iso,
            new,
            unit="abundace ratio",
            comment=" / ".join(stuple(self.comment, comment)),
        )

    def __rtruediv__(self, other):
        """
        return abundance ratio

        TODO - check for X or Y, set unit
        """
        if isinstance(other, AbuSet):
            iso, abu1, abu2 = self._return_matching(
                self.iso, self.X(), other.iso, other.X(), missing2=self._zero_missing
            )
            new = abu2 / abu1
            comment = other.comment
        else:
            try:
                new = other / self.abu
            except:
                return NotImplemented
            comment = str(other)
            iso = self.iso
        return AbuSet(
            iso,
            new,
            unit="abundace ratio",
            comment=" / ".join(stuple(comment, self.comment)),
        )

    def __floordiv__(self, other):
        """
        abundace ratio, try working subset match

        TODO - check for X or Y, unit
        """
        if not isinstance(other, AbuSet):
            return NotImplemented
        iso = []
        new = []
        mis = []
        for i in other.iso:
            a = self[i]
            b = other[i]
            if not np.isnan(a) and b > 0:
                iso.append(i)
                new.append(a / b)
            # elif a == b == 0:
            #     iso.append(i)
            #     new.append(1.)
            # elif a > 0 and  b == 0:
            #     iso.append(i)
            #     new.append(np.inf)
            else:
                mis.append(i)
                break
        if len(mis) > 0:
            iso = []
            new = []
            mis = []
            for i in self.iso:
                a = self[i]
                b = other[i]
                if not np.isnan(a) and b > 0:
                    iso.append(i)
                    new.append(a / b)
                # elif a == b == 0:
                #     iso.append(i)
                #     new.append(1.)
                # elif a > 0 and  b == 0:
                #     iso.append(i)
                #     new.append(np.inf)
                else:
                    mis.append(i)
                    break
        if len(mis) > 0:
            print(f"[floordiv] missing {mis}")
            return NotImplemented
        return AbuSet(
            iso,
            new,
            unit="abundace ratio",
            comment=" // ".join(stuple(self.comment, other.comment)),
        )

    def __add__(self, other):
        """
        return sum of abundances

        TODO - check for X or Y
        """
        if not isinstance(other, AbuSet):
            raise NotImplementedError()
        iso, abu1, abu2 = self._return_matching(
            self.iso, self.X(), other.iso, other.X(), missing=0
        )
        new = abu1 + abu2
        comment = other.comment
        return AbuSet(
            iso, new, unit=self.unit, comment=" + ".join(stuple(self.comment, comment))
        )

    def __sub__(self, other):
        """
        return difference of abundances

        TODO - check for X or Y
        """
        if not isinstance(other, AbuSet):
            return NotImplemented
        iso, abu1, abu2 = self._return_matching(
            self.iso, self.X(), other.iso, other.X(), missing=0
        )
        new = abu1 - abu2
        comment = other.comment
        return AbuSet(
            iso, new, unit=self.unit, comment=" - ".join(stuple(self.comment, comment))
        )

    def __mul__(self, other):
        """
        return product of abundances

        TODO - check for X or Y
        """
        if isinstance(other, AbuSet):
            iso, abu1, abu2 = self._return_matching(
                self.iso, self.X(), other.iso, other.X(), missing=self._zero_missing
            )
            new = abu1 * abu2
            comment = other.comment
        else:
            try:
                new = self.abu * other
            except:
                return NotImplemented
            comment = str(other)
            iso = self.iso
        return AbuSet(
            iso,
            new,
            unit="abundace product",
            comment=" * ".join(stuple(self.comment, comment)),
        )

    __rmul__ = __mul__

    def __pow__(self, other):
        """
        return power of abundances

        TODO - check for X or Y
        """
        if isinstance(other, AbuSet):
            iso, abu1, abu2 = self._return_matching(
                self.iso, self.X(), other.iso, other.X(), missing1=self._zero_missing
            )
            new = abu2**abu1
            comment = other.comment
        else:
            try:
                new = self.abu**other
            except:
                return NotImplemented
            comment = str(other)
            iso = self.iso
        return AbuSet(
            iso,
            new,
            unit="abundace power",
            comment=" ** ".join(stuple(self.comment, comment)),
        )

    def __neg__(self):
        """
        return negative of abundance
        """
        new = -self.abu
        return AbuSet(
            self.iso, new, unit=self.unit, comment=" - ".join(stuple("", self.comment))
        )

    def index(self, iso):
        """
        Return isotope index in iso/abu array

        accepts scalar, np.ndarray, list, tuple, set
        """
        if isinstance(iso, str) or not isinstance(iso, (Iterable, np.ndarray)):
            return np.where(self.iso == iso)[0][0]
        elif isinstance(iso, np.ndarray):
            ii = np.ndarray(iso.shape, dtype=np.int64)
            for i, k in enumerate(iso.flat):
                ii[i] = np.where(self.iso == k)[0][0]
            return ii
        elif isinstance(iso, Iterable):
            ii = []
            for k in iso:
                if isinstance(k, str):
                    i = np.where(self.iso == k)[0][0]
                else:
                    i = self.index(k)
                ii += [i]
            if not isinstance(iso, list):
                ii = type(iso)(ii)
            return ii
        else:
            raise AttributeError("Argument type not supported.")

    def kepler_write_compsurb(self, filename):
        """
        Write link file to set KEPLER surface BURN abundances.
        """
        filename = Path(filename).expanduser()
        with open(filename, "wt") as f:
            f.write("compsurb clear\n")
            for i in self.ions():
                x = self.X(i)
                if x > 0:
                    f.write(f"compsurb {x:12.5e} {i.name():<5s}\n")
            f.write("compsurb show\n")


class IonMatchDB(object):
    """
    DB to hash ion list matches to accelerate abundace operarions

    TODO - use IonList, make unmutable version with hash
    """

    def __init__(self):
        self.db = dict()

    def hash(self, l):
        if isinstance(l, int):
            return int
        if isinstance(l, np.ndarray):
            return hash(tuple(l.tolist()))
        raise Exception(f"Unknown type for hash: {type(l)}")

    def set(self, l1, m1, l2, m2, iso, d1, d2):
        l1 = self.hash(l1)
        l2 = self.hash(l2)
        t1 = m1 is None
        t2 = m2 is None
        ii = (l1, t1, l2, t2)
        if ii in self.db:
            print(f" [{self.__class__.__name__}] replacing entry.")
        self.db[ii] = (iso, d1, d2)

    def get(self, l1, m1, l2, m2):
        l1 = self.hash(l1)
        l2 = self.hash(l2)
        t1 = m1 is None
        t2 = m2 is None
        ii = (l1, t1, l2, t2)
        x = self.db.get(ii, None)
        if x is None:
            ii = (l2, t2, l1, t1)
            x = self.db.get(ii, None)
            if x is not None:
                x = (x[0], x[2], x[1])
        return x

    def __setitem__(self, key, value):
        self.set(*key, *value)

    def __getitem__(self, key):
        return self.get(*key)


ionmatch = IonMatchDB()


class AbuRat:
    """
    Compute log number ratios such as '[X]' or '[X/H]'.

    use as

      r = AbuRat(mydata, refsolar)
      r['C N O/Fe']
      r['C12/C13']
      r['N']

    use spaces as separator between species, though some run-in cases may work;
    isomers need to be separated by spaces
    special species specs such as A=150 work as well
    """

    pattern = re.compile(
        r"|".join(
            (
                r"[ANZ][:=][0-9]{,3}",
                r"[A-Za-z][a-z]?[1-9][0-9]{,2}(?:g|(?:m[0-9]{,3}))",
                r"[A-Z][a-z]?-?[0-9]*",
                r"[a-z]{1,2}-?[0-9]*",
            )
        )
    )

    def __init__(self, val, ref="As09"):
        assert isinstance(val, AbuSet)
        if isinstance(ref, str):
            import abusets

            ref = abusets.SolAbu(ref)
        assert isinstance(ref, AbuSet)
        self.val = val
        self.ref = ref

    def strsplit(self, s):
        """
        Split isotopes by space, and try to separate some run-in cases
        as well.

        Does not work for isomers and 'special' names.
        """
        res = s.split()
        if len(res) == 1:
            # will not work for isomers
            r2 = self.pattern.findall(s)
            if len(r2) > len(res):
                res = r2
        return res

    def __getitem__(self, key):
        """
        convention:
          x = 'X' or Ion()
        support following formats for args
        x
        'X/Y', (x, y) for ratios
        'A B C / X Y Z', ([a, b, c], [x, y, z])
        'A B C', ([a, b, c])

        unfortunately arguments (a, b) will treted as [a / b] not [a b]

        for consistency, whereas tuple length different from 2 could be
        treated correctly - and differently from the ambiguous 2-tuples,
        we just raise an error here.
        """
        if isinstance(key, tuple):
            assert (
                len(key) == 2
            ), """
            supporting only tuples of length 2 that are used for ratios;
            use lists for grouping or __call__ method.
            """
            args = key
        else:
            args = (key,)
        return self.__call__(*args)

    def __call__(self, *args):
        """
        This method is the same as '[..]' but can accept tuples for
        grouping.
        """
        val = []
        ref = []
        assert len(args) in (1, 2)
        if len(args) == 1:
            arg = args[0]
            if isinstance(arg, str):
                parts = arg.split("/")
                assert len(parts) in (1, 2)
                val = self.strsplit(parts[0])
                if len(parts) == 2:
                    ref = self.strsplit(parts[1])
            elif np.iterable(arg):
                val = list(arg)
            else:
                val = [arg]
        else:
            arg = args[0]
            if isinstance(arg, str):
                val = self.strsplit(arg)
            elif np.iterable(arg):
                val = list(arg)
            else:
                val = [arg]
            arg = args[1]
            if isinstance(arg, str):
                ref = self.strsplit(arg)
            elif np.iterable(arg):
                ref = list(arg)
            else:
                ref = [arg]
        assert len(val) > 0
        val = I(list(val))
        ref = I(list(ref))

        res = np.sum(self.val.Y(val)) / np.sum(self.ref.Y(val))
        if len(ref) > 0:
            res /= np.sum(self.val.Y(ref)) / np.sum(self.ref.Y(ref))
        return np.log10(res)
