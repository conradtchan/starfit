"""
Classes for reading and writing STARDB files.
"""

import bz2
import gzip
import hashlib
import lzma
import os
import re
import sys
import textwrap
from collections.abc import Iterable, Mapping
from copy import copy
from enum import IntEnum
from pathlib import Path

import numpy as np
import scipy.sparse

from .abuset import AbuData, IonList
from .human import byte2human
from .isotope import ion as I
from .loader import _loader, loader
from .logged import Logged
from .utils import CachedAttribute, prod, xz_file_size
from .uuidtime import UUID1

# used to replace '^' as anker point
_db_path = "~/starfit_data/db"


class StarDB(AbuData, Logged):
    """
    Class for reading STARDB binary files.

    Compressed files with .gz will be automatically uncompressed.
    This may fail, however, if the file is bigger than 2GB or 4GB.

    Compressed files with .bz2 will be automatically uncompressed.
    This is currently rather inefficient because to determine the
    file size the entire stream need to be read first.

    FIELD_TYPES:
      0 UNDEFINED Undefined
      1 BYTE Byte
      2 INT Integer
      3 LONG Longword integer
      4 FLOAT Floating point
      5 DOUBLE Double-precision floating
      6 COMPLEX Complex floating
      7 STRING String
      8 STRUCT Structure [NOT ALLOWED]
      9 DCOMPLEX Double-precision complex
     10 POINTER Pointer
     11 OBJREF Object reference
     12 UINT Unsigned Integer
     13 ULONG Unsigned Longword Integer
     14 LONG64 64-bit Integer
     15 ULONG64 Unsigned 64-bit Integer
     16 ASCII8 64-it (8 byte) ASCII String [StarDB]
    """

    current_version = 10200
    _extension = "stardb"

    sys_is_le = sys.byteorder == "little"
    native_byteorder = "<" if sys_is_le else ">"
    defaultbyteorder = native_byteorder

    # this should actually become an information class
    abundance_type_names = (
        "isomer",
        "isotope",
        "element",
        "mass number (isobar)",
        "neutron number (isotone)",
    )

    class AbundanceType(IntEnum):
        isomer = 0
        isotope = 1
        element = 2
        isobar = 3
        isotone = 4

    abundance_class_names = (
        "all (raw)",
        "rad (raw + decays)",
        "dec (stable subset of radiso)",
        "mix (mix of different types for special purpose)",
    )

    class AbundanceClass(IntEnum):
        raw = 0
        rad = 1
        dec = 2
        mix = 3

    abundance_unit_names = (
        "mass (g)",
        "mass (solar)",
        "mol",
        "mass fraction (X)",
        "mol fraction (YPS)",
        "log epsion",
        "[ ]",
        "production factor (10^[])",
        "log [Y/Si] + 6",
    )

    class AbundanceUnit(IntEnum):
        g = 0
        solar_mass = 1
        mol = 2
        mass_fraction = 3
        mol_fraction = 4
        log_eps = 5
        bracket = 6
        production_factor = 7
        log_abu = 8

    abundance_total_names = (
        "initial (total, not normalized, fallback is missing mass)",
        "ejecta",
    )

    class AbundanceTotal(IntEnum):
        initial = 0
        ejecta = 1

    abundance_data_names = (
        "all ejecta (SN ejecta + wind)",
        "SN ejecta, no wind",
        "wind only",
        "ejecta including fallback and wind (outside piston)",
    )

    class AbundanceData(IntEnum):
        all_ejecta = 0
        SN_ejecta = 1
        piston = 2

    abundance_sum_names = (
        "mass fraction",
        "number fraction",
    )

    class AbundanceSum(IntEnum):
        mass_fraction = 0
        number_fraction = 1

    # add whether normalized where appropriate ???

    flag_names = (
        "parameter",
        "property",
    )

    class Flags(IntEnum):
        parameter = 0
        property = 1

    type_names = (
        "UNDEFINED",
        "BYTE",
        "INT",
        "LONG",
        "FLOAT",
        "DOUBLE",
        "COMPLEX",
        "STRING",
        "STRUCT",
        "DCOMPLEX",
        "POINTER",
        "OBJREF",
        "UINT",
        "ULONG",
        "LONG64",
        "ULONG64",
        "ASCII64",
    )

    class Type(IntEnum):
        undefined = 0
        byte = 1
        int = 2
        long = 3
        float = 4
        double = 5
        float64 = 5
        complex = 6
        string = 7
        struct = 8
        dcomplex = 9
        pointer = 10
        objref = 11
        unit = 12
        ulong = 13
        long64 = 14
        int64 = 14
        ulong64 = 15
        uint64 = 15
        ascii64 = 16

    class SignatureError(Exception):
        """
        Exception raised when signature could not be read.
        """

        def __init__(self, filename):
            """
            Store file name that causes error.
            """
            self.filename = filename

        def __str__(self):
            """
            Return error message.
            """
            return f"Error reading signature from file {self.filename:s}."

    class VersionError(Exception):
        """
        Exception raised when version mismatch.
        """

        def __init__(self):
            """
            Just set up.
            """

        def __str__(self):
            """
            Return error message.
            """
            return "Version Error."

    class IntegrityError(Exception):
        """
        Exception raised when file integrity seems broken.
        """

        def __init__(self):
            """
            Just set up.
            """

        def __str__(self):
            """
            Return error message.
            """
            return "File Integrity Error."

    class DataError(Exception):
        """
        Exception raised when data seems faulty.
        """

        def __init__(self):
            """
            Just set up.
            """

        def __str__(self):
            """
            Return error message.
            """
            return "Data seems faulty."

    def __init__(self, filename=None, db=None, **kwargs):
        """
        Initialize data fields and open file.

        Optionally the byte order can be specified.
        The default is big endian.
        """

        kw = kwargs.copy()
        kw["filename"] = filename
        kw["db"] = db

        if filename is not None:
            self._from_file(**kw)
        elif db is not None:
            self._from_db(**kw)
        else:
            self._from_data(**kw)

    def abu(self):
        """
        Overwrite inherited routine.

        There is more different kinsd of data than the original
        version can handle.
        """
        return self.data

    def _set_dtypes(self):
        self.swapbyteorder = self.byteorder != self.native_byteorder
        try:
            if self.swapbyteorder:
                self.logger.info("Swapping endian.")
            else:
                self.logger.info("Not swapping endian.")
        except AttributeError:
            pass

        self.dtype_i8 = np.dtype(np.int64).newbyteorder(self.byteorder)
        self.dtype_u8 = np.dtype(np.uint64).newbyteorder(self.byteorder)
        self.dtype_f8 = np.dtype(np.float64).newbyteorder(self.byteorder)
        self.dtype_s8 = np.dtype((np.bytes_, 8)).newbyteorder(self.byteorder)

        self.dtype_i8 = np.dtype(np.int64)
        self.dtype_u8 = np.dtype(np.uint64)
        self.dtype_f8 = np.dtype(np.float64)
        self.dtype_s8 = np.dtype((np.bytes_, 8))

        self.dtypes = np.zeros(len(self.type_names), dtype=object)
        self.dtypes[
            [
                self.Type.float64,
                self.Type.int64,
                self.Type.uint64,
                self.Type.ascii64,
            ]
        ] = [self.dtype_f8, self.dtype_i8, self.dtype_u8, self.dtype_s8]

        self._dtypes = {
            np.uint64: self.Type.uint64,
            np.int64: self.Type.int64,
            np.float64: self.Type.float64,
            np.bytes_: self.Type.ascii64,
        }

    def _from_db(self, db=None, **kwargs):
        """
        Initialize from existing db.

        similar to copy but allow derived class initializtion
        """
        self.byteorder = copy(db.byteorder)
        self._set_dtypes()

        self.version = copy(db.version)
        self.name = copy(db.name)
        self.label = copy(db.label)
        self.comments = db.comments.copy()
        self.ions = db.ions.copy()
        self.data = db.data.copy()
        self.fielddata = db.fielddata.copy()

        self.nstar = db.nstar
        self.nabu = db.nabu
        self.nfield = db.nfield

        self.fieldnames = db.fieldnames.copy()
        self.fieldunits = db.fieldunits.copy()
        self.fieldtypes = db.fieldtypes.copy()
        self.fieldformats = db.fieldformats.copy()
        self.fieldflags = db.fieldflags.copy()

        self.abundance_type = copy(db.abundance_type)
        self.abundance_class = copy(db.abundance_class)
        self.abundance_unit = copy(db.abundance_unit)
        self.abundance_total = copy(db.abundance_total)
        self.abundance_norm = copy(db.abundance_norm)
        self.abundance_data = copy(db.abundance_data)
        self.abundance_sum = copy(db.abundance_sum)

        self.nvalues = db.nvalues.copy()
        self.values = db.values.copy()
        self.indices = db.indices.copy()

        self._from_db_other(db)

    def _from_db_other(self, db):
        """
        Dummy routine to link in data copying by derived DB classes
        """
        pass

    def _from_data(self, **kwargs):
        """
        Create DB from provided data
        """
        silent = kwargs.get("silent", False)
        self.setup_logger(silent=silent)

        self.byteorder = kwargs.get("byteorder", self.defaultbyteorder)
        self._set_dtypes()

        self.version = self.current_version

        self.name = kwargs.get("name", None)
        self.label = kwargs.get("label", None)
        self.comments = kwargs.get("comments", tuple())
        if isinstance(self.comments, str):
            self.comments = (self.comments,)
        self.comments = np.array(self.comments, dtype=object)
        self.comments = np.append(self.comments, f"UUID: {UUID1()}")

        self.ions = kwargs.get("ions", None)
        self.data = kwargs.get("data", None)
        if isinstance(self.data, AbuData):
            if self.ions is None:
                self.ions = self.data.ions.copy()
            self.data = self.data.data.copy()

        self.fielddata = kwargs.get("fielddata", None)

        assert self.data.ndim == 2
        self.nstar = self.data.shape[0]
        self.nabu = self.data.shape[1]
        self.nfield = len(self.fielddata.dtype)

        if not isinstance(self.ions, IonList):
            self.ions = IonList(self.ions)

        assert self.nstar == self.fielddata.shape[0]
        assert self.nabu == len(self.ions)

        self.fieldnames = kwargs.get("fieldnames", None)
        self.fieldunits = kwargs.get("fieldunits", None)
        self.fieldtypes = kwargs.get("fieldtypes", None)
        self.fieldformats = kwargs.get("fieldformats", None)
        self.fieldflags = kwargs.get("fieldflags", None)

        if self.fieldnames is None:
            self.fieldnames = np.array(self.fielddata.dtype.names, dtype=object)
        else:
            self.fieldnames = np.array(self.fieldnames, dtype=object)
        if self.fieldunits is None:
            self.fieldunits = np.array([""] * self.nfield, dtype=object)
        else:
            self.fieldunits = np.array(self.fieldunits, dtype=object)
        if self.fieldflags is None:
            self.fieldflags = np.array([0] * self.nfield, dtype=np.uint64)
        else:
            self.fieldflags = np.array(self.fieldflags, dtype=np.uint64)
        if self.fieldtypes is None:
            fieldtypes = []
            for i in range(self.nfield):
                t = self.fielddata.dtype[i].type
                fieldtypes += [self._dtypes[t]]
            self.fieldtypes = np.array(fieldtypes, dtype=np.uint64)
        else:
            self.fieldtypes = np.array(self.fieldtypes, dtype=np.uint64)
        if self.fieldformats is None:
            self.fieldformats = np.array(["8.2G"] * self.nfield, dtype=object)
        else:
            self.fieldformats = np.array(self.fieldformats, dtype=object)

        assert self.nfield == len(self.fieldnames)
        assert self.nfield == len(self.fieldunits)
        assert self.nfield == len(self.fieldtypes)
        assert self.nfield == len(self.fieldformats)
        assert self.nfield == len(self.fieldflags)

        self.abundance_type = kwargs.get("abundance_type", self.AbundanceType.element)
        self.abundance_class = kwargs.get("abundance_class", self.AbundanceClass.dec)
        self.abundance_unit = kwargs.get(
            "abundance_unit", self.AbundanceUnit.mol_fraction
        )
        self.abundance_total = kwargs.get("abundance_total", self.AbundanceTotal.ejecta)
        self.abundance_norm = kwargs.get("abundance_norm", None)
        self.abundance_data = kwargs.get(
            "abundance_data", self.AbundanceData.all_ejecta
        )
        self.abundance_sum = kwargs.get(
            "abundance_sum", self.AbundanceSum.number_fraction
        )

        if self.abundance_type is None:
            if self.ions[0].is_isomer:
                self.abundance_type = self.AbundanceType.isomer
            elif self.ions[0].is_isotope:
                self.abundance_type = self.AbundanceType.isotope
            elif self.ions[0].is_element:
                self.abundance_type = self.AbundanceType.element
            elif self.ions[0].is_isobar:
                self.abundance_type = self.AbundanceType.isobar
            elif self.ions[0].is_isotone:
                self.abundance_type = self.AbundanceType.isotone
            else:
                raise Exception("unknown isotope type")

        if self.name is None:
            self.name = ""
        if self.label is None:
            self.label = ""
        if self.abundance_norm is None:
            self.abundance_norm = ""

        # compute and add hash
        sha1_comment = kwargs.get("sha1_comment", False)
        if sha1_comment:
            m = hashlib.sha1()
            m.update(self.data.tobytes())
            m.update(self.fielddata.tobytes())
            self.comments = self.comments.append(f"SHA1: {m.hexdigest()}")

        self.compute_fielddata()

        self.print_info()
        self.close_logger(timing="DB generated from data in")

    def write(
        self,
        filename=None,
        silent=False,
        byteorder=None,
    ):
        """
        Write data to file.
        """
        self.setup_logger(silent=silent)
        saved_byteorder = self.byteorder
        if byteorder is not None:
            self.byteorder = byteorder
        self.swapbyteorder = self.byteorder != self.native_byteorder
        self._open_file(filename, mode="write")
        self._write()
        self._close()
        self.byteorder = saved_byteorder
        self.swapbyteorder = self.byteorder != self.native_byteorder
        self.close_logger(timing=f"File {filename!s} written in")

    def _from_file(self, filename=None, silent=False, **kwargs):
        """
        Create DB by loading from file.
        """

        self.setup_logger(silent=silent)
        if (s := str(filename)).startswith("^"):
            s = s[1:]
            if s.startswith("/"):
                s = s[1:]
            filename = Path(_db_path) / s
        self._open_file(filename)
        self.logger.info(f"Loading {self.filename:s}")
        s = f"File size: {byte2human(self.filesize)}"
        if self.compressed:
            s += " (compressed)"
        self.logger.info(s + ".")

        self.byteorder = self._check_signature()
        self._set_dtypes()

        self._load()

        self._close()
        self.close_logger(timing="Data loaded in")

    def _open_file(self, filename, mode="read"):
        """
        Open the file.

        TODO - automatic compression mode if compression == None
        """
        compressions = ["", ".gz", ".xz", ".bz2"]

        self.filename = os.path.expandvars(os.path.expanduser(filename))
        if mode == "read":
            for c in compressions:
                fn = self.filename + c
                if os.path.exists(fn):
                    self.filename = fn
                    break
                else:
                    raise IOError("File not found.")
            access = "rb"
        else:
            access = "wb"
        if self.filename.endswith(".gz"):
            self.compressed = True
            self.compress_mode = "gz"
            self.file = gzip.open(self.filename, access)
            if mode == "read":
                pos = self.file.myfileobj.tell()
                self.file.myfileobj.seek(-4, os.SEEK_END)
                self.filesize = np.ndarray(
                    1, dtype="<u4", buffer=self.file.myfileobj.read(4)
                )[0]
                self.file.myfileobj.seek(pos, os.SEEK_SET)
        elif self.filename.endswith(".bz2"):
            self.compressed = True
            self.compress_mode = "bz2"
            self.file = bz2.BZ2File(self.filename, access, 2**16)
            if mode == "read":
                pos = self.file.tell()
                self.file.seek(0, os.SEEK_END)
                self.filesize = self.file.tell()
                self.file.seek(pos, os.SEEK_SET)
        elif self.filename.endswith(".xz"):
            self.compressed = True
            self.compress_mode = "xz"
            if mode == "read":
                self.filesize = xz_file_size(self.filename)
            self.file = lzma.LZMAFile(self.filename, access)
        else:
            self.file = open(self.filename, access, -1)
            self.stat = os.fstat(self.file.fileno())
            if mode == "read":
                self.filesize = self.stat.st_size
            self.compressed = False
            self.compress_mode = ""

    def _check_signature(self):
        """
        Check file signature and return byte order.
        """
        le_byte_order = "<"
        be_byte_order = ">"
        signature_dtype = np.dtype(le_byte_order + "u8")
        self.signature = np.array(0xADBADBADBABDADBA, dtype=signature_dtype)
        u8 = np.ndarray((), dtype=signature_dtype)
        u8.data.cast("B")[0:8] = self.file.read(8)
        if self.signature == u8:
            return le_byte_order
        u8.byteswap(True)
        if self.signature == u8:
            return be_byte_order
        raise self.SignatureError(self.filename)

    def _write_signature(self):
        """
        Write signature to file.
        """
        # signature_dtype = np.dtype(self.byteorder + 'u8')
        signature_dtype = np.dtype("u8")
        self.signature = np.array(0xADBADBADBABDADBA, dtype=signature_dtype)
        signature = self.signature.copy()
        if self.swapbyteorder:
            signature.byteswap(True)
        self.file.write(signature.data)
        self.filesize += self.signature.data.nbytes

    def _read_uin(self, dim=()):
        """
        Read numpy uint64 array.
        """
        value = np.ndarray(
            dim,
            buffer=self.file.read(8 * int(prod(dim))),
            dtype=self.dtype_u8,
            order="F",
        ).copy()
        if self.swapbyteorder:
            value.byteswap(True)
        return value

    def _write_uin(self, data):
        """
        Write numpy uint64 array.
        """
        if not isinstance(data, np.ndarray):
            data = np.array(data, np.uint64)
        dim = data.shape
        value = np.ndarray(dim, dtype=self.dtype_u8, order="F")
        value[()] = data[()]
        if self.swapbyteorder:
            value.byteswap(True)
        self.file.write(value.data)
        self.filesize += value.data.nbytes

    def _read_int(self, dim=()):
        """
        Read numpy int64 array.
        """
        value = np.ndarray(
            dim,
            buffer=self.file.read(8 * int(prod(dim))),
            dtype=self.dtype_i8,
            order="F",
        )
        if self.swapbyteorder:
            value.byteswap(True)
        return value

    def _write_int(self, data):
        """
        Write numpy int64 array.
        """
        if not isinstance(data, np.ndarray):
            data = np.array(data, np.int64)
        dim = data.shape
        value = np.ndarray(dim, dtype=self.dtype_i8, order="F")
        value[()] = data[()]
        if self.swapbyteorder:
            value.byteswap(True)
        self.file.write(value.data)
        self.filesize += value.data.bytes

    def _read_dbl(self, dim=()):
        """
        Read numpy float64 array.
        """
        value = np.ndarray(
            dim,
            buffer=self.file.read(8 * int(prod(dim))),
            dtype=self.dtype_f8,
            order="F",
        ).copy()
        if self.swapbyteorder:
            value.byteswap(True)
        return value

    def _write_dbl(self, data):
        """
        Write numpy float64 array.
        """
        if not isinstance(data, np.ndarray):
            data = np.array(data, np.float64)
        dim = data.shape
        value = np.ndarray(
            dim, buffer=np.asfortranarray(data).data, dtype=self.dtype_f8, order="C"
        ).copy()
        if self.swapbyteorder:
            value.byteswap(True)
        self.file.write(value.data)
        self.filesize += value.data.nbytes

    def _read_str(self, dim=()):
        """
        Read numpy string array.

        Unfortunately, this is somewhat inefficient in numpy as the
        maximum string size is allocated for all array entries.
        In Python 3.3+, issues with Unicode arise.
        We now just use str object for strings
        """

        n_elements = int(prod(dim))
        strlen = np.ndarray(
            dim, buffer=self.file.read(8 * n_elements), dtype=self.dtype_u8, order="F"
        ).copy()
        if self.swapbyteorder:
            strlen.byteswap(True)
        loadlen = strlen + 7 - np.mod(strlen + 7, 8)
        maxlen = int(strlen.max())
        if maxlen == 0:
            value = np.ndarray(dim, dtype=object)
            value.reshape(-1)[:] = ""
            return value
        buf = self.file.read(int(loadlen.sum()))
        value = np.ndarray(dim, dtype=object)
        # need check byteswap?
        first = np.ndarray(n_elements, dtype=np.uint64)
        first[0] = 0
        first[1:] = loadlen.reshape(-1).cumsum()[:-1]
        last = first + strlen.reshape(-1)
        flat_value = value.reshape(-1)
        for i, ifirst, ilast in zip(range(n_elements), first, last):
            flat_value[i] = buf[ifirst:ilast].decode()
        return value

    def _read_str_arr(self, ndim=1, return_dim=False):
        """
        Read string array, preceded by dimensions
        """
        dim = self._read_uin(ndim)
        value = self._read_str(dim)

        retval = (value,)
        if return_dim:
            retval = retval + (dim,)
        if len(retval) == 1:
            retval = retval[0]
        return retval

    def _write_str(self, data):
        """
        Write numpy string array.

        Unfortunately, this is somewhat inefficient in numpy as the
        maximum string size is allocated for all array entries.
        In Python 3.3+, issues with Unicode arise.
        We now just use str object for strings.
        Data is written as flat data array in Fortran layout
        """

        value = np.asfortranarray(data).reshape(-1)
        n_elements = len(value)
        strlen = np.ndarray(n_elements, dtype=self.dtype_u8, order="F")
        for i, s in enumerate(value):
            strlen[i] = len(s)
        loadlen = strlen + 7 - np.mod(strlen + 7, 8)
        first = np.ndarray(n_elements, dtype=np.uint64)
        first[0] = 0
        first[1:] = loadlen.cumsum()[:-1]
        last = first + strlen

        buf = np.ndarray((loadlen.sum(),), dtype=np.uint8)
        # use space s filler for backward compatibility
        buf[:] = ord(" ".encode())
        for i, ifirst, ilast in zip(range(n_elements), first, last):
            buf.data.cast("B")[ifirst:ilast] = value[i].encode()
        if self.swapbyteorder:
            strlen.byteswap(True)
        self.file.write(strlen.data)
        self.filesize += strlen.data.nbytes
        self.file.write(buf.data)
        self.filesize += buf.data.nbytes

    def _write_str_arr(self, data):
        """
        Write string array, preceded by dimensions
        """
        self._write_uin(data.shape)
        self._write_str(data)

    def _read_stu(self, dim=(), fieldnames=None, fieldtypes=None):
        """
        Read numpy record array.

        Currently only supports 8 byte types float64, int64, uint64.
        """
        nfields = len(fieldnames)
        dtypes = np.choose(np.array(fieldtypes, dtype=np.int64), tuple(self.dtypes))
        dtype = np.dtype({"names": fieldnames, "formats": dtypes})
        value = np.ndarray(
            dim,
            buffer=self.file.read(8 * int(nfields) * int(prod(dim))),
            dtype=dtype,
            order="F",
        ).copy()
        if self.swapbyteorder:
            value.byteswap(True)
        return value

    def _write_stu(
        self,
        data,
    ):
        """
        Write numpy record array.

        Currently only supports 8 byte types float64, int64, uint64.
        """

        for item in data.flatten()[0]:
            assert item.itemsize == 8, "trying to write data of wrong size"
        value = np.asfortranarray(data)
        if self.swapbyteorder:
            value.byteswap(True)
        self.file.write(value.data)
        self.filesize += value.data.nbytes

    def _close(self):
        """Close the file."""
        self.file.close()

    def _read_tail(self):
        """
        read file integrity check.
        """
        self.savedsize = self._read_uin()
        if self.filesize != self.savedsize:
            self.logger.error("file integrity seems broken")
            raise self.IntegrityError()
        self.logger.info("file integrity seems OK")

    def _write_tail(self):
        """
        write file integrity check.
        """

        self._write_uin(self.filesize + 8)

    def _load(self):
        """
        Load the data from file.
        """
        self.version = int(self._read_uin())
        self.name = str(self._read_str())
        if self.version < 10200:
            self.label = ""
        else:
            self.label = str(self._read_str())
        self.ncomment = int(self._read_uin())
        self.comments = self._read_str(self.ncomment)

        self.nstar = int(self._read_uin())
        self.nfield = int(self._read_uin())
        self.nabu = int(self._read_uin())

        if self.version < 10100:
            iabutype = self._read_uin()
            if iabutype != 1:
                self.logger.error("currently only supporting element data (type 1)")
                raise self.VersionError()
            self.abundance_type = 2
            self.abundance_class = 2
            self.abundance_unit = 7
            self.abundance_total = 0
            self.abundance_norm = "Lod03"
            self.abundance_data = 0
            self.abundance_sum = 1
        else:
            self.abundance_type = int(self._read_uin())
            self.abundance_class = int(self._read_uin())
            self.abundance_unit = int(self._read_uin())
            self.abundance_total = int(self._read_uin())
            self.abundance_norm = self._read_str()
            self.abundance_data = int(self._read_uin())
            self.abundance_sum = int(self._read_uin())

        self.fieldnames = self._read_str(self.nfield)
        self.fieldunits = self._read_str(self.nfield)
        self.fieldtypes = self._read_uin(self.nfield)
        fieldformats = self._read_str(self.nfield)
        self.fieldflags = self._read_uin(self.nfield)

        for i in range(self.nfield):
            self.fieldnames[i] = self.fieldnames[i].strip()
            self.fieldunits[i] = self.fieldunits[i].strip()
            fieldformats[i] = fieldformats[i].strip()

        # convert IDL to python formats in version <= 10100
        # will use python formats only in later versions
        self.fieldformats = fieldformats.copy()
        for i in range(self.nfield):
            self.fieldformats[i] = self._idlformat2pyformat(fieldformats[i])

        # l1 = max(len(x) for x in self.fieldnames)
        # l2 = max(len(x) for x in self.fieldunits)
        # l3 = max(len(x) for x in self.type_names)

        abu_Z = self._read_uin(self.nabu)
        if self.version < 10100:
            abu_A = np.zeros(self.nabu, dtype=np.uint64)
            abu_E = np.zeros(self.nabu, dtype=np.uint64)
        else:
            abu_A = self._read_uin(self.nabu)
            abu_E = self._read_uin(self.nabu)

        if self.abundance_type == 0:
            self.ions = np.array(
                [
                    I(Z=int(abu_Z[i]), A=int(abu_A[i]), E=int(abu_E[i]))
                    for i in range(self.nabu)
                ]
            )
        elif self.abundance_type == 1:
            self.ions = np.array(
                [I(Z=int(abu_Z[i]), A=int(abu_A[i])) for i in range(self.nabu)]
            )
        elif self.abundance_type == 2:
            self.ions = np.array([I(Z=int(abu_Z[i])) for i in range(self.nabu)])
        elif self.abundance_type == 3:
            self.ions = np.array([I(A=int(abu_A[i])) for i in range(self.nabu)])
        elif self.abundance_type == 4:
            self.ions = np.array([I(N=int(abu_A[i])) for i in range(self.nabu)])
        else:
            self.logger.error("anundance type not defined.")
            raise self.DataError()

        # load field data
        self.fielddata = self._read_stu(
            self.nstar,
            self.fieldnames,
            self.fieldtypes,
        )

        # load actual data
        self.data = self._read_dbl((self.nabu, self.nstar))
        self.data = np.ascontiguousarray(self.data.transpose())

        if len(np.nonzero(self.data.sum(1) == 0)[0]) > 0:
            self.logger.error("found zero data sets.")
            raise self.DataError()

        self._read_other()

        self._read_tail()

        self.compute_fielddata()
        self.print_info()

    def _read_other(self):
        """
        Dummy routine to link in data writing by derived DB classes
        """
        pass

    def compute_fielddata(self):
        """
        compute field reference values
        """

        nvalues = np.zeros(self.nfield, dtype=np.uint64)
        values = np.ndarray((self.nfield, self.nstar), dtype=np.float64)
        ivalues = np.ndarray((self.nfield, self.nstar), dtype=np.uint64)
        values[:] = np.nan
        for ifield in range(self.nfield):
            # values
            v = self.fielddata[self.fieldnames[ifield]]
            vv, vx = np.unique(v, return_inverse=True)
            nv = len(vv)
            values[ifield, 0:nv] = vv
            nvalues[ifield] = nv
            ivalues[ifield] = vx

        nvalues_max = max(nvalues)
        values = values[:, 0:nvalues_max]

        self.nvalues = nvalues
        self.values = values
        self.indices = ivalues

    def print_info(self):
        """
        Print out DB info
        """
        self.setup_logger()
        self.logger.info(f"Data version: {int(self.version):6d}")
        self.logger.info(f"data set name: {self.name:s}")
        self.logger.info(f"data set label: {self.label:s}")

        self.logger.info("".ljust(58, "="))
        for comment in self.comments:
            self.logger.info(f"COMMENT: {comment:s}")
        self.logger.info("".ljust(58, "="))

        self.logger.info(f"data sets:      {int(self.nstar):6d}")
        self.logger.info(f"abundance sets: {int(self.nabu):6d}")
        self.logger.info("".ljust(58, "-"))
        self.logger.info(
            "abundance type:  {:1d} - {:s}".format(
                int(self.abundance_type), self.abundance_type_names[self.abundance_type]
            )
        )
        self.logger.info(
            "abundance class: {:1d} - {:s}".format(
                int(self.abundance_class),
                self.abundance_class_names[self.abundance_class],
            )
        )
        self.logger.info(
            "abundance unit:  {:1d} - {:s}".format(
                int(self.abundance_unit), self.abundance_unit_names[self.abundance_unit]
            )
        )
        self.logger.info(
            "abundance total: {:1d} - {:s}".format(
                int(self.abundance_total),
                self.abundance_total_names[self.abundance_total],
            )
        )
        s = self.abundance_norm
        if s == "":
            s = "(NONE)"
        self.logger.info(f"abundance norm:      {s}")
        self.logger.info(
            "abundance data:  {:1d} - {:s}".format(
                int(self.abundance_data), self.abundance_data_names[self.abundance_data]
            )
        )
        self.logger.info(
            "abundance sum:   {:1d} - {:s}".format(
                int(self.abundance_sum), self.abundance_sum_names[self.abundance_sum]
            )
        )

        self.logger.info("".ljust(58, "-"))
        self.logger.info(f"{self.nfield:d} data fields: ")

        l1 = max(len(x) for x in self.fieldnames)
        l2 = 0
        l3 = 0
        l4 = max(len(x) for x in self.fieldunits)
        l5 = 0
        re_len = re.compile("^([0-9]+)", flags=re.I)
        for ifield in range(self.nfield):
            # output formatting
            l5 = max(l5, len(self.type_names[self.fieldtypes[ifield]]))
            flen = int(re_len.findall(self.fieldformats[ifield])[0])
            l2 = max(l2, flen)
            l3 = max(l3, len(f"{self.nvalues[ifield]:d}"))

        format = f"{{:{l1:d}s}} [{{:>{l4:d}s}}] ({{:{l5:d}s}}) <{{:s}}>"
        for ifield in range(self.nfield):
            self.logger.info(
                format.format(
                    self.fieldnames[ifield],
                    self.fieldunits[ifield],
                    self.type_names[self.fieldtypes[ifield]],
                    self.flag_names[self.fieldflags[ifield]],
                )
            )
            if (
                len(
                    np.argwhere(
                        self.type_names[self.fieldtypes[ifield]]
                        == np.array(["DOUBLE", "LONG64", "ULONG64"])
                    )
                )
                == 0
            ):
                self.logger.error("data type not yet supported")
                self.logger.error("only supporting 8-byte scalar data types")
                raise self.VersionError()

        self.logger.info("".ljust(58, "-"))
        self.logger.info("ABUNDANCES:")
        s = " ".join(str(i) for i in self.ions)
        s = textwrap.wrap(s, 50)
        for line in s:
            self.logger.info(line)

        # These two b;ock below should become their own function
        xpar = np.argwhere(self.fieldflags == 0)
        if len(xpar) > 0:
            self.logger.info("".ljust(58, "-"))
            self.logger.info("PARAMETER RANGES:")
        for ip in xpar.flat:
            fmax = max(self.values[ip, 0 : self.nvalues[ip]])
            fmin = min(self.values[ip, 0 : self.nvalues[ip]])
            line = (
                self.fieldnames[ip] + ": ",
                ("{:" + self.fieldformats[ip] + "}").format(fmin),
                ("{:" + self.fieldformats[ip] + "}").format(fmax),
                f"{int(self.nvalues[ip]):d}",
            )
            format = (
                "{{:<{:d}s}} {{:>{:d}s}} ... {{:>{:d}s}} ({{:>{:d}s}} values)".format(
                    l1 + 2, l2, l2, l3
                )
            )
            self.logger.info(format.format(*line))

        xprop = np.argwhere(self.fieldflags != 0)
        if len(xprop) > 0:
            self.logger.info("".ljust(58, "-"))
            self.logger.info("PROPERTY RANGES:")
        for ip in xprop.flat:
            fmax = max(self.values[ip, 0 : self.nvalues[ip]])
            fmin = min(self.values[ip, 0 : self.nvalues[ip]])
            line = (
                self.fieldnames[ip] + ": ",
                ("{:" + self.fieldformats[ip] + "}").format(fmin),
                ("{:" + self.fieldformats[ip] + "}").format(fmax),
                f"{int(self.nvalues[ip]):d}",
            )
            format = (
                "{{:<{:d}s}} {{:>{:d}s}} ... {{:>{:d}s}} ({{:>{:d}s}} values)".format(
                    l1 + 2, l2, l2, l3
                )
            )
            self.logger.info(format.format(*line))

        if len(xpar) > 0:
            self.logger.info("".ljust(58, "-"))
            self.logger.info("PARAMETER VALUES:")
        for ip in xpar.flat:
            self.logger.info(self.fieldnames[ip] + ":")
            flen = int(re_len.findall(self.fieldformats[ifield])[0])
            s = ""
            f = " {:" + self.fieldformats[ip] + "}"
            for id in range(self.nvalues[ip]):
                if len(s) >= 50:
                    self.logger.info(s)
                    s = ""
                s += f.format(self.values[ip, id])
            self.logger.info(s)

        maxpropvalues = 100
        if len(xprop) > 0:
            self.logger.info("".ljust(58, "-"))
            self.logger.info("PROPERTY VALUES:")
        for ip in xprop.flat:
            self.logger.info(self.fieldnames[ip] + ":")
            if self.nvalues[ip] > maxpropvalues:
                self.logger.info(f"(more than {maxpropvalues:d} values)")
            else:
                flen = int(re_len.findall(self.fieldformats[ifield])[0])
                s = ""
                f = " {:" + self.fieldformats[ip] + "}"
                for id in range(self.nvalues[ip]):
                    if len(s) >= 50:
                        self.logger.info(s)
                        s = ""
                    s += f.format(self.values[ip, id])
                self.logger.info(s)
        self.logger.info("".ljust(58, "-"))

        # compute hash
        m = hashlib.sha1()
        m.update(self.data.tobytes())
        m.update(self.fielddata.tobytes())

        self.logger.info(f"SHA1: {m.hexdigest()}")
        self.logger.info("".ljust(58, "-"))
        self.close_logger()

    def _write(self):
        """
        Write data to file.
        """

        self.version = self.current_version

        self.filesize = 0
        self._write_signature()
        self._write_uin(self.version)
        self._write_str(self.name)
        self._write_str(self.label)
        self._write_str_arr(self.comments)
        self._write_uin(self.nstar)
        self._write_uin(self.nfield)
        self._write_uin(self.nabu)

        self._write_uin(self.abundance_type)
        self._write_uin(self.abundance_class)
        self._write_uin(self.abundance_unit)
        self._write_uin(self.abundance_total)
        self._write_str(self.abundance_norm)
        self._write_uin(self.abundance_data)
        self._write_uin(self.abundance_sum)

        fieldformats = self.fieldformats.copy()
        for i in range(self.nfield):
            fieldformats[i] = self._pyformat2idlformat(self.fieldformats[i])

        self._write_str(self.fieldnames)
        self._write_str(self.fieldunits)
        self._write_uin(self.fieldtypes)
        self._write_str(fieldformats)
        self._write_uin(self.fieldflags)

        # set A,Z,E from ions
        abu_Z = np.array([ion.Z for ion in self.ions], dtype=np.uint64)
        abu_E = np.array([ion.E for ion in self.ions], dtype=np.uint64)
        if self.abundance_type in (0, 1, 2, 3):
            abu_A = np.array([ion.A for ion in self.ions], dtype=np.uint64)
        elif self.abundance_type == 4:
            abu_A = np.array([ion.N for ion in self.ions], dtype=np.uint64)
        else:
            self.logger.error("anundance type not defined.")
            raise self.DataError()

        self._write_uin(abu_Z)
        self._write_uin(abu_A)
        self._write_uin(abu_E)

        self._write_stu(self.fielddata)

        self._write_dbl(self.data.transpose())

        self._write_other()

        self._write_tail()

    def _write_other(self):
        """
        Dummy routine to link in data reading by derived DB classes
        """
        pass

    def _compute_index(self, array_limit=2**24):
        nindex = np.ndarray((self.nfield, self.nstar), dtype=np.uint64)
        vindex = np.ndarray((self.nfield,), dtype=object)
        imax = np.ndarray((self.nfield,), dtype=np.uint64)
        for ifield in range(self.nfield):
            vindex[ifield] = scipy.sparse.lil_matrix(
                (self.nstar, self.nstar), dtype=np.uint64
            )
            nv = self.nvalues[ifield]
            for iv in range(nv):
                (ii,) = np.where(self.indices[ifield] == iv)
                ni = len(ii)
                nindex[ifield, iv] = ni
                vindex[ifield][iv, :ni] = ii
            nimax = np.max(nindex[ifield, :nv])
            imax[ifield] = nimax
            vindex[ifield] = vindex[ifield][0:nv, 0:nimax]
        nvmax = np.max(self.nvalues)
        nimax = np.max(imax)
        if self.nfield * nvmax * nimax <= array_limit:
            v = np.ndarray((self.nfield, nvmax, nimax), dtype=np.uint64)
            for ifield in range(self.nfield):
                v[ifield, : self.nvalues[ifield], : imax[ifield]] = vindex[
                    ifield
                ].toarray()
            vindex = v
        else:
            print(self.nfield * nvmax * nimax)
        self._nvindex = nindex
        self._vindex = vindex

    @CachedAttribute
    def nvindex(self):
        """
        return number of elements for each value in value index array
        """
        try:
            self._nvindex
        except AttributeError:
            self._compute_index()
        return self._nvindex

    @CachedAttribute
    def vindex(self):
        """
        return value index array

        array of star indices for given parameter value
        """
        try:
            self._vindex
        except AttributeError:
            self._compute_index()
        return self._vindex

    @staticmethod
    def _idlformat2pyformat(format):
        """
        Convert IDL format coding into python format coding.
        """

        fmt_convert = {"I": "D"}
        fmt = format[1:] + fmt_convert.get(format[0], format[0])
        return fmt

    @staticmethod
    def _pyformat2idlformat(format):
        """
        Convert python format coding into IDL format coding.
        """

        fmt_convert = {"D": "I"}
        fmt = fmt_convert.get(format[-1], format[-1]) + format[:-1]
        return fmt

    def get_star_slice(self, *args, **kwargs):
        """
        return slice of star indices for given parameters
        """
        if len(args) == 1 and isinstance(args[0], Mapping):
            fields = args[0]
        else:
            fields = kwargs
        for key in fields:
            assert key in self.fieldnames
        vfilter = dict()
        for i, f in enumerate(self.fieldnames):
            val = fields.get(f, None)
            if val is not None:
                j = np.abs(self.values[i, : self.nvalues[i]] - fields[f]).argmin()
                vfilter[i] = int(j)
        mask = np.tile(True, self.nstar)
        for i, j in vfilter.items():
            imask = self.fielddata[self.fieldnames[i]] == self.values[i, j]
            mask = np.logical_and(mask, imask)
        (indices,) = np.where(mask)
        return indices

    # ------------------------------------------------------
    # data access methods
    # ------------------------------------------------------

    def index_info(self, indices):
        """
        Print field information for list of indices.
        """
        if not isinstance(indices, Iterable):
            indices = (indices,)
        maxlenname = max(len(n) for n in self.fieldnames)
        maxlenvalue = max(int(f[:-1].split(".")[0]) for f in self.fieldformats)
        maxlenunit = max(len(n) for n in self.fieldunits)
        maxlen = maxlenname + maxlenvalue + maxlenunit + 3
        print("-" * maxlen)
        for i in indices:
            print(
                "{name:>{maxlenname}s}: {index:>{maxlenvalue}d}".format(
                    name="index",
                    index=i,
                    maxlenvalue=maxlenvalue,
                    maxlenname=maxlenname,
                )
            )
            for f, n, u in zip(self.fieldformats, self.fieldnames, self.fieldunits):
                print(
                    "{name:>{maxlenname:d}s}: {value:>{maxlenvalue}s} {unit:s}".format(
                        name=n,
                        unit=u,
                        value="{value:{format:s}}".format(
                            value=self.fielddata[i][n], format=f
                        ),
                        maxlenname=maxlenname,
                        maxlenvalue=maxlenvalue,
                    )
                )

            print("-" * maxlen)

    def abundance_info(self, indices):
        """
        Print out abundance info for givien index
        """
        if not isinstance(indices, Iterable):
            indices = (indices,)
        print("=" * (9 + 13 * len(indices)))
        print(("{:>8s}" + " ".join(["{:>12d}"] * len(indices))).format("Ion", *indices))
        print("=" * (9 + 13 * len(indices)))
        for ion, abu in zip(self.ions, self.data[indices, :].transpose()):
            print(("{!s:>8}" + " ".join(["{:12.5e}"] * len(indices))).format(ion, *abu))
        print("=" * (9 + 13 * len(indices)))


load = loader(StarDB, __name__ + ".load")
_load = _loader(StarDB, __name__ + ".load")
loadstardb = load


def test():
    import filecmp

    f0 = os.path.expanduser("~/kepler/znuc/test.stardb")
    f1 = os.path.expanduser("~/test.dat")
    d = StarDB(f0)
    print("=" * 40)
    d.write(f1)
    print("=" * 40)
    x = StarDB(f1)
    print(filecmp.cmp(f0, f1))
    print("=" * 40)
    d = StarDB(db=x)
    d.write(f1)
    x = StarDB(f1)
    print(filecmp.cmp(f0, f1))
