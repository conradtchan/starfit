#!/usr/bin/python3

"""Reads elemental abundance file into an array
File format:

========================================
Version 100200
========================================

Add tetection threshold.  Otherwise identical to and backward
compatible with Format 10100

10200 <--- VERSION NUMBER
SDSS12345-6789 <-- STAR NAME
ApJ 222, 333 (1999) <-- SOURCE/REFERENCE
data automatically generated by Anna Frebel <-- SOME COMMENT
1 1 <--- DATA FORMAT, MODE (DETECTION THRESHOLD PRESENT = 1)
    <--- SYMBOL OF ELEMENT OR VALUE TO WHICH DATA IS NORMALIZED, IF ANY
17  <--- NUMBER OF ELEMENTS
C   -1.23  +0.07  -1.07  +0.05  -0.05
N    1.23  +0.14  -2.14  -0.10  +0.10
O   -1.23  +0.21  -3.21  -0.15  -0.15  <--- MEASUREMNT,
Fe   1.23  +0.28  -4.28  +0.20  +0.20       UNCORRELATED ERROR,
Mg  -1.23  +0.35  -5.35  +0.25  -0.25       DETECTION THRESHOLD,
Si   1.23  +0.42  -6.42  +0.30  +0.30       COVARIANCE VECTORS
Sc   1.23  +0.57  -7.57  -0.35  -0.35
Ca  -1.23  +0.49  -8.49
Ti  -1.23  +0.64  -9.64
V   +1.23  +0.64  -8.64  <--- MEASUREMNT,
Cr   1.23  +0.10  -6.64       UNCORRELATED ERROR,
Mn  -1.23  +0.20  -6.64       DETECTION THRESHOLD
Co  -1.23  -0.30
Ni  -1.23  -0.40
Cu   1.23  -0.50
Sr  -1.23  -0.60 <--- UPPER LIMIT,
Zr   1.23  -0.70      (UNCORRELATED) UNCERTAINTIES
Ba  -1.23  -0.80
Lo09	<--- SOLAR REFERENCE; USE "-" IF NOT PROVIDED

DATA FORMAT
1 log epsilon (PREFERRED DUE TO DIFFERENT SOLAR ABUNDANCES, give H mol/g as norm!)
2 [ ]
3 [X/Y] Y= norm element, provide [Y/H] in column for Y
4 log X/Si + 6 (by number), norm is assumeed log(Si) + 6 (mol/g)
5 log mole fraction (REALLY PREFERRED)
6 [X/H] - provide [H] as norm, otherwise BBN value is assumed
7 X/Si * 1e6 (by number), norm is Si/Si_sun, error is absolute

The first number is the data, the second is the uncorrelated
(random)error in dex.  Use negative random error values to indicate
upper limits.

The 3rd column is the detection threshold, if present.  May be "-" if
missing.

The following columns give first correlated errors.

The covariances can have positive or negative signs, the
uncorrelated (random) error always is inherently non-negative.

The square of the total error is the sum of the squares of all errors.

If no covariances are provide, the entire error is random.

As of the time of this writing, StarFit does not use correlated
upper limits, but that may change in the future.


========================================
Version 10100
========================================

Adds covariant errors, format compatible with Version 10002

10100 <--- VERSION NUMBER
SDSS12345-6789 <-- STAR NAME
ApJ 222, 333 (1999) <-- SOURCE/REFERENCE
data automatically generated by Anna Frebel <-- SOME COMMENT
1 <--- DATA FORMAT (RECOMEND USE LOG EPSILON OR LOG MOL FRACTION)
    <--- NAME OF ELEMENT OR VALUE TO WHICH DATA IS NORMALIZED, IF ANY
17  <--- NUMBER OF ELEMENTS
C   -1.23  +0.07  +0.05  -0.05
N    1.23  +0.14  -0.10  +0.10 <--- MEASUREMNT,
O   -1.23  +0.21  -0.15  -0.15      UNCORRELATED ERROR,
Fe   1.23  +0.28  +0.20  +0.20      COVARIANCE VECTORS
Mg  -1.23  +0.35  +0.25  -0.25
Si   1.23  +0.42  +0.30  +0.30
Ca  -1.23  +0.49  -0.35  -0.35
Sc   1.23  +0.57
Ti  -1.23  +0.64
Cr   1.23  +0.10 <--- MEASUREMENT,
Mn  -1.23  +0.20      UNCORRELATED ERRORS
Ni  -1.23  -0.40
Cu   1.23  -0.50
Sr  -1.23  -0.60 <--- UPPERL LIMIT,
Zr   1.23  -0.70      (UNCORRELATED) UNCERTAINTIES
Ba  -1.23  -0.80
Lo09	<--- SOLAR REFERENCE; USE "-" IF NOT PROVIDED

DATA FORMAT
1 log epsilon (PREFERRED DUE TO DIFFERENT SOLAR ABUNDANCES, give H mol/g as norm!)
2 [ ]
3 [X/Y] Y= norm element, provide [Y/H] in column for Y
4 log X/Si + 6 (by number), norm is assumeed log(Si) + 6 (mol/g)
5 log mole fraction (REALLY PREFERRED)
6 [X/H] - provide [H] as norm, otherwise BBN value is assumed
7 X/Si * 1e6 (by number), norm is Si/Si_sun, error is absolute

The first number is the data, the second is the uncorrelated error in dex.
Use negative random error values to indicate upper limits.

The following columns give first covariant errors.

The covariant errors can have positive or negative signs, the
uncorrelated (random) error always is inherently non-negative.

The square of the total error is the sum of the squares of all errors.

If no covariances are provide, the entire error is random.

As of the time of this writing, StarFit does not use correlated
upper limits, but that may change in the future.


========================================
Version 100002
========================================

10002 <--- VERSION NUMBER (FILE FORMAT)
CS30336-049 <-- STAR NAME
ApJ 222, 333 (1999) <-- SOURCE/REFERENCE
data automatically generated by Anna Frebel <-- SOME COMMENT
1 <--- DATA FORMAT (USE LOG EPSILON OR LOG MOLE FRACTION)
  <--- NORM, IF ANY (ELEMENT OR ABUNDANCE OFFSET)
17  <--- NUMBER OF ELEMENTS
C      -1.23    0.1
N       1.23    0.2
O      -1.23    0.3
Fe      1.23    0.4 <--- MEASUREMENT,
Mg     -1.23    0.5      UNCORRELATED ERRORS
Si      1.23    0.6
Ca     -1.23    0.7
Sc      1.23    0.8
Ti     -1.23    0.9
Cr      1.23   -0.1
Mn     -1.23   -0.2
Co      1.23   -0.3
Ni     -1.23   -0.4 <--- UPPERL LIMIT,
Sr     -1.23   -0.6      (UNCORRELATED) UNCERTAINTIES
Zr      1.23   -0.7
Ba     -1.23   -0.8
Lo09  <--- SOLAR REFERENCE; USE "-" IF NOT PROVIDED

DATA FORMAT
  1 log epsilon (PREFERRED DUE TO DIFFERENT SOLAR ABUNDANCES, give H mol/g as norm!)
  2 [ ]
  3 [X/Y] Y = norm element, provide [Y/H] in column for Y
  4 log X/Si + 6 (by number), norm is assumeed log(Si) + 6 (mol/g)
  5 log mole fraction (REALLY PREFERRED)
  6 [X/H] - provide [H] as norm, otherwise BBN value is assumed
  7 X/Si * 1e6 (by number), norm is Si/Si_sun, error is absolute

The first number is the data, the second is the error in dex.
Use negative values to indicate upper limits; all need to be non-zero.


========================================
Version 100001
========================================

10001 <--- VERSION NUMBER (FILE FORMAT)
CS30336-049 <-- STAR NAME
ApJ 222, 333 (1999) <-- SOURCE/REFERENCE
data automatically generated by Anna Frebel <-- SOME COMMENT
1 <--- DATA FORMAT (USE LOG EPSILON OR LOG MOLE FRACTION)
  <--- NORM, IF ANY (ELEMENT OR ABUNDANCE OFFSET)
17  <--- NUMBER OF ELEMENTS
C      -1.23    0.1
N       1.23    0.2
O      -1.23    0.3
Fe      1.23    0.4 <--- MEASUREMENT,
Mg     -1.23    0.5      UNCORRELATED ERRORS
Si      1.23    0.6
Ca     -1.23    0.7
Sc      1.23    0.8
Ti     -1.23    0.9
Cr      1.23   -0.1
Mn     -1.23   -0.2
Co      1.23   -0.3
Ni     -1.23   -0.4 <--- UPPERL LIMIT,
Sr     -1.23   -0.6      (UNCORRELATED) UNCERTAINTIES
Zr      1.23   -0.7
Ba     -1.23   -0.8

DATA FORMAT
  1 log epsilon (PREFERRED DUE TO DIFFERENT SOLAR ABUNDANCES, give H mol/g as norm!)
  2 [ ]
  3 [X/Y] Y = norm element, provide [Y/H] in column for Y
  4 log X/Si + 6 (by number), norm is assumeed log(Si) + 6 (mol/g)
  5 log mole fraction (REALLY PREFERRED)
  6 [X/H] - provide [H] as norm, otherwise BBN value is assumed
  7 X/Si * 1e6 (by number), norm is Si/Si_sun, error is absolute

The first number is the data, the second is the error in dex.
Use negative values to indicate upper limits; all need to be non-zero.


========================================
Version 100000
========================================

10001 <--- VERSION NUMBER (FILE FORMAT)
CS30336-049, data from a friend  <-- STAR NAME AND COMMENT
1 <--- DATA FORMAT (USE LOG EPSILON OR LOG MOLE FRACTION)
  <--- NORM, IF ANY (ELEMENT OR ABUNDANCE OFFSET)
17  <--- NUMBER OF ELEMENTS
C      -1.23    0.1
N       1.23    0.2
O      -1.23    0.3
Fe      1.23    0.4 <--- MEASUREMENT,
Mg     -1.23    0.5      UNCORRELATED ERRORS
Si      1.23    0.6
Ca     -1.23    0.7
Sc      1.23    0.8
Ti     -1.23    0.9
Cr      1.23   -0.1
Mn     -1.23   -0.2
Co      1.23   -0.3
Ni     -1.23   -0.4 <--- UPPERL LIMIT,
Sr     -1.23   -0.6      (UNCORRELATED) UNCERTAINTIES
Zr      1.23   -0.7
Ba     -1.23   -0.8

DATA FORMAT
  1 log epsilon (PREFERRED DUE TO DIFFERENT SOLAR ABUNDANCES, give H mol/g as norm!)
  2 [ ]
  3 [X/Y] Y = norm element, provide [Y/H] in column for Y
  4 log X/Si + 6 (by number), norm is assumeed log(Si) + 6 (mol/g)
  5 log mole fraction (REALLY PREFERRED)
  6 [X/H] - provide [H] as norm, otherwise BBN value is assumed
  7 X/Si * 1e6 (by number), norm is Si/Si_sun, error is absolute

The first number is the data, the second is the error in dex.
Use negative values to indicate upper limits; all need to be non-zero.
"""
import numpy as np

from . import BBN, REF, SOLAR, STARS
from .autils.abusets import BBNAbu, SolAbu
from .autils.isotope import Ion
from .autils.isotope import ion as I
from .autils.logged import Logged
from .utils import find_data

# low values for detection limits
LOW = np.array([np.nan] + [-199.0] * 6 + [1.0e-199])

DATA_FORMAT = {
    1: "log epsilon",
    2: "[X]",
    3: "[X/X_ref]",
    4: "log(X/Si) + 6",
    5: "log Y",
    6: "[X/H]",
    7: "(X/Si) * 1e6",
}


class StarFileError(Exception):
    pass


# A class for each set of readings
class Star(Logged):
    """All the data read from each file"""

    def __init__(self, filename, silent=False):
        self.silent = silent
        self.setup_logger(silent=silent)

        filename = find_data(STARS, filename)

        # Load BBN data for H reference
        self.BBN_data = BBNAbu(
            name=find_data(REF, BBN),
            silent=silent,
        )
        self._read(filename)
        # Sort star
        self.element_abundances.sort(order="element")
        # Add BBN data
        self._bbn()

        self.filename = filename

    # Write list into numpy recarray
    def list_to_array(self, data_list):

        # Determine shape of data
        data = [d.split() for d in data_list]
        n_data = np.array([len(d) for d in data])
        n_elem = len(n_data)
        if self.version < 10100:
            n_det = 0
            n_cov = 0
            assert np.all(n_data == 3), "error in star data format, require 3 columns"
        elif self.version < 10200:
            n_det = 0
            n_cov = np.max(n_data) - 3
            assert n_cov >= 0, "error in star data format, require at least 3 columns"
            assert np.all(
                (n_data == 3) | (n_data == n_cov + 3)
            ), "require consistent covariance vectors, if present"
        else:
            n_det = int(self.data_mode == 1)
            n_cov = np.max(n_data) - 3 - n_det
            assert n_cov >= 0, "error in star data format, require at least 3 columns"
            assert np.all(
                (n_data == 3) | (n_data == 3 + n_det) | (n_data == n_cov + 3 + n_det)
            ), "require consistent covariance vectors, if present"

        self.data_type = [
            ("element", object),
            ("abundance", np.float64),
            ("error", np.float64),
            ("detection", np.float64),
            ("covariance", np.float64, n_cov),
        ]

        # Initialize zeroed numpy array
        data_array = np.recarray(n_elem, dtype=self.data_type)
        data_array[:] = 0

        # i lines, j elements on each line
        for i, line in enumerate(data):
            for j, item in enumerate(line):
                # Set the corresponding data type
                if j == 0:
                    data_array[i].element = I(item)
                elif j == 1:
                    data_array[i].abundance = float(item)
                elif j == 2:
                    data_array[i].error = float(item)
                    if n_det == 0:
                        data_array[i].detection = LOW[self.data_format]
                    if data_array[i].error == 0.0:
                        raise StarFileError(
                            f"Statistical error can't be zero for element {data_array[i].element!s}"
                        )
                elif j == 2 + n_det:
                    if item in ("-",) or data_array[i].error < 0.0:
                        data_array[i].detection = LOW[self.data_format]
                    else:
                        data_array[i].detection = float(item)
                else:
                    if data_array[i].error < 0.0:
                        data_array[i].covariance[j - 3 - n_det] = 0.0
                    else:
                        data_array[i].covariance[j - 3 - n_det] = float(item)
            if n_det == 1 and len(line) < 4:
                data_array[i].detection = LOW[self.data_format]
        if n_det == 0:
            data_array[:].detection = LOW[self.data_format]

        return data_array

    # Read file
    def _read(self, filename):
        self.setup_logger(silent=self.silent)
        with open(filename, "rt") as f:
            self.logger_file_info(f)
            # Initialize empty array for reading each line
            content = []
            # Initialize the count for number of elements
            n_elements = 9999
            # Number each line, and write into the array, stripping whitespace and carriage returns
            for i, line in enumerate(f):
                content += [line.strip()]
                if i == 0:
                    # Version number
                    self.version = int(content[0])
                    if self.version == 10000:
                        n_elements_line = 4
                    else:
                        n_elements_line = 6
                # Line 6 gives the number of elements in the file
                # Write this into n_elements
                if i == n_elements_line:
                    n_elements = int(line)
                # Read up to the last element, read one more line (solar reference) and then stop
                if i > n_elements_line + n_elements:
                    break
        n = 1

        # Additional lines for the appropriate version numbers
        if self.version < 10001:
            # Star name
            self.name = content[n].split()[0]
            self.source = ""
        else:
            # Star name
            self.name = content[n]
            n += 1
            # Source
            self.source = content[n]
            n += 1
        # Comment
        self.comment = content[n]
        n += 1
        # Format
        items = [int(x) for x in content[n].split()]
        assert (
            len(items) == 1 or self.version >= 10200
        ), f"invalid values {n} lines below start of file"
        self.data_format = items[0]
        if len(items) > 1:
            self.data_mode = items[1]
            assert self.data_mode == 1
        else:
            self.data_mode = 0
        n += 1
        # Normalization element
        if self.data_format == 1:
            # Data format Type 1 specifies H number fraction
            try:
                self.norm_element = float(content[n])
                if self.norm_element <= 0.0:
                    self.norm_element = 10.0**self.norm_element
            except ValueError:
                # BBN H number fraction
                self.norm_element = self.BBN_data.Y("H")
        elif self.data_format == 4:
            # Data format Type 4 specifies log(Si) + 6 in mol/g in this line instead of a normalization element
            try:
                self.norm_element = float(content[n])
            except ValueError:
                # log10(Solar(Lodders2022).Y(Si))
                self.norm_element = 1.4433183352909564
        elif self.data_format == 6:
            # Data format Type 6 specifies log(H/H_sun) in this line instead of a normalization element
            try:
                self.norm_element = float(content[n])
            except ValueError:
                # log10(BBN(Coc19).Y(H) / Solar(Lodders2022).Y(H))
                self.norm_element = 0.02813358159449738
        elif self.data_format == 7:
            # Data format type 7 specifies Si * 1e6 in mol/g in this line instead of a normalization element
            try:
                self.norm_element = float(content[n])
            except ValueError:
                # log10(Solar(Lodders2022).Y(Si))
                self.norm_element = 27.75364315444385
        else:
            try:
                norm = int(content[n])
            except:
                norm = content[n]
            norm = I(norm)
            if not norm.is_element:
                self.norm_element = None
            else:
                self.norm_element = norm

        n += 1
        # Number of elements
        self.n_elements = n_elements
        n += 1
        # Element abundances, the main set of data, converted to format 5
        self.element_abundances = self.list_to_array(content[n : n + n_elements])
        n += n_elements
        # Solar reference
        if self.version >= 10002:
            self.solar_ref = content[n]
            # Default solar ref
            if self.solar_ref in (
                "",
                "-",
            ):
                self.solar_ref = find_data(REF, SOLAR)
        else:
            self.solar_ref = None

        # Convert
        self.input_data_format = self.data_format
        self.element_abundances = self._convert5(
            self.element_abundances, self.data_format, self.norm_element
        )
        # Since the format has been converted, adjust the number accordingly
        self.data_format = 5

        self.close_logger(timing=" Star loaded and converted in")

    def _convert5(self, array, data_format, norm_element):
        """
        Convert abundance array to Format 5
        """
        if data_format == 1:
            # norm_element is H number fraction
            h_ref = self.norm_element
            # array.abundance[:] = np.log10((10 ** (array.abundance - 12)) * h_ref)
            for a in (array.abundance, array.detection):
                a[:] = a - 12 + np.log10(h_ref)
        elif data_format == 2:
            # Load the abundances for the sun, using the appropriate solar reference
            sun = SolAbu(self.solar_ref, silent=True)
            solar = sun.Y(array.element)
            for a in (array.abundance, array.detection):
                a[:] = np.log10((10 ** (a)) * solar)
        elif data_format == 3:
            # Load the Big Bang hydrogen abundance
            h_ref = self.BBN_data.Y("H")
            # Load the abundances for the sun
            sun = SolAbu(self.solar_ref, silent=True)
            solar = sun.Y(array.element)
            # Set the index for the norm element, since this will need to be treated separately
            norm_index = np.where(array.element == norm_element)[0][0]
            # Calculate the linear abundance for the norm element on the star
            norm_abundance = ((10 ** array[norm_index].abundance) * h_ref) * (
                sun.Y(norm_element) / sun.Y("H")
            )
            norm_detection = ((10 ** array[norm_index].detection) * h_ref) * (
                sun.Y(norm_element) / sun.Y("H")
            )
            # Calculate the abundance for the norm element on the sun
            norm_solar = sun.Y(norm_element)
            for a, b in zip(
                (array.abundance, array.detection),
                (norm_abundance, norm_detection),
            ):
                # Calculate abundances for all the other elements
                abundances = 10**a * (norm_abundance) * solar / norm_solar
                # Change the entry for the norm element that we calculated earlier
                abundances[norm_index] = b
                # Log, and then put back into array
                a[:] = np.log10(abundances)
        elif data_format == 4:
            # Since this is format 4, norm_element is actually log(Y(Si)) + 6
            logSi = norm_element
            for a in (array.abundance, array.detection):
                a[:] = a + logSi - 12
        elif data_format == 5:
            # No changes need to be made here
            pass
        elif data_format == 6:
            # Since this is format 6, norm_element is actually log(H/H_sun)
            logH = norm_element
            # Load the abundances for the sun, using the appropriate solar reference
            sun = SolAbu(self.solar_ref, silent=True)
            solar = sun.Y(array.element)
            # array.abundance[:] = np.log10((10 ** (array.abundance + logH)) * solar)
            for a in (array.abundance, array.detection):
                a[:] = a[:] + logH + np.log10(solar)
        elif data_format == 7:
            # sigma is not given logarthmically in dex but absolute
            array.error[:] = array.error / (array.abundance * np.log(10.0))
            array.covariance[:, :] = array.covariance[:, :] / (
                array.abundance.reshape(-1, 1) * np.log(10.0)
            )
            # Since this is format 7, norm_element is actually Y(Si) * 1e6
            Si = norm_element
            for a in (array.abundance, array.detection):
                a[:] = np.log10(a[:] * Si * 1e-12)
        else:
            raise Exception(f"Format type {data_format:d} not found")
        return array

    def get_elements(self):
        """
        return all elements in dataset
        """
        return [
            a.element.Name()
            for a in self.element_abundances
            if a.element not in self.added_elements
        ]

    def get_upper_limits(self):
        return [
            a.element.Name()
            for a in self.element_abundances
            if a.error < 0 and a.element not in self.added_elements
        ]

    def get_measured(self):
        return [
            a.element.Name()
            for a in self.element_abundances
            if a.error > 0 and a.element not in self.added_elements
        ]

    def get_detection_thresholds(self):
        return [a.element.Name() for a in self.element_abundances if a.detection > -80]

    def get_covariances(self):
        return [
            a.element.Name()
            for a in self.element_abundances
            if np.any(a.covariance != 0)
        ]

    def get_n_covariances(self):
        return self.element_abundances.covariance.shape[-1]

    def get_input_data_format(self):
        return DATA_FORMAT[self.input_data_format]

    def get_norm(self):
        if self.norm_element is None:
            return ""
        if isinstance(self.norm_element, Ion):
            return self.norm_element.Name()
        return str(self.norm_element)

    # A nice readable name
    def __str__(self):
        return f"{self.__class__.__name__:s}({self.name:s})"

    __repr__ = __str__

    @staticmethod
    def _rec_insert_1d(array, record, pos):
        """
        Provides obvious functionality missing from NumPy
        """
        assert array.ndim == 1
        assert array.dtype == record.dtype
        new = np.empty_like(array, shape=array.shape[0] + 1)
        new[:pos] = array[:pos]
        new[pos + 1 :] = array[pos:]
        new[pos] = record
        return new

    def _bbn(self):
        """Adds BBN data for H and He to the star"""
        self.setup_logger(silent=self.silent)
        self.logger.info(
            "Adding BBN upper limits of H and He into star data if missing"
        )

        self.added_elements = list()
        for i in range(2):
            element = I(Z=i + 1)
            if self.element_abundances[i].element is not element:
                entry = np.empty_like(self.element_abundances, shape=())[()]
                entry.element = element
                entry.abundance = np.log10(self.BBN_data.Y(entry["element"]))
                entry.error = -0.1  # This is an upper limit, 0.1 is somewhat arbitrary
                entry.detection = LOW[5]
                entry.covariance[:] = 0.0
                self.element_abundances = self._rec_insert_1d(
                    self.element_abundances, entry, i
                )
                self.added_elements.append(element)
        self.n_elements = self.element_abundances.shape[0]


def StarTable(
    filenames=[
        "Caffau.dat",
        "HE0107-5240.dat",
        "HE0557-4840.dat",
        "HE1327-2326.dat",
        "SM0313-6708.dat",
    ],
    starnames=[
        "A",
        "B",
        "C",
        "D",
        "E",
    ],
):
    stars = []
    # Read all the star files into an object, then put into list
    for filename in filenames:
        stars += [Star(filename)]

    # Use only elements that are in at least 1 star
    z = set()
    el = [] * len(stars)
    for i, star in enumerate(stars):
        el += [star.element_abundances["element"].tolist()]
        z |= {element.Z for element in el[i]}
    z = list(z)
    z = [i for i in z if i <= 30]
    z.sort()
    for elrow in el:
        elrow = [e for e in elrow if e.Z <= 30]
    maxcolumns = 7
    tablecount = round((len(z) / maxcolumns) + 0.5)
    print(r"\begin{tabular}{|l|r@{}r|r@{}r|r@{}r|r@{}r|r@{}r|r@{}r|r@{}r|}")
    for block in range(tablecount):
        print(r"\hline")
        print(r"\multicolumn{1}{|c|}{\textbf{Star}}", end="")
        for element in z[block * maxcolumns : block * maxcolumns + maxcolumns]:
            print(
                r" & \multicolumn{2}{c|}{\textbf{"
                + I(int(element)).element_symbol()
                + r"}}",
                end="",
            )
        print(r"\\")
        print(r"\hline")
        # Read those elements from the objects, then sort by charge
        for i, (star, e) in enumerate(zip(stars, el)):
            print(f"{starnames[i]:s}", end="")
            for number in z[block * maxcolumns : block * maxcolumns + maxcolumns]:
                if number in [element.Z for element in e]:
                    row = star.element_abundances[
                        [element.Z for element in e].index(number)
                    ]
                    if row.error < 0:
                        uplimstring = "<"
                    else:
                        uplimstring = " "
                    abustring = r"{:4.2f} \pm {:3.2f}".format(
                        row.abundance, abs(row.error)
                    ).rjust(15)
                    print(r" & ${:s}$ & ${:s}$".format(uplimstring, abustring), end="")
                else:
                    print(r" & $ $ & $             - $", end="")
            print(r"\\")
        print(r"\hline")
        if block == tablecount:
            break
        else:
            print(r"\noalign{\smallskip}")
    print(r"\end{tabular}")


if __name__ == "__main__":
    import sys

    args = sys.argv[1:]
    for i, a in enumerate(args):
        try:
            args[i] = eval(a)
        except:
            pass
    print(Star(*args))
