import time  # always good to have
from pathlib import Path
from socket import gethostname

import numpy as np
from abuset import AbuData, AbuSet
from human.time import time2human
from ionmap import decay
from isotope import ion as I
from stardb import StarDB
from tqdm import tqdm

_translate = {
    "h": "h1",
    "p": "h1",
    "a": "he4",
    "n": "nt1",
}


class YisoNod(object):
    """
    class for Limongi & Chieffi 2012 *.yiso_nod files
    """

    def __init__(self, filename):
        abu = dict()
        with open(filename, "rt") as f:
            next(f)
            for l in f:
                x = l.split()
                i = x[0].lower()
                i = I(_translate.get(i, i)).isotope()
                a = float(x[4])
                try:
                    abu[i] += a
                except:
                    abu[i] = a
        abu = AbuSet(abu, sort=True)
        self.mass = float(Path(filename).name[1:4])
        ejecta = abu.norm()
        self.remnant = self.mass - ejecta
        self.abu = abu.normalized()


class Data(object):
    hostname = gethostname()
    if hostname == "w.2sn.net":
        sdir = Path("/home/alex/starfit_data/contrib/LimongiChieffi/lc12/")
    else:
        raise AttributeError(f"Unknown machine {hostname=}")

    def __init__(self):
        pattern = "*.yiso_nod"
        files = list()
        for filename in (self.sdir).glob(pattern):
            files.append(filename)
        print(f" [{self.__class__.__name__}] Found {len(files)} models.")
        self.files = files

    def load_dumps(self):
        starttime = time.time()
        if not hasattr(self, "dumps"):
            print(f" [{self.__class__.__name__}] loading dumps.")
            dumps = list()
            for filename in tqdm(self.files):
                dumps.append(YisoNod(filename))
            self.dumps = dumps
            print(f" [{self.__class__.__name__}] loaded {len(dumps)} nucleo dumps.")
        print(
            f" [{self.__class__.__name__}] finished in {time2human(time.time() - starttime)}"
        )

    def make_stardb(
        self,
        filename=None,
        mode=(
            "alliso",
            "radiso",
            "deciso",
            "el",
        ),
    ):
        # def make_stardb(self, filename=None, mode='el', ):
        starttime = time.time()
        self.load_dumps()

        if isinstance(mode, (list, tuple)):
            assert filename is None
            dbs = list()
            for m in mode:
                db = self.make_stardb(mode=m)
                dbs.append(db)
            return dbs

        comments = (
            "Zero metallicity yield data set",
            "Limongi & Chieffi, ApJS, 199, 38 (2012).",
        )

        if not (hasattr(self, "abu") and hasattr(self, "fielddata")):
            dtype = np.dtype(
                [
                    ("mass", np.float64),
                    ("remnant", np.float64),
                ]
            )

            fielddata = list()
            abu = list()
            for d in tqdm(self.dumps):
                mass = d.mass
                remnant = d.remnant
                fielddata.append((mass, remnant))
                abu.append(d.abu)

            abu = AbuData.from_abusets(abu)
            fielddata = np.array(fielddata, dtype=dtype)

            ii = np.argsort(fielddata)
            fielddata = fielddata[ii]
            abu.data = abu.data[ii]

            self.abu = abu
            self.fielddata = fielddata

        parms = dict()

        basename = "znuc_lc12"
        if mode == "el":
            data = decay(self.abu, molfrac_out=True, elements=True)
            parms["abundance_type"] = StarDB.AbundanceType.element
            parms["abundance_class"] = StarDB.AbundanceClass.dec
        elif mode == "alliso":
            data = self.abu.as_molfrac()
            parms["abundance_type"] = StarDB.AbundanceType.isotope
            parms["abundance_class"] = StarDB.AbundanceClass.raw
        elif mode == "radiso":
            data = decay(self.abu, molfrac_out=True, decay=True, stable=False)
            parms["abundance_type"] = StarDB.AbundanceType.isotope
            parms["abundance_class"] = StarDB.AbundanceClass.dec
        elif mode == "deciso":
            data = decay(self.abu, molfrac_out=True, decay=True, stable=True)
            parms["abundance_type"] = StarDB.AbundanceType.isotope
            parms["abundance_class"] = StarDB.AbundanceClass.dec
        else:
            raise AttributeError(f"Unknown {mode=}.")

        parms["name"] = f"{basename}.{mode}.y"
        parms["comments"] = comments

        parms["data"] = data
        parms["fielddata"] = self.fielddata

        parms["fieldnames"] = ["mass", "remnant"]
        parms["fieldunits"] = ["solar masses", "M_sun"]
        parms["fieldtypes"] = [StarDB.Type.float64] * 2
        parms["fieldformats"] = ["2G", "6.3F"]
        parms["fieldflags"] = [StarDB.Flags.parameter] + [StarDB.Flags.property]

        parms["abundance_unit"] = StarDB.AbundanceUnit.mol_fraction
        parms["abundance_total"] = StarDB.AbundanceTotal.ejecta
        parms["abundance_norm"] = None
        parms["abundance_data"] = StarDB.AbundanceData.all_ejecta
        parms["abundance_sum"] = StarDB.AbundanceSum.number_fraction

        db = StarDB(**parms)
        if filename is None:
            filename = parms["name"] + ".stardb.xz"
        db.write(filename)

        print(
            f" [{self.__class__.__name__}] finished in {time2human(time.time() - starttime)}"
        )
        return db
