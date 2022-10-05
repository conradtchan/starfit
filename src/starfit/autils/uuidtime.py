"""
Some extra functionallity for UUID-1
"""

import datetime
import time
import uuid

import numpy as np


class UUID(uuid.UUID):
    """Class converting UUID Version 1 time information"""

    time_grain = 10000000
    time_offset = (
        (
            datetime.datetime(1970, 1, 1).toordinal()
            - datetime.datetime(1582, 10, 15).toordinal()
        )
        * time_grain
        * 86400
    )

    def time_uuid2unix(self, time):
        return (time - self.time_offset) / self.time_grain

    def get_time(self):
        time = self.time
        return self.time_uuid2unix(time)

    def ctime(self):
        itime = self.get_time()
        if itime < 0:
            return ""
        else:
            return time.asctime(time.gmtime(itime))

    def ctimex(self):
        s = self.ctime()
        if s != "":
            s = f"({s:s})"
        return s

    def ip(self):
        s = self.hex
        return (
            s[20:22]
            + ":"
            + s[22:24]
            + ":"
            + s[24:26]
            + ":"
            + s[26:28]
            + ":"
            + s[28:30]
            + ":"
            + s[30:32]
        )

    def machine(self):
        return iplib.get(self.ip(), ("unknown", "unknown"))[0]


def UUID1():
    return UUID(uuid.uuid1().hex)


ufunc_UUID = np.frompyfunc(lambda x: UUID(bytes=x), 1, 1)


# import uuidtime
# x = uuidtime.UUID('7c3f0714-bd03-11df-b722-0022154815eb')
# x.ctime()


iplib = {
    "00:22:15:48:15:eb": ("home.2sn.org [MSP]", "AMD Athlon"),  # ?
    "00:21:6a:40:10:26": (
        "w.spa.umn.edu",
        "Intel mobile",
    ),  # ? is this startic or wireless?
    "00:30:48:66:f5:70": ("e.spa.umn.edu", "AMD Opteron"),  # ?
    "60:a4:4c:e6:29:e4": (
        "home.2sn.org [MEL]",
        "Intel(R) Xeon(R) CPU E3-1270 V2 @ 3.50GHz",
    ),  # 00:25:90:58:e8:e5
    "00:24:8c:19:5e:bd": (
        "b.spa.umn.edu",
        "Intel(R) Core(TM) i7 CPU 965 @ 3.20GHz",
    ),  # lk, 00:24:8c:19:5e:bc
    "00:e0:81:79:9f:23": (
        "c.spa.umn.edu",
        "Quad-Core AMD Opteron(tm) Processor 8356",
    ),  # 00:E0:81:79:9F:24
    "00:1b:21:ab:c1:ea": (
        "g.spa.umn.edu",
        "Intel(R) Xeon(R) CPU E31280 @ 3.50GHz",
    ),  # ?
    "00:25:90:06:3f:ee": (
        "i.spa.umn.edu",
        "Intel(R) Xeon(R) CPU X5680 @ 3.33GHz",
    ),  # 00:25:90:06:3F:EF
    "00:25:90:59:5c:1a": (
        "j.spa.umn.edu",
        "Intel(R) Xeon(R) CPU E5-4650 0 @ 2.70GHz",
    ),  # (only one port)
    "00:e0:81:b8:3b:f8": (
        "k.spa.umn.edu",
        "Six-Core AMD Opteron(tm) Processor 8431",
    ),  # 00:E0:81:B8:3B:F6, 00:E0:81:B8:3B:F7
    "00:30:48:fe:48:d7": (
        "l.spa.umn.edu",
        "AMD Opteron(tm) Processor 6176 SE",
    ),  # 00:30:48:FE:48:D6
    "00:25:90:12:9f:b8": (
        "v.spa.umn.edu",
        "Intel(R) Xeon(R) CPU X5680 @ 3.33GHz",
    ),  # 00:25:90:12:9F:B9
    "20:cf:30:f1:7d:c8": (
        "vo.spa.umn.edu",
        "Intel(R) Xeon(R) CPU W3680 @ 3.33GHz",
    ),  # (only one port)
    "20:cf:30:25:1c:dd": (
        "behemoth.spa.umn.edu",
        "Intel(R) Xeon(R) CPU W3680 @ 3.33GHz",
    ),  # 20:CF:30:25:1C:3C
    "48:5b:39:0a:7b:e8": (
        "mc.spa.umn.edu",
        "Intel(R) Xeon(R) CPU W3680 @ 3.33GHz",
    ),  # 48:5b:39:0a:80:23
    "48:5b:39:18:97:18": (
        "hou.spa.umn.edu",
        "Intel(R) Xeon(R) CPU W3680  @ 3.33GHz",
    ),  # 48:5B:39:18:8E:CF
    "00:24:8c:19:5e:be": (
        "wen.spa.umn.edu",
        "Intel(R) Core(TM) i7 CPU 965 @ 3.20GHz",
    ),  # 0:24:8c:19:5e:bf
    "d4:be:d9:65:16:80": (
        "zinc.maths.monash.edu",
        "Intel(R) Core(TM) i7-3820QM CPU @ 2.70GHz",
    ),  # 60:67:20:61:12:30 (wlan)
    "18:03:73:20:10:70": (
        "aurum.maths.monash.edu",
        "Intel(R) Core(TM) i5-3550 CPU @ 3.30GHz",
    ),  # (only one port)
    "00:25:90:91:10:54": (
        "boron.maths.monash.edu",
        "Intel(R) Xeon(R) CPU E5-2690 0 @ 2.90GHz",
    ),  # 00:25:90:91:10:55
    "00:25:90:58:e8:e4": (
        "carbon.sci.monash.edu",
        "Intel(R) Xeon(R) CPU E5-4650 0 @ 2.70GHz",
    ),  # 00:25:90:58:e8:e5
    "18:03:73:1f:f6:1e": (
        "deuterium.maths.monash.edu",
        "Intel(R) Core(TM) i5-3550 CPU @ 3.30GHz",
    ),  # (only one port)
    "90:b1:1c:64:37:35": (
        "erbium.maths.monash.edu",
        "Intel(R) Core(TM) i5-3570 CPU @ 3.40GHz",
    ),  # (only one port)
    "e0:db:55:0b:af:ad": (
        "latte.sci.monash.edu",
        "Intel(R) Xeon(R) CPU E5-2680 0 @ 2.70GHz",
    ),  # (only one port enabled)
    "00:25:90:a2:dd:f0": (
        "chou.sci.monash.edu",
        "Intel(R) Xeon(R) CPU E3-1270 V2 @ 3.50GHz",
    ),  # 00:25:90:a2:dd:f1
}
# e - machine info,
# home.2sn.org[MSP] - machine info
# notebook - machine info
# banzai, medusa, tug's computer
# d, f - no runs on those
