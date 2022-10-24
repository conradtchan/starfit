"""
Some extra functionality for UUID-1
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


def UUID1():
    return UUID(uuid.uuid1().hex)


ufunc_UUID = np.frompyfunc(lambda x: UUID(bytes=x), 1, 1)
