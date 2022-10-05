"""
Provide interface to load data files conveniently
"""

import os
import os.path
from functools import partial
from glob import glob

from .logged import Logged


class Loader(Logged):
    """
    Return load function for class.

    Usage:  (for example)
    in file:
        loadlc = loader(LCData)
    command line
        lcdata = loadlc('filename')
    """

    def __init__(
        self,
        c,
        name=None,
        extension=None,
        graphical=True,
        silent=False,
    ):
        if name is None:
            name = c.__name__
        if extension is None:
            try:
                extension = c._extension
            except AttributeError:
                extension = "*"
        self.c = c
        self.name = name
        self.graphical = graphical
        self.silent = silent
        self.extension = extension

    def __call__(self, *args, **kwargs):
        """
        Check whether file exists and load file.
        Return None if file does not exist.
        """
        silent = kwargs.pop("silent", self.silent)
        self.setup_logger(silent)
        graphical = kwargs.pop("graphical", self.graphical)
        if len(args) > 0:
            try:
                return self.c(*args, **kwargs)
            except IOError:
                filename = args[0]
                # should this be done with a logger?
                if not kwargs.get("silent", silent):
                    self.logger.error(
                        f' [{self.name:s}] ERROR: File "{filename:s}" not found.'
                    )
        # Maybe this should be done as part of the FortranRead class...
        compression_types = ("", ".gz", ".xz", ".bz2")
        files = []
        args = list(args)
        if self.extension != "*":
            for comp in compression_types:
                ext = "." + self.extension + comp
                files += glob("*" + ext)
        files = [f for f in files if not f.startswith("xxx")]
        if len(files) != 1:
            if not graphical:
                return None
            from tkinter.filedialog import askopenfilename

            fd_args = {"title": f"Choose {self.name} file"}
            if len(args) > 0:
                p = os.path.expanduser(os.path.expandvars(args[0]))
                if os.path.isdir(p):
                    fd_args["initialdir"] = p
                elif os.path.isfile(p):
                    fd_args["initialdir"] = os.path.dirname(p)
                elif os.path.isdir(os.path.dirname(p)):
                    fd_args["initialdir"] = os.path.dirname(p)
            if self.extension != "*":
                fd_args["filetypes"] = [
                    (self.name + " files", "." + self.extension + comp)
                    for comp in compression_types
                ]
            filename = askopenfilename(**fd_args)
            args[:1] = [filename]
        else:
            args[:1] = files[0:1]
        try:
            return self.c(*args, **kwargs)
        except IOError:
            self.logger.error(f' [{self.name:s}] ERROR Loading File "{filename:s}".')
            return None


# old interface
def loader(*args, **kwargs):
    return Loader(*args, **kwargs)


_loader = partial(loader, graphical=False, silent=True)
