import os

from numpy.distutils.core import Extension, setup

flags = ["-fPIC", "-Ofast", "-funroll-loops", "-fno-second-underscore", "-w"]

cwd = os.getcwd()
starfit_src = os.path.join(cwd, "src/starfit")

module_solver = Extension(
    "starfit.fitness._solver",
    sources=[
        os.path.join(starfit_src, "fitness", x)
        for x in ["solver.f90", "powell.f90", "norm.f"]
    ],
    extra_f90_compile_args=flags,
)

module_solgen = Extension(
    "starfit.solgen._solgen",
    sources=[os.path.join(starfit_src, "solgen/solgen.f90")],
    extra_f90_compile_args=flags,
)

setup(ext_modules=[module_solver, module_solgen])
