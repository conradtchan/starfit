import os

from numpy.distutils.core import Extension, setup

flags = ["-fPIC", "-Ofast", "-g", "-funroll-loops", "-fno-second-underscore", "-w"]

cwd = os.getcwd()
starfit_src = os.path.join(cwd, "src/starfit")

module_solver = Extension(
    "starfit.fitness._solver",
    sources=(
        [
            os.path.join(starfit_src, "fitness", x)
            for x in [
                "type_def.f90",
                "powell.f90",
                "norm.f90",
                "mleqs.f90",
                "star_data.f90",
                "abu_data.f90",
                "solver.f90",
            ]
        ]
    ),
    extra_f90_compile_args=flags,
    f2py_options=[
        "--f2cmap",
        os.path.join(starfit_src, "fitness", ".f2py_f2cmap"),
    ],
)

module_solgen = Extension(
    "starfit.solgen._solgen",
    sources=[os.path.join(starfit_src, "solgen/solgen.f90")],
    extra_f90_compile_args=flags,
)

setup(ext_modules=[module_solver, module_solgen])
