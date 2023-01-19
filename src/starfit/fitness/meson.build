fortran_sources = [
    'solver.pyf',
    'typedef.f90',
    'utils.f90',
    'powell.f90',
    'norm.f90',
    'linalg.f90',
    'star_data.f90',
    'abu_data.f90',
    'fitting.f90',
    'solver.f90',
]


incdir_numpy = run_command(py3,
    ['-c', 'import os; os.chdir('..'); import numpy; print(numpy.get_include())'],
    check : true
).stdout().strip()

incdir_f2py = run_command(py3,
    ['-c', 'import os; os.chdir('..'); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

module = custom_target('solvermodule.c',
                            input : ['solver.pyf', 'solver.f90'],
                            output : ['solvermodule.c', 'solver-f2pywrappers.f'],
                            command : [ py3, '-m', 'numpy.f2py', '-c', '@INPUT@',
                            '--build-dir', '@OUTDIR@'])

inc_np = include_directories(incdir_numpy, incdir_f2py)

py3.extension_module('fitness',
    fortran_sources,
    module,
    incdir_f2py+'/fortranobject.c',
    include_directories: inc_np,
    dependencies: py3.dependency(),
    install: true,
    subdir: 'starfit',
    )
