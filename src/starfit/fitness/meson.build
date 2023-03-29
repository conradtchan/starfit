name = 'fitness'
path = 'starfit' / name
module_name = '_solver'

python_sources = [
    '__init__.py',
    'solver.py',
]

fortran_library_sources = [
    'typedef.f90',
    'utils.f90',
    'powell.f90',
    'norm.f90',
    'linalg.f90',
    'star_data.f90',
    'abu_data.f90',
    'fitting.f90',
]

fortran_interface_source = 'solver.f90'

fortran_test_sources = [
    'test.f90',
    'tests.f90',
]

py3.install_sources(python_sources, subdir: path)

# Build static library with the actual code
fortran_library = static_library(module_name, fortran_library_sources, dependencies: omp)

# Create the f2py wrapper using only the interface
f2py_wrapper = custom_target(module_name,
    input : [fortran_interface_source, '.f2py_f2cmap'],
    output : [module_name + 'module.c'],
    command : [ py3, '-m', 'numpy.f2py', '@INPUT0@',
    '-m', module_name,
    '--build-dir', '@OUTDIR@',
    '--f2cmap', '@INPUT1@'])

# Build module using the f2py wrapper and static library
py3.extension_module(module_name,
    f2py_wrapper,
    fortran_interface_source,
    incdir_f2py+'/fortranobject.c',
    link_with: fortran_library,
    include_directories: inc_np,
    dependencies: [py3.dependency(), omp],
    install: true,
    subdir: path,
    )

# Build target for directly building the test executable
executable('fitness_test',
    fortran_test_sources,
    link_with: fortran_library,
    build_by_default: false)
