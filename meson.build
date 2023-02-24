project(
  'starfit',
  'fortran',
  'c',
version: run_command('python3', '-m', 'setuptools_scm', check: true).stdout().strip(),
default_options: ['buildtype=debugoptimized', 'optimization=3'],
meson_version: '>=0.64.0'
)

py_mod = import('python')
py3 = py_mod.find_installation('python3', pure: false)

incdir_numpy = run_command(py3,
    ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
    check : true
).stdout().strip()

incdir_f2py = run_command(py3,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

inc_np = include_directories(incdir_numpy, incdir_f2py)

if meson.get_compiler('c').get_id() == 'clang'
  # Meson's OpenMP dependency does not work on mac yet
  omp = declare_dependency(compile_args: ['-Xpreprocessor', '-fopenmp'], link_args: ['-lgomp'])
else
  omp = dependency('openmp')
endif

subdir('src/starfit')
