project(
  'starfit',
  'c',
# Run command to rewrite version number:
# meson rewrite kwargs set project / version $(git describe --tags)
version : '0.0.0',
default_options : ['warning_level=2']
)

add_languages('fortran')

py_mod = import('python')
py3 = py_mod.find_installation('python3')

incdir_numpy = run_command(py3,
    ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
    check : true
).stdout().strip()

incdir_f2py = run_command(py3,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

inc_np = include_directories(incdir_numpy, incdir_f2py)

subdir('src/starfit')
