project(
  'starfit',
  'fortran',
  'c',
version: run_command('python3', '-m', 'setuptools_scm', check: true).stdout().strip(),
default_options: ['buildtype=plain'],
meson_version: '>=0.64.0'
)

fc_flags = ['-Ofast', '-fno-finite-math-only', '-fPIC']
fc = meson.get_compiler('fortran')
foreach f: fc_flags
  if fc.has_argument(f)
    add_project_arguments(f, language: 'fortran')
  endif
endforeach

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

subdir('src/starfit')
