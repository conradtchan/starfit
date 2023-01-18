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
message(py3.path())
message(py3.get_install_dir())

subdir('src/starfit')
