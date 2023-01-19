name = 'solgen'
path = 'starfit' / name
module_name = '_' + name

python_sources = [
    '__init__.py',
]

py3.install_sources(python_sources, subdir: path)

# f2py generates the wrapper automatically from the single
# fortran source file
f2py_wrapper = custom_target(module_name,
    input : ['solgen.f90'],
    output : [module_name + 'module.c'],
    command : [ py3, '-m', 'numpy.f2py', '@INPUT@',
    '-m', module_name, '--build-dir', '@OUTDIR@'])


py3.extension_module(module_name,
    f2py_wrapper,
    'solgen.f90',
    incdir_f2py+'/fortranobject.c',
    include_directories: inc_np,
    dependencies: py3.dependency(),
    install: true,
    subdir: path,
    )
