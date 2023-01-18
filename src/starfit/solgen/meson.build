incdir_numpy = run_command(py3,
    ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
    check : true
).stdout().strip()

incdir_f2py = run_command(py3,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

source = custom_target('solgenmodule.c',
                            input : ['solgen.f90'],
                            output : ['solgenmodule.c', 'solgen-f2pywrappers.f'],
                            command : [ py3, '-m', 'numpy.f2py', '@INPUT@',
                            '-m', 'solgen', '--build-dir', '@OUTDIR@'])

inc_np = include_directories(incdir_numpy, incdir_f2py)

py3.extension_module('solgen',
    'solgen.f90',
    source,
    incdir_f2py+'/fortranobject.c',
    include_directories: inc_np,
    dependencies : py3.dependency(),
    install : true)
