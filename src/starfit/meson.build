path = 'starfit'

python_sources = [
    '__init__.py',
    'direct.py',
    'fit.py',
    'ga.py',
    'multi.py',
    'single.py',
    'star.py',
    'starfit.py',
    'starplot.py',
    'utils.py',
    'autils/__init__.py',
    'autils/abumodel.py',
    'autils/abuset.py',
    'autils/abusets.py',
    'autils/isotope.py',
    'autils/keputils.py',
    'autils/loader.py',
    'autils/logged.py',
    'autils/physconst.py',
    'autils/stardb.py',
    'autils/utils.py',
    'autils/uuidtime.py',
    'autils/human/__init__.py',
    'autils/human/byte.py',
    'autils/human/config.py',
    'autils/human/power.py',
    'autils/human/prefix.py',
    'autils/human/range.py',
    'autils/human/template.py',
    'autils/human/time.py',
    'autils/human/unit.py',
    'autils/human/util.py',
    'autils/human/value.py',
    'autils/human/version.py',
]
py3.install_sources(python_sources, subdir: path, preserve_path: true)

install_subdir('data', install_dir: py3.get_install_dir() / path)

subdir('fitness')
subdir('solgen')
