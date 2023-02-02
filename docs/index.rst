.. starfit documentation master file, created by
   sphinx-quickstart on Thu Feb  2 10:44:49 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to StarFit's documentation!
===================================

Python package for matching stellar abundance measurements against a database of model stellar explosions. Based on the `old IDL code <https://2sn.org/starfit/>`_ by `Alexander Heger <https://2sn.org>`_.

StarFit can match combined abundances of multiple models. For single stars and combinations of multiple stars, a complete search can be found.  For three or more stars, the problem is extremely expensive, so a `Genetic Algorithm <https://en.wikipedia.org/wiki/Genetic_algorithm>`_ has been implemented by Conrad Chan to efficiently find an approximate solution.

An online interface (with a subset of functionality) is available at `starfit.org <https://starfit.org>`_.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation.md
   usage.md
   custom_data.md
   contributions.md


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
