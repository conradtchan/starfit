#!/usr/bin/env bash

# Delete existing dists
rm -rf dist

# Build
python -m build --sdist

# Upload to PyPI
python -m twine upload --repository testpypi dist/*
