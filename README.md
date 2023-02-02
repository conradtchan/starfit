![x](https://github.com/conradtchan/starfit/actions/workflows/test.yml/badge.svg)
![x](https://github.com/conradtchan/starfit/actions/workflows/pre-commit.yml/badge.svg)
![x](https://github.com/conradtchan/starfit/actions/workflows/publish.yml/badge.svg)

![GitHub commits since latest release (by SemVer including pre-releases)](https://img.shields.io/github/commits-since/conradtchan/starfit/latest/main?include_prereleases)
![GitHub Release Date](https://img.shields.io/github/release-date/conradtchan/starfit)

StarFit is a Python package for matching stellar abundance measurements against a database of model stellar explosions. Based on the [old IDL code](https://2sn.org/starfit/) by [Alexander Heger](https://2sn.org).

StarFit can match combined abundances of multiple models. For single stars and combinations of multiple stars, a complete search can be found.  For three or more stars, the problem is extremely expensive, so a [Genetic Algorithm](https://en.wikipedia.org/wiki/Genetic_algorithm) has been implemented by Conrad Chan to efficiently find an approximate solution.

An online interface (with a subset of functionality) is available at [starfit.org](https://starfit.org).
