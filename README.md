# Huoguo

> My (over-)ambitious goal to achieve `Isophote.x` and `galfit` functions in Python without losing the efficiency

## Motivation

* Both 1-D isophotal analysis and 2-D modeling of light distribution of galaxies are useful in dealing with images of galaxies.
* The most commonly used 1-D analysis tool is the [`ellipse`](http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?ellipse.hlp) function in the [`stsdas.analysis.isophote`](http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?isophote.men)
  package used in `IRAF` environment.
    - The source code is basically written in C. It is very efficient, but since it is designed for the outdated `IRAF` environment, it is not very convinient to integrate it into today's image analysis pipeline.
    - The `isophote` package is already available in Python under the [`photutils.isophote` module](https://github.com/astropy/photutils/tree/master/photutils/isophote). But this is a pure-Python implementation, so the efficiency is much slower than the `IRAF` version.  And it also lacks some of the functionalities of the origin `isophote`. e.g., fix the geometry of the isophote, and the force photometry mode.
    - There is a very [good Python wrapper of `ellipse` by Peter Erwin](https://github.com/perwin/ellipsefits).
    - There is a high-order harmonic expansion of `isophote` called [`isofit`](https://github.com/BogdanCiambur/ISOFIT) by Bogdan Ciambur. The [relevant publication is here](https://arxiv.org/pdf/1507.02691.pdf)
* [`GALFIT`](https://users.obs.carnegiescience.edu/peng/work/galfit/galfit.html) is a very popular galaxy image fitting tool designed by Chien Peng. Although right now there are many very good alternative choices (e.g. [`imfit` by Peter Erwin](https://www.mpe.mpg.de/~erwin/code/imfit/), ['ProFit` by ICRAR](https://github.com/ICRAR/ProFit)).
    - However, `GALFIT v3.0` remains one of the most flexible image fitting tool with the capabilities to fit asymmetric, truncated, or evenr coordinate-rotated 2-D components. Such capabilities can be very useful in many occasions.

## Practical Reasons

* `huoguo` is also my approach to learn all the dirty laundries about galaxy photometry in both 1-D and 2-D.
* I also want to learn how to write `CPython` code to enhance the efficiency of Python code.

# Development

* This is one of my (too-)many side projects, so it will probably be updated (very) slowly.
* But do let me know if you are also interested in this topic and want to make contribution.