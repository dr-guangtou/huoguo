# Huoguo

- Customized 1-D photometry for galaxies from the Hyper Suprime-Cam survey and other imaging data sets. 
  - Mainly to support Song Huang's projects about the evolution and galaxy-halo connection of massive galaxies.

## Practical Motivation

- `Huoguo` will not try to re-invent the wheels. It will try to apply and compare a few different available tools or algorithms. 

- `Huoguo` will explore the fitting of 1-D profiles and try to bridge the 1-D and 2-D photometry for interesting science.

- `Huoguo` is also my approach to learn all the dirty laundries behind 1-D galaxy photometry, and see if it will survive in the next 10 years with all the LSST or JWST data.
- I also want to use this opportunity to learn some new tricks to make 1-D photometry more efficient. 
  - e.g., Use `CPython` to speed things up? Or maybe use `jax`? Or maybe rewrite some functions in `Julia`.
## Reference 

- The most commonly used 1-D analysis tool is the [`ellipse`](http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?ellipse.hlp) function in the [`stsdas.analysis.isophote`](http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?isophote.men)
  package used in `IRAF` environment.
    - The source code is basically written in C. It is very efficient, but since it is designed for the outdated `IRAF` environment, it is not very convinient to integrate it into today's image analysis pipeline.
    - The `isophote` package is already available in Python under the [`photutils.isophote` module](https://github.com/astropy/photutils/tree/master/photutils/isophote). But this is a pure-Python implementation, so the efficiency is much slower than the `IRAF` version.  And it also lacks some of the functionalities of the origin `isophote`. e.g., fix the geometry of the isophote, and the force photometry mode.
    - There is a very [good Python wrapper of `ellipse` by Peter Erwin](https://github.com/perwin/ellipsefits).
    - There is a high-order harmonic expansion of `isophote` called [`isofit`](https://github.com/BogdanCiambur/ISOFIT) by Bogdan Ciambur. The [relevant publication is here](https://arxiv.org/pdf/1507.02691.pdf)
    - The [`legacyhalos`](https://github.com/moustakas/legacyhalos) project by John Moustakas is also an important source of reference. 

- The [`AutoProf`](https://github.com/ConnorStoneAstro/AutoProf) (Automatic Isophotal solutions for galaxy images) package developed by Connor Stone and collaborators is a very exciting new tool on the market. 

- In the [`GalfitPyWrap`](https://github.com/Grillard/GalfitPyWrap) repo developed by Gabriel Torrealba contains an interesting implementation of elliptical isophotal fitting algorithm that is worth looking into.

- And [`Pix2Prof`: extracting sequential galaxy data with a language model](https://github.com/Smith42/pix2prof) is an interesting application of machine learning on 1-D photometry.

# Development

- This is one of my (too-)many side projects, so it will probably be updated (very) slowly.
- But do let me know if you are also interested in this topic and want to make contribution.
