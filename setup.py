#!/usr/bin/env python

from distutils.core import *
from distutils import sysconfig
import os.path

# Get numpy include directory (works across versions)
import numpy
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

# CCL extension module
_ccllib = Extension(
            "_ccllib",
               ["pyccl/ccl.i",],
               libraries = ['m', 'gsl', 'gslcblas',  'ccl'], #'gomp',
               include_dirs = [numpy_include, "include/", "class/include", "/opt/local/include"],
               extra_compile_args=['-O4', '-std=c99'],
               swig_opts=['-threads'],
           )

# CCL setup script
setup(  name         = "pyccl",
        description  = "Library of validated cosmological functions.",
        author       = "LSST DESC",
        version      = "0.1",
        packages     = ['pyccl'],
        ext_modules  = [_ccllib,],
)
