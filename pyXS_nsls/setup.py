#!/usr/bin/env python

"""
setup.py file for RQconv module
"""

from distutils.core import setup, Extension
import numpy
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

RQ_module = Extension('_RQconv',
                      sources=['RQconv.i', 'RQconv.c'],
                      include_dirs=[numpy_include],
                      )

setup (name = 'RQconv',
       version = '20091210',
       author      = "Lin Yang",
       description = """C module for scattering pattern conversion""",
       ext_modules = [RQ_module],
       )
