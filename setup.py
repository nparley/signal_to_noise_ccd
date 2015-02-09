"""
Can be used to build the cython code. However, this should be done for you when imported using
pyximport
"""

__author__ = 'Neil Parley'
from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules=cythonize("aperture_sum.pyx"), requires=['Cython', 'flask', 'wtforms']
)