from distutils.core import setup
from Cython.Build import cythonize

#to compile
#python setup_get_non_overlapping_hits.py build_ext --inplace

setup(name='Get non overlapping hits',
      ext_modules=cythonize("get_non_overlapping_hits.pyx"))
