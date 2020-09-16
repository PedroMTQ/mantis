from distutils.core import setup
from Cython.Build import cythonize
from os import path,getcwd
from sys import platform

#to compile
#python setup_get_non_overlapping_hits.py build_ext --inplace

if 'win' in platform.lower():   splitter = '\\'
else:                           splitter = '/'
mantis_folder = path.abspath(path.dirname(__file__)).split(splitter)[0:-1]
mantis_folder = splitter.join(mantis_folder)+splitter
cython_folder=mantis_folder + 'cython_src'+splitter
setup(name='Get non overlapping hits',
      ext_modules=cythonize([cython_folder + "get_non_overlapping_hits.pyx"]))
