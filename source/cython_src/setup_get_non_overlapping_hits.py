from distutils.core import setup
from Cython.Build import cythonize
from os import path,getcwd
from sys import platform

#to compile
#python setup_get_non_overlapping_hits.py build_ext --inplace

if platform.startswith('win'):    SPLITTER = '\\'
else:                                SPLITTER = '/'

MANTIS_FOLDER = path.abspath(path.dirname(__file__)).split(SPLITTER)[0:-1]
MANTIS_FOLDER = SPLITTER.join(MANTIS_FOLDER)+SPLITTER
cython_folder=MANTIS_FOLDER + 'cython_src'+SPLITTER
setup(name='Get non overlapping hits',
      ext_modules=cythonize([cython_folder + "get_non_overlapping_hits.pyx"]))
