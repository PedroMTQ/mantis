import os

from setuptools import setup
from Cython.Build import cythonize
from os import path,getcwd,walk,rename
from shutil import copy

#to compile
#python setup_get_non_overlapping_hits.py build_ext --inplace

def move_o_file():
    current_folder = getcwd()+SPLITTER+'build'+SPLITTER
    for root, dirs, files in walk(current_folder):
        if 'get_non_overlapping_hits.o' in files:
            o_file= path.join(root, 'get_non_overlapping_hits.o')
            try:
                rename(o_file, CYTHON_FOLDER + 'get_non_overlapping_hits.o')
            except:
                copy(o_file, CYTHON_FOLDER + 'get_non_overlapping_hits.o')
                os.remove(o_file)
            return


SPLITTER = '/'
MANTIS_FOLDER = path.abspath(path.dirname(__file__)).split(SPLITTER)[0:-1]
MANTIS_FOLDER = SPLITTER.join(MANTIS_FOLDER)+SPLITTER
CYTHON_FOLDER=MANTIS_FOLDER + 'cython_src'+SPLITTER

setup(name='Get non overlapping hits',
      ext_modules=cythonize([CYTHON_FOLDER + "get_non_overlapping_hits.pyx"]))
      
move_o_file()