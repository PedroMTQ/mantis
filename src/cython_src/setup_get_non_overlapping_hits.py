import os
from os import getcwd, path, rename, walk
from shutil import copy

from Cython.Build import cythonize
from setuptools import setup

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

def move_so_file():
    for root, dirs, files in walk(CYTHON_FOLDER):
        for f in files:
            if f.endswith('.so') and f.startswith('get_non_overlapping_hits.'):
                so_file = path.join(root, f)
                try:
                    rename(so_file, CYTHON_FOLDER + 'get_non_overlapping_hits.so')
                except:
                    copy(so_file, CYTHON_FOLDER + 'get_non_overlapping_hits.so')
                    os.remove(so_file)
                return


SPLITTER = '/'
CYTHON_FOLDER = os.path.abspath(os.path.dirname(__file__))+SPLITTER

setup(name='Get non overlapping hits',
      ext_modules=cythonize([CYTHON_FOLDER + "get_non_overlapping_hits.pyx"]))

move_o_file()
move_so_file()
