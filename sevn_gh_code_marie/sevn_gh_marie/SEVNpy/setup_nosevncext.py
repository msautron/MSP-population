from setuptools import setup, Extension
import sys
import glob
import os
import numpy

""""""""""""""""""""""""""""""""""""
"          SETUP UTILS             "
""""""""""""""""""""""""""""""""""""
def include_all_dirs(root_path="../src"):
    root_path=os.path.abspath(root_path)
    return glob.glob(f"{root_path}/**/*/",recursive=True)

def include_all_sources(root_path="../src"):
    root_path=os.path.abspath(root_path)
    return glob.glob(f"{root_path}/**/*.cpp",recursive=True)
def include_all_header(root_path="../src"):
    root_path=os.path.abspath(root_path)
    return glob.glob(f"{root_path}/**/*.h",recursive=True)

def make_abs(list_path):

    return list(map(os.path.abspath,list_path))

""""""""""""""""""""""""""""""""""""
"          SETUP UTILS  END         "
""""""""""""""""""""""""""""""""""""


# @TODO 1- At a certain point we have to include in the setup.py file a pre-process in which we get the src,
#  include and tables from the SEVN parent directory, copy it there somehere and use this new directory
#  for the installation. In this way we can ship SEVNpy using something like pip or conda.



#Check Python version
if sys.version_info < (3,7):
    sys.exit('Sorry, Python < 3.7 is not supported by SEVNpy')




#To get the version
exec(open('sevnpy/version.py').read())
#To get the long description from README.md
with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name = "SEVNpy",
    version =__version__,
    author = "Giuliano Iorio",
    author_email = "giuliano.iorio.astro@gmail.com",
    description = ("A Python companion module for SEVN"),
    license = "MIT",
    packages=['sevnpy','sevnpy/sevn', 'sevnpy/utility', 'sevnpy/io'],
    long_description=long_description,
    python_requires='>=3.7',
)
