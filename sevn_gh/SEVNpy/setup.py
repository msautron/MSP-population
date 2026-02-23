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



sevn_sources=include_all_sources()
wrap_sources=include_all_sources(root_path="sevnpy/sevn/src")
sevn_includes=include_all_dirs(root_path="../src")+["../include/",]
wrap_includes=["sevnpy/sevn/src/",]
sevn_headers=include_all_header(root_path="../src")+include_all_header(root_path="../include")
wrap_headers=include_all_header(root_path="sevnpy/sevn/src")

all_sources=sevn_sources+wrap_sources
all_includes=sevn_includes+wrap_includes
all_headers=sevn_headers+wrap_headers


inc_dirs=include_all_dirs(root_path="../src")+["../include/",]
module = Extension("sevnpy/sevn/sevnwrap",
                   sources=make_abs(all_sources),
                   include_dirs=make_abs(all_includes+[numpy.get_include(),]),
                   depends=make_abs(all_headers),
                   language="c++")

# @TODO 2- At a certain point we would like to install SEVNpy  directly linking the static or shared SEVN library
# The advantages will be:
#  1 - Much faster compilation (the SEVN libary are already compiled)
#  2- SEVNpy can be totally detached by SEVN (it means is much easier to put it in a repository such as pypy)
#  Below there is an example of using directly the library, but at this stage we have some problems to solve:
# 1- It seems that when import sevnwrap there are missing  names for the static variables, this is strange since
#    when used SEVN as a library with c++ code this is not happening, I should investigate
# 2- The above problem is solved if we add all the sources file in sources, however this generate another problem
#    due to the fact that the loaded utilities are now in the usr local directory and this seems to create problems
#    when deriving the SEVN path. Also in this case this is not happening when using SEVN as a library with C++ codes
#    something to investigate...
# Anyway, at this stage this kind of installation requires a double step from the user, install SEVN + install
# SEVNpy, but the compilation will be much faster
# For this first version we are using this possibilities
"""
module = Extension("sevnpy/sevn/sevnwrap",
                   sources=make_abs(["sevnpy/sevn/src/sevnwrap.cpp",]),
                   include_dirs=["/usr/local/include/sevn",numpy.get_include()],
                   libraries=["sevn_lib_static"],
                   library_dirs=["/usr/local/lib/sevn"],
                   depends=make_abs(inc_dirs+["sevnpy/sevn/src/sevnwrap.h",]),
                   language="c++"
                   )
"""

setup(
    name = "SEVNpy",
    version =__version__,
    author = "Giuliano Iorio",
    author_email = "giuliano.iorio.astro@gmail.com",
    description = ("A Python companion module for SEVN"),
    license = "MIT",
    packages=['sevnpy','sevnpy/sevn', 'sevnpy/utility', 'sevnpy/io'],
    ext_modules= [module],
    long_description=long_description,
    python_requires='>=3.7',
)
