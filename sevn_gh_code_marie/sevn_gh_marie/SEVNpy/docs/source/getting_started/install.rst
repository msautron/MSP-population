************
Installation
************

.. _sevnpy-main-req:

Getting SEVNpy
**************

The module ``sevnpy`` is inside the main SEVN repository at https://gitlab.com/sevncodes/sevn.
The module in included in the folder SEVNpy (https://gitlab.com/sevncodes/sevn/-/tree/SEVN/SEVNpy).
Additional documentation can be found in the SEVN userguide (https://gitlab.com/sevncodes/sevn/-/blob/SEVN/resources/SEVN_userguide.pdf)

Requirements
************

``sevnpy`` has the following strict requirements:

- Python3.7 or later
- numpy
- scipy
- pandas
- setuptools
- typing_extensions
- C++ (std:c++11) to install the SEVN wrappers

The following packages can optionally be used when testing:

-  pytest

.. _installing-sevnpy:

Installing ``sevnpy``
**********************

Using pip
=========

If pip is available is recommended to use it to install sevnpy.
Enter in the sevnpy folder and  run::

    pip install .

Basic Python installation
=========================

If pip is not avaible, enter in the sevnpy folder and just execute the setup.py script run::

    python setup.py install

Installing without SEVN c-extension
===================================
If you don't want to install the SEVN c-extension for any reason, use the setup file setup_npsevncext.py
instead of the classical setupy.py::

    python setup_nosevncext.py install