"""
==============================================================
SEVN constants , (:mod:`sevnpy.sevn.sevnconst`)
==============================================================

This module store the physical constants used in SEVN.
It is currently under development and the current is just a prototype,
use it at your own risk.

Example
-------
To get the value of the const use, constname.value()

>>> from sevnpy.sevn.sevnconsts import G
>>> print(G.value())
>>> 392512559.8496094

To show all the info about the constant just print it

>>> from sevnpy.sevn.sevnconsts import G
>>> print(G())
>>> value: 392512559.8496094
>>> description: gravitational G constant in RSUN^3 YR^-2 MSUN^-1. Taken fro astropy astropy-4.3.1
>>> units: RSUN^3 YR^-2 MSUN^-1

"""

try:
    from . import sevnwrap as sw
except:
    raise ImportError("The sevnwrap is not installed")
from ..utility._privateutility import Singleton


class _SEVNconst:
    _value = None
    _description = None
    _units = None

    @classmethod
    def value(cls):
        return cls._value

    @classmethod
    def description(cls):
        return cls._description

    @classmethod
    def units(cls):
        return cls._units

    def __repr__(self):
        string = ""
        string += f"value: {self._value}\n"
        string += f"description: {self._description}\n"
        string += f"units: {self._units}\n"

        return string


class G(_SEVNconst, metaclass=Singleton):
    _value = sw.sevnconst_G
    _description = "gravitational G constant in RSUN^3 YR^-2 MSUN^-1. Taken from astropy astropy-4.3.1"
    _units = "RSUN^3 YR^-2 MSUN^-1"


class Rsun(_SEVNconst, metaclass=Singleton):
    _value = sw.sevnconst_Rsun
    _description = "Sun radius in cm. Taken from astropy astropy-4.3.1"
    _units = "cm"
