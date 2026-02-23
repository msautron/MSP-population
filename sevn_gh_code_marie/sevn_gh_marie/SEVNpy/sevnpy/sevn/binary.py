"""
==============================================================
Binary , (:mod:`sevnpy.sevn.binary`)
==============================================================

This module contains the class Binary. It is used to initialise binary systems and evolve
them using the SEVN backend.

"""

from __future__ import annotations
import warnings
import numpy as np
import pandas as pd
import copy

from scipy.interpolate import interp1d

try:
    from . import sevnwrap as sw
except:
    raise ImportError("The sevnwrap is not installed")
from ..sevnpy_types import Number, Optional, Union, Dict, Any, Tuple, ListLikeType
from .sevnmanager import SEVNmanager
from .. import utility as ut
from ..io.logreader import readlogstring
from . import star

class Binary:
    """
    The class Binary is used to initialise a binary systems and evolve it using the SEVN backend.

    """

    def __init__(self):
        return