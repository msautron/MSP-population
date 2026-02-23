"""
==============================================================
Remnant , (:mod:`sevnpy.sevn.remnant`)
==============================================================

This module contains methods about the remnants produce from the stellar evolution.

"""

import pandas as pd

from .star import Star
from ..sevnpy_types import Number, Optional, Union


def get_remnant(Mzams: Number,
                Z: Number,
                snmodel: str ="rapid",
                spin: Number =0,
                pureHe: bool = False,
                rseed: Optional[int] = None,
                Mass: Optional[Number] =None,
                MHE: Optional[Number] =None,
                MCO: Optional[Number] =None) -> pd.DataFrame:
    """
    Get the remnant properties from Star with given Mzams, Z (Masses) and for a given SEVN snmodel

    Parameters
    ----------
    Mzams:
        Zero-age main sequence mass of the star [Msun]
    Z:
        Absolute Metallicity of the star
    snmodel:
        SEVN snmodel to be used to create a remnant from the star
    spin:
        Stellar spin (Angular velocity over critical angular velocity)
    pureHe:
        If True consider this star as a pureHe star (Mzams in this case is the mase at the beginning of the core-helium
        burning)
    rseed:
        A random seed to be used for all the random calls during the evolution and the creation of the remnant
        If None a random random seed will be generated
    Mass:
        current total mass of the star [Msun].
        Use this parameter only if it necessary to modify the initial mass given by the star
    MHE:
        current helium core mass of the star [Msun].
        Use this parameter only if it necessary to modify the initial mass given by the star
    MCO:
        current co core mass of the star [Msun].
        Use this parameter only if it necessary to modify the initial mass given by the star

    Returns
    -------
    remnant:
        A Pandas dataframe containing the properties of the remnant

    """
    star_flag = "HE" if pureHe else "H"
    tini = "cheb" if  pureHe else "zams"
    s1 = Star(Mzams=Mzams, Z=Z, spin=spin, tini=tini,
              snmodel=snmodel, star_flag=star_flag,
              rseed=rseed,
              Mass=Mass, MHE=MHE, MCO=MCO)

    return s1.get_remnant()
