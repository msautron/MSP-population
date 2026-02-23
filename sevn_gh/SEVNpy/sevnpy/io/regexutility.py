"""
==============================================================
Regex utility, (:mod:`sevny.io.regexutility`)
==============================================================

This module contains methods and utilities to use regex

"""

from typing import Union, Type, List, Optional, Any, Dict
import pandas as pd
import numpy as np
import re
#from ..utility import utility as ut
from .. import utility as ut


def regex_dataframe_from_string(string: str,
                                matchexpr: str,
                                columns: List[str],
                                columns_type: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Apply a regex match to a string and store the captured values to a pandas DataFrame

    Parameters
    ----------
    string:
        string to read
    matchexpr:
        regex expression to match
    columns:
        name of the capturing groups from the regex match
    columns_type:
        type of the capturing groups from the regex match


    Returns
    -------
    matching_value: pandas DataFrame
        Return a pandas DataFrame storing the values matched by the regex pattern.

    """

    ma = re.findall(matchexpr, string)
    na = np.array(ma)

    if (len(na) == 0):
        _df = pd.DataFrame({name: [] for name in columns})

    else:
        if len(columns) != na.shape[1]: raise ValueError(
            f"The dimension of the columns in input {columns} is not consistent"
            f"with the dimension of the captured items")

        if columns_type is not None:
            if len(columns_type) != len(columns): raise ValueError(
                f"The number of the column types in input {columns_type} is not consistent"
                f"with the number of the columns")

        _df = pd.DataFrame(na, columns=columns)

    if columns_type is not None:
        type_dict = dict(zip(columns, columns_type))
        _df = _df.astype(type_dict)

    return _df


def regex_from_file(filename: str,
                    matchexpr: str,
                    columns: List[str],
                    columns_type: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Apply a regex match to the content of a file and store the captured values to a pandas DataFrame

    Parameters
    ----------
    filename:
        file to read
    matchexpr:
        regex expression to match
    columns:
        name of the capturing groups from the regex match
    columns_type:
        type of the capturing groups from the regex match


    Returns
    -------
    matching_value: pandas DataFrame
        Return a pandas DataFrame storing the values matched by the regex pattern.

    """
    with open(filename, "r") as fo:
        df = regex_dataframe_from_string(string=fo.read(), matchexpr=matchexpr, columns=columns,
                                         columns_type=columns_type)

    return df


capturing = lambda value: f"({value})"
"""Return the string (value) where value is the input"""
notcapturing = lambda value: f"(?:{value})"
"""Return the string (?:value) where value is the input"""

class ReTypeMatch:
    """
    This class is a pure static class used to retrieve the regex matching pattern for various types.
    At the moment the allowed types are:

        - **'int'**, **'type'** or **int**: matching pattern for an integer numer
        - **'id'**: matching pattern for an ID type, i.e. a positive int
        - **'float'** or **float**: matching pattern for a float number
        - **'str'** or **'name'**: generic matching pattern for a string


    Examples
    --------
    ReTypeMatch can be used directly as dictionary providing the type of
    matching pattern we want to retrieve, e.g.

    >>> ReTypeMatch["str"]
       '(?:[0-9|A-Za-z]*\\_)?[0-9]*'
    >>> ReTypeMatch[float]
    '[+|-]?[0-9]+\\.?[0-9]*(?i:e)?[+|-]?[0-9]*|(?i:nan)'

    It is also possibile to guess the matching patter for a given generic input using the static method guess

    >>> ReTypeMatch.guess("hello world")
    '(?:[0-9|A-Za-z]*\\_)?[0-9]*'
    >>> ReTypeMatch.guess("3e10")
    '[+|-]?[0-9]+\\.?[0-9]*(?i:e)?[+|-]?[0-9]*|(?i:nan)'
    >>> ReTypeMatch.guess(13.2)
    '[+|-]?[0-9]+\\.?[0-9]*(?i:e)?[+|-]?[0-9]*|(?i:nan)'
    >>> ReTypeMatch.guess(2)
    '[+|-]?[0-9]+\\.?[0-9]*(?i:e)?[+|-]?[0-9]*|(?i:nan)'

    """

    matchtype = r"[+|-]?\d+"
    matchnum = r"[+|-]?[0-9]+\.?[0-9]*(?i:e)?[+|-]?[0-9]*|(?i:nan)"
    matchid = r"[0-9]+"
    matchname = r"(?:[0-9|A-Za-z]*\_)?[0-9]*"

    transform_dict = {
        "int": matchtype,
        "float": matchnum,
        "id": matchid,
        "type": matchtype,
        "str": matchname,
        "name": matchname,
        int: matchtype,
        float: matchnum,
        str: matchname
    }

    @classmethod
    def guess(cls, guess_value: Union[int, float, str]) -> str:
        """
        Guess the regex pattern to match in the input value

        Parameters
        ----------
        guess_value: int|float|str
            Input value for which we want to guess the related pattern matching

        Returns
        -------
        regex_pattern: str
            The regex pattern matching for the input value

        """
        if isinstance(guess_value, int):
            return cls.transform_dict[int]
        elif isinstance(guess_value, float):
            return cls.transform_dict[float]
        elif isinstance(guess_value, str) and guess_value.isdigit():
            return cls.transform_dict[int]
        elif isinstance(guess_value, str) and ut.str_is_float(guess_value):
            return cls.transform_dict[float]
        elif isinstance(guess_value, str):
            return cls.transform_dict[str]
        else:
            raise ValueError(f"Regex pattern match for input value {guess_value} cannot be guessed")

    def __class_getitem__(cls, key: Union[str, Type[int], Type[float], Type[str]]) -> str:
        try:
            return cls.transform_dict[key]
        except KeyError:
            raise KeyError(f"key \'{key}\' not available.Available keys are {list(cls.transform_dict.keys())}")

        return cls.transform_dict[key]

    @classmethod
    def capturing(cls, key: Union[str, Type[int], Type[float], Type[str]]) -> str:
        """
        Get the regex matching patter considering a capturing group

        Parameters
        ----------
        key: str or type int,str,float
            key for which we want to get the matching pattern

        Returns
        -------
        matching_pattern: str
            return the matching pattern for a capturing group (match_pattern)

        Examples
        --------

        >>> ReTypeMatch.capturing(float)
        >>> '([+|-]?[0-9]+\\.?[0-9]*(?i:e)?[+|-]?[0-9]*|(?i:nan))'

        """
        return capturing(cls[key])

    @classmethod
    def notcapturing(cls, key: Union[str, Type[int], Type[float], Type[str]]) -> str:
        """
        Get the regex matching patter considering a non-capturing group

        Parameters
        ----------
        key: str or type int,str,float
            key for which we want to get the matching pattern

        Returns
        -------
        matchin_pattern: str
            return the matching pattern for a non-capturing group (?:match_pattern)

        Examples
        --------

        >>> ReTypeMatch.capturing(float)
        >>> '(?:[+|-]?[0-9]+\\.?[0-9]*(?i:e)?[+|-]?[0-9]*|(?i:nan))'

        """
        return notcapturing(cls[key])
