"""
==============================================================
Logfile Handling, (:mod:`sevny.io.sevndata`)
==============================================================

This module contains the class used to store the SEVN output


"""

from typing import Dict, Optional, List
import pandas as pd
import copy


class SEVNDataResults:
    """
    Template to to create classes to store results from SEVN runs.
    The data are stored in a pandas DataFrame, the other attributes contain info
    about the data and the columns and the list of the files used to load the data.

    """

    def __init__(self):
        self._data = pd.DataFrame()
        self._data_info = None
        self._columns_info = None
        self._sevn_files = None

    @property
    def data(self) -> pd.DataFrame:
        """
        Copy of the internal DataFrame

        """
        return self._data.copy()

    @property
    def data_info(self) -> str:
        """
        Info about the stored data

        """
        return copy.deepcopy(self._data_info)

    @property
    def columns_info(self) -> Dict:
        """
        Dictionary containg the pairs colum_nam:column_info

        """
        return copy.deepcopy(self._columns_info)

    @property
    def sevn_files(self) -> List:
        """
        List containg the read SEVN output files

        """
        return copy.deepcopy(self._sevn_files)


class SEVNDataLog(SEVNDataResults):
    """
    Specialised :class:`SEVNDataResults` class to handle data from the log reading

    """

    def __init__(self,
                 dataframe: pd.DataFrame,
                 data_info: Optional[str] =None,
                 column_descriptions: Optional[Dict[str,str]]=None,
                 regex_pattern: Optional[str]=None,
                 sevn_files: Optional[List]=None):
        """

        Parameters
        -----------
        dataframe:
            Dataframe containing the log data
        data_info:
            info about the log data
        column_description:
            description of the columns in the dataframe
        regex_pattern:
            regex pattern used to catch values from the logdata
        sevn_files:
            list of the files containing the analysed log data


        """
        super().__init__()

        self._data = dataframe
        self._columns_info = column_descriptions
        self._sevn_files = sevn_files
        self._regex_pattern = regex_pattern
        self._data_info = data_info

    @property
    def regex_pattern(self) -> str:
        """
        regex pattern used in the load data file

        """
        return copy.deepcopy(self._regex_pattern)
