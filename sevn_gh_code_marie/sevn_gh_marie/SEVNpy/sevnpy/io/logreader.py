"""
==============================================================
Logfile Handling, (:mod:`sevny.io.logreader`)
==============================================================

This module contains methods and utilities to read the SEVN logfiles

Methods
-------
    readlogfiles:
        wrapper function to easily read the logfiles
    readlogstr:
        wrapper function to easily read a log message from a string
    logcolumns_info:
        Methods returning info about the columns from the various events in the logfile
    logevents_list:
        Methods returning a list of all the available events

Class
-------
    LogReader:
        this is the basic class used to read the logfile getting the information of a specific event from the logfile.
        It automatically  set the necessary information to read the Header file, but the "schema" of the body needs to be provided.
        The module contains a list of pre-made classes inherited from LogReader to read all the events in the SEVN logfile

Developers: how to add a new specific LogReader class.
______________________________________________________

Let's assume that we want to define a new class that can be used to read the information of a new event called SPAM.
The event only is triggered by the single stellar evolution and produce just two values that we call Atest and Btest,
Aspam is in an integer value while Bspam is a float.
The schema will be

>>> body_schema={"Aspam": (int, "information about Aspam"), "Bspam": (float,"information about Bspam")}

So we will define a new class following the name structure <EVENTNAME>LogReader, since the event is triggered from SSE this class will inherit
from  _SSELogReader (that inherits from _SpecialisedLogReader and hence LogReader) to automatically include the informationa about the header.
In case the log event had been triggered by BSE, the class would have  inherited from _BSELogReader.
Then we have to define just the log_event name and the body schema, and optionally, but is warmly suggested a shor description
of the event in the attribute event_description

>>> class SPAMLogReader(_SpecialisedLogReader):
>>>     log_event = "SPAM"
>>>     event_description = "A simple example"
>>>     body_schema={"Aspam": (int, "information about Aspam"), "Bspam": (float,"information about Bspam")}

If we know that such event need to happen only once for a star or a binary (e.g. the event SN for a star),
it is enough to add as a first parent class  _OnlyOnePerSystemGetLast or  _OnlyOnePerSystemGetFirst,
the first one automatically remove the duplicates base on name,ID and IDfile mantaining the last occurence,
while the last one mantain just the last occurence. So for example

>>> class SPAMLogReader(_OnlyOnePerSystemGetLast,_SpecialisedLogReader):
>>>     log_event = "SPAM"
>>>     event_description = "A simple example"
>>>     body_schema={"Aspam": (int, "information about Aspam"), "Bspam": (float,"information about Bspam")}

Finally we have to include a new item ("SPAM",SPAMLogReader)  in the module attribute _logevents.

At this point the new log events will appear in the available logevents list and  can be read thorugh the wrapper function
readlogfiles.


"""
import warnings
import pandas as pd
import multiprocessing as mp
import glob
from functools import partial
from collections.abc import Iterable
import os
import copy as cp

from typing import Union, List, Type, Tuple, Optional, Dict, Literal

try:
    from typing import TypedDict
except ImportError:
    from typing_extensions import TypedDict

from .. import utility as ut
check_allowed = ut.check_allowed
#from ..utility import utility as ut
#from ..utility.utility import check_allowed
from .logschema import LogSchema
from .sevndata import SEVNDataLog
from . import regexutility as rut

GTUnion = Union[Type[int], Type[float], Type[float]]
FloatInt = Union[str, int]


# @TODO-2 Add all the Log
# @TODO-3 Documentation and examples
# @TODO-4 At certain point we have to include as parameter the SEVN version since the logfiles can change depending on the version


def _regex_from_string(string: str,
                       matchexpr: str,
                       columns: List[str],
                       columns_type: List[str]) -> pd.DataFrame:
    """
    Just a special wrapper of  :func:`~io.regexutility.regex_from_file`

    """
    df = rut.regex_dataframe_from_string(string=string, matchexpr=matchexpr, columns=columns, columns_type=columns_type)

    return df


def _regex_dataframe_from_file_adding_filename(filename: str,
                                               matchexpr: str,
                                               columns: List[str],
                                               columns_type: List[str]) -> Tuple[pd.DataFrame, int]:
    """
    Just a special wrapper of  :func:`~io.regexutility.regex_from_file`
    The only difference is that it returns the DataFrame with an additional column
    'IDfile' containing an 10 digit hascode generated from the read logfile

    """

    df = rut.regex_from_file(filename=filename, matchexpr=matchexpr, columns=columns, columns_type=columns_type)
    IDfile = ut.md5hashcode_from_string(filename, ncut=10)
    df["IDfile"] = IDfile

    return df, IDfile


# Base class
#############
class LogReader:
    """
    Basic class to read a log file

    """

    def __init__(self,
                 log_event: Union[str, List[str], Tuple[str]],
                 log_type: Union[str, List[str], Tuple[str]],
                 body_schema: Union[LogSchema, Dict[str, Tuple[Union[GTUnion, str], str]]],
                 event_description: str = "",
                 names: Optional[Union[FloatInt, List[FloatInt]]] = None,
                 IDs: Optional[Union[FloatInt, List[FloatInt]]] = None):
        """

        Parameters
        ----------
        log_event:
            name of the log event (or log events)
        log_type:
            type of log event, S for SSE event, B for BSE event, S|B for both
        body_schema:
            Schema of the log message
        names:
            if not None, when matching the log message just consider systems with that name(s)
        IDs:
            if not None, when matching the log message just consider systems with that ID(s)

        """

        # Check body schema
        if isinstance(body_schema, LogSchema):
            self._body_schema = body_schema
        elif isinstance(body_schema, dict):
            self._body_schema = LogSchema(body_schema)
        else:
            raise TypeError(f"body schema {body_schema} not allowed. It can be an"
                            f"instance of the class LogSchema or a Dictionary")

        # Log types
        _allowed_log_types = ("S", "B", "S|B", "B|S")
        if isinstance(log_type, str): log_type = [log_type, ]
        for _logt in log_type: check_allowed(_logt, _allowed_log_types)
        self.log_type = "|".join(log_type)

        # Log events
        if isinstance(log_event, str): log_event = [log_event, ]
        for _loge in log_event:
            if not isinstance(_loge, str): raise TypeError(
                f"log event {_loge} not allowed. All log events needs to be strings")
        self.log_event = "|".join(log_event)

        # Description
        self._event_description = event_description

        # Check names
        if names is None:
            _name = "name"
        else:
            try:
                if isinstance(names, str):
                    _name = str(names)
                elif isinstance(names, int):
                    _name = str(names)
                elif isinstance(names, Iterable):
                    _name = "|".join(map(str, names))
                else:
                    raise TypeError()
            except TypeError:
                raise TypeError(f"names {names} not allowed. Names can only be a string, an integer or"
                                f"a collection of strings/integers")

        # Check IDs
        if IDs is None:
            _id = "id"
        else:
            try:
                if isinstance(IDs, str):
                    _id = str(IDs)
                elif isinstance(IDs, int):
                    _id = str(IDs)
                elif isinstance(names, Iterable):
                    _id = "|".join(map(str, IDs))
                else:
                    raise TypeError()
            except TypeError:
                raise TypeError(f"names {IDs} not allowed. Names can only be a string, an integer or"
                                f"a collection of strings/integers")

        # NOTICE: the order is important, for Python>=3.7, the dict are all ordered
        self._header_schema = LogSchema({
            "logtype": (self.log_type, "log type (S for single evolution events, B for binary evolution events)"),
            "name": (_name, "name of the star/binary"),
            "ID": (_id, "id of the star/binary"),
            "event": (self.log_event, "event label"),
            "time": (float, "event time [Myr]")
        })

        # Define the dictionary containing the column schema
        self._column_schema = self._header_schema.column_schema(offset=0)
        self._column_schema.update(self._body_schema.column_schema(offset=len(self._header_schema)))

    # METHODS USEFUL FOR THE USERS
    def readfiles(self, logfiles: Union[str, List[str], Tuple[str]],
                  capturing_cols: Union[str, List[Union[str, int]], Tuple[Union[str, int], ...]] = "default",
                  nproc: int = 1) -> SEVNDataLog:
        """
        Basic method to read the logfile

        Parameters
        ----------
        logfiles:
            Path to the logfile to read or list containing the path of the Logfiles.
            It accepts wildcard *, e.g. logfile_* to list all the logfile in a folder.
        capturing_cols:
            Columns to capture
            If string it has to be a special value, otherwise a list of column names and/or indexes
            The special values are
            - "all": return all the column names
            - "default": return only the column names with a regex patter to catch a str an int or a float
            - "header": return only the colum names of the header
            It is possible to include duplicated values, but they will be neglected in the output
        nproc:
            number of processes to be used for a parallel reading on multiple files.

            .. note::
                The number of processes will be anyway limited to number of logfiles to be read

        Returns
        -------
        SEVNDataLog:
            an instance of the :class:`~sevnpy.io.sevndata.SEVNDataLog` containing the log information

        """

        # Check input
        if not isinstance(nproc, int) or nproc < 0: raise ValueError(
            f"nproc needs to be a positive integer, instead it is {nproc} (type: {type(nproc)}")

        # Get the list of computing names
        _capturing_names = self._generate_capturing_names(capturing_cols)
        # Define the common values
        matchexpr, captured_cols, capturing_types = self.matching_pattern(capturing_names=_capturing_names)

        # Check readfile input, use set to remove duplicates
        _logfiles = set()
        if isinstance(logfiles, str):
            # Assume the logfile could contain wild card * for filename expansion, so use glob
            _logfiles.update(glob.glob(logfiles))
        elif isinstance(logfiles, Iterable):
            for logfile in logfiles:
                if isinstance(logfile, str):
                    _logfiles.update(glob.glob(logfile))
                else:
                    ValueError(f"Input  {logfile} with type {type(logfile)} not allowed as filename/filepath")
        else:
            ValueError(f"Input  {logfiles} with type {type(logfiles)} not allowed as filename/filepath")

        # Finally get a list again considering absolute paths
        _logfiles = list(map(os.path.abspath, _logfiles))

        # Generate partial function
        get_df_from_regex = partial(_regex_dataframe_from_file_adding_filename, matchexpr=matchexpr,
                                    columns=captured_cols, columns_type=capturing_types)

        # Limit the number of processes to the number of files to read.
        # Additional processes are useless and add just additional overheads
        nproc = min(nproc, len(_logfiles))

        # Go
        if nproc == 1:
            dfl_plut_IDfile = list(map(get_df_from_regex, _logfiles))
        elif nproc > 0:
            with mp.Pool(nproc) as pool:
                dfl_plut_IDfile = pool.map(get_df_from_regex, _logfiles)
        else:
            raise ValueError

        # Generate df list and sevn_files dictionary containing pairs of hash_code:filename
        dfl = []
        sevn_files = {}
        # This work only beacuse map function (both from builtin and from miltiprocess) preserv the order
        for (df, IDfile), sevn_file in zip(dfl_plut_IDfile, _logfiles):
            dfl.append(df)
            sevn_files[IDfile] = sevn_file

        df = pd.concat(dfl)

        info_col = {}
        for name in captured_cols:
            if name in self._header_schema.names:
                info_col[name] = self._header_schema.description_schema[name]
            elif name in self._body_schema.names:
                info_col[name] = self._body_schema.description_schema[name]
        # Add the info about the extra column
        info_col["IDfile"] = "Unique ID defining the logfile from which the information in the row has been taken."

        return SEVNDataLog(dataframe=df,
                           data_info=self.event_info,
                           column_descriptions=info_col,
                           regex_pattern=matchexpr,
                           sevn_files=sevn_files)

    def readstring(self, string: str,
                   capturing_cols: Union[str, List[Union[str, int]], Tuple[Union[str, int], ...]] = "default", ) \
            -> pd.DataFrame:
        """
        Basic method to read logs from  a string

        Parameters
        ----------
        string:
            string containing the log
        capturing_cols:
            Columns to capture
            If string it has to be a special value, otherwise a list of column names and/or indexes
            The special values are
            - "all": return all the column names
            - "default": return only the column names with a regex patter to catch a str an int or a float
            - "header": return only the colum names of the header
            It is possible to include duplicated values, but they will be neglected in the output

        Returns
        -------
        LogDataFrame:
            a pandas Dataframe containing the info from the read logfile


        .. note::
            This method returns simply a pandas Dataframe instead of the instance of the class
            :class:`~sevnpy.io.sevndata.SEVNDataLog` returned by the method :func:`~LogReader.readfiles`.
            In order to access the additional information stored in The :class:`~sevnpy.io.sevndata.SEVNDataLog`
            class, save the string to a file and use :func:`~LogReader.readfiles`.

        """

        # Get the list of computing names
        _capturing_names = self._generate_capturing_names(capturing_cols)
        # Define the common values
        matchexpr, captured_cols, capturing_types = self.matching_pattern(capturing_names=_capturing_names)

        dfl = rut.regex_dataframe_from_string(string=string, matchexpr=matchexpr,
                                              columns=captured_cols, columns_type=capturing_types)

        return dfl

    def summary(self) -> Dict[str, Dict[str, str]]:
        """
        Get the total summary of the LogSchema

        Returns
        -------
        log_schema:
            Return a dictionary containing the complete information of the log-schema.
            Each item is a pair name:info_dictionary and info_dictionary is a dictionary structures as follow
            {<item_name>: {"name":<item_name>, "description":<item_description>, "kind":<item_kind>, "type":<item_type>, "pattern":<item_pattern>}}

        """
        return {**self._header_schema.full_schema(), **self._body_schema.full_schema()}

    def columns_summary(self, show_description: bool = False) -> Dict[int, str]:
        """
        Get the column summary including their sorted index

        Parameters
        ----------
        show_description (bool):
            If True, include all the column description in outpt

        Returns
        -------
        column_summary_dictionary:
            - If show_description is False: A dictionary with the pair index:name
            - If show_description is True: A dictionary with the pair index:(name,description)

        """

        if not show_description:
            return cp.deepcopy(self._column_schema)
        else:
            t = {}
            for idx, name in self._column_schema.items():
                t[idx] = (name, self.description[idx])

            return t

    def columnid_to_name(self, ids: Union[int, List[int]]) -> Union[int, List[str]]:

        """
        Get the column name based on their index

        Parameters
        ----------
        ids: str, List[str|int]
            Integer or list of integers

        Returns
        -------
        column_names: str, List[str|int]
            If the input is an integer return just the name of the column, otherwise a list with the column names corresponded to the indexes


        """
        if isinstance(ids, int): return self._column_schema[ids]

        return [self._column_schema[i] for i in ids]

    #######

    ##PROPERTIES
    # return the description
    @property
    def event_info(self):
        return self._event_description

    # return a order list of all the names
    @property
    def names(self):
        return self._header_schema.names + self._body_schema.names

    # return a order list of all the types
    @property
    def types(self):
        return self._header_schema.types + self._body_schema.types

    # return a order list of all the descriptions
    @property
    def description(self):
        return self._header_schema.descriptions + self._body_schema.descriptions

    ###############

    ####REGEX METHODS
    def header_regex(self, capturing_names: Union[Literal["default", "all"], List[str], Tuple[str]] = None) -> Tuple[
        str, Optional[List[str]]]:
        """

        Parameters
        ----------
        capturing_names (str or list): A string or list of capturing names to match. If set to "all",
                all column names will be matched. If set to "default" all column names will be matched,
                except for log_type and event. If a list is provided, only the specified capturing names will be matched.
                If None, none of the columns will be matched.
                Default value is "default".


        Returns
        -------
        header_regex_pattern:
            the regex pattern for the header.

        """
        sep = ";"
        if capturing_names is None:
            return sep.join(self._header_schema.patterns), None
        else:
            if capturing_names == "default":
                reg_pattern, _capturing_names = self._header_schema.regex_pattern(
                    capturing_names=["name", "ID", "time"])
            else:
                reg_pattern, _capturing_names = self._header_schema.regex_pattern(capturing_names)
            return sep.join(list(reg_pattern.values())), _capturing_names

    def body_regex(self, capturing_names: Union[Literal["default", "all"], List[str], Tuple[str]] = None) -> Tuple[
        str, Optional[List[str]]]:
        """

        Parameters
        ----------
        capturing_names (str or list): A string or list of capturing names to match. If set to "all" or "default"
                all column names will be matched.
                If a list is provided, only the specified capturing names will be matched.
                If None, none of the columns will be matched.
                Default value is "default".


        Returns
        -------
        header_regex_pattern:
            the regex pattern for the body.

        """
        sep = ":"
        if capturing_names is None:
            return sep.join(self._body_schema.patterns), None
        else:
            reg_pattern, _capturing_names = self._body_schema.regex_pattern(capturing_names)
            return sep.join(list(reg_pattern.values())), _capturing_names

    def matching_pattern(self,
                         capturing_names: Union[Literal["default", "all", "header_regex"], List[str]] = "default") -> \
            Tuple[Union[str, List[str]], Union[str, List[str]], Union[str, List[str]]]:
        """
        Return the regex matching pattern of capturing names in the header and body of the object.

        Parameters:
            capturing_names (str or list): A string or list of capturing names to match. If set to "all",
                all column names will be matched. If set to "default" all column names will be matched,
                except for log_type and event. If a list is provided, only the specified capturing names will be matched.
                Default value is "default".

        Returns:
            tuple: A tuple containing three elements:
                - A string representing the regex matching pattern.
                - A list of  capturing names from the header and body.
                - A list of types of the matched capturing names from the header and body.

        Raises:
            None.


        """

        if capturing_names == "default" or capturing_names == "all":
            header_capturing_names = body_capturing_names = capturing_names
        elif capturing_names == "header":
            header_capturing_names = "all"
            body_capturing_names = []
        else:
            # Filter header-body capture name
            header_capturing_names = []
            body_capturing_names = []
            for name in capturing_names:
                if name in self._header_schema.names:
                    header_capturing_names.append(name)
                elif name in self._body_schema.names:
                    body_capturing_names.append(name)
                else:
                    warnings.warn(f"capturing name {name} is not present neither in the header not in the body schema")

        match_header, capt_header = self.header_regex(header_capturing_names)
        match_body, capt_body = self.body_regex(body_capturing_names)

        capt_header_type = [self._header_schema.type_schema[name] for name in capt_header]
        capt_body_type = [self._body_schema.type_schema[name] for name in capt_body]

        return f"{match_header};{match_body}", capt_header + capt_body, capt_header_type + capt_body_type

    def _generate_capturing_names(self, capturing_cols: Union[str, List[Union[str, int]]]) -> List[str]:
        """
        Generate a list of names of column to be captured based on a mixed list of names and integers
        or special string values.

        Examples
        --------
        >>> dl=SEVNDataLog(log_event="WD",log_type="S",body_schema={"dummy":(str,"dummy variable")})
        >>> dl._generate_capturing_names("default")
        >>> dl._generate_capturing_names("all")
        >>> dl._generate_capturing_names(["name","ID"])
        >>> dl._generate_capturing_names([1,4])
        >>> dl._generate_capturing_names(["name",1,"id",2,4])


        Parameters
        ----------
        capturing_cols: str, List[str|int]
            If string it has to be a special value, otherwise a list of column names and/or indexes
            The special values are
            - "all": return all the column names
            - "default": return only the column names with a regex patter to catch a str an int or a float
            - "header": return only the colum names of the header
            It is possible to include duplicated values, but they will be neglected in the output

        Returns
        -------
        capturing_names:
            List of str containing the column names. If one of the names does not exist or one of the integer is out of range,
            it returns a ValueError. The returned list will not contain duplicated values even if they are inserted in input


        """
        _capturing_names = ""
        if capturing_cols == "default" or capturing_cols == "all":
            _capturing_names = capturing_cols
        elif capturing_cols == "header":
            _capturing_names = self._header_schema.names
        elif not isinstance(capturing_cols, str):
            _capturing_names = set()
            for element in capturing_cols:
                if isinstance(element, str):
                    _capturing_names.add(element)
                elif isinstance(element, int) and element >= 0 and element < len(self.names):
                    _capturing_names.add(self.columnid_to_name(element))
                elif isinstance(element, int):
                    raise ValueError(
                        f"Integer element {element} not allowed, it is negative or lager than the available names")
                else:
                    raise ValueError(f"Element {element} not allowed")
            _capturing_names = list(_capturing_names)
        else:
            raise ValueError(f"Capturing_cols {capturing_cols} not allowed")

        return _capturing_names

    #########

    ####SPECIAL METHODS
    def __str__(self):

        return self._header_schema.__str__() + self._body_schema.__str__()

    def __repr__(self):

        return self._header_schema.__repr__() + self._body_schema.__repr__()
    ####################


# Auxiliary classes
#############
class _SpecialisedLogReader(LogReader):
    """
    Specialised LogReader  for the SEVN events

    """

    log_event = "UNKOWN"
    log_type = None
    _body_schema = LogSchema({})
    event_description = ""

    def __init__(self,
                 names: Optional[Union[FloatInt, List[FloatInt]]] = None,
                 IDs: Optional[Union[FloatInt, List[FloatInt]]] = None):
        super().__init__(log_event=self.log_event,
                         log_type=self.log_type,
                         event_description=self.event_description,
                         body_schema=self._body_schema,
                         names=names,
                         IDs=IDs)


class _SSELogReader(_SpecialisedLogReader):
    """
    Special addition to the LogReader for the events regarding the SSE.
    It just define the log_type.

    """

    log_type = "S"


class _BSELogReader(_SpecialisedLogReader):
    """
    Special addition to the LogReader for the events regarding the BSE
    It just define the log_type.

    """
    log_type = "B"


class _OnlyOnePerSystemGetLast(_SpecialisedLogReader):
    """
    Special addition to the LogReader for events that need to
    be forced to  have at most  one occurrence for system.
    It overrides some methods and introduces some post-processing methods to remove duplicates
    and take just the last occurrence if there are more than one.


    """

    keep_mode = "last"

    def _prepare_capturing_cols(self, capturing_cols: Union[Literal["default", "all", "header"],
                                                            List[Union[str, int]], Tuple[Union[str, int], ...]]
    = "default") -> Union[str, List[str]]:
        """
        Auxiliary function to be sure to add the column name and ID in the capturing cols

        Parameters
        ----------
        capturing_cols:
            Columns to capture
            If string it has to be a special value, otherwise a list of column names and/or indexes
            The special values are
            - "all": return all the column names
            - "default": return only the column names with a regex patter to catch a str an int or a float
            - "header": return only the colum names of the header
            It is possible to include duplicated values, but they will be neglected in the output


        Returns
        -------
        mod_capturing_cols:
            the capturing_cols modified so that also the name and ID are considered

        """

        if capturing_cols == "default" or capturing_cols == "all" or capturing_cols == "header":
            _capturing_cols = capturing_cols
        elif "name" not in capturing_cols and "ID" not in capturing_cols:
            _capturing_cols = [*capturing_cols, "name", "ID"]
        elif "name" not in capturing_cols:
            _capturing_cols = [*capturing_cols, "name"]
        elif "ID" not in capturing_cols:
            _capturing_cols = [*capturing_cols, "ID"]
        else:
            _capturing_cols = capturing_cols

        return _capturing_cols

    def readfiles(self, logfiles: Union[str, List[str], Tuple[str]],
                  capturing_cols: Union[Literal["default", "all", "header"], List[Union[str, int]], Tuple[
                      Union[str, int], ...]] = "default",
                  nproc: int = 1) -> SEVNDataLog:
        """
        Same as the method in the base class (:func:`~LogReader.readfiles`), but dropping possible duplicates.
        This is used for Log entries that can produce multiple outputs do to check and repeat of SEVN,
        but only one should be considered in ouput (last by defaut, use the class _OnlyOnePerSystemGetFirst
        to maintain the first entry).

        Parameters
        ----------
        logfiles:
            Path to the logfile to read or list containing the path of the Logfiles
        capturing_cols:
            Columns to capture
            If string it has to be a special value, otherwise a list of column names and/or indexes
            The special values are
            - "all": return all the column names
            - "default": return only the column names with a regex patter to catch a str an int or a float
            - "header": return only the colum names of the header
            It is possible to include duplicated values, but they will be neglected in the output
        nproc:
            number of processes to be used for a parallel reading on multiple files.

            .. note::
                The number of processes will be anyway limited to number of logfiles to be read

        Returns
        -------
        SEVNDataLog:
            an instance of the :class:`~sevnpy.io.sevndata.SEVNDataLog` containing the log information

        """

        # Since we have to remove duplicates based on name and ID, add anyway the columns name and ID
        _capturing_cols = self._prepare_capturing_cols(capturing_cols)

        sevno = super().readfiles(logfiles=logfiles, capturing_cols=_capturing_cols, nproc=nproc)
        df = sevno.data
        df = df.drop_duplicates(subset=["name", "ID", "IDfile"], keep=self.keep_mode)

        # Now correct if necessary
        columns_info = sevno.columns_info
        regex_pattern = sevno.regex_pattern

        if capturing_cols != "default" and capturing_cols != "all" and capturing_cols != "header":

            if "name" not in capturing_cols:
                df = df.drop(columns=["name", ])
                columns_info.pop("name")
                regex_pattern_tmp = regex_pattern.split(";")
                # Name is the second element, split and  remove the capturing parenthesis
                regex_pattern_tmp[1] = f"(?:{regex_pattern_tmp[1][1:-1]})"
                regex_pattern = ";".join(regex_pattern_tmp)

            if "ID" not in capturing_cols:
                df = df.drop(columns=["ID", ])
                columns_info.pop("ID")
                regex_pattern_tmp = regex_pattern.split(";")
                # Name is the third element, split and remove the capturing parenthesis
                regex_pattern_tmp[2] = f"(?:{regex_pattern_tmp[2][1:-1]})"
                # recombine
                regex_pattern = ";".join(regex_pattern_tmp)

        return SEVNDataLog(dataframe=df,
                           data_info=sevno.data_info,
                           column_descriptions=columns_info,
                           regex_pattern=regex_pattern,
                           sevn_files=sevno.sevn_files)

    def readstring(self, string: str,
                   capturing_cols: Union[str, List[Union[str, int]], Tuple[Union[str, int], ...]] = "default", ) \
            -> pd.DataFrame:
        """
        Basic method to read logs from  a string

        Parameters
        ----------
        string:
            string containing the log
        capturing_cols:
            Columns to capture
            If string it has to be a special value, otherwise a list of column names and/or indexes
            The special values are
            - "all": return all the column names
            - "default": return only the column names with a regex patter to catch a str an int or a float
            - "header": return only the colum names of the header
            It is possible to include duplicated values, but they will be neglected in the output

        Returns
        -------
        LogDataFrame:
            a pandas Dataframe containing the info from the read logfile without duplicates

        .. note::
            This method returns simply a pandas Dataframe instead of the instance of the class
            :class:`~sevnpy.io.sevndata.SEVNDataLog` returned by the method :func:`~LogReader.readfiles`.
            In order to access the additional information stored in The :class:`~sevnpy.io.sevndata.SEVNDataLog`
            class, save the string to a file and use :func:`~LogReader.readfiles`.


        """

        # Since we have to remove duplicates based on name and ID, add anyway the columns name and ID
        _capturing_cols = self._prepare_capturing_cols(capturing_cols)
        # Get the data
        df = super().readstring(string=string, capturing_cols=_capturing_cols)
        # Remove duplicates
        df = df.drop_duplicates(subset=["name", "ID"], keep=self.keep_mode)

        # Now remove name and ID columns if they were not present among the original capturing_cols
        if capturing_cols != "default" and capturing_cols != "all" and capturing_cols != "header":

            if "name" not in capturing_cols:
                df = df.drop(columns=["name", ])
            if "ID" not in capturing_cols:
                df = df.drop(columns=["ID", ])

        return df


class _OnlyOnePerSystemGetFirst(_SpecialisedLogReader):
    """
    Special addition to the LogReader for events that need to
    be forced to  have at most  one occurrence for system.
    It overrides some methods and introduces some post-processing methods to remove duplicates
    and take just the first occurrence if there are more than one.


    """

    keep_mode = "first"


# Specialised classes
#############
####SUPER IMPORTANT: THE ORDER OF KEY IN THE LOGSCHEMA MATTERS. IT NEEDS TO FOLLOW STRICTLY THE ONE REPORTED IN THE SEVN LOG FILES
class SNLogReader(_OnlyOnePerSystemGetLast, _SSELogReader):
    """
    Specialised :class:`~LogReader`  for the event SN (see SEVN documentation)

    """

    log_event = "SN"
    event_description = "A star turns into a compact remnant"

    _body_schema = LogSchema({
        "Mtot_preSN": (float, "Stellar mass preSN [Msun]"),
        "MHE_preSN": (float, "Stellar He-core mass preSN[Msun]"),
        "MCO_preSN": (float, "Stellar CO-core preSNn [Msun]"),
        "SNMrem": (float, "Remnant mass  [Msun]"),
        "SNRemnantType": ("type", "SEVN Remant_type"),
        "SN_type": ("type", "SEVN SN_type"),
        "Vnk_SN": (float, "Natal kick magnitude [km/s]"),
        "Vfk_SN": (float, "Effective kick magnitude (after possibile corrections) [km/s]"),
        "Vfkx_SN": (float, "x-component of the effective kick [km/s] "
                           "(in binary x is the direction connecting the two stars, in single the direction is arbirary)"),
        "Vfky_SN": (float, "y-component of the effective kick [km/s] "
                           "(in binary y is perpendicular to the line connectinc the two stars"
                           "the postive direction is toward the preSN relative motion, in single the direction is arbirary"),
        "Vfkz_SN": (float, "z-component of the effective kick [km/s]"
                           "(in binary z is the direction normal to the  orbital plane, in single the direction is arbirary)"),

    })

    _SN_type_map = {
        0: "Unkown",
        1: "ECSN",
        2: "CCSN",
        3: "PPISN",
        4: "PISN",
        5: "Ia"
    }

    @classmethod
    def SN_type(cls, SN_type: int) -> str:
        """
        Get the str representation correspondent to the SN_type ID (see SEVN documentation)

        Parameters
        ----------
        SN_type:
            SN type ID (see SEVN documentation)

        Returns
        -------
        SN_type_str:
            str representation of the SN type (see SEVN documentation)


        """

        return cls._SN_type_map[SN_type]

class QHELogReader(_SSELogReader):
    """
    Specialised :class:`~LogReader`  for the event QHE (see SEVN documentation)

    """

    log_event = "QHE"
    event_description = "a stars starts to follow the Quasi-Homogeneous evolution"

    _body_schema = LogSchema({
        "Mass_QHE": (float, "total mass of the Star at the beginning of the QHE [Msun]"),
        "Radius_QHE": (float, "stellar radius at the beginning of the QHE [Rsun]"),
        "dMcumul_RLO_QHE": (float, "amount of mass accumulated during the ongoing RLO (if any) [Msun]"),
    })

class WDLogReader(_OnlyOnePerSystemGetLast, _SSELogReader):
    """
    Specialised :class:`~LogReader`  for the event WD (see SEVN documentation)

    """

    log_event = "WD"
    event_description = "A White Dwarf is formed."

    _body_schema = LogSchema({
        "Mtot_preWD": (float, "Stellar mass preWD formation [Msun]"),
        "MHE_preWD": (float, "Stellar He-core mass preWD formation [Msun]"),
        "MCO_preWD": (float, "Stellar CO-core preWD formation [Msun]"),
        "WDMass": (float, "Final WD mass [Msun]"),
        "WDRemnant_type": ("type", "SEVN Remant_type")
    })


class NSLogReader(_OnlyOnePerSystemGetLast, _SSELogReader):
    """
    Specialised :class:`~LogReader`  for the event NS (see SEVN documentation)

    """
    log_event = "NS"
    event_description = "A NS is formed."

    _body_schema = LogSchema({
        "NSRemnant_type": ("type", "SEVN Remant_type"),
        "NSMass": (float, "NS Mass [Msun]"),
        "NSBmag_natal": (float, "NS Magnetic field at birth [Gauss]"),
        "NSSpin_natal": (float, "NS rotation frequency at birth [1/s]"),
        "NSsina_natal": (float, "sin of the angle between the magnetic axis and the rotation axis at birth")
    })


class HENAKEDLogReader(_SSELogReader):
    """
    Specialised :class:`~LogReader`  for the event HENAKED (see SEVN documentation)

    """
    log_event = "HENAKED"
    event_description = "The envelope of a star is totally stripped and it turns into a pureHe star"

    _body_schema = LogSchema({
        "Mass_preHE": (float, "Stellar mass before becoming HE naked [Msun]"),
        "MHE_preHE": (float, "Stellar HE-core mass before becoming HE naked [Msun]"),
        "MCO_preHE": (float, "Stellar CO-core mass before becoming HE naked [Msun]"),
        "Radius_preHE": (float, "Stellar radius before becoming HE naked [Rsun]"),
        "Radius_postHE": (float, "Stellar radius of the  HE naked star [Rsun]"),
        "Phase_postHE": ("type", "SEVN phase of the HE naked star"),
        "Minterp_preHE": (float, "Mzams of the fake track used for interpolation before becoming HE naked [Msun]"),
        "Minterp_postHE": (float, "Mzams of the  pureHE fake track used for interpolation of the HE naked [Msun]"),
        "Ctrack_outcome": ("type",
                           "Change of track outcome 1:track found, convergence reached, 2:track found, convergence not reached, 3:trak not found")
    })


class CONAKEDLogReader(_SSELogReader):
    """
    Specialised :class:`~LogReader`  for the event CONAKED (see SEVN documentation)

    """

    log_event = "CONAKED"
    event_description = "Bothe Hydrogen and the Helium  envelope of a star is totally stripped and it turns into a nakedCO star"

    _body_schema = LogSchema({
        "Mass_preCO": (float, "Stellar mass before becoming CO naked [Msun]"),
        "MHE_preCO": (float, "Stellar HE-core mass before becoming CO naked [Msun]"),
        "MCO_preCO": (float, "Stellar CO-core mass before becoming CO naked [Msun]"),
        "Radius_preCO": (float, "Stellar radius before becoming CO naked [Rsun]"),
        "RCO_preCO": (float, "Stellar radius of the CO core before becoming CO naked [Rsun]"),
        "Mass_postCO": ("type", "SEVN phase of the HE naked star [Msun]"),
        "Radius_postCO": (float, "Mzams of the fake track used for interpolation before becoming HE naked [Rsun]"),
        "Phase_postCO": (float, "Mzams of the  pureHE fake track used for interpolation of the HE naked"),
    })


class BSNLogReader(_BSELogReader):
    """
    Specialised :class:`~LogReader`  for the event BSN (see SEVN documentation)

    """
    log_event = "BSN"
    event_description = "One of the star in a binary turns into a compact remnant"

    _body_schema = LogSchema({
        "SID_SN": ("id", "Internal ID (0 or 1) of the star exploding as SN"),

        "Mass_preSN_SN": (float, "Stellar mass of the exploding star (ID=SNSID) at the onset of the SN [Msun]"),
        "MHE_preSN_SN": (float, "Stellar He core mass of the exploding star (ID=SNSID) at the onset of the SN [Msun]"),
        "MCO_preSN_SN": (float, "Stellar CO core mass of the exploding star (ID=SNSID) at the onset of the SN [Msun]"),
        "Phase_preSN_SN": ("type", "SEVN stellar phase  of the exploding star (ID=SNSID) at the onset of the SN"),
        "RemnantType_postSN_SN": (
            "type", "SEVN stellar remnant type phase  of the exploding star  (ID=SNSID) after the SN"),

        "SID_Other": ("id", "Internal ID (0 or 1) of the companionf of the star exploding as SN"),
        "Mass_preSN_Other": (float, "Stellar mass of the companion star  (ID=OTSID)  at the onset of the SN [Msun]"),
        "MHE_preSN_Other": (
            float, "Stellar He core mass of the companion star  (ID=OTSID)  at the onset of the SN [Msun]"),
        "MCO_preSN_Other": (
            float, "Stellar CO core mass of the companion star (ID=OTSID)  at the onset of the SN [Msun]"),
        "Phase_preSN_Other": ("type", "SEVN stellar phase  of the companion star (ID=OTSID)  at the onset of the SN"),
        "RemnantType_postSN_Other": (
            "type", "SEVN stellar remnant type phase  of the (ID=OTSID)  companion star after the SN"),

        "Semimajor_preSN": (float, "Semimajor-axis at the onset of the SN [Rsun]"),
        "Eccentricity_preSN": (float, "Eccentricity at the onset of the SN"),
        "Semimajor_postSN": (float, "Semimajor-axis after the SN [Rsun]"),
        "Eccentricity_postSN": (float, "Eccentricity at the onset the SN [Rsun]"),

        "CosNu": (float, "Cosine of the angle between the normals to the preSN and postSN orbital planes"),
        "Vcom": (
            float, "Module of the difference between the binary centre of mass velocity before and after the SN [km/s]")
    })


class RLO_BEGINLogReader(_BSELogReader):
    """
    Specialised :class:`~LogReader`  for the event RLO_BEGIN (see SEVN documentation)

    """
    log_event = "RLO_BEGIN"
    event_description = "A Roche-Lobe Overflow begins"

    _body_schema = LogSchema({
        "SID_RLOBd": ("id", "Internal ID (0 or 1) of the donor star, i.e. the star filling the Roche-Lobe"),
        "Mass_RLOBd": (float, "Stellar mass of the donor star (ID=SID_RLOBd) at the onset of RLO [Msun]"),
        "MHE_RLOBd": (float, "Stellar He core mass of the donor star (ID=SID_RLOBd) at the onset of RLO [Msun]"),
        "MCO_RLOBd": (float, "Stellar CO core mass of the donor star (ID=SID_RLOBd) at the onset of RLO [Msun]"),
        "Phase_RLOBd": ("type", "SEVN stellar phase  of the donor star (ID=SID_RLOBd) at the onset of RLO"),
        "RemnantType_RLOBd": (
            "type", "SEVN stellar remnant type phase  of the donor star (ID=SID_RLOBd) at the onset of RLO"),
        "Radius_RLOBd": (float, "Stellar radius of the donor star (ID=SID_RLOBd) at the onset of RLO [Rsun]"),
        "RL_RLOBd": (float, "Roche-Lobe radius of the donor star (ID=SID_RLOBd) at the onset of RLO [Rsun]"),

        "SID_RLOBa": ("id", "Internal ID (0 or 1) of the accretor star"),
        "Mass_RLOBa": (float, "Stellar mass of the accretor star (ID=SID_RLOBa)  at the onset of RLO [Msun]"),
        "MHE_RLOBa": (float, "Stellar He core mass of the accretor star (ID=SID_RLOBa)  at the onset of RLO [Msun]"),
        "MCO_RLOBa": (float, "Stellar CO core mass of the accretor star (ID=SID_RLOBa)  at the onset of RLO [Msun]"),
        "Phase_RLOBa": ("type", "SEVN stellar phase  of the accretor star (ID=SID_RLOBa)  at the onset of RLO"),
        "RemnantType_RLOBa": ("type", "SEVN stellar remnant of the accretor star (ID=SID_RLOBa)  at the onset of RLO"),
        "Radius_RLOBa": (float, "Stellar radius   of the accretor star star (ID=SID_RLOBa) at the onset of RLO [Rsun]"),
        "RL_RLOBa": (float, "Roche-Lobe radius   of the accretor star (ID=SID_RLOBa)   at the onset of RLO [Rsun]"),

        "q_RLOB": (float, "Mass ratio (donor/accretor) at  the onset of RLO"),
        "qcrit_RLOB": (float, "Critical mass ratio at  the onset of RLO"),

        "Semimajor_RLOB": (float, "Semimajor-axis at the onset of RLO [Rsun]"),
        "Eccentricity_RLOB": (float, "Eccentrcitiy the onset of RLO [Rsun]"),
    })


class RLO_ENDLogReader(_BSELogReader):
    """
    Specialised :class:`~LogReader`  for the event RLO_END (see SEVN documentation)

    """
    log_event = "RLO_END"
    event_description = "Roche-Lobe overflow ends"

    _body_schema = LogSchema({
        "SID_RLOEd": ("id", "Internal ID (0 or 1) of the donor star, i.e. the star filling the Roche-Lobe"),
        "Mass_RLOEd": (float, "Stellar mass of the donor star (ID=SID_RLOEd) at the end of RLO [Msun]"),
        "MHE_RLOEd": (float, "Stellar He core mass of the donor star (ID=SID_RLOEd) at the end of RLO [Msun]"),
        "MCO_RLOEd": (float, "Stellar CO core mass of the donor star (ID=SID_RLOEd) at the end of RLO [Msun]"),
        "Phase_RLOBEd": ("type", "SEVN stellar phase  of the donor star (ID=SID_RLOEd) at the end of RLO"),
        "RemnantType_RLOEd": ("type", "SEVN stellar remnant type  of the donor star (ID=SID_RLOEd) at the end of RLO"),
        "Radius_RLOEd": (float, "Stellar radius of the donor star (ID=SID_RLOEd) at the end of RLO [Rsun]"),
        "RL_RLOEd": (float, "Roche-Lobe radius of the donor star (ID=SID_RLOEd) at the end of RLO [Rsun]"),

        "SID_RLOEa": ("id", "Internal ID (0 or 1) of the accretor star"),
        "Mass_RLOEa": (float, "Stellar mass of the accretor star (ID=SID_RLOEa)  at the end of RLO [Msun]"),
        "MHE_RLOEa": (float, "Stellar He core mass of the accretor star (ID=SID_RLOEa)  at the end of RLO [Msun]"),
        "MCO_RLOEa": (float, "Stellar CO core mass of the accretor star (ID=SID_RLOEa)  at the end of RLO [Msun]"),
        "Phase_RLOEa": ("type", "SEVN stellar phase  of the accretor star (ID=SID_RLOEa)  at the end of RLO"),
        "RemnantType_RLOEa": (
            "type", "SEVN stellar remnant type of the accretor star (ID=SID_RLOEa)  at the end of RLO"),
        "Radius_RLOEa": (float, "Stellar radius   of the accretor star star (ID=SID_RLOEa) at the end of RLO [Rsun]"),
        "RL_RLOEa": (float, "Roche-Lobe radius   of the accretor star (ID=SID_RLOEa)   at the end of RLO [Rsun]"),

        "Mlost_RLOE": (float, "Total mass lost by the donor (ID=SID_RLOEd)  due to the RLO [Msun]"),
        "Maccreted_RLOE": (float, "Total mass accreted by the accretor (ID=SID_RLOEa)  due to the RLO [Msun]"),

        "Semimajor_RLOB": (float, "Semimajor-axis at the end of RLO [Rsun]"),
        "Eccentricity_RLOB": (float, "Eccentrcitiy the end of RLO [Rsun]"),
    })


class CIRCLogReader(_BSELogReader):
    """
    Specialised :class:`~LogReader`  for the event CIRC (see SEVN documentation)

    """
    log_event = "CIRC"
    event_description = "The binary orbit is forced to be circular"

    _body_schema = LogSchema({
        "Semimajor_preCIRC": (float, "Semimajor axis at the onset of the circularisation [Rsun]"),
        "Eccentricity_preCirc": (float, "Eccentricity at the onset of the circularisation"),
        "Semimajor_postCIRC": (float, "Semimajor axis after the circularisation [Rsun]"),
        "Eccentricity_postCirc": (float, "Eccentricity  after the circularisation"),
    })


class CELogReader(_BSELogReader):
    """
    Specialised :class:`~LogReader`  for the event CE (see SEVN documentation)

    """
    log_event = "CE"
    event_description = ""

    _body_schema = LogSchema({

        "SID_CEp": ("id", "Internal ID (0 or 1) of the CE primary star"),
        "Mass_CEp": (float, "Stellar mass of the primary star (ID=SID_CEp) at the onset of the CE [Msun]"),
        "MHE_CEp": (float, "Stellar He core mass of the primary star (ID=SID_CEp) at the onset of the CE [Msun]"),
        "MCO_CEp": (float, "Stellar CO core mass of the primary star (ID=SID_CEp) at the onset of the CE [Msun]"),
        "Phase_CEp": ("type", "SEVN Phase  of the primary star (ID=SID_CEp) at the onset of the CE"),
        "RemnantType_CEp": ("type", "SEVN Remnant type  of the primary star (ID=SID_CEp) at the onset of the CE"),

        "SID_CEs": ("id", "Internal ID (0 or 1) of the CE secondary star"),
        "Mass_CEs": (float, "Stellar mass of the secondary star (ID=SID_CEs) at the onset of the CE [Msun]"),
        "MHE_CEps": (float, "Stellar mass of the secondary star (ID=SID_CEs) at the onset of the CE [Msun]"),
        "MCO_CEs": (float, "Stellar CO core mass of the secondary star (ID=SID_CEs) at the onset of the CE [Msun]"),
        "Phase_CEs": ("type", "SEVN Phase  of the secondary star (ID=SID_CEs) at the onset of the CE"),
        "RemnantType_CEs": ("type", "SEVN Remnant type  of the secondary star (ID=SID_CEs) at the onset of the CE"),

        "Semimajor_preCE": (float, "Semimajor axis at the onset of the CE [Rsun]"),
        "Semimajor_postCE": (float, "Semimajor axis at the end of the CE [Rsun]"),

        "CEoutcome": ("id", "CE outomce: 0: binary survives, 1:stars coalesce")

    })


# @TODO We have to update the SEVN guide since it is still calling the primary as the one that fill the RL
class COLLISIONLogReader(_BSELogReader):
    """
    Specialised :class:`~LogReader`  for the event COLLISION (see SEVN documentation)

    """
    log_event = "COLLISION"
    event_description = "The sum of the stellar radii is smaller than the periastron distance"

    _body_schema = LogSchema({

        "SID_COLp": ("id", "Internal ID (0 or 1) of the primary star"),
        "Mass_COLp": (float, "Stellar mass of the primary star (ID=SID_COLp) at the onset of the Collision [Msun]"),
        "Radius_COLp": (float, "Stellar radius of the primary star (ID=SID_COLp) at the onset of the Collision[Rsun]"),
        "Phase_COLp": ("type", "SEVN Phase  of the primary star (ID=SID_COLp) at the onset of the Collision"),

        "SID_COLs": ("id", "Internal ID (0 or 1) of the secondary star"),
        "Mass_COLs": (float, "Stellar mass of the secondary star(ID=SID_COLs) at the onset of the Collision [Msun]"),
        "Radius_COLs": (float, "Stellar radius of the secondary star(ID=SID_COLs) at the onset of the Collision[Rsun]"),
        "Phase_COLs": ("type", "SEVN Phase  of the secondary star (ID=SID_COLs) at the onset of the Collision"),

        "Semimajor_preCOL": (float, "Semimajor axis at the onset of the Collision [Rsun]"),
        "Eccentricity_preCOL": (float, "Semimajor axis at the end of the Collision [Rsun]"),

        "RL_COLp": (float, "Roche-Lobe radius of the primary star (ID=SID_COLp) at the onset of the Collision[Rsun]"),
        "RL_COLs": (float, "Roche-Lobe radius of the secondary star (ID=SID_COLs) at the onset of the Collision[Rsun]")

    })


class MERGERLogReader(_OnlyOnePerSystemGetLast, _BSELogReader):
    """
    Specialised :class:`~LogReader`  for the event MERGER (see SEVN documentation)

    """
    log_event = "MERGER"
    event_description = "Two stars merge"

    _body_schema = LogSchema({

        "SID_MERGa": ("id", "Internal ID (0 or 1) of the star that survive after the merger (the accretor)"),
        "Mass_MERGa": (float, "Stellar mass of the surviving star (ID=SID_MERGa) at the onset of the merger [Msun]"),
        "MHE_MERGa": (
            float, "Stellar He core mass of the surviving star (ID==SID_MERGa) at the onset of the merger  [Msun]"),
        "MCO_MERGa": (
            float, "Stellar CO core mass the surviving star (ID==SID_MERGa) at the onset of the merger  [Msun]"),
        "Phase_MERGa": ("type", "SEVN Phase  of the surviving star (ID==SID_MERGa) at the onset of the merger "),
        "RemnantType_MERGa": (
            "type", "SEVN Remnant type  of the surviving star(ID==SID_MERGa) at the onset of the merger "),
        "Radius_MERGa": (
            float, "Stellar radius of the surviving star (ID==SID_MERGa) at the onset of the merger  [Rsun]"),

        "SID_MERGd": ("id", "Internal ID (0 or 1) of the star that will be removed after the merger (the donor)"),
        "Mass_MERGd": (
            float, "Stellar mass of the star that will be removed (ID=SID_MERGd) at the onset of the merger  [Msun]"),
        "MHE_MERGd": (
            float, "Stellar He core mass  of the star that will be removed (ID=SID_MERGd) "
                   "at the onset of the merger  [Msun]"),
        "MCO_MERGd": (float,
                      "Stellar CO core mass  of the star that will be removed (ID=SID_MERGd) "
                      "at the onset of the merger  [Msun]"),
        "Phase_MERGd": (
            "type", "SEVN Phase   of the star that will be removed (ID=SID_MERGd) at the onset of the merger "),
        "RemnantType_MERGd": (
            "type", "SEVN Remnant type   of the star that will be removed (ID=SID_MERGd) at the onset of the merger "),
        "Radius_MERGd": (
            float, "Stellar radius of the star that will be removed (ID=SID_MERGd) at the onset of the merger  [Rsun]"),

        "Mass_postMERG": (float, "Total mass of the final merger product [Msun]"),

        "Semimajor_preMERG": (float, "Semimajor axis at the onset of the Merger [Rsun]"),
        "Eccentricity_preMERG": (float, "Eccentricity   at the onset of the Merger  [Rsun]"),

    })


class SWALLOWEDLogReader(_BSELogReader):
    """
    Specialised :class:`~LogReader`  for the event SWALLOWED (see SEVN documentation)

    """
    log_event = "SWALLOWED"
    event_description = "During an interaction one of the star is swallowed and totally destroyed by the companion"

    _body_schema = LogSchema({

        "SID_SWd": ("id", "ID of the destroyed star"),
        "Mass_SWd": (float, "Stellar mass of the destroyed star [Msun]"),
        "MHE_SWd": (float, "Stellar He core mass of the destroyed star [Msun]"),
        "MCO_SWd": (float, "Stellar CO core mass of the destroyed star [Msun]"),
        "Phase_SWd": ("type", "SEVN Phase  of the destroyed star (ID==SID_SWd)"),
        "RemnantType_SWd": ("type", "SEVN Remnant type of the destroyed star (ID==SID_SWd)"),
        "SID_SWa": ("id", "ID of the destroyed star"),
        "Mass_SWa": (float, "Stellar mass of the other (swallowing) star [Msun]"),
        "MHE_SWa": (float, "Stellar He core mass of the other (swallowing) star [Msun]"),
        "MCO_SWa": (float, "Stellar CO core mass of the other (swallowing) star [Msun]"),
        "Phase_SWa": ("type", "SEVN Phase  of the other (swallowing) star (ID==SID_SWa)"),
        "RemnantType_SWa": ("type", "SEVN Remnant type of the other (swallowing) star (ID==SID_SWa)"),
        "Macc_SWa": (float, "Total accreted mass by the swallowing star (ID==SID_SWa) from the"
                            "swallowed destroyed star"),

    })


# Useful functions and attribute
#################
_logevents = {
    "SN": SNLogReader,
    "QHE": QHELogReader,
    "WD": WDLogReader,
    "NS": NSLogReader,
    "HENAKED": HENAKEDLogReader,
    "CONAKED": CONAKEDLogReader,
    "BSN": BSNLogReader,
    "RLO_BEGIN": RLO_BEGINLogReader,
    "RLO_END": RLO_ENDLogReader,
    "CIRC": CIRCLogReader,
    "CE": CELogReader,
    "COLLISION": COLLISIONLogReader,
    "MERGER": MERGERLogReader,
    "SWALLOWED": SWALLOWEDLogReader,
}


def logevents_list() -> List[str]:
    """
    Return the list of the available log events

    Returns
    -------
    logevent_list: List
        The list of the available log events

    """

    return list(_logevents.keys())


def logcolumns_info(event: str, show_description: bool = False) -> Union[Dict[int, str], Dict[int, Tuple[str, str]]]:
    """
    Return the info about the values stored in a given event

    Parameters
    ----------
    event:
        Name of the event
    show_description:
        If True show also the description of that column, if False only the data type

    Returns
    -------
    logcolum_info: Dictionary
        a dictionary containg  the pair property_name:type of property_name:(type,description)
        if show_description is True

    Raises
    -------
    KeyError:
        If the event in input is not present in the available event list

    """

    if event not in _logevents: raise KeyError(
        f"Event {event} not in the event list. Avaialable events are {events_list()}")

    return _logevents[event]().columns_summary(show_description=show_description)


def readlogfiles(logfiles: Union[str, List[str], Tuple[str]],
                 events: Union[str, Union[List[str], Tuple[str]]],
                 capturing_cols: Union[str, List[Union[str, int]], Tuple[str]] = "default",
                 nproc: int = 1,
                 names: Optional[Union[FloatInt, List[FloatInt]]] = None,
                 IDs: Optional[Union[FloatInt, List[FloatInt]]] = None) -> Union[SEVNDataLog, Dict[str, SEVNDataLog]]:
    """
    Read  a series of logfiles retrieving the information for the events in input

    Parameters
    ----------
    logfiles:
        Path to the logfile or a list of logfiles paths.
        It accepts wildcard *, e.g. logfile_* to list all the logfile in a folder.
    events:
        label of the events to catch, call :func:`~logevents_list` to get the list of available events.
        It could be both a single string with an event label, or a list of label strings
    capturing_cols:
            Columns to capture
            If string it has to be a special value, otherwise a list of column names and/or indexes
            The special values are:

                - "all": return all the column names
                - "default": return only the column names with a regex patter to catch a str an int or a float
                - "header": return only the colum names of the header

            It is possible to include duplicated values, but they will be neglected in the output
    nproc:
            number of processes to be used for a parallel reading on multiple files.

            .. note::
                The number of processes will be anyway limited to number of logfiles to be read
    names:
        List of names to consider in the logfile. If None do not filter the names.
    IDs:
        List of IDs to consider in the logfile. If None do not filter the IDs.

    Returns
    -------
    log_info : SEVNDataLog or Dictionary
        - If the *events* input is a string returns a single instance of :class:`~sevnpy.io.sevndata.SEVNDataLog`
        - If the *events* input is a list of event labels, return a dictionary containing pairs
          event_label::class:`~sevnpy.io.sevndata.SEVNDataLog`


    Raises
    -------
    KeyError:
        If the event in input is not present in the available event list



    Examples
    --------
    Read a single logfile looking for the events SN and RLO_BEGIN

    >>> readlogfiles("logfile_0.dat",events=["SN","RLO_BEGIN"])

    Read all the logfiles in the folder using wildcards

    >>> readlogfiles("logfile_*.dat",events=["SN","RLO_BEGIN"])

    """

    if isinstance(events, str): events = [events, ]

    # Check that all the events are legit before the proper reading
    for event in events:
        if event not in _logevents: raise KeyError(
            f"Event {event} not in the event list. Available events are {logevents_list()}")

    # Reading
    output = {event: _logevents[event](names=names, IDs=IDs).readfiles(logfiles,
                                                                       capturing_cols=capturing_cols,
                                                                       nproc=nproc)
              for event in events}

    # If the output contains just one object just return the object
    if len(output) == 1:
        return output[events[0]]

    return output


def readlogstring(string: Union[str, List[str], Tuple[str]],
                  events: Union[str, Union[List[str], Tuple[str]]],
                  capturing_cols: Union[str, List[Union[str, int]], Tuple[str]] = "default",
                  names: Optional[Union[FloatInt, List[FloatInt]]] = None,
                  IDs: Optional[Union[FloatInt, List[FloatInt]]] = None) -> Union[
    pd.DataFrame, Dict[str, pd.DataFrame]]:
    """
    Read  a log string  retrieving the information for the events in input

    Parameters
    ----------
    string:
        string containing the log
    capturing_cols:
        Columns to capture
        If string it has to be a special value, otherwise a list of column names and/or indexes
        The special values are
        - "all": return all the column names
        - "default": return only the column names with a regex patter to catch a str an int or a float
        - "header": return only the colum names of the header
        It is possible to include duplicated values, but they will be neglected in the output

    Returns
    -------
    log_info : DataFrame or Dictionary
        - If the *events* input is a string return a pandas DataFrame with the log info
        - If the *events* input is a list of event labels, return a dictionary of pairs
          labe_name:DataFrame with the log info

    .. note::
        This method returns simply a pandas Dataframe instead of the instance of the class
        :class:`~sevnpy.io.sevndata.SEVNDataLog` returned by the method :func:`~LogReader.readfiles`.
        In order to access the additional information stored in The :class:`~sevnpy.io.sevndata.SEVNDataLog`
        class, save the string to a file and use :func:`~LogReader.readfiles`.


    Raises
    -------
    KeyError:
        If the event in input is not present in the available event list

    """

    if isinstance(events, str): events = [events, ]

    # Check that all the events are legit before the proper reading
    for event in events:
        if event not in _logevents: raise KeyError(
            f"Event {event} not in the event list. Avaialable events are {logevents_list()}")

    # Reading
    output = {event: _logevents[event](names=names, IDs=IDs).readstring(string=string,
                                                                        capturing_cols=capturing_cols)
              for event in events}

    # If the output contains just one object just return the object
    if len(output) == 1:
        return output[events[0]]

    return output
