"""
==============================================================
SEVNmanager , (:mod:`sevnpy.sevn.sevnmanager`)
==============================================================

This module contains the singleton static class SEVNmanager that is used to handle all the necessary
connections with the SEVN C++ backend

"""

import numpy as np
import copy

try:
    from . import sevnwrap as sw
except:
    raise ImportError("The sevnwrap is not installed")

from ..utility._privateutility import RestrictiveSingleton
from .. import utility as ut
from ..sevnpy_types import ListLikeType, StandardValue, Dict, List,Tuple,Set, Any, Optional, Union


class SEVNmanager(metaclass=RestrictiveSingleton):
    """
    This class is a Singleton (but not a classical one, read below) almost pure static container,
    and it controls all the connections with the
    SEVN C++ back-end. It  is possible to create only one instance of this class, but anyway this is useful only for
    initialising a context manager, in all the other cases only the static methods are used.
    Therefore, there will be only one common SEVNmanager for the whole life of  a process (multiple processes will have
    multiple SEVNmanagers), however the same SEVNmanager can be initialised (and closed) more tha once with different
    parameters. Each initialisation produce a new load of the tables (even if th tables option have not changed) and
    update an internal sequential ID (that can be retrieved with
    the static method :func:`~SEVNmanager.get_ID`. An internal dictionary store
    the parameters used for each ID (can be retrieved with
    the static method :func:`~SEVNmanager.get_sevnParams_history`).
    The SEVN parameters can be retrived with the static methods :func:`~SEVNmanager.get_sevnParamsDefault`
    The static method :func:`~SEVNmanager.sevnparam_describe`can be used to get a short description
    of the SEVN parameters. The static method sevn_path can be used to get the SEVN path where the libraries are
    located. The static method  :func:`~SEVNmanager.tables_info` return the info of the loaded tables
    (only avaible after at least one initilisation).


    .. note::
        The SEVN parameters in the SEVNmanager does not include all the IO parameters from SEVN since some of them
        are just used in the main SEVN excutables (e.g. io_logfile, scol, bcol,... ecc). The full list of the SEVN
        parameters can be retrieved with the static method  :func:`~SEVNmanager.get_sevnparams_raw`


    **How to initiliase a SEVN session**

    The SEVNmanager is necessary to use the SEVN backend and can be used in two ways.

    **Option-1: Init/Close**
        In this case the static method  :func:`~SEVNmanager.init`
        shoud be called at the beinning of the part of the code
        using the SEVN methods, and it needs to be closed at the end.

        >>> SEVNmanager.init() #default parameters
        >>> ..... rest of the code .....
        >>> SEVNmanager.close() #default parameters

        Calling :func:`~SEVNmanager.init`  without parameters will initiliase
        the SEVN parameters with the default value, if one
        wants to change some (of all the parameters) this can be given to init using a Dictionary with pairs
        param_name:value, e.g. if we want to change the tables and the common envelope parameter alpha.

        >>> new_param = {"tables":"abs_path_to_tables", "ce_alpha":5}
        >>> SEVNmanager.init(new_param) #custom parameters for tables and ce_alpha, default for all the others
        >>> ..... rest of the code .....
        >>> SEVNmanager.close() #default parameters

        In this case all the other parameters not included in new_param will have  their default value

        Warning
        ---------
        After calling :func:`~SEVNmanager.close`, the SEVN parameters will be restored
        to their default values

    **Option-2: Context manager/Close**
        Alternitavely to che Init/Close case, a SEVNmanager instance can be used to create a context manager (actually
        this is the only use of a SEVNmanager instance).
        All the code that use the SEVN functions needs to be included within the context manager.
        At the exit of the context manager the SEVNmanager instance will be automatically closed and cleaned, so the
        user does not need to worry about it.

        Example:

        >>> with SEVNmanager() as sm:
        >>>     ..... rest of the code .....

        Same as the first option, the empy initilisation will use the default parameters, if we want to change some of
        them we can use a dictionary, e.g.

        >>> new_param = {"tables":"abs_path_to_tables", "ce_alpha":5}
        >>> with SEVNmanager(new_param) as sm:
        >>>     ..... rest of the code .....

    .. warning::
        In order to have more readable codes, please  use just one of the two options and do not mix them in a single
        script/process

    .. note::
        SEVNmanager can work without never instantiating a class just using static methods.
        The users should never directly create a SEVNmanager instance, they are only created (and destroyed)
        when SEVNmanager is used in a contex manager (see below).
        Anyway, when an instance is created, SEVNmanager  is also automatically initialised, e.g.

        >>> sm = SEVNmanager() #Create a SEVNmanager instance and initialise with default SEVN parameters

        or

        >>> sm = SEVNmanager({"ce_alpha":2}) #Create a SEVNmanager instance and initialise it with defatul SEVN parameters except for ce_alpha

        When an instance is created, it will added to a private dictionary and any other try to instanciate another
        object will result in a RuntimeError. When the instance is closed, the private dictionary will be cleared
        and a new instance can be generated.

        TODO: This is not true anymore, if a new SEVNmanager is called the old one is automatically closed

        >>> sm=SEVNmanager()
        >>> sm2=SEVNmanager() #!RuntimeError, sm has not been closed

        >>> sm=SEVNmanager()
        >>> sm.close()
        >>> sm2=SEVNmanager() #!OK, sm has  been closed

    ------


    """

    _static_open_flag = False
    _static_close_flag = True
    _param_to_not_include = ("initerror_stop", "io_logfile", "Z", "bcol", "scol", "ibmode",
                             "list", "o", "omode", "snmode", "spin","myself","io_literal_phases",
                             "max_z", "min_z", "max_z_he", "min_z_he", "max_zams", "min_zams",
                             "max_zams_he", "min_zams_he", "nthreads","tf","tini","dtout")

    _ID = 0
    _sevnParams = ut.copy_dict_and_exclude(sw.sevnParams, exclude_keys=_param_to_not_include)
    _sevnParamsDefault = ut.copy_dict_and_exclude(sw.sevnParams, exclude_keys=_param_to_not_include)
    _sevnParamsDescription =ut. copy_dict_and_exclude(sw.sevnio_param()[1], exclude_keys=_param_to_not_include)
    _sevnParams_history = {_ID:_sevnParams}

    def __init__(self,params: Optional[Dict[str,StandardValue]] =None):
        """
        Initialise a SEVNmanager instance

        Parameters
        ----------
        params:
            Dictionary containing the SEVN parameters, the parameters that are not included in the dictionary
            are set to their default value
            (use the static method :func:`~SEVNmanager.get_sevnParamsDefault` to
            check what are the default parameters). If None use the default values for all the parameters

        Warning
        -----
        Never directly create a SEVNmanager instance (see general description of the class for additional info)

        """
        SEVNmanager.init(params)


    @staticmethod
    def get_ID() -> int:
        """
        Get the current SEVNmanager ID. Each time the SEVNmanager is initiliased the ID
        is increased by 1. To get the history of all the parameters set for a given ID
        use the static method :func:`~SEVNmanager.get_sevnParams_history`

        Returns
        -------
        ID: int
            The current SEVNmanager ID

        """

        return copy.deepcopy(SEVNmanager._ID)

    @staticmethod
    def get_sevnParams() -> Dict[str,StandardValue]:
        """
        Return A copy of the Dictionary containing the pair param_name:value where value is the current
        value stored in teh SEVNmanager.

        Returns
        -------
        sevnParams : Dictionary
            Dictionary containing the current SEVNparameters

        Note
        ------
        The Dictionary is a copy of the internal one, changes
        on the dictionary will not affect the internal one, so that the current status will be always preserved
        until a new initialisation is called.

        """

        return copy.deepcopy(SEVNmanager._sevnParams)

    @staticmethod
    def get_sevnParamsDefault() -> Dict[str,StandardValue]:
        """
        Return A copy of the default SEVN parameters

        Returns
        -------
        sevnParams: Dictionary
            Dictionary containing the current SEVNparameters

        Note
        ------
        The Dictionary is a copy of the internal one, changes
        on the dictionary will not affect the internal one.

        """

        return copy.deepcopy(SEVNmanager._sevnParamsDefault)

    @staticmethod
    def get_sevnParams_history() -> Dict[int,Dict[str,StandardValue]]:
        """
        For each intialisation the SEVNmanager intenral ID is increased by 1.
        An internal dictionary stores all the SEVN parameters used for the initialisation corresponded to the ID.
        This method return the dictionary storing all those parameters

        Returns
        -------
        sevnParam_history: Dictionary
            Dictionary containing the parid ID:sevnParams, where sevnParams is the
            set of SEVN parameters for the SEVNmanager initialisation correspondent to ID

        Note
        ------
        The Dictionary is a copy of the internal one, changes
        on the dictionary will not affect the internal one.

        """

        return copy.deepcopy(SEVNmanager._sevnParams_history)

    @staticmethod
    def init(params: Optional[Dict[str,StandardValue]] =None):
        """
        Initialise the SEVNmanager with a given set of parameters.
        This call is necessary to proper set the connection with the SEVN C++ background.
        After the call all the method calling SEVN C++ methods will use the given set of SEVN parameters.
        If a new init is called, the old one is automatically closed using the method
        :func:`~SEVNmanager.close`

        Parameters
        ----------
        params:
            Dictionary containing the SEVN parameters, the parameters that are not included in the dictionary
            are set to their default value
            (use the static method :func:`~SEVNmanager.get_sevnParamsDefault` to
            check what are the default parameters). If None use the default values for all the parameters

        Example
        ------

        >>> from sevnpy.sevn import SEVNmanager
        >>> SEVNmanager.init() #Initialise with the default parameters
        >>> SEVNmanager.init({"ce_alpha":2}) #Initialise with a new set of parameters, close is automatically called

        It is equivalent to

        >>> from sevnpy.sevn import SEVNmanager
        >>> SEVNmanager.init() #Initialise with the default parameters
        >>> SEVNmanager.close()
        >>> SEVNmanager.init({"ce_alpha":2}) #Initialise with a new set of parameters

        If a unknown SEVN parameters is  included in the parameter list, the SEVNmanager will not be initialised and a
        proper error is raised informing the user of the unknown parameters

        """
        if params is not None:
            for param in params:
                # Before to update the parameter dictionary check if there are unknown parameters
                if not param in SEVNmanager._sevnParams:
                    raise ValueError(f"input parameter {param} is not a valid SEVN parameter. Use "
                                     f"SEVNmanager.get_sevnParams().keys() to obtain a list of valise SEVN parameters")
                SEVNmanager._sevnParams[param]=params[param]

        sw.sevnio_initialise(SEVNmanager._sevnParams)
        SEVNmanager._switch_open()
        SEVNmanager._ID+=1 #Each initialisation increment the ID, so that we know how many initiliasation have been done
        SEVNmanager._sevnParams_history[SEVNmanager._ID] = copy.deepcopy(SEVNmanager._sevnParams)

    @staticmethod
    def close():
        """
        Close the current SEVN session

        """

        sw.sevnio_finalise()
        SEVNmanager._sevnParams.update(SEVNmanager._sevnParamsDefault)
        SEVNmanager._switch_close()
        SEVNmanager._clear()

    @staticmethod
    def sevn_path() -> str:
        """
        Return the path to the SEVN folder from which the SEVN C++ extensions have been installed

        Returns
        -------
        SEVNpath: str
            Path of the SEVN folder

        """
        return copy.deepcopy(sw.sevnParams["myself"])

    @staticmethod
    def tables_info() -> Dict[str,StandardValue]:
        """
        Return a dictionary containing the information about the loaded SEVN tables

        Returns
        -------
        table_info: Dictionary
            Dictionary containing the info about the loaded tables, the key are:

            - tables: path to the tables loaded for the H-star
            - tables_HE: path to the tables loaded for the He-star
            - max_z: Highest metallicity in the H-star tables
            - min_z: Lowest metallicity in the H-star tables
            - max_z_he:  Highest metallicity in the He-star tables
            - min_z_he: Lowest metallicity in the He-star tables
            - max_zams: Highest zero age main sequence mass  in the H-star tables
            - min_zams: Lowest zero age main sequence mass  in the H-star tables
            - max_zams_he: Highest zero age main sequence mass  in the He-star tables
            - min_zams_he: Lowest zero age main sequence mass  in the He-star tables

        Warning
        -------
        The method raises a RuntimeError if the tables have not been loaded

        """

        if SEVNmanager._static_close_flag: raise RuntimeError("Cannot retrieve table info because SEVNmanager has not been initiliased")
        _key_tables = ("tables","tables_HE","max_z", "min_z", "max_z_he", "min_z_he", "max_zams", "min_zams","max_zams_he", "min_zams_he")
        _raw_param=SEVNmanager.get_sevnparams_raw()
        return {key:_raw_param[key]for key in _key_tables}

    @staticmethod
    def sevnparam_describe(param_name: Optional[Union[str,ListLikeType]] = None) -> Union[str,Dict[str,str]]:
        """
        Get the short description of a SEVN parameter or multiple SEVN parameters

        Parameters
        ----------
        params:
            Three possibilities:

            - A: a string with the SEVN parameter name
            - B: a collection of string referring to SEVN parameters
            - C: None

        Returns
        -------
        sevnparam_description: str or Dictionary
            Three possibilities depending on the input:

            - A: return the description of the parameters with the given name
            - B: return a Dictionary with pair param_name:param_description containing
            - C: return the Dictionary containing all the parameters


        .. note::
            The SEVN parameters in the SEVNmanager does not include all the IO parameters from SEVN since some of them
            are just used in the main SEVN excutables (e.g. io_logfile, scol, bcol,... ecc). The full list of the SEVN
            parameters can be retrieved with the static method  :func:`~SEVNmanager.get_sevnparams_raw`

        """

        if param_name is None:
            return copy.deepcopy(SEVNmanager._sevnParamsDescription)
        elif isinstance(param_name,list) or isinstance(param_name,tuple) or isinstance(param_name,np.ndarray):
            return {key:copy.deepcopy(SEVNmanager._sevnParamsDescription[key]) for key in param_name}

        return copy.deepcopy(SEVNmanager._sevnParamsDescription[param_name])

    @staticmethod
    def get_sevnparams_raw() -> Dict[str,StandardValue]:
        """
        SEVNpy does not use all the SEVN parameters, some of them are just useful in the context of the C++
        SEVN executables.
        The methods :func:`~SEVNmanager.get_sevnParams` and :func:`~SEVNmanager.get_sevnParamsDefault` just
        returns the SEVN parameters used in SEVNpy.
        This method instead return the full list of SEVN parameters with their default values.

        Returns
        -------
        sevnParams_all: Dictionary
            Dictionary containing the pairs param_name:param_default_value for all the SEVN parameters

        Note
        ------
        The Dictionary is a copy of the internal one, changes
        on the dictionary will not affect the internal one.

        """

        return copy.deepcopy(sw.sevnio_param()[0])

    @staticmethod
    def _switch_open():
        #Set the open and close flag
        SEVNmanager._static_open_flag = True
        SEVNmanager._static_close_flag = False

    @staticmethod
    def _switch_close():
        #Set the open and close flag
        SEVNmanager._static_open_flag = False
        SEVNmanager._static_close_flag = True

    @staticmethod
    def _clear():
        #Reset the internal dictiomary generated from the metaclass
        RestrictiveSingleton._instances.clear()

    @staticmethod
    def check_initiliased() -> bool:
        """
        Check if the SEVNmanager is currently initialised

        Returns
        -------
        intialised_flag: bool
            True if SEVNmanager is currently initialised, False otherwise

        Examples
        --------

        >>> from sevnpy.sevn import SEVNmanager
        >>> SEVNmanager.check_initiliased() #False, SEVNmanager it not initialised
        >>> SEVNmanager.init()
        >>> SEVNmanager.check_initiliased() #True, SEVNmanager is initialised
        >>> SEVNmanager.close()
        >>> SEVNmanager.check_initiliased() #False, SEVNmanager it not initialised
        >>> with SEVNmanager() as sm:
        >>>     sm.check_initiliased() #True, SEVNmanager is initialised within the context manager
        >>>     #SEVNmanager.check_initiliased() #Equivalent to the above call
        >>> SEVNmanager.check_initiliased() #False, SEVNmanager has been closed when exiting from the context manager

        """
        return SEVNmanager._static_open_flag

    def __enter__(self):
        # Just check that is not closed, it should not be, given it is initilised when the instance is created
        if not SEVNmanager.check_initiliased():
            raise RuntimeError(
                f"The instance of the class {self.__class__.__name__} is not opened, "
                f"open it before using a context manager")
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        SEVNmanager.close()
