"""
==============================================================
Star , (:mod:`sevnpy.sevn.star`)
==============================================================

This module contains the class Star. It is used to initialise stars and evolve
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


def _check_masses(Mass: Optional[float] = None, MHE: Optional[float] = None, MCO: Optional[float] = None):
    """
    Auxiliary function for the class Star to check the input modifications to the masses

    Raises
    ------
    ValueError
        If the Mass hierarchy is not satisfied (Mass, MHE, MCO), a mass type cannot be not None if the precedent
        mass type is None. The error is raised also if the mass is larger than the preceded mass type

    """
    # Check hierarchy: MHE and MCO not None, only if Mass is not None. MCO not None, only if MHE not None
    if Mass is None and (MHE is not None or MCO is not None):
        raise ValueError("In Star initialisation the input Mass is None, but MHE or MCO is not None,"
                         " this is not allowed")
    if MHE is None and MCO is not None:
        raise ValueError("In Star initialisation the input MHE is None, but MCO is not None,"
                         " this is not allowed")
    # Check problems
    if Mass is not None and MHE is not None:
        if MHE > Mass:
            raise ValueError(f"In Star initialisation the input Mass ({Mass}) is larger "
                             f"than helium core mass ({MHE}), this is not allowed")
    if MCO is not None:
        if MCO > MHE:
            raise ValueError(f"In Star initialisation the input MHE ({MHE}) is larger "
                             f"than CO core mass ({MCO}), this is not allowed")

    return True


class Star:
    """
    The class Star is used to initialise a star and evolve it using the SEVN backend.

    Initilisation
    --------------
    The basic parameter used to initialise the star are the zero-age-main-seqeunce mass Mzams and Z.
    These two values identify a unique interpolating track (see the SEVN userguide) for the evolution.
    The user need also to define the initial age of the star (set by defaul to "zams").
    In addition, the used can set the initial stellar rotation, the snmodel used to create the stellar remnant.
    The user can also decide to modify the stadandard property of the  interpolating track setting the current
    value of the total mass and mass of the cores.
    In addition to a normal H-star, the parameter star_flag can be used to intialise a pureHe star (in this case
    the Mzams value is the mass of the star at the beginning of the core helium burning), or a remnant such
    as a black hole, a neutron star or a whith dwarf (in this case Mzams represents the mass of the remnant).
    See the :func:`~Star.__init__` method for additional details.

    The Star can also be initialised from the current properties of another star using the
    class method :func:`~Star.from_star`. Similarly a star can create a new Star instance from the
    current properties using the method :func:`~Star.to_star`

    Any instance of the class Star needs to be created within an initialised SEVNmanager,
    otherwise an error is raised, example:

    >>> from sevnpy.sevn import Star, SEVNmanager
    >>> s=Star(10,0.02) # !ERROR, SEVNmanager is not initialised yet
    >>> SEVNmanager.init()
    >>> s=Star(10,0.02) # !OK, SEVNmanager is initiliased

    When a Star instance is created it get the SEVNparameters from the current SEVNmanager initialisation.
    The Star evolution methods can be call  only within the same SEVNmanager session, is the SEVNmanager
    is re-intialised any further evolution will raise an error, e.g.

    >>> from sevnpy.sevn import Star, SEVNmanager
    >>> # OK!
    >>> SEVNmanager.init()
    >>> s=Star(10,0.02)
    >>> s.evolve()
    >>> SEVNmanager.close()
    >>> # Error!
    >>> SEVNmanager.init({"ce_alpha:2"})
    >>> s.evolve() #Error, the SEVNmanager has been reinitialised with respect to the Star initialisation



    Get the stellar properties
    --------------------------
    The stellar properties are stored in an internal evolution table reporting the star properties at each
    time step. After the initilisation the evolution table contains one single row storing the initial
    stellar properties.
    The methods :func:`~Star.getp` and :func:`~Star.getp_array` can be used to retrieve the stellar properties
    as a pandas DataFrame and a numpy array, respectively.
    As a shortchut the :func:`~Star.getp` method can be called directly using the [] operator on the class instance,
    e.g.

    >>> s=Star(10,0.02)
    >>> s["Mass"] #return the Mass column from the current evolved Dataframe
    >>> s["Worldtime","Mass","Radius"] #return the Worldtime,Mass and Radius columns from the  current evolved Dataframe


    Single stellar evolution
    ------------------------
    The main purpose of the class Star is to faciliate the call of the SEVN evolution backend.
    this is done through the methods :func:`~Star.evolve` and :func:`~Star.evolve_for`:

        - :func:`~Star.evolve`: Evolve the Star from the initial age (the one set at the initilisation), up to
            a given time in input. During the evolution, the internal evolution table  is filled with
            properties at each timestep. Each time  the method is called, the evolve_table will be overwritten.
        - :func:`~Star.evolve_for`: Evolve the Star for an amount of time in input starting from the current
            properties of the star, e.g. the last row in the evolution table. This method do not overwrite the
            evolution table rather it appends the new value to it.

    When the star is evolved, the log messages (see SEVN documentation) are also internally store.
    To retrieve it use the property :func:`~Star.log`. Same as for the evolution table, each
    call to :func:`~Star.evolve` ovewrites the log, while calls to :func:`~Star.evolve_for` append
    their log messages.


    Remnant
    --------
    It is possibile to retreive directly the property of the remnant generated by the Star using the method
    :func:`~Star.get_remnant`. Internally a complete evolution will be performed, but whitout updating or overwriting
    the evolution table, therefore the Star status will not be changed. The method returns a DataFrame
    with all the properties of the remnant.



    """

    _static_ID = 0  # static counter to generate IDs

    def __init__(self,
                 Mzams: Number,
                 Z: Number,
                 spin: Number = 0,
                 tini: Union[Number, str] = "zams",
                 snmodel: str = "rapid_gauNS",
                 star_flag: str = "H",
                 rseed: Optional[int] = None,
                 ID: Optional[int] = None,
                 Mass: Optional[Number] = None,
                 MHE: Optional[Number] = None,
                 MCO: Optional[Number] = None,
                 Radius: Optional[Number] = None
                 ):
        """
        Initialise a Star instance. During the initialisation, the evolve table will be filled with one single row
        containing the initial Stellar properties.

        Parameters
        ----------
        Mzams:
            Zero-age-main-sequence mass of the star (or better of the interpolating track) [Msun]
        Z:
            Metallicity of the star (or better of the interpolating track)
        spin:
            Initial stellar spin, i.e. the ratio between the rotational angular velocity and the critical angular
            velocity
        tini:
            Initial age of the star to initialise, the options are:

            - age number, a float in Myr
            - phase initialisation, initialise the star at a given phase, using the string:

                - *zams*: Zero age main sequence
                - *tams*: Terminal main sequence
                - *shb*: shell H burning
                - *cheb*: core helium burning
                - *tcheb*: terminal core helium burning
                - *sheb*: shel He burning

            - percentage phase initialisation, using the string *%<P>:<IDphase>* where *<P>* is the percentage of the
              phase, and *<IDphase>* is the integer ID depending on the phase:

                - 1: *zams*
                - 2: *tams*
                - 3: *shb*
                - 4: *cheb*
                - 5: *tcheb*
                - 6: *sheb*

                so, for example to initialise the Star at 48% of the cheb phase: *%48:4*

        snmodel:
            SEVN snmodel to use to transform a Star to a remnant (see the SEVN userguide)
        star_flag:
            String that  defines the type of Star:

            - *H*: initialise an Hydrogen star
            - *HE*: initialise a pureHe star (use the pureHe trakcs, see the SEVN userguide)
            - *HEWD*: initialise a Helium White Dwarf remnant with mass equal to Mzams
            - *COWD*: initialise a Carbon-Oxygen White Dwarf remnant with mass equal to Mzams
            - *ONEWD*: initialise a Oxygen-Neon White Dwarf remnant with mass equal to Mzams
            - *NS*: initialise a Neutron Star with mass equal to Mzams
            - *BH*: initialise a Black Hole with mass equal to Mzams

        rseed:
            Random seed to be used in the stellar evolution, if None a random seen will be automatically generated
        ID:
            ID of the Star, if None a ID will be automatically assigned
        Mass:
            Use this value to modify the initial Mass of the star [Msun].
             If None the Mass will be the one from the interpolating track.
             If star_flag="HE", to modify the total mass use Mass instead of MHE
        MHE:
            Use this value to modify the initial helium core mass (MHE) of the star [Msun].
            If None MHE will be the one from the interpolating track.
            MHE can be not None only if Mass is not None. It cannot be larger than Mass, if equal to Mass
            a pureHe star will be initialised (forcing the star_model to be "HE").
        MCO:
            Use this value to modify the initial carbon-oxygen core mass (MCO) of the star [Msun].
            If None MCO will be the one from the interpolating track.
            MCO can be not None only if MHE is not None. It cannot be larger than MHE.
        Radius:
            Use this value to modify the initial value of the stellar radius [Rsun].



        .. note::
            If you want to initialise a pureHe star with a custom initial mass
            use the star_flag="HE" and set the initial mass with the parameter Mass (not MHE)
            E.g. to initialise a pureHE star following the interpolating track with Zams mass 30 Msun,
            Z =0.01 but with initial Mass 25 Msun (instead of 30 Msun).

            >>> s=Star(Mzams=30, Z=0.01, tini="cheb", star_flag="HE", Mass=30)

        Raises
        ------
        ValueError
            If the  hierarchy  of the mass in input is not satisfied (Mass, MHE, MCO),
            a mass type cannot be not None if the precedent mass type is None.
            The error is raised also if the mass is larger than the preceded mass type

        """
        # Check_masses
        _check_masses(Mass, MHE, MCO)

        # Check star_flag
        if Mass is not None and MHE is not None:
            if ut.check_equality(Mass, MHE) and star_flag != "HE":
                star_flag = "HE"
                warnings.warn(f"In star initialisation the input Mass ({Mass}) is equal to "
                              f"the helium core mass ({MHE}), the Star will be initialised as pureHE star")
            elif star_flag == "HE" and Mass > MHE:
                star_flag = "H"
                warnings.warn(f"In star initialisation the star_flag is HE, but input Mass ({Mass}) is larger than  "
                              f"the helium core mass ({MHE}), the Star will be initialised as a H star")
        elif Mass is not None and star_flag == "HE":
            MHE = Mass

        self._Mzams = Mzams
        self._Z = Z
        self._spin = spin
        self._tini = tini
        self._massini = Mass
        self._mheini = MHE
        self._mcoini = MCO
        self._radiusini = Radius
        self._snmodel = snmodel
        self._star_flag = star_flag
        self._evolve_counter = 0
        self._properties_interpolators = {"last_counter": 0}

        if rseed is None:
            self._rseed = np.random.randint(1, int(1E15))
        else:
            self._rseed = rseed

        self._SEVNmanager_ID = None
        self._used_sevnParams = None

        self._evolve_dataframe = pd.DataFrame()

        self._tlife = None
        self._name = None
        self._log = None

        # Set ID
        if ID is None:
            self._ID = Star._static_ID
            Star._static_ID += 1
        else:
            self._ID = ID

        # Initialise the star to the starting value
        self._goto(t=tini, Mass=self._massini, MHE=self._mheini, MCO=self._mcoini, Radius=self._radiusini)
        self._tlife_at_initialisation = self._tlife  # Same as tini but never change (unless a reinitialisation is called)

    @classmethod
    def from_star(cls, star: Star, ID: Optional[int] = None, rseed: Optional[Union[int, str]] = "same") -> Star:
        """
        Construct a new Star instance from the current state (i.e. the last points in the current evolution dataframe)
        of a Star in input

        Parameters
        ----------
        star:
            An instance of class :class:`Star`. The new instance will be created using the current properties of star,
            i.e., the properties in the last row of the evolution dataframe.
        ID:
            The ID of the new star, if None automatically generate one from a static ID counter
        rseed:
            The random seed of the new star, if equal to the string *same* get the random seed from star

        Returns
        -------
        new_star : Star
            a new instance of the class :class:`Star`

        """
        # Get star properties
        zams, spin, tini, tnow, mass, mhe, mco, rad, remtype = star.getp_array(
            properties=["Zams", "Spin", "Localtime", "Worldtime", "Mass", "MHE", "MCO", "Radius", "RemnantType"],
            mode="last")

        if rseed is None or isinstance(rseed, int):
            pass
        elif rseed == "same":
            rseed = star.rseed
        else:
            raise ValueError(f"Rseed {rseed} not allowed")

        return Star(Mzams=zams, Z=star.Zmet,
                    spin=spin, tini=tini,
                    snmodel=star.snmodel,
                    star_flag=star.star_flag,
                    rseed=rseed,
                    ID=ID, Mass=mass,
                    MHE=mhe, MCO=mco,
                    Radius=rad
                    )

    def to_star(self, ID: Optional[int] = None, rseed: Optional[Union[int, str]] = "same") -> Star:
        """
        Create a new Star instance from the current state (i.e. the last points in the current evolution dataframe)

        Parameters
        ----------
        star:
            An instance of class :class:`Star`. The new instance will be created using the current properties of star,
            i.e., the properties in the last row of the evolution dataframe.
        ID:
            The ID of the new star, if None automatically generate one from a static ID counter
        rseed:
            The random seed of the new star, if equal to the string *same* get the random seed from star

        Returns
        -------
        new_star : Star
            a new instance of the class :class:`Star`

        """

        return self.from_star(self, ID, rseed)

    def _evolve_basic(self,
                      tstart: Union[str, Number] = "zams",
                      tend: Union[str, Number] = "end",
                      just_init: bool = False,
                      Mass: Optional[Number] = None,
                      MHE: Optional[Number] = None,
                      MCO: Optional[Number] = None,
                      Radius: Optional[Number] = None,
                      Mzams: Optional[Number] = None,
                      spin: Optional[Number] = None) -> Tuple[Dict, Dict]:
        """
        This function represent the interface to call  the function evolve_star from the SEVN C++ wrapper.
        Every evolve-like call, even the one used just to initialise the star need to pass through this
        function since in addition to the evolve call it sets some other class quantities.
        The initialised variable Z cannot be change

        Parameters
        ----------
        tstart:
            Initial age of the star to initialise, the options are:

            - age number, a float in Myr
            - phase initialisation, initialise the star at a given phase, using the string:

                - *zams*: Zero age main sequence
                - *tams*: Terminal main sequence
                - *shb*: shell H burning
                - *cheb*: core helium burning
                - *tcheb*: terminal core helium burning
                - *sheb*: shel He burning

            - percentage phase initialisation, using the string *%<P>:<IDphase>* where *<P>* is the percentage of the
              phase, and *<IDphase>* is the integer ID depending on the phase:

                - 1: *zams*
                - 2: *tams*
                - 3: *shb*
                - 4: *cheb*
                - 5: *tcheb*
                - 6: *sheb*

                so, for example to initialise the Star at 48% of the cheb phase: *%48:4*
        tend:
            Stopping time of the simulation, can be a number in Myr or the word end.
            If end, the evolution will stop when a remnant is generated

        just_init:
            If True the star will be just initialised and the evolution is not called, i.e. the
            results will just store the initialisation values. In this case, the parameter tend will
            be not considered.

        Mass:
            Use this value to modify the initial Mass of the star [Msun].
             If None the Mass will be the one from the interpolating track.
        MHE:
            Use this value to modify the initial helium core mass (MHE) of the star [Msun].
            If None MHE will be the one from the interpolating track.
            MHE can be not None only if Mass is not None. It cannot be larger than Mass, if equal to Mass
            a pureHe star will be initialised (forcing the star_model to be "HE").
        MCO:
            Use this value to modify the initial carbon-oxygen core mass (MCO) of the star [Msun].
            If None MCO will be the one from the interpolating track.
            MCO can be not None only if MHE is not None. It cannot be larger than MHE.
        Radius:
            Use this value to modify the initial value of the stellar radius [Rsun].
        Mzams:
            Zero-age-main-sequence mass  of the interpolating track [Msun]. I Mzams is None use the one
            defined at the star initialisation
        spin:
            Initial stellar spin, i.e. the ratio between the rotational angular velocity and the critical angular
            velocity. If None use the one defined at the star initialisation

        Returns
        -------
        evolution_results: Dictionary
            Dictionary containing the evolution results. Each key is a star property and each item is a
            numpy array storing the given properties during the stellar evolution

        additiona_evolution_info: Dictionary
            Dictionary storing additional evolution info such as the name assigned to the star, the logfile,
            etc..

        """

        # Check_masses
        _check_masses(Mass, MHE, MCO)

        star_flag = self._star_flag
        if Mass is not None and MHE is not None:
            if ut.check_equality(Mass, MHE) and star_flag != "HE":
                star_flag = "HE"
                warnings.warn(f"In Star::_evolve_basic  the input Mass ({Mass}) is equal to "
                              f"the helium core mass ({MHE}), the Star will be initialised as pureHE star")

        if Mzams is None: Mzams = self._Mzams
        if spin is None: spin = self._spin

        self._check_SEVNmanager_synchronisation()
        results, extra_info = sw.evolve_star(Mzams=Mzams,
                                             Z=self._Z,
                                             spin=spin,
                                             tstart=tstart,
                                             tend=tend,
                                             star_flag=star_flag,
                                             snmodel=self._snmodel,
                                             rseed=self._rseed,
                                             Mass=Mass,
                                             MHE=MHE,
                                             MCO=MCO,
                                             Radius=Radius,
                                             just_init=just_init,
                                             )
        self._after_evolve_duties()  # Increment counters and update other private variables

        return results, extra_info

    def _goto(self, t: Union[str, Number],
              Mass: Optional[Number] = None,
              MHE: Optional[Number] = None,
              MCO: Optional[Number] = None,
              Radius: Optional[Number] = None):
        """
        Auxiliary function to initialise the property of the Star

        """

        results, extra_info = self._evolve_basic(tstart=t,
                                                 Mass=Mass,
                                                 MHE=MHE,
                                                 MCO=MCO,
                                                 Radius=Radius,
                                                 just_init=True)
        # Update and overwrite quantities
        self._evolve_dataframe = pd.DataFrame(results)
        self._name = extra_info["name"]
        self._log = extra_info["Log"]
        self._tlife = extra_info["tlife"]

    def evolve(self, tend: Union[str, Number] = "end"):
        """
        Main evolve function. The star is evolved from their initial age (set at the initialisation)
        to the value in input.

        .. warning::
            Each time the evolve function is called, the past evolution of the Star is overwritten.
            If you want to continue the evolution from the current Star state use the method
            :func:`~Star.evolve_for`

        Parameters
        ----------
        tend:
            Stopping time of the simulation, can be a number in Myr or the word *end*.
            If *end*, the evolution will stop when a remnant is generated

        """
        results, extra_info = self._evolve_basic(tstart=self._tini,
                                                 tend=tend,
                                                 just_init=False,
                                                 Mass=self._massini,
                                                 MHE=self._mheini,
                                                 MCO=self._mcoini,
                                                 Radius=self._radiusini,
                                                 )
        # Update and overwrite quantities
        self._evolve_dataframe = pd.DataFrame(results)
        self._name = extra_info["name"]
        self._log = extra_info["Log"]
        self._tlife = extra_info["tlife"]

    def evolve_for(self, dt: Number):
        """
        Evolve the star from the current status for a given time interval.

        Parameters
        ----------
        dt:
            Evolution time interval [Myr]

        Examples
        --------
        It is possible to use evolve_for to perform a step-by-step stellar evolution directly in Python.
        For example assume that we want to evolve a star  for 10 Myr and get the properties each 1 Myr

        >>> from sevnpy.sevn import SEVNmanager, Star
        >>> import pandas as pd
        >>> SEVNmanager.init()
        >>> s1=Star(10,0.02)
        >>> res_df = s1.getp(mode="last") #Store initial properties
        >>> t=0
        >>> dt=1
        >>> while t<10:
        >>>     s1.evolve_for(dt)
        >>>     res_df = pd.concat(res_df, s1.getp(mode="last"))
        >>>     t+=dt
        >>> SEVNmanager.close()

        """
        # Get current properties
        zams, spin, tini, tnow, mass, mhe, mco, rad, phase = self.getp_array(
            properties=["Zams", "Spin", "Localtime", "Worldtime", "Mass", "MHE", "MCO", "Radius", "Phase"],
            mode="last")

        # Set star flag to deal with remnants
        old_star_flag = self._star_flag
        self._star_flag = self.star_flag

        # Deal with remnant
        if int(phase) == 7:
            zams = mass
            mass = None
            mco = None
            mhe = None
            rad = None

        results, extra_info = self._evolve_basic(tstart=tini,
                                                 tend=dt,
                                                 just_init=False,
                                                 Mass=mass,
                                                 MHE=mhe,
                                                 MCO=mco,
                                                 Radius=rad,
                                                 Mzams=zams,
                                                 spin=spin
                                                 )
        # Reset starflag
        self._star_flag = old_star_flag

        # Update quantities (not overwrite, except for tlife)
        current_results = pd.DataFrame(results).iloc[1:]
        current_results["Worldtime"] = current_results["Worldtime"] + tnow
        if int(phase) == 7:
            current_results["Zams"] = zams
            current_results["Localtime"] = current_results["Worldtime"]
        self._evolve_dataframe = pd.concat([self._evolve_dataframe, current_results])
        self._log += extra_info["Log"]
        self._tlife = extra_info["tlife"]

    def get_remnant(self) -> pd.DataFrame:
        """
        Get the property of the remnant. This method retrieve the properties of the remnant produced
        by evolving the star from the properties set at the star initialization, therefore it neglect the
        current status of the star.

        Returns
        -------
        remnant: pandas DataFrame
            A Pandas DataFrame containing the remnant properties


        .. note::
            Even if the method will internally evolve the star, the results of this evolution
            will not be stored in the class attribute, but only returned by the method

        """

        results, extra_info = self._evolve_basic(tstart=self._tini,
                                                 tend="end",
                                                 just_init=False,
                                                 Mass=self._massini,
                                                 MHE=self._mheini,
                                                 MCO=self._mcoini,
                                                 Radius=self._radiusini,
                                                 )

        return pd.DataFrame(results).tail(1)

    def look_at_track(self, t: Union[str, Number]) -> pd.DataFrame:
        """
        The method returns the properties of the *interpolating tracks* followed by the star
        at a given time t. Since this will just consider the interpolating track, the initialisation properties
        Mass, MHE and MCO will not be considered.

        Parameters
        ----------
        t:
            Time in Myr at which look for the properties of the interpolating track

        Returns
        -------
        star_properties:  pandas DataFrame
            A Pandas DataFrame containing the  properties of the interpolating tracks at time t.
            If the time t is larger than the stellar lifetime, the properties of the generated remnant
            is returned (see :func:`~Star.get_remnant`)


        .. warning::
            The time inserted in this method refers to the Localtime, while  the time used in
            :func:`~Star.evolver` and :func:`~Star.evolve_for` refers to the Worldtime

        """
        if t >= self._tlife_at_initialisation:
            warnings.warn("The input t in look_at_tracks is larger than the life time, a remnant is returned")
            # Save the initial mass values
            _Mass, _MHE, _MCO = self._massini, self._mheini, self._mcoini
            # Set the None to call remnant
            self._massini = None
            self._mheini = None
            self._mcoini = None
            ret = self.get_remnant()
            # Restore values
            self._massini = _Mass
            self._mheini = _MHE
            self._mcoini = _MCO

            ret["Worldtime"] = t
            return ret
        else:
            results, _ = sw.evolve_star(Mzams=self._Mzams,
                                        Z=self._Z,
                                        spin=self._spin,
                                        tstart=t,
                                        tend="end",
                                        star_flag=self._star_flag,
                                        snmodel=self._snmodel,
                                        rseed=self._rseed,
                                        Mass=None,
                                        MHE=None,
                                        MCO=None,
                                        Radius=None,
                                        just_init=True,
                                        )
            ret = pd.DataFrame(results)
            ret["Worldtime"] = t
            return ret

    def reinit(self,
               spin: Optional[float] = None,
               tini: Optional[Union[str, Number]] = None,
               snmodel: Optional[str] = None,
               Mass: Optional[Number] = None,
               MHE: Optional[Number] = None,
               MCO: Optional[Number] = None,
               Radius: Optional[Number] = None):
        """
        Use this method to re-initialise the star changing the initial properties, except for Mzams and Z.
        If you need to change Mzams and Z create a new star.

        Parameters
        -----------
        spin:
            Initial stellar spin, i.e. the ratio between the rotational angular velocity and the critical angular
            velocity
        tini:
            Initial age of the star to initialise, the options are:

            - age number, a float in Myr
            - phase initialisation, initialise the star at a given phase, using the string:

                - *zams*: Zero age main sequence
                - *tams*: Terminal main sequence
                - *shb*: shell H burning
                - *cheb*: core helium burning
                - *tcheb*: terminal core helium burning
                - *sheb*: shel He burning

            - percentage phase initialisation, using the string *%<P>:<IDphase>* where *<P>* is the percentage of the
              phase, and *<IDphase>* is the integer ID depending on the phase:

                - 1: *zams*
                - 2: *tams*
                - 3: *shb*
                - 4: *cheb*
                - 5: *tcheb*
                - 6: *sheb*

                so, for example to initialise the Star at 48% of the cheb phase: *%48:4*

        snmodel:
            SEVN snmodel to use to transform a Star to a remnant (see the SEVN userguide)

        star_flag:
            String that  defines the type of Star:

            - *H*: initialise an Hydrogen star
            - *HE*: initialise a pureHe star (use the pureHe trakcs, see the SEVN userguide)
            - *HEWD*: initialise a Helium White Dwarf remnant with mass equal to Mzams
            - *COWD*: initialise a Carbon-Oxygen White Dwarf remnant with mass equal to Mzams
            - *ONEWD*: initialise a Oxygen-Neon White Dwarf remnant with mass equal to Mzams
            - *NS*: initialise a Neutron Star with mass equal to Mzams
            - *BH*: initialise a Black Hole with mass equal to Mzams

        rseed:
            Random seed to be used in the stellar evolution, if None a random seen will be automatically generated

        ID:
            ID of the Star, if None a ID will be automatically assigned

        Mass:
            Use this value to modify the initial Mass of the star [Msun].
             If None the Mass will be the one from the interpolating track.
             If star_flag="HE", to modify the total mass use Mass instead of MHE

        MHE:
            Use this value to modify the initial helium core mass (MHE) of the star [Msun].
            If None MHE will be the one from the interpolating track.
            MHE can be not None only if Mass is not None. It cannot be larger than Mass, if equal to Mass
            a pureHe star will be initialised (forcing the star_model to be "HE").

        MCO:
            Use this value to modify the initial carbon-oxygen core mass (MCO) of the star [Msun].
            If None MCO will be the one from the interpolating track.
            MCO can be not None only if MHE is not None. It cannot be larger than MHE.

        Radius:
            Use this value to modify the initial value of the stellar radius [Rsun].

        """

        self.__init__(Mzams=self._Mzams,
                      Z=self._Z,
                      spin=self._spin if spin is None else spin,
                      tini=self._tini if tini is None else tini,
                      snmodel=self._snmodel if snmodel is None else snmodel,
                      star_flag=self._star_flag,
                      rseed=self._rseed,
                      ID=self._ID,
                      Mass=self._massini if Mass is None else Mass,
                      MHE=self._mheini if MHE is None else MHE,
                      MCO=self._mcoini if MCO is None else MCO,
                      Radius=self._radiusini if Radius is None else Radius
                      )

    def getp(self, properties: Optional[Union[str, ListLikeType]] = None,
             mode: str = "all", t: Optional[Number, ListLikeType] = None) -> pd.DataFrame:
        """
        Return  the stellar properties as a Pandas DataFrame.

        Parameters
        ----------
        properties:
            single property name or list of property names (see SEVN documentation for the property names).
            If None return all the available properties.

        mode:
            Set the return type:

                - *all*: return all the timesteps for the stellar evolution
                - *last*: return just the current Stellar properties (last row form the evolution tables)
                - *first*: return the initial stellar properties, i.e. the properties set at the initialisation

        t:
            If not None, overwrite the mode property and return the properties at time specified by the input value(s)
            In order to return the values at time  not stored in the evolution table, the proeprties are interpolated
            linearly

        Returns
        -------
        star_properties: pandas DataFrame
            A pandas dataframe containing the values of the properties chosen in input

        Raises
        ------
        RuntimeError
            If the  evolution table is empty


        Examples
        ---------

        >>> s = Star(10,0.02)
        >>> s.evolve()
        >>> dfevolve  = s.getp(mode="all") #Get the complete evolution table
        >>> dfevolve  = s.getp(mode="last") #Get all the current properties
        >>> t = np.linspace(0.1,10)
        >>> dfevolve  = s.getp(properties=["Worldtime","Mass","Radius"], t=t) #Get the evolution of time mass and radius
        >>> #at given t


        .. warning::
            Use the t option wisely. The interpolation  consists on a simple 1D interpolation, so it is not comparable
            with the detailed interpolation implemente in SEVN. THerefore, although the SEVN adaptive timestep catches
            all the important changes during the evolution do not bindly rely on the results obtain from this method
            when t is used (if t=None, the results will be robust since they are the one derived by SEVN directly),

        """

        if self._evolve_dataframe.empty:
            raise RuntimeError("Cannot retrieve stellar property with get without "
                               "at least one evolution call")

        if properties is None:
            properties = self._evolve_dataframe.columns
        properties = np.atleast_1d(properties)

        if t is not None:
            t = np.atleast_1d(t)
            self._set_property_interpolator()
            df = pd.DataFrame({prop: self._properties_interpolators[prop](t)
                               for prop in properties})
            return df

        elif mode == "all":
            return self._evolve_dataframe[properties]
        elif mode == "last":
            return self._evolve_dataframe[properties].tail(1)
        elif mode == "first":
            return self._evolve_dataframe[properties].head(1)

    def getp_array(self, properties: Optional[Union[str, ListLikeType]] = None,
                   mode: str = "all", t: Optional[Number, ListLikeType] = None) -> np.ndarray[np.float]:
        """
        Return  the stellar properties as a numpy array.

        Parameters
        ----------
        properties:
            single property name or list of property names (see SEVN documentation for the property names).
            If None return all the available properties.

        mode:
            Set the return type:

                - *all*: return all the timesteps for the stellar evolution
                - *last*: return just the current Stellar properties (last row form the evolution tables)
                - *first*: return the initial stellar properties, i.e. the properties set at the initialisation

        t:
            If not None, overwrite the mode property and return the properties at time specified by the input value(s)
            In order to return the values at time  not stored in the evolution table, the proeprties are interpolated
            linearly

        Returns
        -------
        star_properties: numpy array
            A numpy array  containing the values of the properties chosen in input. The shape depends on the input data
            and mode option:

                - If mode="all" or t is not None, 2D array with  the  is equal to NxC, where N is the number of row
                  in the evolved dataframe and C are the number of  properties in input.
                - if mode="last" or mode="first", 1D array with length=C, where C are the  number of properties in input.




        Raises
        ------
        RuntimeError
            If the  evolution table is empty


        Examples
        ---------

        >>> s = Star(10,0.02)
        >>> s.evolve()
        >>> dfevolve  = s.getp_array(mode="all") #Get the complete evolution table
        >>> dfevolve  = s.getp_array(mode="last") #Get all the current properties
        >>> t = np.linspace(0.1,10)
        >>> dfevolve  = s.getp_array(properties=["Worldtime","Mass","Radius"], t=t) #Get the evolution of time mass and radius
        >>> #at given t


        .. warning::
            Use the t option wisely. The interpolation  consists on a simple 1D interpolation, so it is not comparable
            with the detailed interpolation implemente in SEVN. THerefore, although the SEVN adaptive timestep catches
            all the important changes during the evolution do not bindly rely on the results obtain from this method
            when t is used (if t=None, the results will be robust since they are the one derived by SEVN directly),

        """

        retarray = self.getp(properties=properties, mode=mode, t=t).values

        if t is not None:
            return retarray
        elif mode=="all":
            return retarray

        #if  retarray.shape[1] == 1: retarray = retarray[:,0]
        if len(retarray) == 1 \
                or (len(retarray.shape) == 2 and retarray.shape[1] == 1): retarray = retarray[0]
        #if len(retarray) == 1:  retarray = retarray[0]

        return retarray

    @property
    def pnames(self) -> List:
        """
        List of names of available properties

        """
        return self.getp().keys()

    @property
    def Zmet(self) -> float:
        """
        Stellar metallicity in input

        """
        return self._Z

    @property
    def snmodel(self) -> str:
        """
        SEVN model used to pass from a Star to a remnant (see the SEVN documentation)

        """
        return self._snmodel

    @property
    def tlife(self) -> float:
        """
        Star lifetime based on the Mzams and Z in input [Myr]

        """
        return self._tlife

    @property
    def rseed(self) -> int:
        """
        Random seed used in the evolution

        """
        return self._rseed

    @property
    def name(self) -> float:
        """
        Unique identifier of this object.

        """
        if self._name is None: warnings.warn("Star name is set just after an evolution call")
        return self._name

    @property
    def ID(self) -> int:
        """
        ID identifier of this object.

        """
        return self._ID

    @property
    def Mzams(self) -> float:
        """
        Zero-age-main-sequence mass used to initialise the star [Msun]

        """
        return self._Mzams

    @property
    def star_flag(self) -> str:
        """
        Stellar star flag based on the current properties (e.g. the last entry of the evolved table).

        """

        # Get current status
        mass, mhe, remtype = self.getp_array(properties=["Mass", "MHE", "RemnantType"], mode="last")

        # Check if it is a remnant or a pureHE and set star_flag
        star_flag = "H"
        if remtype == 1:
            star_flag = "HEWD"
        elif remtype == 2:
            star_flag = "COWD"
        elif remtype == 3:
            star_flag = "ONEWD"
        elif remtype == 4:
            star_flag = "NSEC"
        elif remtype == 5:
            star_flag = "NS"
        elif remtype == 6:
            star_flag = "BH"
        elif ut.check_equality(mass, mhe):
            star_flag = "HE"

        return star_flag  # but never change the initial self._star_flag

    @property
    def log(self) -> str:
        """
        Log created during the stellar evolution

        """
        return copy.deepcopy(self._log)

    @property
    def used_sevnParams(self) -> Dict[str:Union[float, str, bool]]:
        """
        Dictionary containing the SEVN parameters used during the stellar evolution

        """
        return copy.deepcopy(self._used_sevnParams)

    @property
    def SEVNmanager_ID(self) -> int:
        """
        Return the ID of the SEVNmanager associated with this Star at the Star initialisation

        """

        return self._SEVNmanager_ID

    @property
    def tables_info(self) -> Dict[str:Union[str, float]]:
        """
        Get a dictionary with the info  of the loaded stellar tables (paths and minimum/maximum Zams - Z values)

        """
        if SEVNmanager.check_initiliased():
            return SEVNmanager.tables_info()
        else:
            warnings.warn("The tables info are filled just after the SEVNmanager initialisation")

        return {}

    @property
    def evolve_table(self) -> pd.DataFrame:
        """
        Return the pandas dataframe containing the table with the stellar properties and their evolution

        """
        if self._evolve_dataframe.empty: warnings.warn("Star evolved dataframe is empty before the evolution calls")
        return copy.deepcopy(self._evolve_dataframe)

    def get_SN_kick(self) -> Dict[str:Any]:
        """
        Get the info about the SN kick velocities.

        Returns
        -------
        SNkick_properties: Dictionary
            A Dictionary containing the following keys:

                - SNtime: the time of the SN explosion (Worldtime)
                - Vkick: a numpy array containing the three Cartesian component of the kick (the axes orientation
                  are arbitrary)

            .. warning::
                If a SN never exploded in the current evolution table, SNtime will be equal to nan and the three components
                of the Vkick will be nan too

        """
        df_log = readlogstring(self.log, events=["SN"])
        SNtime = np.nan
        Vkick = np.array([np.nan, np.nan, np.nan])
        if not df_log.empty:
            SNtime = df_log["time"][0]
            Vkick = np.array([df_log.Vfkx_SN[0], df_log.Vfky_SN[0], df_log.Vfkz_SN[0]])

        return {"SNtime": SNtime, "Vkick": Vkick}

    def _set_property_interpolator(self):
        """
        Set the linear interpolator for all the properties.

        Raises
        -------
        RuntimeError:
            If the current evolve table contains just one row

        """
        if self._properties_interpolators["last_counter"] == self._evolve_counter:
            return
        else:
            self._properties_interpolators["last_counter"] = self._evolve_counter

        if len(self._evolve_dataframe) < 2:
            raise RuntimeError("Star is trying to set the property interpolators, but the evolve table has only "
                               "one value")

        for prop in self._evolve_dataframe.columns:
            prop_int_list = ["Phase", "RemnantType", "PhaseBSE", "Event"]

            if prop == "Worldtime":
                self._properties_interpolators[prop] = lambda x: x
            else:
                if prop in prop_int_list:
                    kind = "previous"
                else:
                    kind = "linear"
                self._properties_interpolators[prop] = interp1d(self._evolve_dataframe["Worldtime"],
                                                                self._evolve_dataframe[prop],
                                                                bounds_error=False, kind=kind)

    def _check_SEVNmanager_synchronisation(self):
        """
        Check if the SEVNmanager is initialised and if the star is running within the context of the same
        SEVNmanager (this assure that there are not unexpected change of SEVN parameters)

        Raises
        -------
        RuntimeError
            If the SEVNmanager is not initialised or if the SEVmanagerID is not the one saved at the star initilisation

        """

        if not SEVNmanager.check_initiliased():
            raise RuntimeError("Cannot evolve  Star, the SEVNmanager has not been initialised")

        elif self._SEVNmanager_ID is not None and self._SEVNmanager_ID != SEVNmanager.get_ID():
            raise RuntimeError("Cannot  evolve Star, the SEVNmanager has been reinitialised")

    def _after_evolve_duties(self):
        """
        Some assignment to perform after an evolution call, i.e. increases the evolve_counter and save
        the SEVNmanager ID and the SEVNmanager parameter related to the ID.

        """
        self._evolve_counter += 1
        self._SEVNmanager_ID = SEVNmanager.get_ID()
        self._used_sevnParams = SEVNmanager.get_sevnParams()

    def __getitem__(self, item):
        """
        Same as :func:`~Star.getp` but with mode set to "all" and t is None
        """
        return self.getp(properties=item,mode="all",t=None)

    def __repr__(self):
        ret = "*** Init value ***"
        ret += f"\nMzams={self._Mzams} Msun Z={self._Z}"
        ret += f"\nSpin={self._spin} Age={self._tini}"
        ret += f"\nSN model={self._snmodel}"
        ret += "\n*** Current properties ***"
        properties = self.getp(mode="last")
        for p in properties.columns:
            ret += f"\n{p}={properties[p].values[0]:.3g}"
        ret += "\n*** Tables  ***"
        Htable = self.tables_info["tables"]
        HEtable = self.tables_info["tables_HE"]
        ret += f"\nH-tables: {Htable}"
        ret += f"\nHe-tables: {HEtable}"

        return ret
