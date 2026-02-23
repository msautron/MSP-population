**********
Quickstart
**********

.. _sevnpy-sevnsession:

The SEVN session
****************

All the functions and classes inside the :mod:`~sevnpy.sevn` submodule have to be run within a SEVN session.
The creation and the properties  of a SEVN session are handled by the static class :class:`~sevnpy.sevn.sevnmanager.SEVNmanager`

The SEVN session is initialised with the static method :func:`~sevnpy.sevn.sevnmanager.SEVNmanager.init` and closed
with the method :func:`~sevnpy.sevn.sevnmanager.SEVNmanager.close`.
All the classes and functions from the :mod:`~sevnpy.sevn` submodule needs to be used within an active SEVN session

.. code-block:: python

    from sevnpy.sevn import SEVNmanager
    SEVNmanager.init() #init the SEVN session and load the tables
    #Use functions and classes from the sevn modules
    .....
    ....
    SEVNmanager.close() #close the Session

The :func:`~sevnpy.sevn.sevnmanager.SEVNmanager.init`  call is also used to set the SEVN runtime parameters.
They can be put in a dictionary (parameter_name:value) and added as input paramter of the :func:`~sevnpy.sevn.sevnmanager.SEVNmanager.init`
static method, e.g.


.. code-block:: python

    from sevnpy.sevn import SEVNmanager
    SEVNmanager.init({"ce_alpha":2, "w_alpha":1, "w_beta":0.1}) #init the SEVN session and load the tables
    #Use functions and classes from the sevn modules
    .....
    ....
    SEVNmanager.close() #close the Session


In this case the parameters *ce_alpha*, *w_alpha* and *w_beta* are set with the values in the dictionary,
all the others maintain their default value.

For additional info about the :class:`~sevnpy.sevn.sevnmanager.SEVNmanager` read the API documentation.

.. _sevnpy-sevnsse:

Single Stellar evolution (SSE)
******************************

Create a star
=============

In order to use the SEVN backend to evolve stars, uses the class :class:`~sevnpy.sevn.star.Star`.
First of all a star object need to created

.. code-block:: python

    from sevnpy.sevn import SEVNmanager, Star
    SEVNmanager.init() #init the SEVN session and load the tables
    # Create a Star with initial zero age main sequence mass of 10 Msun, metallicity Z=0.02, initialise it at zams age
    # and use the rapid model for the remnant creation
    s = Star(Mzams=10, Z=0.02, tini="zams", snmodel="rapid")
    SEVNmanager.close() #close the Session

Evolve a star
=============

There are two methods that can be used to evolve the star:

    - :func:`~sevnpy.sevn.star.Star.evolve`:
        evolve the star up to a time tend in Myr or do not use any input parameter
        to evolve the star until the remnant formation.  Each time this method is called the past evolution of the star
        is overwritten.

        .. code-block:: python

            from sevnpy.sevn import SEVNmanager, Star
            SEVNmanager.init() #init the SEVN session and load the tables
            # Create a Star with initial zero age main sequence mass of 10 Msun, metallicity Z=0.02, initialise it at zams age
            # and use the rapid model for the remnant creation
            s = Star(Mzams=10, Z=0.02, tini="zams", snmodel="rapid")
            # Evolve the Star until the remnant formation
            s.evolve()
            # Evolve the star  until 2 Myr
            s.evolve(2)
            # NOTICE: after this call the first evolution until the remnant formation is overwritten
            SEVNmanager.close() #close the Session

    - :func:`~sevnpy.sevn.star.Star.evolve_for`:
        evolve for a given time interval in Myr. Each time this method is called
        the evolution status of the star is updated.

        .. code-block:: python

            from sevnpy.sevn import SEVNmanager, Star
            SEVNmanager.init() #init the SEVN session and load the tables
            # Create a Star with initial zero age main sequence mass of 10 Msun, metallicity Z=0.02, initialise it at zams age
            # and use the rapid model for the remnant creation
            s = Star(Mzams=10, Z=0.02, tini="zams", snmodel="rapid")
            # Evolve the Star for 1 Myr
            s.evolve_for(1)
            # Evolve the star  for 2 Myr
            s.evolve_for(2)
            # NOTICE: after this call the evolution status is updated, so now the current status of the Star is
            # the one after 3 Myr (1 + 2), so it is equivalent of a s.evolve(3) call
            SEVNmanager.close() #close the Session

Retrieve the star properties
============================

To retrieve the stellar properties use :func:`~sevnpy.sevn.star.Star.getp` and :func:`~sevnpy.sevn.star.Star.getp_array`.
The first one returns a pandas Dataframe containing also the properties names (e.g. the columns names), the second one
returns a numpy array.

.. code-block:: python

    from sevnpy.sevn import SEVNmanager, Star
    SEVNmanager.init() #init the SEVN session and load the tables
    # Create a Star with initial zero age main sequence mass of 10 Msun, metallicity Z=0.02, initialise it at zams age
    # and use the rapid model for the remnant creation
    s = Star(Mzams=10, Z=0.02, tini="zams", snmodel="rapid")
    #Get the initial properties of the star
    df = s.getp() # It returns a Pandas Dataframe containg all the stored properties
    print(df.columns)
    # Index(['Worldtime', 'Localtime', 'Mass', 'Radius', 'Inertia', 'Luminosity',
    #         'Temperature', 'MHE', 'MCO', 'RHE', 'RCO', 'Phase', 'RemnantType',
    #         'PhaseBSE', 'Spin', 'Ebind', 'Zams', 'Event', 'dMRLOdt', 'dMaccwinddt',
    #         'Plife'],
    #        dtype='object')
    print(df)
    #    Worldtime  Localtime      Mass    Radius    Inertia   Luminosity   Temperature  MHE  MCO  RHE  ...  Phase  RemnantType  PhaseBSE  Spin  Ebind  Zams  Event  dMRLOdt  dMaccwinddt  Plife
    #     0        0.0   0.548376  9.999887  3.850697  14.827698  5708.495393  25567.445359  0.0  0.0  0.0  ...    1.0          0.0       1.0   0.0    0.0  10.0   -1.0      0.0          0.0    0.0
    #  [1 rows x 21 columns]
    s.evolve_for(2) # Evolve the star for 2 Myr
    df = s.getp(["Worldtime","Mass","Radius","Phase"]) # get only some properties
    print(df)
    #       Worldtime      Mass    Radius  Phase
    #  0  0.000000e+00  9.999887  3.850697    1.0
    #  1  8.415000e-08  9.999887  3.850697    0.0
    #  2  8.490000e-08  9.999887  3.850697    0.0
    #  3  8.510000e-08  9.999887  3.850697    1.0
    #  4  1.925391e+00  9.999467  4.048132    1.0
    #  5  2.000000e+00  9.999450  4.055983    1.0
    # It is also possible to directly use the operator [, (equivalent to getp
    # but forcing mode=all and t=None)
    df = s["Worldtime","Mass","Radius","Phase"]
    print(df)
    #       Worldtime      Mass    Radius  Phase
    #  0  0.000000e+00  9.999887  3.850697    1.0
    #  1  8.415000e-08  9.999887  3.850697    0.0
    #  2  8.490000e-08  9.999887  3.850697    0.0
    #  3  8.510000e-08  9.999887  3.850697    1.0
    #  4  1.925391e+00  9.999467  4.048132    1.0
    #  5  2.000000e+00  9.999450  4.055983    1.0
    dfarr = s.getp_array(["Worldtime","Mass","Radius","Phase"]) # the same but return directly a numpy array
    print(dfarr)
    # array([[0.00000000e+00, 9.99988717e+00, 3.85069672e+00, 1.00000000e+00],
    #       [8.41500000e-08, 9.99988717e+00, 3.85069673e+00, 0.00000000e+00],
    #       [8.49000000e-08, 9.99988717e+00, 3.85069673e+00, 0.00000000e+00],
    #       [8.51000000e-08, 9.99988717e+00, 3.85069673e+00, 1.00000000e+00],
    #       [1.92539112e+00, 9.99946677e+00, 4.04813204e+00, 1.00000000e+00],
    #       [2.00000000e+00, 9.99945048e+00, 4.05598312e+00, 1.00000000e+00]])
    SEVNmanager.close() #close the Session

More
====

See the :class:`~sevnpy.sevn.star.Star` for additional info about the usage the Star class
and its methods.


Binary Stellar evolution (BSE)
******************************

At the moment, SEVNpy includes the binary stellar evolution wrapper *evolve_binary* in the
module *sevn.sevnwrap*. The complete python wrapper (similar to :class:`~sevnpy.sevn.star.Star`)
will be released soon.
The evolve_binary method gets in input the following parameters:

    - Semimajor: initial semimajor axis in Rsun (Required)
    - Eccentricity: initial eccentricity (Required)
    - Mzams_0: Zams mass of the first star in Msun (Required)
    - Z_0: absolute metallicity of the first  star (Required)
    - Mzams_1: Zams mass of the second star in Msun (Required)
    - Z_1: absolute metallicity of the second star (Required)
    - spin_0: initial spin (angular velocity over critical angular velocity) of the first star (Optional, default=0)
    - tstart_0: initial age of the first star in Myr or the label of a specific phase (Optional, default="zams")
    - spin_1: initial spin (angular velocity over critical angular velocity) of the second star (Optional, default=0)
    - tstart_1: initial age of the second star in Myr or the label of a specific phase (Optional, default="zams")
    - tend: ending time of the evolution, can be a number in Myr, or the keyword "end" to stop when both stars are
      remnant, or the keyword "broken" to stop when the binary is broken. (Optional, default="end")
    - snmodel: SN model to use to form a remnant (Optional, default="rapid")
    - star_flag_0: use this flag to create special stars, "HE" to create a pureHE, "BH" for a BH, "NS" for a NS
      "HEWD", "COWD", "ONEWD" for a WD. The default value "H" refers to a classic hydrogen star.
    - star_flag_1: same as above but for the second star. (Optional, default="H")
    - rseed: Set the random seed, if 0 a "random" random seed will be assigned. (Optional, default=0)
    - Mass_0,MHE_0,MCO_0,Radius_0: use this values to force the value of the mass, core masses (in Msun)
      and stellar radius (Rsun) of the first star to be of a given value. (Optionals, default=None)
    - Mass_1,MHE_1,MCO_1,Radius_1: same as above but for the second star  (Optionals, default=None).
    - just_init: If True the binary is just initialised and not evolved (Optional, default=False)

    The function returns a  dictionary containing the evolution of all the binary and star parameters and a dictionary
    containing extra information from the evolution (e.g. the lof information).

    >>> from sevnpy.sevn import SEVNmanager, sevnwrap
    >>> SEVNmanager.init()
    >>> evolvedict,extrainfo=sevnwrap.evolve_binary(1000,0.2,10,0.02,8,0.02,tend=10)
    >>> SEVNmanager.close()


Examples
******************************
The folder examples in SEVNpy (https://gitlab.com/sevncodes/sevn/-/tree/SEVN/SEVNpy/examples)
contains Python scripts that use SEVNpy. Give them a look!
