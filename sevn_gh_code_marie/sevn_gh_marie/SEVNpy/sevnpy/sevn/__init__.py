"""
==============================================================
SEVN , (:mod:`sevnpy.sevn`)
==============================================================

This module contains all the tools necessary to access to the SEVN backend.
In particular, it contains the class SEVNmanager to initialise and handle the connection
with the SEVN backend, and the class Star and Binary to access the SEVN SSE and BSE backend.
ALl the functions and classes contained in this module needs to be used within a SEVNmanager session.
To initialise a SEVNmanager session use :func:`~SEVNmanager.init`, to close it :func:`~SEVNmanager.close`,
e.g.

>>> from sevnpy.sevn import SEVNmanage
>>> SEVNmanager.init() # initialisation
>>> ...... # Use tools from the sevn module
>>> SEVNmanager.close() # intialisation

During the initilisation it is possibile to set the SEVN parameters (otherwise the default ones will be used),
e.g.

>>> from sevnpy.sevn import SEVNmanage
>>> SEVNmanager.init({"ce_alpha:2",}) # Initialisation using all the default SEVN paramer except for ce_alpha
>>> ...... # Use tools from the sevn module
>>> SEVNmanager.close() # intialisation

If a unknown SEVN parameters is  included in the parameter list, the SEVNmanager will not be initialised and a
proper error is raised informing the user of the unknown parameters

"""

from .sevnmanager import SEVNmanager
from .star import Star
