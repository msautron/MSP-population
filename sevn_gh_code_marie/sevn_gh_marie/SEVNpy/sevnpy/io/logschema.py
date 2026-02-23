"""
==============================================================
Logfile Handling, (:mod:`sevny.io.logschema`)
==============================================================

This module contains the class LogSchema that is used to bring at a code-level the SEVN-like log structure and the interface with regex pattern.
The SEVN-like log structure is defined as a series of  items of different  types.
The class LogSchema define each item with five attributes:

    - name
    - kind
    - type
    - pattern
    - description


"""

import warnings
from typing import Dict, Union, List, Type, Tuple, Literal

try:
    from typing import TypedDict
except ImportError:
    from typing_extensions import TypedDict

from .regexutility import ReTypeMatch, capturing, notcapturing
from .. import utility as ut

GTUnion = Union[Type[int], Type[float], Type[float]]
IPDict = Dict[str, Union[Union[str, GTUnion], Tuple[Union[str, GTUnion], str]]]


class LogSchema:
    """
    This class is used to define and handle the structure of a given log file and produce the correspondent matching pattern
    There are four important keyword used in the class:

        - **name:** is the label used to identify the data stored in a given position of the log structure, this can be
                    used as the name of the column on a dataframe.
        - **description:** short description of the data stored at a given name
        - **kind:** is the type of the data identified by the name, it could be "int","id","float","name",int,str,float
                    or a string. If it is a string (but not "int","id","float","name","type"), it is called *named_kind*
                    and the kind specifies exactly the searching pattern otherwise the searching pattern will defined by the value stored in kind.
                    For example if kind="S", the searching pattern will contain exactly "S",  but if kind="id" the
                    searching pattern will be "[0-9]+"
        - **pattern:** is the regex searching pattern associated with name and kind
        - **type:** data type associated with each name, it depends on the kind value:

                    - "id", "type", "int" or int -> int
                    - "float" or float -> float
                    - "name" -> str
                        In all the other cases, i.e. when kind is a string (but not "int","id","float","name","type"),
                        the matching pattern will be inferred:
                    - If the string can be transformed to an integer -> int
                    - If the string can be transformed to a float -> float
                    - Otherwise -> str




    So, each element in a LogSchema is defined by the name, kind and pattern. Name and kind are defined by the user
    at the class instantiation  or using the method add_item (see example below).  The pattern and the type
    is instead generated  by the class based on the kind value (see above).
    Actually there is a fifth element that is the regex_pattern that is estimated only  when the method
    regex_pattern is called and it  includes the pattern + the extra paranthesis needed to create a capturing or
    non capturing group based on method input


    Examples
    ----------
    The class is used to define an obejct containing the information about the structure of a SEVN-like logfile,
    for example assume that we have the following log structure
    B;<name>;<id>;CIRC;<tiem>;<semimajor_axis_ini>:<eccentricity_ini>:<semimajor_axis_post>:eccentricity_post>
    e.g.: B;857175750378006;0;CIRC;1.874849e+01;38.2411:0.000633313:38.2411:0

    We can use initialize a LogSchema for the header like

    >>> header=LogSchema({"logtype":("B",""), "name":("name","unique identifier"), "ID":("id",""),"event":("CIRC",""),"Worldtime":(float,"time in Myr")})

    While for the  body we can use another LogSchema, let'start from a empty initilisation

    >>> body=LogSchema()
    >>> body.add_item("semimajor_axis_ini",float,"pre circularisation semimajor axis Rsun")
    >>> body.add_item("eccentriciy_ini",float,"pre circularisation eccentricity")
    >>> body.add_item("semimajor_axis_post",float,"post circularisation semimajor axis Rsun")
    >>> body.add_item("eccentriciy_post",float,"pre circularisation eccentricity")

    If we want to check the different schema

    >>> body.kind_schema
    >>> {'semimajor_axis_ini': <class 'float'>, 'eccentriciy_ini': <class 'float'>, 'semimajor_axis_post': <class 'float'>, 'eccentriciy_post': <class 'float'>}
    >>> body.type_schema
    >>> {'semimajor_axis_ini': <class 'float'>, 'eccentriciy_ini': <class 'float'>, 'semimajor_axis_post': <class 'float'>, 'eccentriciy_post': <class 'float'>}
    >>> body.pattern_schema
    >>> {'semimajor_axis_ini': '[+|-]?[0-9]+\\.?[0-9]*(?i:e)?[+|-]?[0-9]*|(?i:nan)', 'eccentriciy_ini': '[+|-]?[0-9]+\\.?[0-9]*(?i:e)?[+|-]?[0-9]*|(?i:nan)', 'semimajor_axis_post': '[+|-]?[0-9]+\\.?[0-9]*(?i:e)?[+|-]?[0-9]*|(?i:nan)', 'eccentriciy_post': '[+|-]?[0-9]+\\.?[0-9]*(?i:e)?[+|-]?[0-9]*|(?i:nan)'}

    Finally, if we want to get the regex pattern considering only the semimajor axis properteis as capturing items

    >>> body.regex_pattern(capturing_names=('semimajor_axis_ini','semimajor_axis_post'))
    >>> {'semimajor_axis_ini': '([+|-]?[0-9]+\\.?[0-9]*(?i:e)?[+|-]?[0-9]*|(?i:nan))', 'eccentriciy_ini': '(?:[+|-]?[0-9]+\\.?[0-9]*(?i:e)?[+|-]?[0-9]*|(?i:nan))', 'semimajor_axis_post': '([+|-]?[0-9]+\\.?[0-9]*(?i:e)?[+|-]?[0-9]*|(?i:nan))', 'eccentriciy_post': '(?:[+|-]?[0-9]+\\.?[0-9]*(?i:e)?[+|-]?[0-9]*|(?i:nan))'}

    """

    def __init__(self,
                 schema: IPDict = {}):
        """

        Parameters
        ----------
        schema:
            dictionary containing the schema with the triple name:kind:description
        """
        # Auxiliary dictionary to map from kind to type
        self._kind_to_type_map = {
            "int": int,
            "id": int,
            "type": int,
            int: int,
            float: float,
            "float": float,
            "name": str,
            str: str,
        }

        # Check input
        self._check_input_dict(schema)

        # Fundametal schema
        self._schema = schema

        # These are derived schemas
        self._kind_schema = None
        self._type_schema = None
        self._pattern_schema = None
        self._description_schema = None


    class _ItemDict(TypedDict):
        """
        Auxiliary TypedDict containing the full keywords for a given name

        """
        name: str
        description: str
        kind: Union[GTUnion, str]
        type: GTUnion
        pattern: str

    def _check_input(self, name: str, kind: Union[str, GTUnion]):
        # Check if the kind is in the allowed kinds
        if isinstance(kind, str) == False and kind not in self._kind_to_type_map:
            raise ValueError(f"kind {kind} not allowed. Possible values are {list(self._kind_to_type_map.keys())}")
        else:
            pass

    def _check_input_dict(self, input_dict: IPDict):
        """
        Check if the input dictionary can be used as a item for the LogSchema.
        If it cannot, the method raises an error.
        Parameters
        ----------
        input_dict:
            input dictionary


        """
        # First check if the dictionary is in the right format
        for key in input_dict:
            if not isinstance(key, str): raise ValueError(f"input key {key} is not a string this is not allowed")
            # second try if the value is a len 2 collection, in case just add an empty description
            if not len(input_dict[key]) == 2: raise ValueError(
                f"Item of key {key} is not a collection with two elements")

        # check if all the kind in the input schema are present in the allowed kinds.
        # Otherwise, return a meaningful error
        not_allowed_names = []
        not_allowed_kinds = set()
        for name, (kind, _) in input_dict.items():
            if isinstance(kind, str) == False and kind not in self._kind_to_type_map:
                not_allowed_names.append(name)
                not_allowed_kinds.add(kind)

        if len(not_allowed_names) > 0:
            err_mess = f"The following names {not_allowed_names} contain not allowed kinds:{not_allowed_kinds}."
            err_mess += f"The allowed links are: {list(self._kind_to_type_map.keys())}"
            raise ValueError(err_mess)

    def update(self, input_dict: IPDict):
        # Check input
        self._check_input_dict(input_dict)
        # Update
        self._schema.update(input_dict)
        # Reset derivative schema so that they are updated at the first new call
        self._reset_derivate_schema()

    def add_item(self, name: str, kind: Union[str, GTUnion], description: str = ""):
        """
        Add an item to the Schema

        Parameters
        ----------
        name:
            name of the new item
        kind:
            kind of the new item
        description:
            a short description of the item

        """
        # Check input
        self._check_input(name, kind)

        self._schema.update({name: (kind, description)})
        # Reset derivative schema so that they are updated at the first new call
        self._reset_derivate_schema()

    def _kind_to_type(self, kind: Union[str, int, float]) -> Union[Type[int], Type[float], Type[str]]:
        """
        Get a type from the kind
        Parameters
        ----------
        kind:
            kind param
        Returns
        -------
        type: int|float|str
            the type
        """
        if kind in self._kind_to_type_map:
            return self._kind_to_type_map[kind]
        elif ut.str_is_int(kind):
            return int
        elif ut.str_is_float(kind):
            return float
        else:
            return str

    def _kind_to_pattern(self, kind: Union[str, int, float]) -> str:
        """
        Get a pattern  from the type
        If kind is a string (but not "int","id","float","name"), the kind specify exactly the
        searching pattern otherwise the searching pattern will be defined by the value stored in kind.

        Parameters
        ----------
        kind:
            the kind

        Returns
        -------
        pattern: str
            the regex pattern
        """
        if kind in self._kind_to_type_map:
            return ReTypeMatch[kind]
        else:
            return kind

    @property
    def schema(self) -> IPDict:
        """
        The schema of the log reader

        """
        return self._schema

    @property
    def kind_schema(self) -> Dict:
        """
        dictionary containing the pair name:kind

        """
        if self._kind_schema is None:
            self._kind_schema = {name: kind for name, (kind, _) in self._schema.items()}
        return self._kind_schema

    @property
    def description_schema(self) -> Dict:
        """
d       Dctionary containing the pair name:description

        """
        if self._description_schema is None:
            self._description_schema = {name: desc for name, (_, desc) in
                                        self._schema.items()}

        return self._description_schema

    @property
    def type_schema(self):
        # Check if empty
        if self._type_schema is None:
            self._type_schema = {name: self._kind_to_type(kind) for name, kind in self.kind_schema.items()}
        return self._type_schema

    @property
    def pattern_schema(self):
        """
        Dictionary containing the pair name:type

        """
        # Check if empty
        if self._pattern_schema is None:
            self._pattern_schema_dict = {name: self._kind_to_pattern(kind)
                                         for name, kind in self.kind_schema.items()}

        return self._pattern_schema_dict

    def default_capturing_names(self) -> List[str]:
        """
        Get the default capturing name, i.e. the item with a kind that is a string but not "name","id","str","float","int".
        Returns
        -------
        capturing_names_temp: List
            A list with the default capturing names (already sorted)
        """
        # Lambda function to check if is a "named" kind, i.e. a kind that is a string but not "name","id","int","float"
        _check_named_kind = lambda _kind: _kind == self._kind_to_pattern(_kind)

        # Cycle over the kind schema (name:kind) to check if it is a named item or not.
        # If it is not a named item include it in the capturing names
        _capturing_names_temp = []
        for name, kind in self.kind_schema.items():
            if not _check_named_kind(kind): _capturing_names_temp.append(name)

        return _capturing_names_temp

    def regex_pattern(self,
                      capturing_names: Union[Literal["default", "all"], List[str], Tuple[str, ...]] = "default") -> \
    Tuple[Dict[str, str], List[str]]:
        """

        Parameters
        ----------
        capturing_names: "default","all", iterable
            A list of names to be included as capturing members.
            If the string "default" is used all the items in the  schema will be captured except for the named item, i.e.
            the item with a kind that is a string but not "name","id","str","float","int".
            If the string "all" is used all the  items  in the schema will be captured

        Returns
        -------
        regex_pattern: Dictionary
            a string containing the pair name:regex_pattern
        capturing_names: List
            A sorted list of the capturing names

        """

        if capturing_names == "default":
            capturing_names = self.default_capturing_names()
        elif capturing_names == "all":
            # Include all the names in the capturing names
            capturing_names = self.names

        elif not isinstance(capturing_names, str):
            # Use the user provided list of names
            # First check if it is actually a collection
            try:
                len(capturing_names)
            except TypeError:
                raise TypeError(f"The input parameter capturing_names need to be a Collectable, or the string \'all\'"
                                f"or the string \'default\', it is instead {capturing_names}")

            # Chek if a name in the capturing names are not present in the schema and raise a warning.
            _capturing_names = []
            for name in capturing_names:
                if name not in self.names:
                    warnings.warn(f"Name {name} is not present in the schema names. "
                                  f"This capturing name will not be included")
                else:
                    _capturing_names.append(name)
            capturing_names = _capturing_names

        _searching_pattern = {}
        _sorted_capturing_names = []
        for name, pattern in self.pattern_schema.items():
            if name in capturing_names:
                _searching_pattern[name] = capturing(pattern)
                _sorted_capturing_names.append(name)
            else:
                _searching_pattern[name] = notcapturing(pattern)

        return (_searching_pattern, _sorted_capturing_names)

    def reset(self):
        """
        Reset the schema

        """
        self._schema.clear()

    def _reset_derivate_schema(self):
        """
        Reset the other dictionary derived from the main schema

        """
        self._pattern_schema.clear()
        self._type_schema.clear()
        self._kind_schema.clear()
        self._description_schema.clear()

    def pop(self, name: str):
        """
        Remove an item of given name form the schema

        Parameters
        ----------
        name:
            Name of the item in the schema to remove

        """
        # Remove the item with a given name
        self._schema.pop(name)


    @property
    def names(self) -> List[str]:
        """
        List of al the names

        """
        # Return all the names
        return list(self._schema.keys())

    @property
    def kinds(self) -> List[Union[str, GTUnion]]:
        """
        List of all the kinds

        """
        # Return all the kinds
        return list(self.kind_schema.values())

    @property
    def types(self) -> List[GTUnion]:
        """
        List of all the types

        """
        # Return all the types
        return list(self.type_schema.values())

    @property
    def patterns(self) -> List[str]:
        """
        List of all the patterns

        """
        # Return all the patterns
        return list(self.pattern_schema.values())

    @property
    def descriptions(self) -> List[str]:
        """
        List of all the descriptions

        """

        return list(self.description_schema.values())

    def column_schema(self, offset: int = 0) -> Dict[int, str]:
        """
        Return a dictionary to map the index of the column in the logschema to the correspondent name
        Parameters
        ----------
        offset: int
            Offset to use to the define the columns

        Returns
        -------
        column_schema: Dictionary
            A Dictionary containg index:name pairs

        """
        return {idx + offset: name for idx, name in enumerate(self.names)}

    def full_schema(self) -> Dict[str, _ItemDict]:
        """
        Create and return a dictionary containing all the information about the schema

        Returns
        -------
        full_schema: Dictionary
            A dictionary in which each item is a pair name:dictionary structured as follow
            {<item_name>: {"name":<item_name>, "description":<item_description>, "kind":<item_kind>, "type":<item_type>, "pattern":<item_pattern>}}

        """
        return {name: self.__getitem__(name) for name in self.names}

    def __getitem__(self, locator: Union[str, int]) -> _ItemDict:

        if isinstance(locator, int):
            try:
                locator = self.column_schema()[locator]
            except KeyError:
                raise KeyError(
                    f"Column index {locator} do not have a correspondence in the LogSchema {self.column_schema()}")

        try:
            return {"name": locator,
                    "kind": self.kind_schema[locator],
                    "description": self.description_schema[locator],
                    "type": self.type_schema[locator],
                    "pattern": self.pattern_schema[locator]}
        except KeyError:
            raise KeyError(f"Name \'{locator}\' not available.Available names are {self.names}")

    def __len__(self):
        return len(self.names)

    def __repr__(self) -> str:
        return str(self.full_schema())

    def __str__(self) -> str:

        fs = ""
        for key1, int_dic in self.full_schema().items():
            fs += f"\n{key1}:\n   "
            for key2, value in int_dic.items():
                fs += f"{key2}={value}, "

        return fs
