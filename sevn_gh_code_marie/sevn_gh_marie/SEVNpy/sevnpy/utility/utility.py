import copy
import hashlib
import sys
from ..sevnpy_types import ListLikeType, StandardValue, Dict, List, Tuple, Set, Any, Optional, Union, Number


def copy_dict_and_exclude(dict: Dict, exclude_keys: Optional[ListLikeType] = None) -> Dict:
    """
    Produce a hard copy of a dictionary excluding some keys

    Parameters
    ----------
    dict: dictionary
        Input value for which we want to guess the related pattern matching
    exclude_keys: None or array_like
        collection of keys to exclude from the output dictionary, if None just copy the dictionary

    Returns
    -------
    newdict: dictionary
        A copy of the input dictionary with the keys in exclude_keys removed
    """
    if exclude_keys is None: return copy.deepcopy(dict)

    newdict = copy.deepcopy(dict)
    for key in exclude_keys:
        if key in newdict: newdict.pop(key)

    return newdict


def str_is_float(param:str) -> bool:
    """This function check if the input parameter can be transformed to a float number

    Parameters
    ----------
    param: str
        The string to check

    Returns
    -------
    result: bool
        True if the input string can be transformed to a float, False otherwise

    """

    try:
        float(param)
        return True
    except:
        return False

def str_is_int(param:str) -> bool:
    """This function check if the input parameter can be transformed to a intger number

    Parameters
    ----------
    param: str
        The string to check

    Returns
    -------
    result: bool
        True if the input string can be transformed to a int, False otherwise

    """

    try:
        int(param)
        return True
    except:
        return False

def check_allowed(input: Union[str,float,int],
                  allowed_list:Optional[tuple] = None) -> Union[str,float,int]:
    """
    Check if the input value is in the allowed_list

    Parameters
    ----------
    input: str,float,int
        input value to check
    allowed_list: tuple or None
        tuple of allowed values

    Return
    -------
    result: str|float|int
        Return the input value if it is in the allowed_list or if allowed_list if None otherwise raise  a Value Error


    Raises
    ------
    ValueError
        if the input value is not in the allowed_list

    Examples
    -------
    >>> a="S"
    >>> b=check_allowed("S", allowed_list=("S","B"))
    'S'

    """

    if allowed_list is not None:
        if not input in allowed_list:
            raise ValueError(f"input value {input} not in the allowed list: {allowed_list}")


    return input

def md5hashcode_from_string(string: str,
                            ncut: Optional[int] = None) -> str:
    """
    Generate a md5 hashcode from a string

    Parameters
    ----------
    string: str
        string for which we want to generate an ash code
    ncut: int,None
        if not None this indicates the maximum number or character of the hashcode we want to return (startingt form the first)
        Default=None

    Returns
    -------
    hash_code: str
        The generated hashcode (it can be truncated if ncut is not None)

    """

    hash_object = hashlib.md5(string.encode())
    hash_code   = hash_object.hexdigest()

    if ncut is not None: hash_code=hash_code[:min(ncut,len(hash_code))]

    return hash_code

def check_equality(val1: Number, val2: Number) -> bool:
    """
    Check if two numbers can be considered equal taking into account the machine precision

    Parameters
    ----------
    val1:
        First value
    val2:
        Second value

    Returns
    -------
    result: bool
        If the two values are equal within the machine precision

    """

    eps = sys.float_info.epsilon
    maxXYOne = max(1.0, max(abs(val1),abs(val2)))
    return abs(val1 - val2) <= eps*maxXYOne

