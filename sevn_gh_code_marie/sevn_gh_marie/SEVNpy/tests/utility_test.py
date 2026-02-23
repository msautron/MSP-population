import pytest
import sevnpy.utility as ut
import numpy as np


def test_copy_dict_and_exclude():
    keys = np.arange(np.random.randint(10, 20))
    dict_original = {key: np.random.randint(0, 1) for key in keys}

    exclude_cols = np.random.choice(keys, replace=False, size=np.random.randint(1, len(keys)))

    dict_new = ut.copy_dict_and_exclude(dict_original, exclude_keys=exclude_cols)

    assert dict_original != dict_new

    assert len(dict_new) == (len(dict_original) - len(exclude_cols))
