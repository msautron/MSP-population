import pytest
from sevnpy.sevn import SEVNmanager


def test_ID():
    assert SEVNmanager.get_ID() == 0

    SEVNmanager.init()
    assert SEVNmanager.get_ID() == 1
    SEVNmanager.close()

    SEVNmanager.init()
    assert SEVNmanager.get_ID() == 2
    SEVNmanager.close()

    SEVNmanager.init()
    assert SEVNmanager.get_ID() == 3
    SEVNmanager.close()

    with SEVNmanager() as a:
        assert SEVNmanager.get_ID() == 4


def test_initialisation_warning():
    SEVNmanager.init()
    # Calling two init without a close should  raise a warning
    with pytest.warns():
        SEVNmanager.init()

def test_initilisation_wrong_parameter():
    # Initilisation with an unkown parameter must raise an InputError
    with pytest.raises(ValueError):
        SEVNmanager.init({"unkown_parameter":3})


def test_initialise_finalise():
    try:
        SEVNmanager.init()
    except Exception as err:
        pytest.fail(f"SEVNmanager initilisation failed with message {err}")

    try:
        SEVNmanager.close()
    except Exception as err:
        pytest.fail(f"SEVNmanager finalisation failed with message {err}")


def test_double_context_manager():
    with SEVNmanager() as sm:
        assert sm.check_initiliased() == True

    assert SEVNmanager.check_initiliased() == False

    with SEVNmanager() as sm:
        assert sm.check_initiliased() == True
