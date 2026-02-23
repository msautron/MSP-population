import pytest


def test_compilation():
    try:
        from sevnpy.sevn import sevnwrap
    except ImportError:
        pytest.fail(f"sevnwrap has not been compiled")


def test_initialisation():
    from sevnpy.sevn import sevnwrap

    try:
        sevnwrap.sevnio_initialise()
    except Exception as err:
        pytest.fail(f"SEVNmanager initilisation failed with message {err}")

    try:
        sevnwrap.sevnio_finalise()
    except Exception as err:
        pytest.fail(f"SEVNmanager initilisation failed with message {err}")


def test_star_evolve():
    from sevnpy.sevn import sevnwrap

    sevnwrap.sevnio_initialise()

    df, _ = sevnwrap.evolve_star(10, 0.02, tstart="zams", tend=2)

    sevnwrap.sevnio_finalise()

    assert df["Worldtime"][-1] == 2

def test_binary_evolve():
    from sevnpy.sevn import sevnwrap

    sevnwrap.sevnio_initialise()

    df, _ = sevnwrap.evolve_binary(1000,0.2,20,0.02,10,0.01,tend=2)

    sevnwrap.sevnio_finalise()

    assert df["Worldtime"][-1] == 2



