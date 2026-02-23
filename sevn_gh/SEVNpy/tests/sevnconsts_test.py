import pytest
from sevnpy.sevn import sevnconst as sc


def test_G():
    assert sc.G.value() == pytest.approx(3.925125598496094e8, 0.01)


def test_Rsuncgs():
    assert sc.Rsun.value() == pytest.approx(6.95700e10, 0.01)
