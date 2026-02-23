# Signifies our desired python version
PYTHON = python
PIP = pip

install:
	${PIP} install .

install_direct:
	${PYTHON} setup.py install

build:
	${PYTHON} setup.py build

build_ext:
	${PYTHON} setup.py build_ext --inplace

tests: install
	pytest tests

doc:
	cd docs && source sphynx_doc_command.sh

clean:
	rm -rf build
	rm -rf dist
	rm -rf *.so
	rm -rf SEVNpy.egg-info
	rm -rf .pytest_cache
	rm -rf docs/source/api/_autosummary