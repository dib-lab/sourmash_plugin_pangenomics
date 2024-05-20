.PHONY: dist test_workflow

all: test_workflow

test: 
	python -m pytest

install-dev:
	python -m pip install -e .

install:
	python -m pip install .

test_workflow:
	cd test_workflow && make

cleanrun:
	cd test_workflow && make cleanall

dist:
	python -m build
