.PHONY: dist test_workflow clean cleanall

all: test_workflow

test: 
	python -m pytest

install-dev:
	python -m pip install -e .

install:
	python -m pip install .

test_workflow:
	cd test_workflow && make

clean: cd test_workflow && make clean

cleanall:
	cd test_workflow && make cleanall

dist:
	python -m build
