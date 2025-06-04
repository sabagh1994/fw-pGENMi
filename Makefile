.PHONY: venv mamba bedtools
SHELL := /bin/bash
PROJBASE := $(shell dirname $(abspath $(lastword $(MAKEFILE_LIST))))
MPIPRECMD := $(shell command -v mpirun >/dev/null 2>&1 && echo "mpirun -n 10")
PROJNAME := fwpgenmi
MAMBABASE := ${PROJBASE}/mamba
PREFENV := venv

#################################################
######## 	 Bedtools Installation       ########
#################################################

# bedtools is used in the preprocessing of TFBS and epigenetic marks
# https://bedtools.readthedocs.io/en/latest/content/installation.html
# bedtools 2.28.0 used in this project, can be accessed from,
# https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools

bedtools:
	mkdir -p software
	wget --directory-prefix=software/ https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz
	cd software
	tar -zxvf bedtools-2.28.0.tar.gz
	cd bedtools2
	make

##########################################################
#####################      Venv     ######################
##########################################################
venv:
	python -m venv venv
	source ./activate venv && python -m pip install --upgrade pip
	source ./activate venv && python -m pip install gdown
	source ./activate venv && python -m pip install -r requirements.txt
	source ./activate venv && python -m pip install -e .
	rm -rf *.egg-info

##########################################################
#####################      Mamba     #####################
##########################################################
mamba:
	mkdir -p ${MAMBABASE}
	cd ${MAMBABASE}; \
	curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba;
	set -e; \
	export MAMBA_ROOT_PREFIX=${MAMBABASE}; \
	eval "$$(${MAMBABASE}/bin/micromamba shell hook -s posix)"; \
	micromamba activate; \
	micromamba create  -y -n ${PROJNAME} python=3.11 -c conda-forge; \
	micromamba activate ${PROJNAME}; \
	micromamba install -y -c conda-forge openssh; \
	export TMPDIR=${PROJBASE}/pip; mkdir -p $${TMPDIR}; \
	python -m pip install --upgrade pip; \
	python -m pip install jupyter; \
	python -m pip install gdown; \
	python -m pip install -r requirements.txt; \
	python -m pip install -e .; \
	rm -rf *.egg-info; \
	rm -r $${TMPDIR};
