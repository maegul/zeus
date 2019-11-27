#! /usr/bin/env python3
import os
from sys import argv

if 'CONDA_DEFAULT_ENV' in os.environ:
	conda_env = os.environ['CONDA_DEFAULT_ENV'] 
else: 
	conda_env = None

assert conda_env == 'ephys', (
	f'Conda env must be "ephys", currently: {conda_env}'
	)

from pyth.zeus import athena, hermes



print(argv)
nb_type, exp, unit = argv[1:]

run = argv[4] if len(argv) == 5 else None

proj = athena.load('Athena_lenTun_371c9ed00b2d7a.pkl')

hermes.mk_new_nb(proj, exp, unit, run)

