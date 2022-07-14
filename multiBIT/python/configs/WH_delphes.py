#!/usr/bin/env python

from TMB.multiBIT.configs.WH_delphes_bkgs import *

training_samples = [WH]

for derivative in bit_derivatives:
    bit_cfg[derivative]['n_trees'] = 120
    bit_cfg[derivative]['min_size'] = 30
