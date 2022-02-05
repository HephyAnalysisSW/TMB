#!/usr/bin/env python

from TMB.BIT.configs.ZH_delphes_bkgs import *

training_samples = [ZH]

for derivative in bit_derivatives:
    bit_cfg[derivative]['n_trees'] = 120
    bit_cfg[derivative]['min_size'] = 30
