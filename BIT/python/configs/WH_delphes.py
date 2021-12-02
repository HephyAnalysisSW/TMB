#!/usr/bin/env python

from TMB.BIT.configs.WH_delphes_bkgs import *

training_samples = [WH]

for derivative in bit_derivatives:
    bit_cfg[derivative]['n_trees'] = 120
