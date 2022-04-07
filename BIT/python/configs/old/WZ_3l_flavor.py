#!/usr/bin/env python

from TMB.BIT.configs.ttZ_3l_flavor import *

# training selection
from tWZ.Tools.cutInterpreter import cutInterpreter
selectionString = cutInterpreter.cutString( 'trilepT-onZ1' )

from tWZ.samples.nanoTuples_Autumn18_nanoAODv6_private_SMEFTsim_fast_postProcessed import *
training_samples = [ WZ ]

bit_derivatives  = [ ('cHq1Re11',), ('cHq1Re22',), ('cHq1Re33',), ('cHq1Re11','cHq1Re11'), ('cHq1Re22','cHq1Re22'), ('cHq1Re33','cHq1Re33'),
                     ('cHq3Re11',), ('cHq3Re22',), ('cHq3Re33',), ('cHq3Re11','cHq3Re11'), ('cHq3Re22','cHq3Re22'), ('cHq3Re33','cHq3Re33')]

assert len(training_samples)==len(set([s.name for s in training_samples])), "training_samples names are not unique!"
