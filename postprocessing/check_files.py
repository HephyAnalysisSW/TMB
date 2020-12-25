import os

# dir for postprocessed jobs
gen_dir = "/afs/hephy.at/data/llechner01/TTXPheno/skims/gen/v3/"
#gen_dir = "/afs/hephy.at/data/rschoefbeck01/TTXPheno/skims/gen/v3/"

# number of files expected (splitted files)
nJobs = 200
order = 2

# define samples to check if all jobs have finished postprocessing
samples = [
#    'fwlite_ttZ_ll_LO_order3_8weights',
#    'fwlite_ttW_LO_order3_8weights',
#    'fwlite_ttgamma_LO_order3_8weights', 

#    'fwlite_ttZ_ll_LO_order2_15weights',
#    'fwlite_ttW_LO_order2_15weights',
#    'fwlite_ttgamma_LO_order2_15weights', 
#
#    'fwlite_ttZ_ll_LO_order2_15weights_ref',
#    'fwlite_ttW_LO_order2_15weights_ref',
#    'fwlite_ttgamma_LO_order2_15weights_ref', 

#    'fwlite_tt_LO_order2_15weights_ref',
#    'fwlite_wZ_LO_order2_15weights',

#    'fwlite_tt_lep_LO_order2_15weights_ref',
#    'fwlite_WZ_lep_LO_order2_15weights',
]


# check all sample dirs
for sample in samples:
    dir = os.path.join(gen_dir, sample)

    # raise error if sample directory not found
    if not os.path.isdir(dir):
        raise ValueError('Directory not found: %s' %sample)

    # get indices of all files in sample directory
    fileIndices = [int(file.rstrip('.root').split('_')[-1]) for file in os.listdir(dir + '/') if file != 'delphes' and file.split('.')[-1] != 'log' and file.split('.')[-1] != 'rw_log']
    fileIndices = sorted(fileIndices)

    # list of expected indices
    indices = list(range(nJobs))

    # list of not processed files
    nonProcessedFiles = [item for item in indices if item not in fileIndices]

    print(sample)

    # run missing files locally or on batch
    for job in nonProcessedFiles:
        print(job)
#        os.system('python genPostProcessing.py --nJobs %s --job %s --overwrite --logLevel DEBUG --sample %s --addReweights --interpolationOrder %s --delphes > logs/%s_%s.log 2>&1' %(nJobs, job, sample, order, sample, job))
#        os.system('submitBatch.py --dpm "./genPostProcessing.py --nJobs %s --job %s --overwrite --logLevel DEBUG --sample %s --addReweights --interpolationOrder %s --delphes"' %(nJobs, job, sample, order, sample, job))

