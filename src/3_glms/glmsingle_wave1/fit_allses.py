import os
import glob
import numpy as np
import nibabel as nib
import pandas as pd
import matplotlib.pyplot as plt
from glmsingle.design.make_design_matrix import make_design
from glmsingle.glmsingle import GLM_single
import time
import h5py

## constants ----

stimdur = 1  ## stimulus duration (s)
tr = 1.2  ## sec
subjlist = "wave1_unrel"
glmname = "glmsingle_wave1"
sessions = ["baseline", "proactive", "reactive"]
sess = ["Bas", "Pro", "Rea"]
variables = ["target", "distractor", "congruency"]
runs = ["run1", "run2"]
runencs = ["1_AP", "2_PA"]

proj_path = os.getcwd()
path_bluearc = "/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/fMRIPrep_AFNI_ANALYSIS/"
data_path = os.path.join(path_bluearc, '{}', "INPUT_DATA", "Stroop", "{}")  ## subj, session
design_path = os.path.join(proj_path, 'out', 'glms', "{}", "wave1", "RESULTS", "Stroop", glmname)  ## subj
subjs = open(os.path.join("out", "subjlist_" + subjlist + ".txt"), "r").read().split('\n')

## run ----



## loop over variables
## loop over subjects

variable = variables[0]
subj = subjs[0]



## get data and design matrices

data = []
design = []
for session, ses in zip(sessions, sess):
    for i, (run, runenc) in enumerate(zip(runs, runencs)):
        ## read data:
        fname_l = os.path.join(data_path.format(subj, session), "lpi_scale_tfMRI_Stroop" + ses + runenc + "_L.func.gii")
        fname_r = os.path.join(data_path.format(subj, session), "lpi_scale_tfMRI_Stroop" + ses + runenc + "_R.func.gii")
        gii_l = np.array(nib.load(fname_l).agg_data())
        gii_r = np.array(nib.load(fname_r).agg_data())
        y = np.concatenate((gii_l, gii_r), axis = 1).transpose()
        data.append(y)
        ## read design:
        filename = os.path.join(design_path.format(subj), "design_" + variable + ".h5")
        with h5py.File(filename, "r") as f:
            x = np.array(f[session])[i, :, :].transpose()
            design.append(x)


outdir = os.path.join(design_path.format(subj), 'out')
opt = dict()
opt['wantfileoutputs'] = [1, 1, 1, 1]
opt['wantmemoryoutputs'] = [1, 1, 1, 1]
opt['wantlss'] = 0

glmsingle_obj = GLM_single(opt)
start_time = time.time()
results = GLM_single().fit(
    design[0:3],
    data[0:3],
    stimdur,
    tr,
    outputdir = outdir)
elapsed_time = time.time() - start_time
print(
    'elapsedtime: ',
    f'{time.strftime("%H:%M:%S", time.gmtime(elapsed_time))}'
)
