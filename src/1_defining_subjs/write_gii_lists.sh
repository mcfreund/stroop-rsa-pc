#!/usr/bin/env bash

out_dir=/data/nil-external/ccp/freund/stroop-rsa-pc/in/subjs

## DMCC2
cd /data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/HCP_SUBJECTS_BACKUPS/fMRIPrep_AFNI_ANALYSIS
find -wholename '*/INPUT_DATA/Stroop/*func.gii' > $out_dir/giis_dmcc2.txt

## DMCC3
cd /data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/DMCC_Phase3
find -wholename '*/INPUT_DATA/Stroop/*func.gii' > $out_dir/giis_dmcc3.txt

## DMCC4
cd /data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS/DMCC_Phase4
find -wholename '*/INPUT_DATA/Stroop/*func.gii' > $out_dir/giis_dmcc4.txt
