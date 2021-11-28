#!/usr/bin/env bash

proj_dir=/data/nil-external/ccp/freund/stroop-rsa-pc
glm="condition_1rpm"  ## name of subdirectory in scripts; will be added to subdirectory name in out
subjects_file=${proj_dir}/out/glms/allsubjs.txt
sessions=(baseline proactive reactive)
tasks=(Stroop)
waves=(1 2)
do_single_subj=false  ## for dev/debugging


source $proj_dir/src/3_glms/${glm}/build.sh