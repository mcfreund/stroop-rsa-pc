#!/usr/bin/env bash


## setup ----

cd $proj_dir
source src/string_funs.sh
function deconvolve {
    ## build xmat
    /usr/local/pkg/afni_18/3dDeconvolve \
    -local_times \
    -force_TR 1.2 \
    -x1D_stop \
    -allzero_OK \
    -input ${image} \
    -polort A \
    -float \
    -censor ${in_dir}/movregs_FD_mask_run${run}.txt \
    -num_stimts 1 \
    -stim_times_IM 1 ${in_dir}/${subject}_wave${wave}_Stroop_${session}_alltrials_shift600_run${run}.txt 'BLOCK(1,1)' -stim_label 1 alltrials \
    -ortvec ${in_dir}/Movement_Regressors_${task}${sess}${run_suffix}.1D movregs \
    -x1D ${out_dir}/X.xmat_run${run}.1D \
    -xjpeg ${out_dir}/X_run${run}.jpg \
    -nobucket
}


## for interactive / dev

if [ $do_single_subj = false ]
then
    mapfile -t subjects < $subjects_file    
else  ## for dev
    unset subjects
    subjects=130518
    subject_i=0
    wave=1
    task_i=0
    session_i=0
    run_i=0
fi



## execute ----

for subject_i in ${!subjects[@]}; do

    subject=${subjects[$subject_i]}
    echo ${subject}

    for wave_i in ${!waves[@]}; do

        for task_i in ${!tasks[@]}; do

            for session_i in ${!sessions[@]}; do

                for run_i in ${!runs[@]}; do
                    
                    ## get strings
                    
                    wave=${waves[$wave_i]}
                    task=${tasks[$task_i]}
                    session=${sessions[$session_i]}
                    run=${runs[$run_i]}
                    sess=$(get_sess $session)
                    run_suffix=$(get_run_suffix $run)
                    wave_bluearc=$(get_wave_bluearc $wave)
                    
                    ## define paths
                    
                    in_dir=${xmat_dir}/${subject}/wave${wave}/INPUT_DATA/${task}/${session}
                    out_dir=${xmat_dir}/${subject}/wave${wave}/RESULTS/${task}/${session}_${glm}
                    image=${bluearc_dir}/${wave_bluearc}/fMRIPrep_AFNI_ANALYSIS/${subject}/INPUT_DATA/${task}/${session}/lpi_scale_tfMRI_${task}${sess}${run_suffix}_L.func.gii  ## either hemi OK
                    
                    ## check if input directories exist, skip if not
                    
                    if [ ! -d ${in_dir} ]
                    then
                        continue
                    fi

                    ## build xmat

                    mkdir -p ${out_dir}
                    cd ${out_dir}
                    deconvolve < /dev/null > ${out_dir}/runtime_3dDeconvolve.log 2>&1 &
                    cd ${proj_dir}

                done

            done

        done

    done

    wait  ## run subjs serially, all else parallel

done



