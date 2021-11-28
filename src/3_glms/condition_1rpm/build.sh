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
    -num_stimts 32 \
    -stim_times 1 ${in_dir}/${subject}_wave${wave}_Stroop_${session}_blueBLUE_shift600_run${run}.txt 'BLOCK(1,1)' -stim_label 1 blueBLUE \
    -stim_times 2 ${in_dir}/${subject}_wave${wave}_Stroop_${session}_bluePURPLE_shift600_run${run}.txt 'BLOCK(1,1)' -stim_label 2 bluePURPLE \
    -stim_times 3 ${in_dir}/${subject}_wave${wave}_Stroop_${session}_blueRED_shift600_run${run}.txt 'BLOCK(1,1)' -stim_label 3 blueRED \
    -stim_times 4 ${in_dir}/${subject}_wave${wave}_Stroop_${session}_blueWHITE_shift600_run${run}.txt 'BLOCK(1,1)' -stim_label 4 blueWHITE \
    -stim_times 5 ${in_dir}/${subject}_wave${wave}_Stroop_${session}_purpleBLUE_shift600_run${run}.txt 'BLOCK(1,1)' -stim_label 5 purpleBLUE \
    -stim_times 6 ${in_dir}/${subject}_wave${wave}_Stroop_${session}_purplePURPLE_shift600_run${run}.txt 'BLOCK(1,1)' -stim_label 6 purplePURPLE \
    -stim_times 7 ${in_dir}/${subject}_wave${wave}_Stroop_${session}_purpleRED_shift600_run${run}.txt 'BLOCK(1,1)' -stim_label 7 purpleRED \
    -stim_times 8 ${in_dir}/${subject}_wave${wave}_Stroop_${session}_purpleWHITE_shift600_run${run}.txt 'BLOCK(1,1)' -stim_label 8 purpleWHITE \
    -stim_times 9 ${in_dir}/${subject}_wave${wave}_Stroop_${session}_redBLUE_shift600_run${run}.txt 'BLOCK(1,1)' -stim_label 9 redBLUE \
    -stim_times 10 ${in_dir}/${subject}_wave${wave}_Stroop_${session}_redPURPLE_shift600_run${run}.txt 'BLOCK(1,1)' -stim_label 10 redPURPLE \
    -stim_times 11 ${in_dir}/${subject}_wave${wave}_Stroop_${session}_redRED_shift600_run${run}.txt 'BLOCK(1,1)' -stim_label 11 redRED \
    -stim_times 12 ${in_dir}/${subject}_wave${wave}_Stroop_${session}_redWHITE_shift600_run${run}.txt 'BLOCK(1,1)' -stim_label 12 redWHITE \
    -stim_times 13 ${in_dir}/${subject}_wave${wave}_Stroop_${session}_whiteBLUE_shift600_run${run}.txt 'BLOCK(1,1)' -stim_label 13 whiteBLUE \
    -stim_times 14 ${in_dir}/${subject}_wave${wave}_Stroop_${session}_whitePURPLE_shift600_run${run}.txt 'BLOCK(1,1)' -stim_label 14 whitePURPLE \
    -stim_times 15 ${in_dir}/${subject}_wave${wave}_Stroop_${session}_whiteRED_shift600_run${run}.txt 'BLOCK(1,1)' -stim_label 15 whiteRED \
    -stim_times 16 ${in_dir}/${subject}_wave${wave}_Stroop_${session}_whiteWHITE_shift600_run${run}.txt 'BLOCK(1,1)' -stim_label 16 whiteWHITE \
    -stim_times 17 ${in_dir}/${subject}_wave${wave}_Stroop_${session}_blackBLACK_shift600_run${run}.txt 'BLOCK(1,1)' -stim_label 17 blackBLACK \
    -stim_times 18 ${in_dir}/${subject}_wave${wave}_Stroop_${session}_blackGREEN_shift600_run${run}.txt 'BLOCK(1,1)' -stim_label 18 blackGREEN \
    -stim_times 19 ${in_dir}/${subject}_wave${wave}_Stroop_${session}_blackPINK_shift600_run${run}.txt 'BLOCK(1,1)' -stim_label 19 blackPINK \
    -stim_times 20 ${in_dir}/${subject}_wave${wave}_Stroop_${session}_blackYELLOW_shift600_run${run}.txt 'BLOCK(1,1)' -stim_label 20 blackYELLOW \
    -stim_times 21 ${in_dir}/${subject}_wave${wave}_Stroop_${session}_greenBLACK_shift600_run${run}.txt 'BLOCK(1,1)' -stim_label 21 greenBLACK \
    -stim_times 22 ${in_dir}/${subject}_wave${wave}_Stroop_${session}_greenGREEN_shift600_run${run}.txt 'BLOCK(1,1)' -stim_label 22 greenGREEN \
    -stim_times 23 ${in_dir}/${subject}_wave${wave}_Stroop_${session}_greenPINK_shift600_run${run}.txt 'BLOCK(1,1)' -stim_label 23 greenPINK \
    -stim_times 24 ${in_dir}/${subject}_wave${wave}_Stroop_${session}_greenYELLOW_shift600_run${run}.txt 'BLOCK(1,1)' -stim_label 24 greenYELLOW \
    -stim_times 25 ${in_dir}/${subject}_wave${wave}_Stroop_${session}_pinkBLACK_shift600_run${run}.txt 'BLOCK(1,1)' -stim_label 25 pinkBLACK \
    -stim_times 26 ${in_dir}/${subject}_wave${wave}_Stroop_${session}_pinkGREEN_shift600_run${run}.txt 'BLOCK(1,1)' -stim_label 26 pinkGREEN \
    -stim_times 27 ${in_dir}/${subject}_wave${wave}_Stroop_${session}_pinkPINK_shift600_run${run}.txt 'BLOCK(1,1)' -stim_label 27 pinkPINK \
    -stim_times 28 ${in_dir}/${subject}_wave${wave}_Stroop_${session}_pinkYELLOW_shift600_run${run}.txt 'BLOCK(1,1)' -stim_label 28 pinkYELLOW \
    -stim_times 29 ${in_dir}/${subject}_wave${wave}_Stroop_${session}_yellowBLACK_shift600_run${run}.txt 'BLOCK(1,1)' -stim_label 29 yellowBLACK \
    -stim_times 30 ${in_dir}/${subject}_wave${wave}_Stroop_${session}_yellowGREEN_shift600_run${run}.txt 'BLOCK(1,1)' -stim_label 30 yellowGREEN \
    -stim_times 31 ${in_dir}/${subject}_wave${wave}_Stroop_${session}_yellowPINK_shift600_run${run}.txt 'BLOCK(1,1)' -stim_label 31 yellowPINK \
    -stim_times 32 ${in_dir}/${subject}_wave${wave}_Stroop_${session}_yellowYELLOW_shift600_run${run}.txt 'BLOCK(1,1)' -stim_label 32 yellowYELLOW \
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



