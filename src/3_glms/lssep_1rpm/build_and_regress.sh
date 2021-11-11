#!/usr/bin/env bash



## setup ----

conda deactivate  ## make sure python3 env is not attached

cd $proj_dir
source src/string_funs.sh
function lssfit {
    ## build xmat
    /usr/local/pkg/afni_18/3dLSS \
    -matrix ${in_dir}/X.xmat_run${run}.1D \
    -input ${image} \
    -save1D pinv_run${run}_${hemi}.1D \
    -prefix STATS_${subject}_run${run}_${hemi} \
    -verb
}


## for interactive / dev

if [ $do_single_subj = false ]
then
    mapfile -t subjects < $subjects_file    
else  ## for dev
    unset subjects
    subjects=130518
    subject_i=0
    subject=130518
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
                    
                    for hemi_i in ${!hemis[@]}; do

                        ## get strings
                        
                        wave=${waves[$wave_i]}
                        task=${tasks[$task_i]}
                        session=${sessions[$session_i]}
                        run=${runs[$run_i]}
                        hemi=${hemis[$hemi_i]}
                        sess=$(get_sess $session)
                        run_suffix=$(get_run_suffix $run)
                        wave_bluearc=$(get_wave_bluearc $wave)
                        
                        ## define paths
                        
                        in_dir=${xmat_dir}/${subject}/wave${wave}/RESULTS/${task}/${session}_${xmat}
                        out_dir=${xmat_dir}/${subject}/wave${wave}/RESULTS/${task}/${session}_${glm}
                        image=${bluearc_dir}/${wave_bluearc}/fMRIPrep_AFNI_ANALYSIS/${subject}/INPUT_DATA/${task}/${session}/lpi_scale_tfMRI_${task}${sess}${run_suffix}_${hemi}.func.gii
                        
                        ## check if input directories exist, skip if not
                        
                        if [ ! -d ${in_dir} ]
                        then
                            continue
                        fi

                        ## regress
                        
                        mkdir -p ${out_dir}
                        cd ${out_dir}
                        lssfit < /dev/null > ${out_dir}/runtime_3dLSS_run${run}_${hemi}.log 2>&1 &
                        cd ${proj_dir}

                    done
                    
                done

            done

        done

        wait  ## run subj*waves serially, all else parallel

    done

done



