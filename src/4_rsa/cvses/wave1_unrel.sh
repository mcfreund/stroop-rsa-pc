atlas_names=("glasser2016" "schaefer2018_17_200")
#roi_cols=("network" "parcel")
roi_cols="network"
sessions=("baseline" "proactive" "reactive")

for session in ${sessions[@]}
do
    echo $session ================================================================================
    for atlas_name in ${atlas_names[@]}
    do

        echo $atlas_name ----------------------------------------
        
        for roi_col in ${roi_cols[@]}
        do
            echo $roi_col ~~~~~~~~~

            Rscript ./src/4_rsa/parcellate_giftis.r \
                --atlas_name $atlas_name \
                --roi_col $roi_col \
                --space "fsaverage5" \
                --glmname "lsall_1rpm" \
                --subjlist "wave1_unrel" \
                --waves "wave1" \
                --sessions $session \
                --n_cores 10 \
                --overwrite "TRUE"\
                --delete_files "TRUE"
                
        done

    done

done

