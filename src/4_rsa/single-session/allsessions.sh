## looping variables
atlas_names=("glasser2016" "schaefer2018_17_200")
ttype_subsets=("bias" "pc50")
sessions=("baseline" "proactive" "reactive")
waves=("wave1" "wave2")
roi_cols=("parcel" "superparcel")

## constants
measure="crcor"
glmname="lsall_1rpm"
prewh="none"
n_cores=24
overwrite="FALSE"
delete_files="FALSE"
n_resamples=10000
space="fsaverage5"


## parcellate ----

for atlas_name in ${atlas_names[@]}
do

    echo --------------------- $atlas_name ---------------------

    for roi_col in ${roi_cols[@]}
    do

        if [[ ($atlas_name == "schaefer2018_17_200" && $roi_col == "superparcel") ]]
        then 
            continue
        fi

        echo --------------------- $roi_col:

        for wave in ${waves[@]}
        do
            for session in ${sessions[@]}
            do

            subjlist=${session}_$wave
            
            Rscript ./src/4_rsa/parcellate_giftis.r \
                --glmname $glmname \
                --atlas_name $atlas_name \
                --space $space \
                --roi_col $roi_col \
                --subjlist $subjlist \
                --waves $wave \
                --sessions $session \
                --n_cores $n_cores \
                --overwrite $overwrite \
                --delete_files $delete_files

            done
        done
    done    
done


## within-session RSA ----

## estimate distances:

for atlas_name in ${atlas_names[@]}
do

    echo --------------------- $atlas_name ---------------------

    for roi_col in ${roi_cols[@]}
    do

        if [[ ($atlas_name == "schaefer2018_17_200" && $roi_col == "superparcel") ]]
        then 
            continue
        fi

        echo --------------------- $roi_col:

        for wave in ${waves[@]}
        do
            for session in ${sessions[@]}
            do

                subjlist=${session}_$wave

                for ttype_subset in ${ttype_subsets[@]}
                do
                    
                    echo -------- $wave $session $ttype_subset:

                    Rscript ./src/4_rsa/estimate_distances.r \
                        --glmname $glmname \
                        --atlas_name $atlas_name \
                        --space $space \
                        --roi_col $roi_col \
                        --subjlist $subjlist \
                        --waves $wave \
                        --sessions $session \
                        --measure $measure \
                        --prewh $prewh \
                        --ttype_subset $ttype_subset \
                        --n_cores $n_cores \
                        --n_resamples $n_resamples \
                        --overwrite $overwrite
                    
                    Rscript ./src/4_rsa/regress_distances.r \
                        --glmname $glmname \
                        --atlas_name $atlas_name \
                        --space $space \
                        --roi_col $roi_col \
                        --subjlist $subjlist \
                        --waves $wave \
                        --sessions $session \
                        --measure $measure \
                        --prewh $prewh \
                        --ttype_subset $ttype_subset \
                        --suffix "__seswave-"$session"_"$wave

                done
            done
        done
    done    
done



# measure="crcor"
# wave="wave2"
# atlas_name="schaefer2018_17_200"
# ttype_subsets=("bias" "pc50")
# sessions=("proactive" "reactive")
# roi_col="parcel"
# space="fsaverage5"
# prewh="none"
# for session in ${sessions[@]}
# do
#     subjlist=${session}_$wave
#     for ttype_subset in ${ttype_subsets[@]}
#     do
#         echo -------- $wave $session $ttype_subset:  
#         Rscript ./src/4_rsa/regress_distances.r \
#             --glmname $glmname \
#             --atlas_name $atlas_name \
#             --space $space \
#             --roi_col $roi_col \
#             --subjlist $subjlist \
#             --waves $wave \
#             --sessions $session \
#             --measure $measure \
#             --prewh $prewh \
#             --ttype_subset $ttype_subset \
#             --suffix "__seswave-"$session"_"$wave
#     done
# done