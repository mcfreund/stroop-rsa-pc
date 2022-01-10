## looping variables
atlas_names=("glasser2016" "schaefer2018_17_200")
ttype_subsets=("bias" "pc50" "all")
measures=("crcor")
subjlists=("mi1" "mc1")
sessions=("proactive" "baseline")

## constants
wave=("wave1")
glmname="lsall_1rpm"
roi_col="parcel"
prewh="none"
n_cores=28
overwrite="TRUE"
n_resamples=10000
space="fsaverage5"


for i in ${!subjlists[@]}
do
    
    subjlist=${subjlists[$i]}
    session=${sessions[$i]}

    for ttype_subset in ${ttype_subsets[@]}
    do

        if [[ ($ttype_subset == "bias" || $ttype_subset == "pc50") && $session == "proactive" ]]
        then 
            continue
        fi

        echo $session $ttype_subset ===============================================================

        for atlas_name in ${atlas_names[@]}
        do
            
            echo $atlas_name $(date)  ------------------

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

