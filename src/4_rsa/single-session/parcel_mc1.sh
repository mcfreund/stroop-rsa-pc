atlas_names=("glasser2016" "schaefer2018_17_200")
ttype_subsets=("bias" "pc50")

measures=("crcor")
subjlist=("mc1")
wave=("wave1")
session=("baseline")
glmname="lsall_1rpm"
roi_col="parcel"
prewh="none"
n_cores=28
overwrite="TRUE"
n_resamples=10000
space="fsaverage5"

for atlas_name in ${atlas_names[@]}
do

    echo $atlas_name
    
    Rscript ./src/4_rsa/parcellate_giftis.r \
        --glmname $glmname \
        --atlas_name $atlas_name \
        --space $space \
        --roi_col $roi_col \
        --subjlist $subjlist \
        --waves $wave \
        --sessions $session \
        --n_cores $n_cores \
        --overwrite $overwrite

done


for measure in ${measures[@]}
do

    echo $measure

    for atlas_name in ${atlas_names[@]}
    do

        echo $atlas_name

        for ttype_subset in ${ttype_subsets[@]}
        do

            echo estimating $ttype_subset $(date)
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
            
            echo regressing $ttype_subset $(date)
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
