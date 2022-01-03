atlas_names=("glasser2016" "schaefer2018_17_200")
ttype_subsets=("bias" "pc50")
measures=("cveuc" "crcor")

subjlist=("mi1")
wave=("wave1")
session=("proactive")
glmname="lsall_1rpm"
roi_col="parcel"
prewh="none"
n_cores=18
overwrite="TRUE"
n_resamples=10000
space="fsaverage5"

## excluded due to missing (2021-12-17):
# mi1 DMCC9478705


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

echo $(date)
for measure in ${measures[@]}
do

    echo $measure

    for atlas_name in ${atlas_names[@]}
    do

        echo $atlas_name

        for ttype_subset in ${ttype_subsets[@]}
        do

            echo $ttype_subset

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
echo $(date)