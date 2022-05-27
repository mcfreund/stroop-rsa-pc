subjlist=("wave1_unrel_pilot")
variables=("target_denoise" "distractor_denoise" "congruency_denoise")
glmtypes=("hrf" "hrf_rr")


session=("baseline")
ttype_subsets=("bias" "pc50")
glmname="glmsingle_wave1"
atlas_name="glasser2016"
roi_col="parcel"
measure="crcor"
wave="wave1"
prewh="none"
n_cores=28
overwrite="TRUE"
n_resamples=10000
space="fsaverage5"


for ttype_subset in ${ttype_subsets[@]}
do
    for variable in ${variables[@]}
    do
        for glmtype in ${glmtypes[@]}
        do

            glmname_i=$glmname"_"$variable"_"$glmtype

            echo estimating $ttype_subset $(date)
            Rscript ./src/4_rsa/estimate_distances.r \
                --glmname $glmname_i \
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
                --glmname $glmname_i \
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
