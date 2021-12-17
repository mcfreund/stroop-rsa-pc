## 1. initial run.

## cross-run correlation, no prewhitening

subjlists=("mc1" "mi1" "mc2" "mi2")
waves=("wave1" "wave1" "wave2" "wave2")
sessions=("baseline" "proactive" "baseline" "proactive")
ttype_subsets=("bias" "pc50")
glmname="lsall_1rpm"
roiset="Schaefer2018Parcel200"
measure="crcor"
prewh="none"
n_cores=18
overwrite="FALSE"
n_resamples=10000

for i in ${!subjlists[@]}
do

    session=${sessions[$i]}
    wave=${waves[$i]}
    subjlist=${subjlists[$i]}
    suffix=__seswave-$session_$wave

    echo $session $wave $subjlist $suffix
    
    Rscript ./src/4_rsa/parcellate_giftis.r \
        --glmname $glmname \
        --roiset $roiset \
        --subjlist $subjlist \
        --waves $wave \
        --sessions $sessions \
        --n_cores $n_cores \
        --overwrite $overwrite

    for ttype_subset in ${ttype_subsets[@]}
    do

        echo $ttype_subset

        Rscript ./src/4_rsa/estimate_distances.r \
            --glmname $glmname \
            --roiset $roiset \
            --subjlist $subjlist \
            --waves $wave \
            --sessions $sessions \
            --measure $measure \
            --prewh $prewh \
            --ttype_subset $ttype_subset \
            --n_cores $n_cores \
            --n_resamples $n_resamples
            
        Rscript ./src/4_rsa/regress_distances.r \
            --glmname $glmname \
            --roiset $roiset \
            --subjlist $subjlist \
            --waves $wave \
            --sessions $sessions \
            --measure $measure \
            --prewh $prewh \
            --ttype_subset $ttype_subset \
            --suffix $suffix
            
    done

done
