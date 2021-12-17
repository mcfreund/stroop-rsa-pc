## MIMC analysis: proactive wave 1, baseline wave 2 ----


## 1. initial runthru:

subjlist="mimc"
waves=("wave1" "wave2")
sessions=("proactive" "baseline")
glmname="lsall_1rpm"
roiset="Schaefer2018Network"
measure="cveuc"
ttype_subsets=("bias" "pc50")
prewhs=("none" "obsbias" "obspc50")

for ttype_subset in ${ttype_subsets[@]}
do
    echo starting $ttype_subset
    for seswav_i in ${!sessions[@]}
    do
        for prewh in ${prewhs[@]}
        do
            Rscript ./src/4_rsa/estimate_distances.r \
                --glmname $glmname \
                --roiset $roiset \
                --subjlist $subjlist \
                --waves ${waves[$seswav_i]} \
                --sessions ${sessions[$seswav_i]} \
                --measure $measure \
                --prewh $prewh \
                --ttype_subset $ttype_subset \
                --n_cores 20 \
                --overwrite "TRUE"
                
            Rscript ./src/4_rsa/regress_distances.r \
                --glmname $glmname \
                --roiset $roiset \
                --subjlist $subjlist \
                --waves ${waves[$seswav_i]} \
                --sessions ${sessions[$seswav_i]} \
                --measure $measure \
                --prewh $prewh \
                --ttype_subset $ttype_subset \
                --suffix __seswave-${sessions[$seswav_i]}_${waves[$seswav_i]}
        done
    done
done