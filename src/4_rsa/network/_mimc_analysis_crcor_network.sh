## MIMC analysis: proactive wave 1, baseline wave 2 ----


## 1. initial runthru:
## manually remove: DMCC9478705, DMCC3963378 (wave2)

subjlist="mimc"
waves=("wave1" "wave2")
sessions=("proactive" "baseline")
glmname="lsall_1rpm"
roiset="Schaefer2018Network"
measure="crcor"
ttype_subsets=("bias" "pc50")

prewh="none"

for seswav_i in ${!sessions[@]}
do
    Rscript ./src/4_rsa/parcellate_giftis.r \
        --glmname $glmname \
        --roiset $roiset \
        --subjlist $subjlist \
        --waves ${waves[$seswav_i]} \
        --sessions ${sessions[$seswav_i]} \
        --n_cores 30 \
        --overwrite "FALSE"
done

for ttype_subset in ${ttype_subsets[@]}
do
    for seswav_i in ${!sessions[@]}
    do
        echo estimate: ${sessions[$seswav_i]} ${waves[$seswav_i]} $(date)
        Rscript ./src/4_rsa/estimate_distances.r \
            --glmname $glmname \
            --roiset $roiset \
            --subjlist $subjlist \
            --waves ${waves[$seswav_i]} \
            --sessions ${sessions[$seswav_i]} \
            --measure $measure \
            --prewh $prewh \
            --ttype_subset $ttype_subset \
            --n_cores 16 \
            --n_resamples 10000 \
            --overwrite "FALSE"

        echo regress: ${sessions[$seswav_i]} ${waves[$seswav_i]} $(date)    
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


## 2. prewhiten coefs
unset prewh
prewhs=("obspc50" "obsbias")
for prewh in ${prewhs[@]}
do
    echo start: $prewh "baseline wave2" $(date)
    Rscript ./src/4_rsa/prewhiten_coefs.r \
        --glmname $glmname \
        --roiset $roiset \
        --subjlist $subjlist \
        --waves "wave2" \
        --sessions "baseline" \
        --prewh $prewh \
        --n_cores 12 \
        --overwrite "FALSE"
    echo stop: $prewh "baseline wave2" $(date)
done

for prewh in ${prewhs[@]}
do
    echo start: $prewh "proactive wave1" $(date)
    Rscript ./src/4_rsa/prewhiten_coefs.r \
        --glmname $glmname \
        --roiset $roiset \
        --subjlist $subjlist \
        --waves "wave1" \
        --sessions "proactive" \
        --prewh $prewh \
        --n_cores 12 \
        --overwrite "FALSE"
    echo stop: $prewh "proactive wave1" $(date)
done


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
                --n_cores 12 \
                --overwrite "FALSE" \
                --n_resamples 10000
                
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
