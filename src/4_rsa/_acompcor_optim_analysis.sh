## proactive/reactive wave1 analysis for alt glms ----

## TEST 1. different first-level time-series models.

## on first run, had to remove DMCC5820265 and DMCC9478705 from subject list (manually)

subjlists=("all_retest" "ispc_retest")
waves="wave1"
sessions=("proactive" "reactive")
glmname="lsall_acompcor06_1rpm"
roiset="Schaefer2018Dev"
prewh="none"
measures=("cveuc" "crcor")
ttype_subset="bias"

## 1. parcellate giftis

for session_i in ${!subjlists[@]}
do
    subjlist=${subjlists[$session_i]}
    session=${sessions[$session_i]}

    Rscript ./src/4_rsa/parcellate_giftis.r \
        --glmname $glmname \
        --roiset $roiset \
        --subjlist $subjlist \
        --waves $waves \
        --sessions $session \
        --n_cores 26
done


## 2-3 estimate distances, regress distances

for session_i in ${!subjlists[@]}
do
    subjlist=${subjlists[$session_i]}
    session=${sessions[$session_i]}

    for measure in ${measures[@]}
    do
        echo estimating distances
        Rscript ./src/4_rsa/estimate_distances.r \
            --glmname $glmname \
            --roiset $roiset \
            --subjlist $subjlist \
            --waves $waves \
            --sessions $session \
            --measure $measure \
            --prewh $prewh \
            --ttype_subset $ttype_subset \
            --n_cores 26
        
        echo regressing distances
        Rscript ./src/4_rsa/regress_distances.r \
            --glmname $glmname \
            --roiset $roiset \
            --subjlist $subjlist \
            --waves $waves \
            --sessions $session \
            --measure $measure \
            --prewh $prewh \
            --ttype_subset $ttype_subset
        done

    done
done