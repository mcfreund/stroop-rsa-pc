## mimc (proactive + baseline wave1) analysis for alt glms ----


## 1. initial runthru:

## - no prewhitening
## - get timing

## on first run, had to remove  from subject list (manually):
## - DMCC9478705

subjlist="mcmi"
wave="wave1"
sessions="proactive|baseline"
glmname="lsall_1rpm"
roiset="Schaefer2018Network"
measure="crcor"
ttype_subsets=("bias" "pc50")

prewh="none"

Rscript ./src/4_rsa/parcellate_giftis.r \
    --glmname $glmname \
    --roiset $roiset \
    --subjlist $subjlist \
    --waves $wave \
    --sessions $sessions \
    --n_cores 26

## started at 2:15p, 2 dec
for ttype_subset in ${ttype_subsets[@]}
do
    Rscript ./src/4_rsa/estimate_distances.r \
        --glmname $glmname \
        --roiset $roiset \
        --subjlist $subjlist \
        --waves $wave \
        --sessions $sessions \
        --measure $measure \
        --prewh $prewh \
        --ttype_subset $ttype_subset \
        --n_cores 30 \
        --n_resamples 10000
        
    Rscript ./src/4_rsa/regress_distances.r \
        --glmname $glmname \
        --roiset $roiset \
        --subjlist $subjlist \
        --waves $wave \
        --sessions $sessions \
        --measure $measure \
        --prewh $prewh \
        --ttype_subset $ttype_subset
done



