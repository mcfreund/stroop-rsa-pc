## MCMI analysis: baseline wave1, proactive wave 1 ----


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
measure="cveuc"
ttype_subsets=("bias" "pc50")

for ttype_subset in ${ttype_subsets[@]}
do

    ## use patterns with prewhitening method that matches ttype_subset:
    if [[ $ttype_subset = "bias" ]]
    then
        prewh="obsbias"
    else
        prewh="obspc50"
    fi
    
    Rscript ./src/4_rsa/estimate_distances.r \
        --glmname $glmname \
        --roiset $roiset \
        --subjlist $subjlist \
        --waves $wave \
        --sessions $sessions \
        --measure $measure \
        --prewh $prewh \
        --ttype_subset $ttype_subset \
        --n_cores 10 \
        --overwrite "TRUE"
        
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