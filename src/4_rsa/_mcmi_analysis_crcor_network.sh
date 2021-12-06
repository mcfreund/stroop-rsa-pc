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
measure="crcor"
ttype_subsets=("bias" "pc50")

prewh="none"

Rscript ./src/4_rsa/parcellate_giftis.r \
    --glmname $glmname \
    --roiset $roiset \
    --subjlist $subjlist \
    --waves $wave \
    --sessions $sessions \
    --n_cores 26 \
    --overwrite "FALSE"

## NB: runtime one iteration (1E4 resamples) ~ 1.5 hr
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




## 2. prewhiten coefs

unset prewh
prewhs=("obsresampbias" "obsresamppc50")

for prewh in ${prewhs[@]}
do
    Rscript ./src/4_rsa/prewhiten_coefs.r \
        --glmname $glmname \
        --roiset $roiset \
        --subjlist $subjlist \
        --waves $wave \
        --sessions $sessions \
        --prewh $prewh \
        --n_cores 25 \
        --overwrite "FALSE" \
        --n_resamples 100
done


## 3. estimate distances among prewhitened coefs

for ttype_subset in ${ttype_subsets[@]}
do

    ## use patterns with prewhitening method that matches ttype_subset:

    if [[ttype_subset = "bias"]]
    then
        prewh="obsresampbias"
    elif [[ttype_subset = "pc50"]]
        prewh="obsresamppc50"
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