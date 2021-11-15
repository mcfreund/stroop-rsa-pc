## reactive test/retest analysis ----

## TEST 1. different first-level time-series models.

## on first run, had to remove DMCC5820265 and DMCC9478705 from subject list (manually)

subjlist="ispc_retest"
waves="wave1|wave2"  ## for loopin over within R
sessions="reactive"
glmnames=("lsall_1rpm" "lssep_1rpm")
roiset="Schaefer2018Dev"
prewh="none"
measures=("cveuc" "crcor")
ttype_subset="bias"

## 1. parcellate giftis

for glmname in ${glmnames[@]}
do
    Rscript ./src/4_rsa/parcellate_giftis.r \
        --glmname $glmname \
        --roiset $roiset \
        --subjlist $subjlist \
        --waves $waves \
        --sessions $sessions \
        --n_cores 26
done


## 2-3 estimate distances, regress distances

for glmname in ${glmnames[@]}
do
    for measure in ${measures[@]}
    do
        echo estimating distances
        Rscript ./src/4_rsa/estimate_distances.r \
            --glmname $glmname \
            --roiset $roiset \
            --subjlist $subjlist \
            --waves $waves \
            --sessions $sessions \
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
            --sessions $sessions \
            --measure $measure \
            --prewh $prewh \
            --ttype_subset $ttype_subset

    done
done






## TEST 2. 

## 2. prewhiten coefs

for prewh in ${prewhs[@]}
do
    Rscript ./src/4_rsa/prewhiten_coefs.r \
        --glmname $glmname \
        --roiset $roiset \
        --subjlist $subjlist \
        --waves $waves \
        --sessions $sessions \
        --prewh $prewh \
        --n_cores 26
done




## Misc: assess test-retest reliabilities of mean patterns ----

## removed DMCC5820265, DMCC9478705, DMCC3963378 (baseline wave2 run 2)


# Rscript ./src/4_rsa/parcellate_giftis.r \
#     --glmname "lsall_1rpm" \
#     --roiset "Schaefer2018Dev" \
#     --subjlist "out/subjlist_all_retest.txt" \
#     --waves "wave1 wave2" \
#     --sessions "reactive proactive baseline"


# Rscript ./src/4_rsa/pattern_reliability.r \
#     --glmname "lsall_1rpm" \
#     --roiset "Schaefer2018Dev" \
#     --subjlist "out/subjlist_all_retest.txt" \
#     --waves "wave1 wave2" \
#     --sessions "reactive proactive baseline" \
#     --outfile_prefix "trr-meanpatt__subjlist-ispc_retest"


# Rscript ./src/4_rsa/pattern_reliability.r \
#     --glmname "lsall_1rpm" \
#     --roiset "Schaefer2018Dev" \
#     --subjlist "out/subjlist_all_retest.txt" \
#     --prewh "none" \
#     --waves "wave1 wave2" \
#     --sessions "reactive proactive baseline" \
#     --outfile_prefix "trr-hilo__subjlist-ispc_retest"
