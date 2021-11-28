## reactive test/retest analysis ----

## TEST 1. different first-level time-series models.

## on first run, had to remove DMCC5820265 and DMCC9478705 from subject list (manually)

subjlist="all_retest"
waves="wave1|wave2"  ## for loopin over within R
sessions="proactive"
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


## 1.1 Add in condition_1rpm models

glmname="condition_1rpm"
measure="cveuc"
echo parcellating coefficients
Rscript ./src/4_rsa/parcellate_giftis.r \
    --glmname $glmname \
    --roiset $roiset \
    --subjlist $subjlist \
    --waves $waves \
    --sessions $sessions \
    --n_cores 4

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
    --n_cores 10

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



## 2. prewhiten coefs

#prewh="obsall"
prewh="obsresampbias"
glmname="lsall_1rpm"
expected_min=3
overwrite="FALSE"
Rscript ./src/4_rsa/prewhiten_coefs.r \
    --glmname $glmname \
    --roiset $roiset \
    --subjlist $subjlist \
    --waves $waves \
    --sessions $sessions \
    --prewh $prewh \
    --n_cores 16 \
    --expected_min $expected_min \
    --overwrite $overwrite

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
