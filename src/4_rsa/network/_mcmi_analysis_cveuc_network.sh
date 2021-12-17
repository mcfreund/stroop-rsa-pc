## MCMI analysis: baseline wave1, proactive wave 1 ----

## had to remove  from subject list (manually):
## - DMCC9478705

subjlist="mcmi"
wave="wave1"
sessions="proactive|baseline"
glmname="lsall_1rpm"
roiset="Schaefer2018Network"
measure="cveuc"
ttype_subsets=("bias" "pc50")
prewhs=("none" "obsbias" "obspc50")

for ttype_subset in ${ttype_subsets[@]}
do
    echo $ttype_subset ----
    for prewh in ${prewhs[@]}
    do
        echo $prewh ----
        Rscript ./src/4_rsa/estimate_distances.r \
            --glmname $glmname \
            --roiset $roiset \
            --subjlist $subjlist \
            --waves $wave \
            --sessions $sessions \
            --measure $measure \
            --prewh $prewh \
            --ttype_subset $ttype_subset \
            --n_cores 12 \
            --overwrite "FALSE"
            
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
done