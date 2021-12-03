## jneurosci reanalysis (Glasser2016 proactive bias wave1) analysis ----


## 1. initial runthru:

## on first run, had to remove  from subject list (manually):


subjlist="jneurosci"
wave="wave1"
session="proactive"
glmnames=("lsall_1rpm" "condition_1rpm")
roiset="Glasser2016Parcel"
measures=("crcor" "cveuc")
ttype_subset="bias"
prewh="none"


for glmname in ${glmnames[@]}
do
    Rscript ./src/4_rsa/parcellate_giftis.r \
        --glmname $glmname \
        --roiset $roiset \
        --subjlist $subjlist \
        --waves $wave \
        --sessions $session \
        --n_cores 26
done

for measure in ${measures[@]}
do

    if [[measure = "crcor"]]
    then
        glmname="lsall_1rpm"
    elif [[measure = "cveuc"]]
    then
        glmname="condition_1rpm"
    fi

    Rscript ./src/4_rsa/estimate_distances.r \
        --glmname $glmname \
        --roiset $roiset \
        --subjlist $subjlist \
        --waves $wave \
        --sessions $session \
        --measure $measure \
        --prewh $prewh \
        --ttype_subset $ttype_subset \
        --n_cores 24 \
        --n_resamples 10000
        
    Rscript ./src/4_rsa/regress_distances.r \
        --glmname $glmname \
        --roiset $roiset \
        --subjlist $subjlist \
        --waves $wave \
        --sessions $session \
        --measure $measure \
        --prewh $prewh \
        --ttype_subset $ttype_subset

done



