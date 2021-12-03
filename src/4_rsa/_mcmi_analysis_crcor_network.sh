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
