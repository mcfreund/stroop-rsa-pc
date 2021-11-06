

## on first run, had to remove DMCC5820265 and DMCC9478705 from subject list (manually)

waves=("wave1" "wave2")
sessions=("reactive")
Rscript ./src/4_rsa/1_parcellate_giftis.r \
    --glmname "lsall_1rpm" \
    --roiset "Schaefer2018_control" \
    --subjlist "out/subjlist_ispc_retest.txt" \
    --waves "wave1 wave2" \
    --sessions "reactive"

