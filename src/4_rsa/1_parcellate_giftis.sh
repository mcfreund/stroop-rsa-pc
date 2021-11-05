
waves=("wave1" "wave2")
sessions=("reactive")
echo ${waves[*]}
Rscript ./src/4_rsa/1_parcellate_giftis.r \
    --glmname "lsall_1rpm" \
    --roiset "Schaefer2018_control" \
    --subjlist "out/subjlist_ispc_retest.txt" \
    --waves "wave1 wave2" \
    --sessions "reactive"