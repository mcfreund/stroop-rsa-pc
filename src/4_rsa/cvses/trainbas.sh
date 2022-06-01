## looping variables
atlas_names=("glasser2016" "schaefer2018_17_200")
waves=("wave1" "wave2")
roi_cols=("superparcel" "parcel")

## constants
glmname="lsall_1rpm"
n_cores=24
decoder="svm"

## parcellate ----

## no preproc

for atlas_name in ${atlas_names[@]}
do
    echo ------------------------------------------ $atlas_name ------------------------------------------
    for roi_col in ${roi_cols[@]}
    do
        if [[ ($atlas_name == "schaefer2018_17_200" && $roi_col == "superparcel") ]]
        then 
            continue
        fi
        echo --------------------- $roi_col:
        for wave_i in ${!waves[@]}
        do
            subjlist=${waves[$wave_i]}
            echo --------- $subjlist:

            Rscript ./src/4_rsa/cvses/decode_trainbas.R $atlas_name $roi_col $subjlist \
                --glmname $glmname \
                --dowave $((wave_i+1)) \
                --n_cores $n_cores \
                --decoder $decoder \
                --demean --divnorm

        done
    done    
done




## with preproc


for atlas_name in ${atlas_names[@]}
do
    echo ------------------------------------------ $atlas_name ------------------------------------------
    for roi_col in ${roi_cols[@]}
    do
        if [[ ($atlas_name == "schaefer2018_17_200" && $roi_col == "superparcel") ]]
        then 
            continue
        fi
        echo --------------------- $roi_col:
        for wave_i in ${!waves[@]}
        do
            subjlist=${waves[$wave_i]}
            echo --------- $subjlist:

            Rscript ./src/4_rsa/cvses/decode_trainbas.R $atlas_name $roi_col $subjlist \
                --glmname $glmname \
                --dowave $((wave_i+1)) \
                --n_cores $n_cores \
                --decoder $decoder \
                --center --detrend --degree 1 --demean --divnorm

        done
    done    
done

