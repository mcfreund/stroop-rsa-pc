#!/usr/bin/env bash

## string funs

function get_wave_bluearc { ## get wave dir on bluearc

    if [ $1 = 1 ]
    then
        echo HCP_SUBJECTS_BACKUPS
    elif [ $1 = 2 ]
    then
        echo DMCC_Phase3
    elif  [ $1 = 3 ]
    then
        echo DMCC_Phase4
    fi

}

function get_sess { ## get sess short

    if [ $1 = baseline ]
    then
        echo Bas
    elif [ $1 = proactive ]
    then
        echo Pro
    elif  [ $1 = reactive ]
    then
        echo Rea
    fi

}


function get_run_suffix { ## get run suffix
	
    if [ $1 = 1 ]
    then
        echo 1_AP
    elif [ $1 = 2 ]
    then
        echo 2_PA
    fi

}


# strings

hemis=(L R)
runs=(1 2)
xmat_dir=${proj_dir}/out/glms
bluearc_dir=/data/nil-bluearc/ccp-hcp/DMCC_ALL_BACKUPS
