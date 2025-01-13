#! /usr/bin/env bash

#   Require at least one argument, the list name
if [ "$#" -ne 2 ]; then
    echo "[INFO] Usage: $0 [user] [runlist] "
    exit 1
fi

#   Defines where the lists are stored
_LIST_DIRECTORY="./lists/"

#   Creates log for download file
touch "download_data.log"

#   Creates directory for data
mkdir -p ./Data/

LIST=$2
for RUN in $(cat ${_LIST_DIRECTORY}${LIST}); do
    echo "[INFO] Start download of run : ${RUN}"
    mkdir -p ./Data/${RUN}/
    scp -r $1@lxplus:/eos/project/s/sipm4all/data/2023-testbeam/physics-good/${RUN}/process-data-v1.0/*.root ./Data/${RUN}/
    if [ "$?" -ne "0" ]; then
        echo "[ERROR] Failed to download ${RUN}"
    fi
    echo "[INFO] Done downloading run : ${RUN}"
done

wait