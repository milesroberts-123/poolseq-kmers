#!/bin/bash

envname=$1

pwd

if [ -d /.singularity.d ]; then
    echo Singularity detected! Using conda env in container
    conda init
    source ~/.bashrc
    echo $envname
    conda activate $(echo $envname)
fi
