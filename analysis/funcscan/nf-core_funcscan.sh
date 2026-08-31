#!/bin/bash

NXF_VER=25.10.4 nextflow run ./pipeline/3_0_0/main.nf \
    -resume -profile apptainer -c custom.config -params-file params.yaml