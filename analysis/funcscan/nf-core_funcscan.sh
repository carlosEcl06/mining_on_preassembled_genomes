#!/bin/bash

nextflow run /path/to/funcscan/main.nf \
    -resume -profile <docker/apptainer> -c custom.config -params-file params.yaml