#!/bin/bash

## Command to run nf-core/funcscan on pre-assembled genomes with docker
## ADJUST FOR APPTAINER

nextflow run /home/carlos_ecl06/projects/whole_metagenome_plus_barcoding/analysis/wgs-shotgun/funcscan/3_0_0/main.nf \
    -resume -profile docker -c custom.config -params-file params.yaml