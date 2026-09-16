#!/usr/bin/env nextflow
// Smoke test — NOT part of the real pipeline. Validates two unknowns flagged
// in the plan before committing the real 554-way Snippy fan-out to Nextflow:
//   1. Does the PBS executor actually dispatch to a compute node (not pne2)?
//   2. Can a compute node pull + run an Apptainer container (registry egress
//      was only ever confirmed from the pne2 login node, not compute nodes)?

process helloPBS {
    executor 'pbs'
    queue 'workq'
    cpus 1
    memory '1 GB'
    time '5m'

    output:
    stdout

    script:
    """
    echo "[plain] host=\$(hostname)"
    echo "[plain] date=\$(date -u)"
    echo "[plain] no container, bare PBS executor"
    """
}

process helloApptainer {
    executor 'pbs'
    queue 'workq'
    cpus 1
    memory '1 GB'
    time '10m'
    container 'docker://ubuntu:22.04'  // alpine lacks /bin/bash, which Nextflow's container wrapper needs — confirmed the hard way

    output:
    stdout

    script:
    """
    echo "[container] host=\$(hostname)"
    echo "[container] os=\$(cat /etc/os-release | head -1)"
    echo "[container] date=\$(date -u)"
    """
}

workflow {
    helloPBS().view { "PLAIN RESULT:\n${it}" }
    helloApptainer().view { "APPTAINER RESULT:\n${it}" }
}
