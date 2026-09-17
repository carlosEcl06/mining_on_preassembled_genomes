#!/usr/bin/env nextflow
// Real pipeline, stage 1: Snippy (--ctgs, per genome) -> snippy-core.
//
// nf-core/modules' snippy_run does NOT support contig input (only
// --R1/--R2/--se reads — confirmed by reading its source directly,
// correcting the plan's original assumption) so SNIPPY_CTGS below is
// hand-written, mirroring that module's container/output conventions so it
// stays a drop-in match for SNIPPY_CORE (which IS generic: per-sample VCF +
// aligned FASTA + reference, indifferent to how each was produced).
// Container pin (4.6.0--hdfd78af_2) taken directly from the nf-core module
// source itself — a live, already-vetted pin, not guessed.

// Plain FASTA, not the annotated .gbk: Snippy auto-builds a SnpEff database
// from GenBank references, which failed here (CDS/protein check files not
// produced — a real bug, not a config issue). We don't need SnpEff's
// variant-effect annotation for this pipeline anyway (only snps.vcf /
// snps.aligned.fa feed snippy-core) — a plain FASTA has no gene models, so
// Snippy skips the SnpEff build step entirely.
// Reference switched from CP079087 (HpGP resequencing) to NC_000915 (the
// classic Tomb et al. 1997 RefSeq assembly): GrafGen's ancestry marker panel
// is coordinate-anchored to NC_000915, confirmed directly against its own
// shipped example VCF (CHROM=NC_000915) and reference dataframe -- running
// GrafGen against a CP079087-called VCF produced a degenerate, biologically
// implausible result (all 552 genomes collapsed to one Refpop, E_percent=0
// uniformly) because the two accessions are different sequencing efforts of
// strain 26695 with non-corresponding coordinates, not the same numbering.
// One reference is used for both the tree and the ancestry panel, so this
// swap re-runs the whole fan-out rather than maintaining two alignments.
params.reference   = "${projectDir}/../rawdata/reference/NC_000915.fasta"
params.seqkit_dir  = "/data2/projects/LGMB-009/NGS/analysis/funcscan/results/bgc/seqkit"
params.results_dir = "${projectDir}/../results/snippy"
params.exclude     = ['GCA_024409305.1', 'GCA_015904695.1']  // <1.4 Mb, see plan
params.snippy_sif  = 'docker://quay.io/biocontainers/snippy:4.6.0--hdfd78af_2'
// pne5 is missing /usr/bin/apptainer entirely (confirmed directly: file does
// not exist there, on all other 5 free nodes it does) — a real gap in that
// node's setup, not fixable without root. PBS here doesn't support host!=
// exclusion syntax ("Illegal attribute or resource value"), so instead every
// process explicitly round-robins across the 5 known-good nodes via
// clusterOptions (validated on a standalone 10-task test: all 10 landed on
// pne3/4/6/7/10, never pne5, before trusting this on the real 552-genome run).
params.good_nodes  = ['pne3', 'pne4', 'pne6', 'pne7', 'pne10']
// For smoke-testing on a handful of genomes, override --seqkit_dir to point
// at a directory with just a few symlinked *_long.fasta files, rather than
// slicing the list here — every attempt at in-workflow list slicing
// (.take(), .subList(), with or without `def`) failed with a bizarre
// "Missing process or function X(...)" DSL error not worth chasing further.

process SNIPPY_CTGS {
    tag "${sample_id}"
    executor 'pbs'
    clusterOptions { "-q workq -l select=1:ncpus=2:mem=4gb:host=${params.good_nodes[task.index % params.good_nodes.size()]}" }
    time '30m'
    container params.snippy_sif
    errorStrategy 'retry'  // dynamic closure form crashed while handling an input-staging-time
    maxRetries 2           // failure (task context not yet initialized) — static form is fine
    publishDir "${params.results_dir}", mode: 'copy'  // filenames already embed sample_id, no need for a subdir

    input:
    tuple val(sample_id), path(contigs)
    path reference

    output:
    tuple val(sample_id), path("${sample_id}.vcf"), path("${sample_id}.aligned.fa"), emit: core_input
    path "${sample_id}.log", emit: log
    path "${sample_id}.txt", emit: stats

    script:
    """
    snippy --ctgs ${contigs} --ref ${reference} --outdir out --cpus 2 --force
    cp out/snps.vcf ${sample_id}.vcf
    cp out/snps.aligned.fa ${sample_id}.aligned.fa
    cp out/snps.log ${sample_id}.log
    cp out/snps.txt ${sample_id}.txt
    """
}

process SNIPPY_CORE {
    executor 'pbs'
    clusterOptions "-q workq -l select=1:ncpus=4:mem=8gb:host=${params.good_nodes[0]}"  // single job, pin to the first good node explicitly
    time '2h'
    container params.snippy_sif
    publishDir "${params.results_dir}/core", mode: 'copy'

    input:
    val sample_data  // list of [sample_id, vcf, fasta] tuples
    path reference

    output:
    path "core.aln", emit: alignment
    path "core.vcf", emit: vcf
    path "core.txt", emit: stats
    path "core.tab", emit: tab

    script:
    def links = sample_data.collect { sid, vcf, fasta ->
        "mkdir -p samples/${sid} && ln -s ${vcf.toRealPath()} samples/${sid}/snps.vcf && ln -s ${fasta.toRealPath()} samples/${sid}/snps.aligned.fa"
    }.join('\n    ')
    """
    ${links}
    snippy-core --ref ${reference} \$(ls -d samples/*/ | sed 's:/\$::') --prefix core
    """
}

workflow {
    reference_ch = Channel.fromPath(params.reference, checkIfExists: true)

    // Plain Groovy (not channel operators) for the glob/filter/limit step —
    // sidesteps Nextflow's internal collection-type quirks (ArrayBag etc.)
    // entirely; only the final list touches a Channel.
    def all_files = new File(params.seqkit_dir).listFiles({ f -> f.name.endsWith('_long.fasta') } as FileFilter) as List
    def tuples = all_files
        .collect { f -> [f.name.replace('_long.fasta', ''), f.toPath()] }  // toPath(): java.io.File isn't a valid Nextflow `path` input type on its own
        .findAll { sid, f -> !(sid in params.exclude) }
        .sort { it[0] }

    genomes_ch = Channel.fromList(tuples)

    snippy_out = SNIPPY_CTGS(genomes_ch, reference_ch.first())

    core_out = SNIPPY_CORE(
        snippy_out.core_input.toList(),  // .collect() flattened the tuples into one long scalar list; .toList() preserves each [sid,vcf,fasta] record
        reference_ch.first()
    )

    core_out.stats.view { "snippy-core done: ${it}" }
}
