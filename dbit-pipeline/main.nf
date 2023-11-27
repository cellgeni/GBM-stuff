#!/usr/bin/env nextflow

def errorMessage() {
    log.info"""
    ===========
    Input error
    ===========
    You failed to provide the --samplefile parameter
    Please provide this parameter as follows:
      --samplefile /full/path/to/sample/file
    The pipeline has exited with error status 1.
    """.stripIndent()
    exit 1
}

process filter_primer {

  label 'filter'

  publishDir "${params.outdir}/QC", mode: 'copy', pattern: "*._stats.primer.txt"

  input:
  tuple val(sampleid), path(fastq_R1), path(fastq_R2)

  output:
  tuple val(sampleid), path("*_raw_qc_primer_R1.fastq.gz"), path("*_raw_qc_primer_R2.fastq.gz"), emit: fq_tuple
  path("*_stats.primer.txt")
  
  shell:
  '''
  !{projectDir}/bin/bbmap/bbduk.sh \
    in1=!{fastq_R1} \
    in2=!{fastq_R2} \
    literal=CAAGCGTTGGCTTCTCGCATCT \
    outm1="!{sampleid}_raw_qc_primer_R1.fastq.gz" \
    outm2="!{sampleid}_raw_qc_primer_R2.fastq.gz" \
    stats="!{sampleid}_stats.primer.txt" \
    k=22 mm=f rcomp=f restrictleft=30 skipr1=t \
    hdist=2 \
    threads=!{task.cpus}
  '''
}

process filter_L1 {

  label 'filter'
 
  publishDir "${params.outdir}/QC", mode: 'copy', pattern: "*._stats.linker1.txt"

  input:
  tuple val(sampleid), path(fastq_R1), path(fastq_R2)

  output:
  tuple val(sampleid), path("*_raw_qc_linker1_R1.fastq.gz"), path("*_raw_qc_linker1_R2.fastq.gz"), emit: fq_tuple
  path("*_stats.linker1.txt")
  
  shell:
  '''
  !{projectDir}/bin/bbmap/bbduk.sh \
    in1=!{fastq_R1} \
    in2=!{fastq_R2} \
    literal=GTGGCCGATGTTTCGCATCGGCGTACGACT \
    outm1="!{sampleid}_raw_qc_linker1_R1.fastq.gz" \
    outm2="!{sampleid}_raw_qc_linker1_R2.fastq.gz" \
    stats="!{sampleid}_stats.linker1.txt" \
    k=30 mm=f rcomp=f restrictleft=108 skipr1=t \
    hdist=3 \
    threads=!{task.cpus}
  '''
}

process filter_L2 {

  label 'filter'

  publishDir "${params.outdir}/QC", mode: 'copy', pattern: "*._stats.linker2.txt"

  input:
  tuple val(sampleid), path(fastq_R1), path(fastq_R2)

  output:
  tuple val(sampleid), path("*_raw_qc_R1.fastq.gz"), emit: r1_tuple
  tuple val(sampleid), path("*_raw_qc_R2.fastq.gz"), emit: r2_tuple
  path("*_stats.linker2.txt")

  shell:
  '''
  !{projectDir}/bin/bbmap/bbduk.sh \
    in1=!{fastq_R1} \
    in2=!{fastq_R2} \
    literal=ATCCACGTGCTTGAGAGGCCAGAGCATTCG \
    outm1="!{sampleid}_raw_qc_R1.fastq.gz" \
    outm2="!{sampleid}_raw_qc_R2.fastq.gz" \
    stats="!{sampleid}_stats.linker2.txt" \
    k=30 mm=f rcomp=f restrictleft=70 skipr1=t \
    hdist=3 \
    threads=!{task.cpus}
  '''
}

process additional_filter {
  //inspecting fastqs with pm19 looked like the first 31 bases were technical and not biological so needed to remove this
  label 'filter'

  input:
  tuple val(sampleid), path(fastq_R1)

  output:
  tuple val(sampleid), path("*_S1_L001_R3_001.fastq.gz"), emit: r3_tuple

  shell:
  '''
  !{projectDir}/bin/bbmap/bbduk.sh \
    in=!{fastq_R1} \
    forcetrimleft=32 \
    out="!{sampleid}_S1_L001_R3_001.fastq.gz" \
    threads=!{task.cpus}
  '''
}

process bc_process {

  input:
  tuple val(sampleid), path(fastq_R2)

  output:
  tuple val(sampleid), path("*_S1_L001_R1_001.fastq.gz"), path("*_S1_L001_R2_001.fastq.gz"), emit: fq_tuple
  
  shell:
  '''
  python !{projectDir}/bin/BC_process.py --input !{fastq_R2} --output_R1 "!{sampleid}_S1_L001_R1_001.fastq" --output_R2 "!{sampleid}_S1_L001_R2_001.fastq"
  gzip "!{sampleid}_S1_L001_R1_001.fastq"
  gzip "!{sampleid}_S1_L001_R2_001.fastq"
  '''
}


process cell_ranger {

  publishDir "${params.outdir}/${sampleid}", mode: 'copy', pattern: '*'

  input:
  tuple val(sampleid), path(fastq_R3, stageAs: 'fastqs/*'), path(fastq_R1, stageAs: 'fastqs/*'), path(fastq_R2, stageAs: 'fastqs/*')
  //tuple val(sampleid), path(fastq_R3, stageAs: {filename -> filename.toString().substring(0, filename.toString().indexOf('raw')) + "S1_L001_R3_001.fastq.gz"}) , path(fastq_R1), path(fastq_R2)

  //--reference="/nfs/srpipe_references/downloaded_from_10X/refdata-cellranger-arc-GRCh38-2020-A" \
  shell:
  '''
  /lustre/scratch127/cellgen/cellgeni/tickets/tic-2598/cellranger-atac-2.0.0/bin/cellranger-atac count \
    --id=!{sampleid} \
    --reference="/nfs/srpipe_references/downloaded_from_10X/refdata-cellranger-arc-GRCh38-2020-A-2.0.0" \
    --fastqs="fastqs" \
    --sample=!{sampleid} \
    --localcores=!{task.cpus} \
    --localmem="128000"
  '''
}

workflow {
  
  main:
    ch_sample_list = params.samplefile != null ? Channel.fromPath(params.samplefile) : errorMessage()
    ch_sample_list | flatMap{ it.readLines() } | map { it -> [ it.split()[0], it.split()[1], it.split()[2] ] } | filter_primer
    filter_L1(filter_primer.out.fq_tuple)
    filter_L2(filter_L1.out.fq_tuple)
    bc_process(filter_L2.out.r2_tuple)
    //filter_L2.out.r2_tuple | map { val, filename -> file_prefix=filename.toString().substring(0, filename.toString().indexOf('raw')); println "${file_prefix}S1_L001_R3_001.fastq.gz" }
    //filter_L2.out.r2_tuple | map { val, filename -> println filename.toString().substring(0, filename.toString().indexOf('raw')) + "S1_L001_R3_001.fastq.gz" }
    additional_filter(filter_L2.out.r1_tuple)
    additional_filter.out.r3_tuple.join(bc_process.out.fq_tuple) | cell_ranger 
}
