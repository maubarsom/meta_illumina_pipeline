#!/usr/bin/env nextflow

params.in_folder= '/media/virushd/preproc_libs/illumina'
params.ref_fasta = 'ref.fa'
params.ref_name = 'ref'

ref_name = params.ref_name

base_folder = Channel.value( params.in_folder )
references = Channel.fromPath(params.ref_fasta)

process create_index{
	module 'bwa/0.7.15'

    input:
    file references

    output:
    file 'ref*' into bwa_indexes

    script:
    """
    bwa index -p ref $references
    """
}

process find_fastq{  
    input:
    val base_folder

    output:
    stdout m into fastq_files

    script:
    """
    find ${base_folder} -name "*.fq.gz" | sort 
    """
}

// Split into two channels
fastq_pe = Channel.create()
fastq_se = Channel.create()

fastq_files.map{it.trim()}.splitText().map{file(it.trim())}.choice(fastq_pe,fastq_se) { f -> f.name.contains('pe.fq.gz') ? 0 : 1 }

mapping_sets = fastq_pe.combine( bwa_indexes.map({[it,]}) )

mapping_cpus=8

process map_pe_reads_to_reference{
	echo true
	cpus mapping_cpus

	module 'samtools/1.3'
	module 'bwa/0.7.15'

    input:
    set file(my_pe),file(bwa_idx) from mapping_sets

    output:
    val sample_id into mapped_samples
    file "${sample_id}_${ref_name}.bam" into bam_files
    file "${sample_id}_${ref_name}.bam.bai" into bai_files

    script:
    sample_id = (my_pe.name=~/(P[0-9]+_[0-9]+)_pe.fq.gz/)[0][1]
    """
	echo "Processing ${sample_id}"
    bwa mem -t ${mapping_cpus} -T 30 -M -p ref ${my_pe} | samtools view -hSb -F256 -f2 - | samtools sort - > ${sample_id}_${ref_name}.bam
    samtools index ${sample_id}_${ref_name}.bam
    """
}

/**
Join output of mapping with sample id
and 
duplicate channel for the two different processes
**/
mapped_samples.merge(bam_files,bai_files){a,b,c->[a,b,c]}.into{bam_to_flagstat;bam_to_depth}

process bam_flagstats{
	module 'samtools/1.3'
	cpus 1

    input:
    set val(sample_id),'mapped.bam','mapped.bam.bai' from bam_to_flagstat

    output:
    val sample_id into flgstat_samples
    file "${sample_id}_${ref_name}.bam.flagstat" into flgstat_files

    script:
    """
    samtools flagstat mapped.bam > ${sample_id}_${ref_name}.bam.flagstat
    """
}

process bam_depth{
	module 'samtools/1.3'
	cpus 1

    input:
    set val(sample_id),'mapped.bam','mapped.bam.bai' from bam_to_depth

    output:
    val sample_id into sam_depth_samples
    file "${sample_id}_${ref_name}.bam.depth.gz" into depth_files

    script:
    """
    samtools depth mapped.bam | gzip > ${sample_id}_${ref_name}.bam.depth.gz
    """
}
