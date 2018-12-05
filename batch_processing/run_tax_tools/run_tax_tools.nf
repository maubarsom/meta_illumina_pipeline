#!/usr/bin/env nextflow

params.in_folder= '/media/virushd/preproc_libs/illumina'

base_folder = Channel.value( params.in_folder )

process find_fastq{
    input:
    val base_folder

    output:
    stdout m into find_output

    script:
    """
    find ${base_folder} -name "*.fq.gz" | sort
    """ }

fastq_files = find_output.map{it.trim()}.splitText().map{file(it.trim())}

// Create a separate channel for each tool to run
fastq_files.filter(~/.*_pe.fq.gz/).map{ [(it.name=~/(P[0-9]+_[0-9]+)_pe.fq.gz/)[0][1],it]}.into{metaphlan2_in;kraken_in;kaiju_in}

mapping_cpus=8
mpa2_pkl="\${mpa2_dir}/db_v20/mpa_v20_m200.pkl"
mpa2_bowtie2db="\${mpa2_dir}/db_v20/mpa_v20_m200"

process metaphlan2{
    cpus mapping_cpus

    errorStrategy 'ignore'

    module 'bowtie2/2.3.0'
    module 'metaphlan2/2.6.0'

    tag {sample_id}

    publishDir 'metaphlan2'

    input:
    set val(sample_id),file(fastq_pe) from metaphlan2_in

    output:
    set val(sample_id),"${sample_id}_mpa2.txt" into mpa2_out

    script:
    """
    metaphlan2.py --mpa_pkl ${mpa2_pkl} --bowtie2db ${mpa2_bowtie2db} --bowtie2out ${sample_id}.bowtie2.bz2 --nproc ${mapping_cpus} --input_type multifastq --sample_id_key ${sample_id} ${fastq_pe} ${sample_id}_mpa2.txt
    """
}

export2graphlan_bin='${mpa2_dir}/utils/export2graphlan/export2graphlan.py'
graphlan_dir='/labcommon/tools/graphlan'

process graphlan{
	cpus 1

	errorStrategy 'ignore'

	tag {sample_id}
	
	module 'metaphlan2/2.6.0'

	publishDir 'metaphlan2/graphlan'

	input:
	set val(sample_id),file(mpa2_file) from mpa2_out

	output:
	set val(sample_id),"${sample_id}_mpa2*"

	script:
    """
	${export2graphlan_bin} -i ${mpa2_file} --tree ${sample_id}_mpa2.tree --annotation ${sample_id}_mpa2.annot --abundance_threshold 0.5
	${graphlan_dir}/graphlan_annotate.py --annot ${sample_id}_mpa2.annot ${sample_id}_mpa2.tree ${sample_id}_mpa2_graphlan.xml
	${graphlan_dir}/graphlan.py --dpi 300 ${sample_id}_mpa2_graphlan.xml ${sample_id}_mpa2.png 
	"""
}

