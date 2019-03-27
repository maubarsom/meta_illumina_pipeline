#!/usr/bin/env nextflow

/*
How to run:
nextflow -C preprocessing.nf.config run preprocessing.nf --fastq_files preprocessing -profile hamlet
or
nextflow -C preprocessing.nf.config run preprocessing.nf --fastq_files preprocessing -profile bianca
*/

params.fastq_dir='reads/'
// Assume all fastq files are in the read folder
fastq_files = Channel.fromFilePairs("${params.fastq_dir}/*.fastq.gz"){ (it.name =~ /P[0-9]{3,5}_[0-9]{3,5}/)[0]}
// Uncomment this is fastq files are nested deep
//fastq_files = Channel.fromFilePairs("${params.fastq_dir}/**/*.fastq.gz"){ (it.name =~ /P[0-9]{3,5}_[0-9]{3,5}/)[0]}

fastq_files.into{qf_raw_fastqc_in;
                 qf_trimgalore_in}


process qf_raw_fastqc{
  tag {"${sample_id}"}
  publishDir "preprocessing/${sample_id}/1_raw", mode: 'link'

  input:
  set sample_id,file(reads) from qf_raw_fastqc_in

  output:
  file '*.html' into raw_fastqc_out

  script:
  """
  fastqc --noextract -k 10 -t ${task.cpus} *.fastq.gz
  """
}

process qf_trimgalore{
  tag {"${sample_id}"}
  publishDir "preprocessing/${sample_id}", mode: 'copy', pattern: "2_trimgalore"

  input:
  set sample_id,file(reads) from qf_trimgalore_in

  output:
  set sample_id,'*_val_{1,2}.fq.gz' into trimgalore_pe_out
  set sample_id,'*_unpaired_{1,2}.fq.gz' into trimgalore_unpaired_out
  file '2_trimgalore' into trimgalore_stats_out

  script:
  """
  mkdir 2_trimgalore
  trim_galore -q 20 --fastqc --fastqc_args '-k 10 -t ${task.cpus}' --illumina --paired --gzip \
  	--stringency 5 --length 60 --output_dir 2_trimgalore --trim1 \
  	--retain_unpaired -r1 85 -r2 85 ${reads[0]} ${reads[1]}

  # Make fastq files accessible for staging
  mv 2_trimgalore/*_{val,unpaired}_{1,2}.fq.gz ./
  """
}

process qf_concat_unpaired_reads{
  tag {"${sample_id}"}

  input:
  set sample_id,'unpaired*.fq.gz' from trimgalore_unpaired_out

  output:
  set sample_id,'unpaired.fq.gz' into concat_unpaired_reads_out

  script:
  """
  cat *.fq.gz > unpaired.fq.gz
  """

}

process qf_remove_sispa_adapters_pe{
  tag {"${sample_id}"}
  publishDir "preprocessing/${sample_id}/3_sispa", mode: 'link', pattern: "*.log"

  input:
  set sample_id,'r*.fq.gz' from trimgalore_pe_out

  output:
  set sample_id,'r{1,2}_nosispa.fq.gz' into cutadapt_pe_out
  file "${sample_id}_pe_sispa_cutadapt.log" into cutadapt_pe_log

  script:
  """
  cutadapt --cut=3 -g ^GCCGGAGCTCTGCAGATATC -g ^GGAGCTCTGCAGATATC \
    -G ^GCCGGAGCTCTGCAGATATC -G ^GGAGCTCTGCAGATATC \
    --no-indels --error-rate=0.1 -f 'fastq' \
    --pair-filter=both -m 65 \
    -o r1_nosispa.fq.gz -p r2_nosispa.fq.gz r1.fq.gz r2.fq.gz | tee ${sample_id}_pe_sispa_cutadapt.log
  """
}

process qf_remove_sispa_adapters_se{
  tag {"${sample_id}"}
  publishDir "preprocessing/${sample_id}/3_sispa", mode: 'link', pattern: "*.log"

  input:
  set sample_id,'unpaired*.fq.gz' from concat_unpaired_reads_out

  output:
  set sample_id,'unpaired_nosispa.fq.gz' into cutadapt_unpaired_out
  file "${sample_id}_unpaired_sispa_cutadapt.log" into cutadapt_unpaired_log

  script:
  """
  cutadapt --cut=3 -g ^GCCGGAGCTCTGCAGATATC -g ^GGAGCTCTGCAGATATC \
    --no-indels --error-rate=0.1 -m 65 \
    -o unpaired_nosispa.fq.gz concat_unpaired.fq.gz | tee ${sample_id}_unpaired_sispa_cutadapt.log
  """
}

/**
TODO: This should remove both PhiX and Human
**/
process hostrm_map_to_grch38_pe{
  tag {"${sample_id}"}
  publishDir "preprocessing/${sample_id}/4_hostrm", mode: 'link', pattern: "*.log"

  input:
  set sample_id, 'r*.fq.gz' from cutadapt_pe_out

  output:
  set sample_id, val('pe'), 'pe.sam' into hostrm_map_to_grch38_pe_out
  file "${sample_id}_pe_bowtie2.log" into bowtie2_pe_log

  script:
  """
  bowtie2 --local --very-sensitive-local -t -p ${task.cpus} -x ${params.hostrm_bowtie2_idx} -1 r1.fq.gz -2 r2.fq.gz > pe.sam 2> ${sample_id}_pe_bowtie2.log
  """
}

process hostrm_map_to_grch38_unpaired{
  tag "${sample_id}"
  publishDir "preprocessing/${sample_id}/4_hostrm", mode: 'link', pattern: "*.log"

  input:
  set sample_id, 'unpaired.fq.gz' from cutadapt_unpaired_out

  output:
  set sample_id, val('unpaired'), 'unpaired.sam' into hostrm_map_to_grch38_unpaired_out
  file "${sample_id}_unpaired_bowtie2.log" into bowtie2_se_log

  script:
  """
  bowtie2 --local --very-sensitive-local -t -p ${task.cpus} -x ${params.hostrm_bowtie2_idx} -U unpaired.fq.gz > unpaired.sam 2> ${sample_id}_unpaired_bowtie2.log
  """
}

/*
Split the sam files into two channels, one to convert sam to fastq, and another to
calculate mapping stats
*/

hostrm_map_to_grch38_pe_out.into{ pe_stats_in ;
                                  sam_pe_to_fastq_in }

hostrm_map_to_grch38_unpaired_out.into{ unpaired_stats_in ;
                                        sam_unpaired_to_fastq_in }

//Combine both pe and unpaired sams into single channel for hostrm_mapping_stats
mapping_stats_in = pe_stats_in.mix(unpaired_stats_in)

process hostrm_mapping_stats{
  tag "${sample_id}_${read_type}"
  publishDir "preprocessing/${sample_id}/4_hostrm", mode:'copy'


  input:
  set sample_id, read_type, 'mapped.sam' from mapping_stats_in

  output:
  file "${sample_id}_${read_type}.flagstat" into hostrm_mapping_stats_out

  script:
  """
  samtools view -hSb mapped.sam | samtools flagstat - > ${sample_id}_${read_type}.flagstat
  """
}

process hostrm_sam_pe_to_fastq{
  tag {"${sample_id}"}
  publishDir "preprocessing/${sample_id}", mode:'link'

  input:
  set sample_id, read_type, 'pe.sam' from sam_pe_to_fastq_in

  output:
  set sample_id, read_type, "${sample_id}_*.fq.gz" into sam_pe_to_fastq_out

  script:
  """
  samtools view -hSb -f12 -F256 pe.sam | samtools fastq -1 ${sample_id}_1.fq.gz -2 ${sample_id}_2.fq.gz -
  """
}

process hostrm_sam_unpaired_to_fastq{
  tag{"${sample_id}"}
  publishDir "preprocessing/${sample_id}", mode:'link'

  input:
  set sample_id, read_type, 'unpaired.sam' from sam_unpaired_to_fastq_in

  output:
  set sample_id, read_type, "${sample_id}_unpaired.fq.gz" into sam_unpaired_to_fastq_out

  script:
  """
  samtools view -hSb -f4 -F256 unpaired.sam | samtools fastq -0 ${sample_id}_unpaired.fq.gz -
  """
}
