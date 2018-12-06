#!/usr/bin/env nextflow

params.fastq_dir='/users/maubar/fsbio/humanrm'
fastq_files = Channel.fromFilePairs("${params.fastq_dir}/**/*_{1,2,se}.fastq.gz",size:3)

/**
ASSEMBLY Module

TODO: Decide if run megahit in meta-sensitive mode or in default?
TODO: Run also SPAdes standard pipeline or only MetaSpades?
**/

process asm_megahit{
  tag { "${sample_id}" }
  publishDir "results/${sample_id}/megahit", mode: 'copy', pattern: "1_assembly"

  input:
  set sample_id, reads from asm_megahit_in

  output:
  set sample_id,val('megahit'), "contigs.fa" optional true into asm_megahit_out
  file "1_assembly" optional true into asm_megahit_dir_out

  script:
  """
  megahit -t ${task.cpus} --presets meta-sensitive -1 ${reads[0]} -2 ${reads[1]} -r ${reads[2]} --cpu-only  -o 1_assembly
  if [ -s 1_assembly/final.contigs.fa ]; then ln 1_assembly/final.contigs.fa contigs.fa; fi
  """
}

/* NOTE: MetaSPAdes does not support incorporating the single reads */
process asm_metaspades{
  tag { "${sample_id}" }
  publishDir "results/${sample_id}/metaspades"

  input:
  set sample_id, reads from asm_megahit_in

  output:
  set sample_id,val('metaspades'),"contigs.fa" optional true into asm_spades_out
  file "1_assembly" optional true into asm_spades_dir_out

  script:
  """
  spades.py --meta -t ${task.cpus} -1 ${reads[0]} -2 ${reads[1]} -o 1_assembly --tmp-dir \${TMP_DIR} -m 124
  if [ -s 1_assembly/contigs.fasta ]; then ln 1_assembly/contigs.fasta contigs.fa; fi
  """
}

process asm_filter_contigs{
  tag { "${sample_id}/${assembler}" }
  publishDir "results/${sample_id}/${assembler}/2_filt_contigs"

  input:
  set sample_id,assembler,"contigs.fa" from assemblies

  output:
  set sample_id,assembler,"contigs_filt.fa" into filtered_contigs

  script:
  """
  seqtk seq -L ${params.min_ctg_size} contigs.fa > contigs_filt.fa
  """
}

/**
NOTE: bbwrap/bbmap parameters
  * kfilter is the minimum length of consecutive exact matches for an alignment (min kmer size of dBg assembler)
  * maxindel limits the indel size
  * subfilter limits the number of mismatches for an alignment
*/
process asm_map_reads_to_contigs{
  tag { "${sample_id}/${assembler}" }

  publishDir "results/${sample_id}/${assembler}/2_filt_contigs"

  input:
  set sample_id, assembler,"contigs_filt.fa" from map_reads_in

  output:
  set sample_id, assembler,"reads_to_contigs.sam.gz"

  script:
  """
  bbwrap.sh ref=contigs_filt.fa in=${reads[0]},${reads[2]} in2=${reads[1]},null out=reads_to_contigs.sam.gz usejni=t kfilter=22 subfilter=15 maxindel=80
  """
}

process asm_mapping_stats{
  input:
  set sample_id, assembler,"reads_to_contigs.sam.gz"

  output:
  set sample_id, assembler,"samtools_flagstat.txt"

  script:
  """
  gunzip -c ${input_bam} | samtools flagstat - > samtools_flagstat.txt
  """
}

process asm_per_ctg_coverage{
  input:
  set sample_id, assembler,"reads_to_contigs.sam.gz"

  output:
  set sample_id, assembler,"reads_to_contigs.cov.txt"

  script:
  """
  pileup.sh in=reads_to_contigs.sam.gz out=reads_to_contigs.cov.txt
  """
}
/**
TAX ASSIGNMENT - READS
**/

process tax_reads_metaphlan2{
  script:
  """
  metaphlan2.py --mpa_pkl ${mpa2_pkl} --bowtie2db ${mpa2_bowtie2db} \
  		--bowtie2out /dev/null --nproc ${task.cpus} --input_type multifastq \
  		--sample_id_key ${sample_id} ${reads[0]} ${reads[1]}
  """
}

/**
TODO: Choose between  NCBI RefSeq or  IMG/VR database
**/
process tax_reads_FastViromeExplorer{
  script:
  """
  java -cp bin FastViromeExplorer -1 ${reads[0]} -2 ${reads[1]} -i ${FastViromeExplorer_index} -o $outputDirectory
  """
}

process tax_reads_kraken2{

}

process tax_reads_kraken2_report{

}

/**
TAX ASSIGNMENT - CONTIGS
**/

process tax_contigs_kraken2{

}

process tax_contigs_diamond{
  script:
  """
  """
}

/*
NOTE: https://github.com/EnvGen/toolbox/tree/master/scripts/assign_taxonomy_from_blast
*/
process tax_diamond_lca{

}

process tax_contigs_virfinder{
  script:
  """
  #!/usr/bin/env Rscript
  library(VirFinder)
  predResult <- VF.pred("${lolol}")
  write.tsv(predResult,file="jojoj.tsv")
  """
}

process tax_contigs_virsorter{
  script:
  """
  """
}

process tax_contigs_FragGeneScan{
  script:
  """
  $(FGS_BIN) -genome=$^ -out=$(basename $@) -complete=0 -train=illumina_10
  """
}

process tax_orfs_hmmscan{
  input:
  set sample_id,contigs_id from dummy_in
  each pfam from pfams

  script:
  """
  hmmscan --cpu ${task.cpus} $(hmmscan_opts) -o $@ $| $<
  """
}
