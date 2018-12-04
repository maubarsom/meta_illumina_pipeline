#!/usr/bin/env nextflow

params.fastq_dir='/users/maubar/fsbio/humanrm'
fastq_files = Channel.fromFilePairs("${params.fastq_dir}/**/*_{1,2}.fastq.gz")

/**
ASSEMBLY
**/

process asm_megahit{
  tag { "${sample_id}" }
  publishDir "results/${sample_id}/megahit/contigs"

  input:
  set sample_id, reads from asm_megahit_in

  output:
  set sample_id,'megahit',file "megahit/final.contigs.fa" into asm_megahit_out

  script:
  """
  megahit -m 5e10 -l $$(( 2*$(READ_LEN) )) --k-step 4 --k-max 81 -1 ${reads[0]} -2 ${reads[1]} --cpu-only -t ${task.cpus} -o megahit
  """

}

process asm_metaspades{
  tag { "${sample_id}" }
  publishDir "results/${sample_id}/spades/contigs"

  input:
  set sample_id, reads from asm_megahit_in

  output:
  set sample_id,'metaspades',file "spades/contigs.fasta" optional true into asm_spades_out

  script:
  """
  spades.py -t ${task.cpus} --pe1-12 $< --pe1-s $(word 2,$^) --s1 $(word 3,$^) -o metaspades --tmp-dir \${TMP_DIR} -m 64
  """
}

// process asm_iva{}

params.min_ctg_size=500

process asm_filter_contigs{
  tag { "${sample_id}/${assembler}" }

  input:
  set sample_id,assembler,"contigs.fa" from assemblies

  output:
  set sample_id,assembler,"contigs_min${params.min_ctg_size}bp.fa" into assemblies

  script:
  """
  seqtk seq -L ${params.min_ctg_size} contigs.fa > contigs_min${params.min_ctg_size}bp.fa
  """
}

process asm_contig_index{
  tag { "${sample_id}/${assembler}" }

}

/**
NOTE: Map all as single ends maybe(?)
NOTE: Use BURST maybe or BBmap ?
*/
process asm_map_reads_to_contigs{
  input:
  set sample_id, assembler,

  output:
  set sample_id, assembler,"${sample_id}_${assembler}.bam"

  script:
  """
  bwa mem -t ${task.cpus} -T 30 -M -p ${bwa_idx} ${reads} | samtools view -hSb -o $@ -
  """
}

process asm_mapping_stats{

  script:
  """
  samtools flagstat ${input_bam} > ${sample_id}_${assembler}_flagstat.txt
  """
}

process asm_mapping_depth{
  script:
  """
  samtools sort $< | $(SAMTOOLS_BIN) depth - | gzip > $@
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


/**
TAX ASSIGNMENT - CONTIGS
**/

process tax_contigs_kraken2{

}

//process tax_contigs_mmseq2{ }

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
