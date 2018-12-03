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
  $(MEGAHIT_BIN) -m 5e10 -l $$(( 2*$(READ_LEN) )) --k-step 4 --k-max 81 -1 ${reads[0]} -2 ${reads[1]} --cpu-only -t ${task.cpus} -o megahit
  """

}

process asm_metaspades{
  tag { "${sample_id}" }

  input:
  set sample_id, reads from asm_megahit_in

  output:
  set sample_id,'metaspades',file "spades/contigs.fasta" optional true into asm_spades_out

  script:
  """
  $(SPADES_BIN) -t ${task.cpus} --pe1-12 $< --pe1-s $(word 2,$^) --s1 $(word 3,$^) -o metaspades --tmp-dir \${TMP_DIR} -m 64
  """
}

// process asm_iva{}

params.min_ctg_size=500

process asm_filter_contigs{
  tag { "${sample_id}/${assembler}" }
  publishDir

  input:
  set sample_id,assembler,"contigs.fa" from assemblies

  output:
  set sample_id,assembler,"contigs_min${params.min_ctg_size}bp.fa" into assemblies

  script:
  """
  $(SEQTK_BIN) seq -L params.min_ctg_size contigs.fa > contigs_min${params.min_ctg_size}bp.fa
  """
}


/**
NOTE: Map all as single ends maybe(?)
NOTE: Use BURST maybe or BBmap ?
*/
process asm_map_reads_to_contigs{
  script:
  """
  $(BWA_BIN) mem -t $(threads) -T 30 -M -p $(basename $(word 2,$^)) $< | $(SAMTOOLS_BIN) view -hSb -o $@ -
  """
}

/**
TAX ASSIGNMENT - READS
**/

process tax_reads_metaphlan2{
  script:
  """
  ifndef TMP_DIR
  $(mpa_bin) --mpa_pkl $(mpa_pkl) --bowtie2db $(mpa_bowtie2db) \
  		--bowtie2out $(TMP_DIR )/nohuman.bowtie2.bz2 --nproc $(threads) --input_type multifastq \
  		--sample_id_key $(sample_name) $< $(word 1,$@)
  """
}

/**
TODO: Choose between  NCBI RefSeq or  IMG/VR database
**/
process tax_reads_FastViromeExplorer{
  script:
  """
  java -cp bin FastViromeExplorer -1 $read1File -2 $read2File -i /path-to-index-file/ncbi-virus-kallisto-index-k31.idx -o $outputDirectory
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
  hmmscan --cpu $(threads) $(hmmscan_opts) -o $@ $| $<
  """
}
