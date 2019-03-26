#!/usr/bin/env nextflow

/*
How to run:
nextflow -C discovery.nf.config run discovery.nf --fastq_files preprocessing -profile hamlet
*/

params.fastq_dir='preprocessing'
fastq_files = Channel.fromFilePairs("${params.fastq_dir}/**/*_{1,2,unpaired}.fq.gz",size:3)


fastq_files.into{
  asm_megahit_in;
  asm_metaspades_in;
  reads_to_map;
  tax_reads_metaphlan2_in;
  tax_reads_kraken2_in;
  tax_reads_FastViromeExplorer_in
 }

/**
ASSEMBLY Module

TODO: Decide if run megahit in meta-sensitive mode or in default?
TODO: Run also SPAdes standard pipeline or only MetaSpades?
**/

process asm_megahit{
  tag { "${sample_id}" }
  publishDir "${params.publish_base_dir}/${sample_id}/megahit", mode: 'copy', pattern: "1_assembly"

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

/*
NOTE: MetaSPAdes does not support incorporating the single reads
TODO: figure out how to use a scratch folder? --tmp-dir opt
*/
process asm_metaspades{
  tag { "${sample_id}" }
  publishDir "${params.publish_base_dir}/${sample_id}/metaspades", mode: 'copy', pattern: "1_assembly"

  input:
  //set sample_id, reads from asm_metaspades_in
  set sample_id, reads from Channel.empty()

  output:
  set sample_id,val('metaspades'),"contigs.fa" optional true into asm_metaspades_out
  file "1_assembly" optional true into asm_metaspades_dir_out

  script:
  """
  spades.py --meta -t ${task.cpus} -1 ${reads[0]} -2 ${reads[1]} -o 1_assembly -m ${params.max_spades_mem}
  if [ -s 1_assembly/contigs.fasta ]; then ln 1_assembly/contigs.fasta contigs.fa; fi
  """
}

/*
NOTE: Combine contigs into a single channels
*/
all_assemblies = asm_megahit_out.mix(asm_metaspades_out)

process asm_filter_contigs{
  tag { "${sample_id}/${assembler}" }
  publishDir "${params.publish_base_dir}/${sample_id}/${assembler}/2_filt_contigs", mode:'link'

  input:
  set sample_id,assembler,"contigs.fa" from all_assemblies

  output:
  set sample_id,assembler,"contigs_filt.fa" into asm_filter_contigs_out

  script:
  """
  seqtk seq -L ${params.min_ctg_size} contigs.fa > contigs_filt.fa
  """
}

/*
Clone the contigs into different channels
*/
asm_filter_contigs_out.into{
  contigs_to_be_mapped;
  tax_contigs_diamond_in;
  tax_contigs_kraken2_in;
  tax_contigs_virfinder_in;
  tax_contigs_fgs_in;
  tax_contigs_virsorter_in
  /*TODO: COMPLETE THIS LIST */
}

/*
Combine reads and contigs for mapping
*/
asm_map_reads_to_contigs_in = reads_to_map.cross(contigs_to_be_mapped).map{ it[0]+it[1][1,2] }

/**
NOTE: bbwrap/bbmap parameters
  * kfilter is the minimum length of consecutive exact matches for an alignment (min kmer size of dBg assembler)
  * maxindel limits the indel size
  * subfilter limits the number of mismatches for an alignment
*/
process asm_map_reads_to_contigs{
  tag { "${sample_id}/${assembler}" }

  publishDir "${params.publish_base_dir}/${sample_id}/${assembler}/2_filt_contigs", mode:'link'

  input:
  set sample_id, 'reads_*.fq.gz',assembler,"contigs_filt.fa" from asm_map_reads_to_contigs_in

  output:
  set sample_id, assembler,"reads_to_contigs.sam.gz" into asm_map_reads_to_contigs_out

  script:
  """
  bbwrap.sh ref=contigs_filt.fa in=reads_1.fq.gz,reads_3.fq.gz in2=reads_2.fq.gz,null out=reads_to_contigs.sam.gz usejni=t kfilter=22 subfilter=15 maxindel=80
  """
}

/* Send the mapped reads to
1) mapping stats with samtools flagstat
2) Calculate coverage with bbtools pileup.sh
*/
asm_map_reads_to_contigs_out.into{ asm_mapping_stats_in;
                                   asm_per_ctg_coverage_in }


process asm_mapping_stats{
  tag {"${sample_id}/${assembler}"}

  publishDir "${params.publish_base_dir}/${sample_id}/${assembler}/2_filt_contigs", mode:'link'

  input:
  set sample_id, assembler,'reads_to_contigs.sam.gz' from asm_mapping_stats_in

  output:
  set sample_id, assembler, "${sample_id}_${assembler}_flagstat.txt" into asm_mapping_stats_out

  script:
  """
  samtools flagstat reads_to_contigs.sam.gz > ${sample_id}_${assembler}_flagstat.txt
  """
}

process asm_per_ctg_coverage{
  tag {"${sample_id}/${assembler}"}

  publishDir "${params.publish_base_dir}/${sample_id}/${assembler}/2_filt_contigs", mode:'link'

  input:
  set sample_id, assembler,"reads_to_contigs.sam.gz" from asm_per_ctg_coverage_in

  output:
  set sample_id, assembler,"reads_to_contigs.cov.txt" into asm_per_ctg_coverage_out

  script:
  """
  pileup.sh in=reads_to_contigs.sam.gz out=reads_to_contigs.cov.txt
  """
}
/**
TAX ASSIGNMENT - READS
**/

process tax_reads_metaphlan2{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/reads", mode:'link'

  input:
  set sample_id, 'reads_*.fq.gz' from tax_reads_metaphlan2_in

  output:
  set sample_id, "${sample_id}_metaphlan2.tsv" into tax_reads_metaphlan2_out

  script:
  """
  metaphlan2.py --nproc ${task.cpus} --input_type multifastq --sample_id_key '#clade' \
  		--sample_id '${sample_id}' <(zcat reads_*.fq.gz) ${sample_id}_metaphlan2.tsv
  """
}

process tax_reads_kraken2{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/reads", mode:'link'

  input:
  set sample_id, 'reads_*.fq.gz' from tax_reads_kraken2_in

  output:
  set sample_id, "${sample_id}_kraken2.txt","${sample_id}_kraken2_report.txt" into tax_reads_kraken2_out

  script:
  """
  kraken2 --db ${params.kraken2_db} --threads ${task.cpus} --output ${sample_id}_kraken2.txt \
    --report ${sample_id}_kraken2_report.txt --gzip-compressed reads_*.fq.gz
  """
}

/**
TODO: Choose between  NCBI RefSeq or  IMG/VR database
**/
process tax_reads_FastViromeExplorer{
  tag {"${sample_id}"}

  publishDir "${params.publish_base_dir}/${sample_id}/reads", mode:'link'

  input:
  set sample_id, 'reads_*.fq.gz' from tax_reads_FastViromeExplorer_in

  output:
  file "${sample_id}_fastviromeexplorer.sam" into fve_out_1
  file "${sample_id}_fastviromeexplorer_abundance.tsv" into fve_out_2

  script:
  """
  java -cp \${FVE_PATH}/bin FastViromeExplorer -i ${params.FVE_index} -l ${params.FVE_viruslist} -1 reads_1.fq.gz -2 reads_2.fq.gz  -o ./
  mv FastViromeExplorer-reads-mapped-sorted.sam ${sample_id}_fastviromeexplorer.sam
  mv FastViromeExplorer-final-sorted-abundance.tsv ${sample_id}_fastviromeexplorer_abundance.tsv
  """
}

/**
TAX ASSIGNMENT - CONTIGS
**/

process tax_contigs_kraken2{
  tag {"${sample_id}_${assembler}"}

  publishDir "${params.publish_base_dir}/${sample_id}/${assembler}", mode:'link'

  input:
  set sample_id, assembler, 'contigs.fa' from tax_contigs_kraken2_in

  output:
  set sample_id, assembler,"${sample_id}_${assembler}_kraken2.txt","${sample_id}_${assembler}_kraken2_report.txt" into tax_contigs_kraken2_out
  file "${sample_id}_${assembler}_kraken2_unmapped.fa" into tax_contigs_kraken2_unmapped

  script:
  """
  kraken2 --db ${params.kraken2_db} --threads ${task.cpus} --use-names \
    --unclassified-out ${sample_id}_${assembler}_kraken2_unmapped.fa \
    --output ${sample_id}_${assembler}_kraken2.txt \
    --report ${sample_id}_${assembler}_kraken2_report.txt contigs.fa
  """
}


process tax_contigs_diamond{
  tag {"${sample_id}_${assembler}"}

  publishDir "${params.publish_base_dir}/${sample_id}/${assembler}", mode:'link'

  input:
  set sample_id, assembler, 'contigs.fa' from tax_contigs_diamond_in

  output:
  set sample_id, assembler,"${sample_id}_${assembler}_diamond.tsv" into tax_contigs_diamond_out
  file "${sample_id}_${assembler}_diamond_unmapped.fa" into diamond_unmapped_out
  file 'diamond_blast_cols.txt' into tax_diamond_blast_cols

  script:
  blast_cols='qseqid sseqid qstart qend sstart send evalue score length pident nident mismatch positive ppos gapopen gaps staxids'
  """
  diamond blastx --sensitive --masking 1 -p ${task.cpus} --db ${params.diamond_db} \
    --taxonmap ${params.diamond_taxonmap} --taxonnodes ${params.diamond_taxonnodes} \
    --query contigs.fa --out ${sample_id}_${assembler}_diamond.tsv \
    --outfmt 6 ${blast_cols} \
    --un ${sample_id}_${assembler}_diamond_unmapped.fa
  echo '${blast_cols}' > diamond_blast_cols.txt 
  """
}

process tax_contigs_diamond_view{
  tag {"${sample_id}_${assembler}"}

  publishDir "${params.publish_base_dir}/${sample_id}/${assembler}", mode:'link'

  input:
  set sample_id, assembler, 'diamond.daa' from Channel.empty()

  output:
  set sample_id, assembler,"${sample_id}_${assembler}_diamond.daa" into tax_contigs_diamond_view_out

  script:
  """
  diamond view -p ${task.cpus} --db ${params.diamond_db} \
    --outfmt 6 qseqid sseqid qstart qend sstart send evalue score length pident nident mismatch positive ppos gapopen gaps staxids \
    --daa diamond.daa \
    --out ${sample_id}_${assembler}_diamond.tsv
  """
}


/*
NOTE: https://github.com/EnvGen/toolbox/tree/master/scripts/assign_taxonomy_from_blast
*/
//process tax_diamond_lca{}


process tax_contigs_virfinder{
  tag {"${sample_id}_${assembler}"}

  publishDir "${params.publish_base_dir}/${sample_id}/${assembler}", mode:'link'

  input:
  set sample_id, assembler, 'contigs.fa' from tax_contigs_virfinder_in

  output:
  set sample_id, assembler,"${sample_id}_${assembler}_virfinder.csv" into tax_contigs_virfinder_out

  script:
  """
  #!/usr/bin/env Rscript
  library(VirFinder)
  predResult <- VF.pred("contigs.fa")
  write.csv(predResult,file="${sample_id}_${assembler}_virfinder.csv")
  """
}


process tax_contigs_FragGeneScan{
  tag {"${sample_id}_${assembler}"}

  publishDir "${params.publish_base_dir}/${sample_id}/${assembler}/orfs", mode:'link'

  input:
  set sample_id, assembler, 'contigs.fa' from tax_contigs_fgs_in

  output:
  set sample_id, assembler,"${sample_id}_${assembler}_fgs_orfs.*" into tax_contigs_fgs_out

  script:
  """
  run_FragGeneScan.pl -thread ${task.cpus} -complete 1 -train=illumina_10 \
    -genome=contigs.fa -out=${sample_id}_${assembler}_fgs_orfs
  """
}

/*
process tax_orfs_hmmscan{
  input:
  set sample_id,contigs_id from dummy_in
  each pfam from pfams

  script:
  """
  hmmscan --cpu ${task.cpus} $(hmmscan_opts) -o $@ $| $<
  """
}
*/

//process tax_contigs_virsorter{}
