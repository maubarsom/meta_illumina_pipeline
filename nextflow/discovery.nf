/**
ASSEMBLY STEPS
**/

process asm_megahit{


  output:
  file "megahit/final.contigs.fa" into asm_megahit_out

  script:
  """
  $(MEGAHIT_BIN) -m 5e10 -l $$(( 2*$(READ_LEN) )) --k-step 4 --k-max 81 --12 $< -r $(word 2,$^),$(word 3,$^) --cpu-only -t $(threads) -o megahit
  """

}

process asm_metaspades{
  output:
  file "spades/contigs.fasta" into asm_spades_out

  script:
  """
  $(SPADES_BIN) -t $(threads) --pe1-12 $< --pe1-s $(word 2,$^) --s1 $(word 3,$^) -o $(dir $@) --tmp-dir $(TMP_DIR) -m 64
  """
}

process asm_iva{

}

process asm_filter_contigs{
  script:
  """
  $(SEQTK_BIN) seq -L 500 $< | awk -vPREFIX=$*  'BEGIN {counter=1; split(PREFIX,fields,"_asm_"); ASM=fields[2];} /^>/ {print ">" ASM "_" counter ;counter+=1;} ! /^>/{print $0;}' > $@
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
Taxonomical assignment
**/

process tax_reads_metaphlan{
  script:
  """
  ifndef TMP_DIR
  $(mpa_bin) --mpa_pkl $(mpa_pkl) --bowtie2db $(mpa_bowtie2db) \
  		--bowtie2out $(TMP_DIR )/nohuman.bowtie2.bz2 --nproc $(threads) --input_type multifastq \
  		--sample_id_key $(sample_name) --biom $*_nohuman_metaphlan.biom $< $(word 1,$@)
  """
}

process tax_contigs_mmseq2{
  script:
  """
  """
}

process tax_diamond{

  script:
  """
  """
}
