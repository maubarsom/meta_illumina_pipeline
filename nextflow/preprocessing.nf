
process trimgalore_qf{

  script:
  """
  trim_galore -q 20 --fastqc --fastqc_args "-k 10 -t 16" --illumina --paired --gzip \
  	--stringency 5 --length 60 --output_dir ${OUT_DIR} --trim1 \
  	--retain_unpaired -r1 85 -r2 85 $(ls ${READS_DIR}/*.fastq.gz | sort)
  """
}

process remove_sispa_adapters_pe{
  script:
  """
  cutadapt --cut=3 -g ^GCCGGAGCTCTGCAGATATC -g ^GGAGCTCTGCAGATATC --no-indels --error-rate=0.1 -f 'fastq' \
    -o R1_PIPE ${IN_DIR}/*_val_1.fq.gz 2>&1 > log/${SAMPLE_ID}_overhangs_R1.log &

  cutadapt --cut=3 -g ^GCCGGAGCTCTGCAGATATC -g ^GGAGCTCTGCAGATATC --no-indels --error-rate=0.1 -f 'fastq' \
    -o R2_PIPE ${IN_DIR}/*_val_2.fq.gz 2>&1 > log/${SAMPLE_ID}_overhangs_R2.log &
  """
}

process remove_sispa_adapters_se{

  script:
  """
  #Merge unpaired into a single file
  cat ${IN_DIR}/*_unpaired_2.fq.gz >> ${IN_DIR}/*_unpaired_1.fq.gz && rm ${IN_DIR}/*_unpaired_2.fq.gz

  #Requires cutadapt
  cutadapt --cut=3 -g ^GCCGGAGCTCTGCAGATATC -g ^GGAGCTCTGCAGATATC --no-indels --error-rate=0.1 \
  	-o SE_PIPE ${IN_DIR}/*_unpaired_1.fq.gz 2>&1 > log/${SAMPLE_ID}_overhangs_unpaired.log &
  """
}


/**
TODO: This should remove both PhiX and Human
**/
process remove_human_pe{
  script:
  """
  bowtie2 --local --very-sensitive-local -t -p 12 -x ${BOWTIE_DB} -1 R1_PIPE -2 R2_PIPE | tee >(samtools view -hSb - | samtools flagstat - > log/${SAMPLE_ID}_bowtie_pe.flagstat ) | samtools view -hSb -f12 -F256 - | samtools fastq - | gzip > ${OUT_FILE}
  """
}

process remove_human_se{
  script:
  """
  bowtie2 --local --very-sensitive-local -t -p 8 -x ${BOWTIE_DB} -U SE_PIPE | tee >(samtools view -hSb - | samtools flagstat - > log/${SAMPLE_ID}_bowtie_unpaired.flagstat ) | samtools view -hSb -f4 -F256 - | samtools fastq - | gzip > ${OUT_FILE}
  """
}
