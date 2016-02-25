
  ```mkdir ~/ngs_course/rnaseq/sailfish
  cd ~/ngs_course/rnaseq/sailfish
  
  # sailfish index sailfish index -p <num of cores> -k <kmer size> -t <fasta of gene sequences> -o <folder name>
  
  export PATH=/groups/bcbio/bcbio/anaconda/bin:/opt/bcbio/local/bin:$PATH
  
  sailfish quant -i /groups/hbctraining/sailfish-run/sailfish.ensembl2.idx/ -l SR -r ngs_course/rnaseq/data/untrimmed_fastq/Mov10_oe_1.subset.fq --useVBOpt -o Mov10_oe_1.subset.sailfish
  
  less Mov10_oe_1.subset.sailfish/quant.sf```
  
  
