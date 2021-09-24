# CircOrf

CircOrf is a Python3-base pipeline for the prediction and visualization of circRNA and ORF.

# Requirement
## Data:

- genome fasta file
- genome gtf file
- ribo-seq sra file
- rRNA fasta file
- circRNA fasta file
- circRNA gtf file
- riboseq_adapters fasta file

## Software:

- bedtools (v.2.26.0+)
- bowtie (v.1.2.2+)
- STAR (v.2.7.1+)
- Python3 (v.3.6.5+)
- R language (v.3.4.4+)
- ViennaRNA (v 2.4.18+)

## python3 package:

- Biopython (v.1.72+)
- Pandas (v.0.23.3+)

## R package:

- UpsetR
- edgeR
- ggplot2

# Quick Start
Usually you can download the package from github simply,and then:
```
git clone git://github.com/Lilab-SNNU/CircOrf.git
cd CircOrf
tar -zxvf requiredSoft.tar.gz
cd requiredSoft
chmod 777 STAR
```


# Usage

Attention: Before you begin to use this package, you need to make sure that you have install the required software and add them to the environment variables.


1. Fill the config file (config.yaml), input absolute path of each required file.

2. If you want to run CircOrf by one script:

 ```
   sh one_step_script.sh config.yaml
 ```

3. If you want to run CircOrf step by step:


  - Make index genomes

  ```
   python make_index_genomes.py -y config.yaml
  ```
  
  - Filter and obtain the reads which belongs to circRNA

  ```
   python filter_coding_reads.py -y config.yaml
  ```
  
  - Filter and obtain the coding circRNAs

  ```
   python filter_coding_circ.py -y config.yaml
  ```
  
  - Identify the ORF from coding circRNAs by some features
  
  ```
  python filter_coding_orf.py -y config.yaml
  ```
  
  - Visualize the IRES region secondary structure, the circ annoated graph, the express_analysis of the circ and the distribution of the predicted ORFs
  
  ```
  python visual_circ_orf.py -y config.yaml
  ```

### How to fill in the config.yaml file?

When opening the config file in text format, there are some lines that need to be filled in, they are:

 - genome_name: It is your program name, so you can fill in any value here in text form (English letters).

 - genome_fasta: Fill in the absolute path of the species genome related to circRNA(not the relative path!).
   
 - genome_gtf: Fill in the absolute path of the species genome annoated file related to circRNA(not the relative path!).
   
 - raw_reads: Fill in the absolute path of the Ribo-Seq data(One or two tissue samples data in the same species) related to your interest species(not the relative path!).
   
 - ribosome_fasta: Fill in the absolute path of the rRNA data related to your interset species(not the relative path!).
 
 - circrnas: Fill in the absolute path of the candidate circRNA(not the relative path!).You need to ensure that the id format of your circRNA is similar to '>hsa_circ_0000003|chr1:1423242-1459777+|NM_031921|ATAD3B'.
 
 - circrna_gtf: Fill in the absolute path of the candidate circRNA annoated file(not the relative path!).
 
 - riboseq_adapters: The absolute path of the adapters file related to the Ribo-Seq data(not the relative path!).
 
 - result_file_location: Fill in the absolute path of a folder to save the final results(not the relative path!).
 
 - tmp_file_location: Fill in the absolute path of a folder to save the temporary files(not the relative path!).
 
 - coverage_counts: Set the minimum counts of the reads that mapping to the junctions.
 
 - thread: Fill the number of threads running.
 
 - lenth_need: Fill the IRES region lenth of the ORF upstream that you need.
 
 - ires_score: Set the minimum ires score that you think.
 
 - orf_score: Set the minimum orf score that you think.
 
 - trimmomatic_jar:The absolute path to the trimmomatic_jar file, we have provided the trimmomatic_jar file in CircOrf package, you just need to fill in the absolute path of this file on your computer.
   
 - reads_type: The type of sequencing data(single or pair).
 
 - ribotype: The type of Ribo-seq data(sra or fastq.gz).
  
  
## Citation

CircOrf: A useful tool for Identification and Visualization of circRNA with coding potential

## Contact us

If you encounter any problems while using CircOrf, please send an email to mr1997@snnu.edu.cn, and we will resolve it as soon as possible.
