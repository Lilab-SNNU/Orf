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

- bedtools (v 2.27.1)
- bowtie (v 1.2.3)
- STAR (v 2.7.8)
- Python3 (v 3.7.6)
- R language (v 3.6.3)
- ViennaRNA (v 2.4.18)
- cutadapt (v 1.18)

## python3 package:

- Biopython (v 1.78)
- Pandas (v 1.0.1)

## R package:

- UpSetR(v 1.4.0)
- edgeR(v 3.28.1)
- ggplot2(v 3.3.5)

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
  
  - Filter and obtain the reads which belongs to circRNA

  ```
   python filter_coding_reads_circ.py -y config.yaml
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

 - transcript_fasta: Fill in the absolute path of the species transcript genome related to circRNA(not the relative path!).
   
 - raw_reads: Fill in the absolute path of the Ribo-Seq data(One or two tissue samples data in the same species) related to your interest species(not the relative path!).
   
 - ribosome_fasta: Fill in the absolute path of the rRNA data related to your interset species(not the relative path!).
 
 - circrnas: Fill in the absolute path of the candidate circRNA(not the relative path!).You need to ensure that the id format of your circRNA is the same as the circrna_gtf file.
 
 - circrna_gtf: Fill in the absolute path of the candidate circRNA annoated file(not the relative path!). This gtf file should be like this type(make by yourself), samples are optional, you can input NaN if you don't need.

'''
circRNA_id	strand	circRNA_lenth	gene_id	samples
TC-hsa-RERE_0160	-	185	ENSG00000142599.17	NaN
TC-hsa-UBE2J2_0010	-	607	ENSG00000160087.20	NaN
'''
 
 - result_file_location: Fill in the absolute path of a folder to save the final results(not the relative path!).
 
 - tmp_file_location: Fill in the absolute path of a folder to save the temporary files(not the relative path!).
 
 - thread: Fill the number of threads running.
 
 - trimmomatic_jar: Fill in the absolute path, we have provided the trimmomatic_jar file in CircOrf package.

 - pos_fa/neg_fa: Fill in the absolute path of the fasta file for positive samples/negative samples(the length of sequences should be no more than 101bp)
 
 - ribotype: The type of Ribo-seq data(sra or fastq).
  
  
## Citation

CircOrf: A useful tool for Identification and Visualization of circRNA with coding potential

## Contact us

If you encounter any problems while using CircOrf, please send an email to mr1997@snnu.edu.cn, and we will resolve it as soon as possible.
