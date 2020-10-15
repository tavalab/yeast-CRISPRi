# yeast-CRISPRi
scripts and supporting files for processing yeast CRISPRi library.

"An inducible CRISPR-interference library for genetic interrogation of Saccharomyces cerevisiae biology"
Amir Momen-Roknabadi, Panos Oikonomou, Maxwell Zegans, Saeed Tavazoie
https://doi.org/10.1101/2020.03.05.978619

command: python process_01_Sequencing_yeast_CRISPRi.py

- sequencing files *.fastq.gz should be placed in the folder data/ in the same directory as the script.
- exemplar fastq.gz are provided with this distribution 
- if new sequencing files are provided they should be appropriately listed in process_01_Sequencing_yeast_CRISPRi.py: in_prefixes = ["Idx3_TTAGGCAT", "Idx1_ATCACGAT"]
- sequencing files should follow the naming scheme AllLanes_Idx1_ATCACGAT_R1_merged.fastq.gz, 
  where ATCACGAT is the illumina index used that should be present in the first column of data/barcodes/barcode_map.tsv
- all internal indexes listed in data/barcodes/myBarcodes.fa should be present in the second column of data/barcodes/barcode_map.tsv
- the output of this script is a count matrix file: all_counts.tsv

command: python process_02_Calculate_stats_yeast_CRISPRi.py

- This script reads the matrix containing read counts for every sample, normalizes
- deplection scores are calculated between sample pairs listed in the dictionary sampleComparisons
- calculates the depletion scores for every gene from its corresponding gRNAs
- calculates the values for the Synthetic scrambled genes
- saves all the values for both the genes and the Synthetic scrambled genes in dir_SAVE
- if DO_PSEUDOCORRECTION=1, it calculates z-scores for the genes based on the distribution of the Synthetic scrambled genes
- average depletion scores are saved in Data_Enrichment_Score.txt
- all_counts.example.1.tsv is provided as an example of a count matrix, should be replaced with the output of process_01_Sequencing_yeast_CRISPRi.py

File: data/library/crspri_gRNA_properties.tsv

- contains all the attributes of the sgRNAs for the Random Forest model


