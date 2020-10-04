# yeast-CRISPRi
scripts and supporting files for processing yeast CRISPRi library.

"An inducible CRISPR-interference library for genetic interrogation of Saccharomyces cerevisiae biology"
Amir Momen-Roknabadi, Panos Oikonomou, Maxwell Zegans, Saeed Tavazoie
https://doi.org/10.1101/2020.03.05.978619


- command: python process_01_Sequencing_yeast_CRISPRi.py
- sequencing files *.fastq.gz should be placed in the folder data/ in the same directory as the script.
- exemplar fastq.gz are provided with this distribution 
- if new sequencing files are provided they should be appropriately listed in process_01_Sequencing_yeast_CRISPRi.py: in_prefixes = ["Idx3_TTAGGCAT", "Idx1_ATCACGAT"]
- sequencing files should follow the naming scheme AllLanes_Idx1_ATCACGAT_R1_merged.fastq.gz, 
  where ATCACGAT is the illumina index used that should be present in the first column of data/barcodes/barcode_map.tsv
- all internal indexes listed in data/barcodes/myBarcodes.fa should be present in the second column of data/barcodes/barcode_map.tsv
- the output of this script is a count matrix file: all_counts.tsv
