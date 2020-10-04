##      Panos Oikonomou, Columbia University
##      Saeed Tavazoie Lab
##
## 	(c) 2020 The Trustees of Columbia University in the City of New York. 
##	This work may be reproduced and distributed for academic 
##	non-commercial purposes only without further authorization, 
##	but rightsholder otherwise reserves all rights.
##
##	
##	requirements: python 2.7, cutadapt v1.8.3, bowtie2 v2.2.6, samtools v1.2, bedtools v2.17.0
##	

import subprocess
import shlex
import tempfile
import os
import re
import glob 
import string
import numpy as np 

from multiprocessing import Pool

from Bio import SeqIO
import csv

DO_TRIM_ADAPTERS	= False;
DO_ALIGN_SEQUENCES	= True;
DO_PROCESS_SAMFILES	= True;
DO_MERGE_FILES		= True;

READDIR='data/'
SAVEDIR='data/out/'
TEMPDIR='data/out/intermediate/'
TRIMDIR='data/out/fastq_trimmed/'
CNTDIR ='data/out/counts/'


# adapter sequences to clip
PHRED_BASE=33
ADAP_SEQ_TAIL 	 = "GTTTTAGAGCTAGAAATAGC"
ADAP_SEQ_HEAD	 = "GAAACTCTGGGAGCTGCGATTGGCAG"

BARCODE_SEQ_FILE 	= "data/barcodes/myBarcodes.fa"
BARCODE_SEQ_FILE1 	= "data/barcodes/myBarcodes.anchored.fa"
BARCODE_MAP	 	= "data/barcodes/barcode_map.tsv"
SGRNA_INDEX	 	= "data/library/crspri_grna-index";

MIN_LENGTH	 = 19;

# target sequences
COMMON_PREFIX		= "AllLanes_"
COMMON_SUFFIX_R1	= "_R1_merged"
COMMON_SUFFIX_R2 	= "_R2_merged"
COMMON_SUFFIX 		= ".fastq.gz" 

in_prefixes = ["Idx3_TTAGGCAT", "Idx1_ATCACGAT"]
out_prefixes = [];

if not os.path.exists(SAVEDIR):
    os.makedirs(SAVEDIR)
if not os.path.exists(TEMPDIR):
    os.makedirs(TEMPDIR)
if not os.path.exists(TRIMDIR):
    os.makedirs(TRIMDIR)
if not os.path.exists(CNTDIR):
    os.makedirs(CNTDIR)

SAMPLE_MAP = {};
with open(BARCODE_MAP,'r') as tsvin:
    tsvin = csv.reader(tsvin, delimiter='\t')
    for row in tsvin:
     illumIndx = row[0]
     myIndx = row[1]
     popLabel = row[2]
     print illumIndx + "\t" + myIndx + "\t" + popLabel;
     SAMPLE_MAP[(illumIndx, myIndx)] = popLabel;      

fasta_sequences = SeqIO.parse(open(BARCODE_SEQ_FILE),'fasta');
myIndex_dict={}
for fasta in fasta_sequences:
   name, sequence = fasta.id, fasta.seq.tostring()
   myIndex_dict[name] = sequence;



#################################################################################################################################################################################################
#################################################################################################################################################################################################

def get_prefix(reNAME, infilename):
  # return the prfix of the read information in a given filename
  mymatch = reNAME.match(infilename)
  if mymatch is None:
    return None
  else:
    return mymatch.groups()[0]

#################################################################################################################################################################################################

def return_input_file(inprefix):
  return os.path.join(READDIR,COMMON_PREFIX + inprefix + COMMON_SUFFIX_R1 + COMMON_SUFFIX);

def return_temp_file(inprefix):
  return os.path.join(TEMPDIR,COMMON_PREFIX + inprefix + COMMON_SUFFIX_R1)

def return_copy_file(sample_name, illumIdx, myIdx):
  return os.path.join(TRIMDIR, sample_name + "_illumIdx_" + illumIdx + "_myIdx_" + myIdx + ".fully_trimmed" + COMMON_SUFFIX);

def return_count_file(sample_name, illumIdx, myIdx):
  return os.path.join(CNTDIR, sample_name + "_illumIdx_" + illumIdx + "_myIdx_" + myIdx + ".cnt");


#################################################################################################################################################################################################

def preprocess_gz_file(inprefix):
  # do some initial preprocessing of a gz file, including trimming and quality score filtering
  # give the sample-specific identifier and output prefix that should be used

  infile_fwd 	= return_input_file(inprefix);
  tmpfile_fwd 	= return_temp_file(inprefix);

  #  do some quality trimming and write a processed file
  #  cutadapt -m 20 -q 20 -n 2 --match-read-wildcards -a TGAAGAGCGAGCGGATACAG  -o out.cut5padapt.All_Lanes_Idx1_CGATGTAT_R1_merged.fastq.gz

  # first:  clip the 3' universal adapter sequence
  cutadapt_cmd_1 = "cutadapt --quality-base=%i -q 12  -a %s -n 1 -o %s.cutadap_3p.fastq.gz %s > %s.cutadap_3p.log 2> %s.cutadap_3p.err" % (PHRED_BASE, ADAP_SEQ_TAIL, tmpfile_fwd, infile_fwd, tmpfile_fwd, tmpfile_fwd);
  # second: clip the 5' custom index sequence and demultiplex
  cutadapt_cmd_2 = "cutadapt --quality-base=%i -q 20 -g file:%s -O 3 --no-indels -e 0.1 --no-trim --untrimmed-o %s.myBarcodeUntrimmed.fastq.gz -o %s.trimmed-{name}.myBarcode.fastq.gz %s.cutadap_3p.fastq.gz  > %s.cutadap_myI.log 2> %s.cutadap_myI.err" % (PHRED_BASE, BARCODE_SEQ_FILE1, tmpfile_fwd, tmpfile_fwd, tmpfile_fwd, tmpfile_fwd, tmpfile_fwd);
  # third:  clip the 5' universal adapter sequence
  cutadapt_cmd_3 = "ls %s.trimmed-*.myBarcode.fastq.gz | while read line; do (nline=${line/myBarcode.fastq.gz/myBarcode.cut_5p.fastq.gz}; cutadapt --quality-base=%i -m %i -g %s -o $nline $line); done > %s.cutadap_5p.log 2> %s.cutadap_5p.err" % (tmpfile_fwd, PHRED_BASE, MIN_LENGTH, ADAP_SEQ_HEAD, tmpfile_fwd, tmpfile_fwd);


  # here we do the forward and reverse reads separately; for the forward read, *only* keep reads with the correct adapter
  print "running... " +  cutadapt_cmd_1;
  print;
  subprocess.call(cutadapt_cmd_1,shell=True, executable="/bin/bash")

  print "running... " +  cutadapt_cmd_2;
  print;
  subprocess.call(cutadapt_cmd_2,shell=True, executable="/bin/bash")

  print "running... " +  cutadapt_cmd_3;
  print;
  subprocess.call(cutadapt_cmd_3,shell=True, executable="/bin/bash")


  demultiplexed_files = glob.glob("%s.trimmed-*.myBarcode.cut_5p.fastq.gz" % tmpfile_fwd)
  print "%s.trimmed-*.myBarcode.cut_5p.fastq.gz" % tmpfile_fwd;
  print;

  fname_re_myIndx    = re.compile(r"%s.trimmed-(.*).myBarcode.cut_5p.fastq.gz" % tmpfile_fwd)  

  tmpBASE = os.path.join(TEMPDIR, COMMON_PREFIX)
  fname_re_illumIndx = re.compile(r"%sIdx\d_(.*)%s.trimmed-(.*).myBarcode.cut_5p.fastq.gz" % (tmpBASE, COMMON_SUFFIX_R1) )  
  print "%sIdx\d_(.*)%s.trimmed-(.*).myBarcode.cut_5p.fastq.gz" % (tmpBASE, COMMON_SUFFIX_R1);
  print;

  prefixes = set()
  for filename in demultiplexed_files: 
   myIndx    = myIndex_dict[("%s"   % get_prefix(fname_re_myIndx, filename))]
   illumIndx = ("%s"   % get_prefix(fname_re_illumIndx, filename))
   copyfile = return_copy_file(SAMPLE_MAP[(illumIndx, myIndx)], illumIndx, myIndx);
   print filename+ "\tx" + myIndx + "x\tx" + illumIndx + "x\t That's that\t" + copyfile;   
   mv_cmd = "cp %s %s" %(filename, copyfile);
   subprocess.call(mv_cmd,shell=True)

  print "\nDone with preprocessing iteration!";
  print;

#################################################################################################################################################################################################

def run_bowtie(trimmedfile):
  rp_suffix = ".fully_trimmed" + COMMON_SUFFIX;
  samoutput = string.replace(trimmedfile, rp_suffix, ".sam");
  logoutput = string.replace(trimmedfile, rp_suffix, ".sam.log");
  erroutput = string.replace(trimmedfile, rp_suffix, ".sam.err");

  cmdline = 'bowtie2 --very-sensitive-local  -N 1 --norc -x %s -U %s -S %s > %s 2> %s' % (SGRNA_INDEX,  trimmedfile, samoutput, logoutput, erroutput);
  print cmdline
  print;
  subprocess.call(cmdline, shell=True)
  print cmdline



def run_postprocess_sam(trimmedfile):
  rp_suffix  = ".fully_trimmed" + COMMON_SUFFIX;
  samoutput  = string.replace(trimmedfile, rp_suffix, ".sam");
  bamoutput  = string.replace(trimmedfile, rp_suffix, ".bam");
  srtoutput  = string.replace(trimmedfile, rp_suffix, ".sorted");
  srtoutput2 = string.replace(trimmedfile, rp_suffix, ".sorted.bam");
  bedoutput  = string.replace(trimmedfile, rp_suffix, ".bed");
  cntoutput  = string.replace(trimmedfile, rp_suffix, ".cnt");

  cmd_1 = 'samtools view -F 0x04 -Sb %s > %s' % (samoutput, bamoutput);
  cmd_2 = 'samtools sort %s %s' % (bamoutput, srtoutput);
  cmd_3 = 'bamToBed -i %s > %s' % (srtoutput2, bedoutput);
  cmd_4 = 'cut -f 1 %s | sort | uniq -c | awk \'{print $2 \"\\t\" $1}\' | sort -k2nr > %s'  % (bedoutput, cntoutput);


  print cmd_1
  print;
  subprocess.call(cmd_1, shell=True)
  print cmd_2
  print;
  subprocess.call(cmd_2, shell=True)
  print cmd_3
  print;
  subprocess.call(cmd_3, shell=True)
  print cmd_4
  print;
  subprocess.call(cmd_4, shell=True)

  cp_cmd = "cp %s %s" %(cntoutput, CNTDIR);
  print cp_cmd
  print; 
  subprocess.call(cp_cmd,shell=True)



#################################################################################################################################################################################################

def read_counts_file(cntfile):
  utrToCnt = {};
  with open(cntfile,'r') as tsvin:
    tsvin = csv.reader(tsvin, delimiter='\t')
    #headers = tsvin.next()
    for row in tsvin:
      utr = row[0];
      cnt = int(row[1]);
      utrToCnt[utr] = cnt;
  return utrToCnt;

#################################################################################################################################################################################################
#################################################################################################################################################################################################












if (__name__ == "__main__"):

  ##############################
  ## preprocess each data set
  ##############################


  if(DO_TRIM_ADAPTERS):
   pool1=Pool(processes=24)
   for i in range(len(in_prefixes)):
     inpr = in_prefixes[i]
     print "beginning %s" % (inpr)
     apply(preprocess_gz_file, [inpr])
   pool1.close()
   pool1.join()

  files_for_alignemnt = set();
  for key, sample_name in SAMPLE_MAP.items():
   illumIndx = key[0];
   myIndx    = key[1];
   fcp = return_copy_file(sample_name, illumIndx, myIndx);
   if os.path.exists(fcp):
    files_for_alignemnt.add(fcp);
    print fcp;
 
  for afile in files_for_alignemnt:
   if(DO_ALIGN_SEQUENCES):
    print "\n\nDO_ALIGN_SEQUENCES\n\n";
    run_bowtie(afile);
   if(DO_PROCESS_SAMFILES):
    print "\n\nDO_PROCESS_SAMFILES\n\n";
    run_postprocess_sam(afile);

  if(DO_MERGE_FILES):
   files_for_counts = set();
   utrs_in_assay    = set();
   counts_dict = {};
   print "\n\nDO_MERGE_FILES\n\n";
   for key, sample_name in SAMPLE_MAP.items():
    illumIndx = key[0];
    myIndx    = key[1];
    sample    = SAMPLE_MAP[(illumIndx, myIndx)];

    fcp = return_count_file(sample_name, illumIndx, myIndx);
    if os.path.exists(fcp):
     files_for_counts.add(fcp);
     c_dict = {};
     print "read file\t" +fcp;
     c_dict = read_counts_file(fcp);
     print "process sample\t" + sample;
     counts_dict[sample] = c_dict; 
     for utr, cnts in c_dict.items():
       utrs_in_assay.add(utr);

   writefile  =  'all_counts.tsv';
   writefile2 =  'all_count_frequencies.tsv';
   sum_c_dict = {};
   sum_r_dict = {};
   with open(writefile,  'w') as outfile:
    with open(writefile2, 'w') as outfile2:
     outfile.write("gRNAname");
     ckeys = counts_dict.keys();
     ckeys.sort();
     for sample_key in ckeys:
      outfile.write("\t" + sample_key);
      sum_c_dict[sample_key] = sum(counts_dict[sample_key].values())/1000000.0
      print sum_c_dict[sample_key];
     outfile.write("\n");
     for utr in utrs_in_assay:
      outfile.write(utr);
      outfile2.write(utr);
      for sample_key in ckeys:
       if utr in counts_dict[sample_key]:
        outfile.write("\t" + np.str(counts_dict[sample_key][utr]));
        outfile2.write("\t" + np.str('{0:.3f}'.format(np.log10(counts_dict[sample_key][utr]/sum_c_dict[sample_key]))));
       else:
        outfile.write("\t0");
        outfile2.write("\t-2.00");
      outfile.write("\n");
      outfile2.write("\n");
    outfile2.close(); 
   outfile.close(); 

  exit();


