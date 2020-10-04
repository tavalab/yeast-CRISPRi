##      Panos Oikonomou, Columbia University
##      Saeed Tavazoie Lab
##
## 	(c) 2020 The Trustees of Columbia University in the City of New York. 
##	This work may be reproduced and distributed for academic 
##	non-commercial purposes only without further authorization, 
##	but rightsholder otherwise reserves all rights.
##
##	
##	requirements: python 2.7
##	

##	This script reads the matrix containing read counts for every sample, normalizes
##	deplection scores are calculated between sample pairs listed in the dictionary sampleComparisons
##	calculates the depletion scores for every gene from its corresponding gRNAs
##	calculates the values for the Synthetic scrambled genes
##	saves all the values for both the genes and the Synthetic scrambled genes in dir_SAVE
##	if DO_PSEUDOCORRECTION=1, it calculates z-scores for the genes based on the distribution of the Synthetic scrambled genes
##	average depletion scores are saved in Data_Enrichment_Score.txt

import pandas as pd
import string
import numpy as np
import re
import scipy.stats as ss
import random
import operator


############################################################################################################################################

THR			= np.log2(1e-5);	# threshold; only consider sgRNAs that have reads above this threshold in either 
						# the sample itself or the control soample it is compared to (see sampleComparisons below)
NR_RANDOMIZATIONS	= 200;			# generate NR_RANDOMIZATIONS*NR_OF_GENES randomly shuffled control genes
DO_PSEUDOCORRECTION	= 1;			# 1 for z-scoring to the distribution of the Synthetic scrambled genes; 0 for not
dir_SAVE		= "data/out/example/" 	# location to save output
sampleComparisons 	= {'CRISPR-SMPL-1':'CRISPR-CNTRL-1', 'CRISPR-SMPL-2':'CRISPR-CNTRL-2', 'CRISPR-SMPL-3':'CRISPR-CNTRL-3'};
count_FILE		= 'all_counts.example.1.tsv'	# this is the output of script process_01_Sequencing_yeast_CRISPRi.py

############################################################################################################################################

def gRNA_measure(gRNA_list):
	if len(gRNA_list)<4:
		return np.mean(gRNA_list);
	else:
		top = np.mean(sorted(gRNA_list)[-3:]);
		bot = np.mean(sorted(gRNA_list)[:3]);
		if np.abs(top)>np.abs(bot):
			return top;
		else:
			return bot;

def write_data_file(dict, writefile):
	with open(writefile, 'w') as outfile:
		outfile.write("ID\tValue\n");
		for gn in sorted(dict.items(), key=operator.itemgetter(1), reverse=True):
			str =  "%s\t%f\n" %(gn[0], dict[gn[0]])
			outfile.write(str);
		outfile.close();

def write_list_file(llist, writefile):
	with open(writefile, 'w') as outfile:
		for nr in llist:
			str =  "%f\n" %nr
			outfile.write(str);
		outfile.close();

############################################################################################################################################

##	1. Read the count matrix
dtcCounts	= pd.read_csv(count_FILE, header=0, index_col=0, delimiter="\t");
dtcCounts	= dtcCounts.fillna(0)

##	2. Calculate the total reads per sample
print("Calculate total reads per sample!");
print(dtcCounts.shape)
USE_KEYS = dtcCounts.keys();
samplesToTotals = dtcCounts[USE_KEYS].sum();
print(samplesToTotals);

##	3. read the sgRNA properties
seqFile		= 'data/library/crspri_gRNA_properties.tsv';
dtc		= pd.read_csv(seqFile, header=0, index_col=None, delimiter="\t");
dtc		= dtc.loc[[ix for ix in dtc.index if (dtc['new index'].loc[ix] in dtcCounts.index)]]
print dtc

##	4. store both counts and properties on the same dataframe
for col in dtcCounts.columns:
	dtc[col] = dtc['new index'].apply(lambda x: dtcCounts[col].loc[x]);

##	5. Calculate log normalized frequencies
print("Calculate Log Normalized Frequencies");
for sample in USE_KEYS:
	dtc[sample + '-log2'] = np.log2(dtc[sample]+1) - np.log2(samplesToTotals[sample])
dtc['ind']  = dtc['new index']

##	6. Calculate enrichments/deletions for sample comparisons
geneToValueList = {};
for sample2 in sampleComparisons:
	sample1 = sampleComparisons[sample2];
	dtc[sample2 + '-' + sample1] = dtc[sample2 + '-log2'] - dtc[sample1 + '-log2']
	cond  = ((dtc[sample2 + '-log2']>THR) | (dtc[sample1 + '-log2']>THR))
	tagR  = sample2 + '-' + sample1;
	geneToValueList[tagR] = {k:g[tagR].values for k,g in dtc[cond][['gene', tagR]].groupby('gene')}
depletionTags	 = [tag for tag in geneToValueList];
geneToValueLists = [geneToValueList[tag] for tag in geneToValueList]

##	7. Report stats of randomly shuffled sgRNAs
kk = 0;
all_genes = set();
for geneToValueList in geneToValueLists:
	print("Replicate %i: Number of random gRNAs, Mean, Std: %i\t%f\t%f" % ((kk+1), geneToValueLists[kk]['randomgrnas'].shape[0], np.mean(geneToValueLists[kk]['randomgrnas']), np.std(geneToValueLists[kk]['randomgrnas'])))
	all_genes =  all_genes.union(set(geneToValueLists[kk].keys()))
	kk = kk +1;

## 	8. Calculate values for Synthetic scrambled genes
rand_array   = []
rand_arrays = []
lens = []
for geneToValueList in geneToValueLists:
	rand_arrays.append([])
	lens.append(len(geneToValueList['randomgrnas']))

for i in range(1, NR_RANDOMIZATIONS):
	print i
	for g in (all_genes):
		cc = 0;
		avg = 0;
		kk = 0;
		for geneToValueList in geneToValueLists:
			if g in geneToValueList:
				smp = random.sample(range(0, lens[kk]), len(geneToValueList[g]))
				gv = gRNA_measure(geneToValueList['randomgrnas'][smp]);
				avg += gv;
				rand_arrays[kk].append(gv)
				cc = cc +1
			kk = kk+1
		if cc>0:
			avg = avg/cc;
			rand_array.append(avg)

if (DO_PSEUDOCORRECTION):
	mmx  = np.mean(rand_array)
	ssx  = np.std(rand_array) 
	mmxs = []
	ssxs = []
	for rand_arrayX in rand_arrays:
		mmxs.append(np.mean(rand_arrayX))
		ssxs.append(np.std(rand_arrayX))
else:
	mmx   = 0;
	ssx   = 1;
	mmxs = []
	ssxs = []
	for rand_arrayX in rand_arrays:
		mmxs.append(0)
		ssxs.append(1)

for ii in range(len(rand_array)):
	rand_array[ii] = (rand_array[ii]-mmx)/ssx
kk = 0;
for rand_arrayX in rand_arrays:
	for ii in range(len(rand_arrayX)):
		rand_arrayX[ii] = (rand_arrayX[ii]-mmxs[kk])/ssxs[kk]
	kk = kk+1

kk = 0;
for rand_arrayX in rand_arrays:
	print("Synthetic scrambled genes Rep%i: Mean, Std, 5th Perc, 10th Perc, length" % (kk+1))
	print(np.percentile(rand_arrayX, 95))
	print(np.percentile(rand_arrayX, 90))
	print("%.4f\t%.4f" %(np.mean(rand_arrayX), np.std(rand_arrayX)))
	print(np.percentile(rand_arrayX, 10))
	print(np.percentile(rand_arrayX, 5))
	print(len(rand_arrayX))
	print("\n")
	kk = kk+1;


kk = 0;
for rand_arrayX in rand_arrays:
	write_list_file(rand_arrayX, dir_SAVE + "/Synthetic_Scrambled_Genes_values.PseudoCorrection_YN_%i.Rep%i.txt" % (DO_PSEUDOCORRECTION, kk))
	kk = kk+1;


## 	9. Map gRNAs to genes and calculate enrichment values at the gene level
all_dct     = {}
all_dcts    = []
all_array   = []
all_arrays  = []
for geneToValueList in geneToValueLists:
	all_dcts.append({})
	all_arrays.append([])

for g in all_genes:
	if g != 'randomgrnas':
		cc = 0;
		avg = 0;
		kk = 0;
		for geneToValueList in geneToValueLists:
			if g in geneToValueList:
				gv = (gRNA_measure(geneToValueList[g])-mmxs[kk])/ssxs[kk];
				avg += gv
				all_arrays[kk].append(gv)
				all_dcts[kk][g] = gv;
				cc = cc +1
			kk = kk+1;
		if cc>0:
			if (DO_PSEUDOCORRECTION):
				avg = avg/np.sqrt(cc);
			else:
				avg = avg/cc;
			all_array.append(avg)
			all_dct[g] = avg;


DATA_LIST = []
kk = 0;
for geneToValueList in geneToValueLists:
	print("Yeast Genes Rep%i: Mean, Std, 5th Perc, 10th Perc, length" % (kk+1))
	print(np.percentile(all_arrays[kk], 95))
	print(np.percentile(all_arrays[kk], 90))
	print("%.4f\t%.4f" %(np.mean(all_arrays[kk]), np.std(all_arrays[kk])))
	print(np.percentile(all_arrays[kk], 10))
	print(np.percentile(all_arrays[kk], 5))
	print(len(all_arrays[kk]))
	print("\n")
	kk = kk+1;
	dtc['Gene Value Rep %i' % kk] = dtc['gene'].map(all_dcts[kk-1])
	DATA_LIST.append('Gene Value Rep %i' % kk)


##	10. Save results
write_file = dir_SAVE + "/Data_Enrichment_Score.txt"
write_data_file(all_dct, write_file);
kk = 0;
for geneToValueList in geneToValueLists:
	write_file = dir_SAVE + "/Data_Enrichment_Score.Rep_%i.DO_PSEUDOCORRECTION_YN_%i.Avg.txt" % (kk+1, DO_PSEUDOCORRECTION)
	write_data_file(all_dcts[kk], write_file);
	kk = kk +1;

