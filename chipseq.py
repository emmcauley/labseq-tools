import sys
import warnings
import xlrd
import csv
import pandas as pd
import numpy as np
from pprint import PrettyPrinter

if len(sys.argv) < 4:
	sys.stderr.write("Usage: %s <filename> <cell-line#1> <cell-line#2> [gene-name ...]\n" % (sys.argv[0]))
	sys.exit(1)

#filename from Tzuni, cell line 1 and cell line 2
filename = sys.argv[1]
cl1 = sys.argv[2]
cl2 = sys.argv[3]

pp=PrettyPrinter()

#turn off warnings
warnings.simplefilter(action="ignore",category=UserWarning)

#read in CSV file to a dataframe, based on given file name and 2 cell lines you are comparing
nvs = pd.read_csv(sys.argv[1], sep='\t')
prostate = pd.read_csv('CWRR1-PrEC70-same.tsv', sep='\t')

prostateset = set(prostate.gene_name)
#print(prostateset)
stem = pd.read_csv('NCCIT-WA01-same.tsv', sep='\t')
stemset = set(stem.gene_name)
#dictionary of results
precFrames = {}

#sort for protein-coding genes uniquely bound by Sox2 in between 2 cell lines.
precFrames['different'] = nvs[nvs.gene_biotype=='protein_coding'][(nvs['Sox2_'+cl1+'Sox2_unique']>=1) | (nvs['Sox2_'+cl2+'Sox2_unique']>=1) | (nvs['Sox2_'+cl1+'Sox2_close']>=1) | (nvs['Sox2_'+cl2+'Sox2_close']>=1)]

#sort for protein-coding genes that are bound by Sox2 in the same position. genes must have >= K4 marks. K27 coverage fraction is not currently addressed. 
precFrames['same'] = nvs[nvs.gene_biotype=='protein_coding'][(nvs['Sox2_'+cl1+'Sox2_same']>=1) & (nvs['Sox2_'+cl2+'Sox2_same']>=1)]


#for every gene that met the filtering criteria, sort them based on <...>, tells you which files are written.
for key in precFrames:
	outfile = "%s-%s-%s.tsv"%(cl1,cl2,key)
	precFrames[key].sort_values(['Expr_PPDE', cl1+'K4', cl2+'K4'], ascending=[False, False, True]).to_csv(outfile, sep='\t')
	#precFrames[key].sort_values(['Expr_' + cl1 + '_Mean', 'Expr_' + cl2 + '_Mean'], ascending=[False,True]).to_csv(outfile, sep='\t')

	sys.stderr.write("%s results written to file %s\n" % (key,outfile) )
sys.stderr.write("\n")	

#if you passed the code a list of specific genes to cross-reference quickly --> loop through all result files, if the code matches a substring of argv3 then store in a dict of results, which then prints out to user

for gene in  sys.argv[4:]:
	results = {}
	for key in precFrames:
		results[key] = precFrames[key][precFrames[key].gene_name.str.contains(gene)]

		print("-- %s results  for %s --" % (key,gene))
		if (len(results[key]) ==0 ):
			print("NONE")
		else:
			print(results[key][ ["gene_name", "Expr_PostrFC","Expr_PPDE", cl1+'K4', cl2+'K4'] ].to_string(index=False))
