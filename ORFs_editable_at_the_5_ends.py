# -*- coding: utf-8 -*-
import sys
import re
import subprocess
from Bio import SeqIO
import numpy as np

argvs = sys.argv
argc = len(argvs)

if argc < 2:
	print("Usage: python3 {} [reference_seq.fasta]".format(argvs[0]))
else:
	reference_seq_name = argvs[1]

	ORF_list = []
	for record in SeqIO.parse(reference_seq_name, 'fasta'):
		ORF_name = record.id
		ORF_1000_seq = record.seq

		SpORF = np.zeros(len(ORF_1000_seq))
		SaORF = np.zeros(len(ORF_1000_seq))
		AsORF = np.zeros(len(ORF_1000_seq))

		Sp_PAM = 0
		Sa_PAM = 0
		As_PAM = 0

		ORF_data = []
		for i in range(len(ORF_1000_seq)):
			#SpCas9
			if i+1 < len(ORF_1000_seq):
				#SpCas9 Forward PAM
				if ORF_1000_seq[i]=='G' and ORF_1000_seq[i+1]=='G':
					Sp_PAM = Sp_PAM + 1
					if i-12 >= 0:
						SpORF[i-12] = True
					if i-11 >= 0:
						SpORF[i-11] = True
					if i-10 >= 0:
						SpORF[i-10] = True
					if i-9 >= 0:
						SpORF[i-9] = True
					if i-8 >= 0:
						SpORF[i-8] = True
					if i-7 >= 0:
						SpORF[i-7] = True
					if i-6 >= 0:
						SpORF[i-6] = True
					if i-5 >= 0:
						SpORF[i-5] = True
					if i-4 >= 0:
						SpORF[i-4] = True
					if i-3 >= 0:
						SpORF[i-3] = True
					if i-2 >= 0:
						SpORF[i-2] = True
					#if i-1 >= 0:
					#	SpORF[i-1] = True
					if i >= 0:
						SpORF[i] = True
					if i+1 < len(ORF_1000_seq):
						SpORF[i+1] = True
				#SpCas9 Reverse PAM
				if ORF_1000_seq[i]=='C' and ORF_1000_seq[i+1]=='C':
					Sp_PAM = Sp_PAM + 1
					if i < len(ORF_1000_seq):
						SpORF[i] = True
					if i+1 < len(ORF_1000_seq):
						SpORF[i+1] = True
					if i+2 < len(ORF_1000_seq):
						SpORF[i+2] = True
					if i+3 < len(ORF_1000_seq):
						SpORF[i+3] = True
					if i+4 < len(ORF_1000_seq):
						SpORF[i+4] = True
					if i+5 < len(ORF_1000_seq):
						SpORF[i+5] = True
					if i+6 < len(ORF_1000_seq):
						SpORF[i+6] = True
					if i+7 < len(ORF_1000_seq):
						SpORF[i+7] = True
					if i+8 < len(ORF_1000_seq):
						SpORF[i+8] = True
					if i+9 < len(ORF_1000_seq):
						SpORF[i+9] = True
					if i+10 < len(ORF_1000_seq):
						SpORF[i+10] = True
					if i+11 < len(ORF_1000_seq):
						SpORF[i+11] = True
					if i+12 < len(ORF_1000_seq):
						SpORF[i+12] = True
					if i+13 < len(ORF_1000_seq):
						SpORF[i+13] = True
			#SaCas9
			if i+3 < len(ORF_1000_seq):
				#SaCas9 Forward PAM
				if ORF_1000_seq[i]=='G' and (ORF_1000_seq[i+1]=='A' or ORF_1000_seq[i+1]=='G') and (ORF_1000_seq[i+2]=='A' or ORF_1000_seq[i+2]=='G') and ORF_1000_seq[i+3]=='T':
					Sa_PAM = Sa_PAM + 1
					if i-13 >= 0:
						SaORF[i-13] = True
					if i-12 >= 0:
						SaORF[i-12] = True
					if i-11 >= 0:
						SaORF[i-11] = True
					if i-10 >= 0:
						SaORF[i-10] = True
					if i-9 >= 0:
						SaORF[i-9] = True
					if i-8 >= 0:
						SaORF[i-8] = True
					if i-7 >= 0:
						SaORF[i-7] = True
					if i-6 >= 0:
						SaORF[i-6] = True
					if i-5 >= 0:
						SaORF[i-5] = True
					if i-4 >= 0:
						SaORF[i-4] = True
					if i-3 >= 0:
						SaORF[i-3] = True
					#if i-2 >= 0:
					#	SaORF[i-2] = True
					#if i-1 >= 0:
					#	SaORF[i-1] = True
					if i >= 0:
						SaORF[i] = True
					#if i+1 < len(ORF_1000_seq):
					#	SaORF[i+1] = True
					#if i+2 < len(ORF_1000_seq):
					#	SaORF[i+2] = True
					if i+3 < len(ORF_1000_seq):
						SaORF[i+3] = True

				#SaCas9 Reverse PAM
				if ORF_1000_seq[i]=='A' and (ORF_1000_seq[i+1]=='C' or ORF_1000_seq[i+1]=='T') and (ORF_1000_seq[i+2]=='C' or ORF_1000_seq[i+2]=='T') and ORF_1000_seq[i+3]=='C':
					Sa_PAM = Sa_PAM + 1
					if i < len(ORF_1000_seq):
						SaORF[i] = True
					#if i+1 < len(ORF_1000_seq):
					#	SaORF[i+1] = True
					#if i+2 < len(ORF_1000_seq):
					#	SaORF[i+2] = True
					if i+3 < len(ORF_1000_seq):
						SaORF[i+3] = True
					#if i+4 < len(ORF_1000_seq):
					#	SaORF[i+4] = True
					#if i+5 < len(ORF_1000_seq):
					#	SaORF[i+5] = True
					if i+6 < len(ORF_1000_seq):
						SaORF[i+6] = True
					if i+7 < len(ORF_1000_seq):
						SaORF[i+7] = True
					if i+8 < len(ORF_1000_seq):
						SaORF[i+8] = True
					if i+9 < len(ORF_1000_seq):
						SaORF[i+9] = True
					if i+10 < len(ORF_1000_seq):
						SaORF[i+10] = True
					if i+11 < len(ORF_1000_seq):
						SaORF[i+11] = True
					if i+12 < len(ORF_1000_seq):
						SaORF[i+12] = True
					if i+13 < len(ORF_1000_seq):
						SaORF[i+13] = True
					if i+14 < len(ORF_1000_seq):
						SaORF[i+14] = True
					if i+15 < len(ORF_1000_seq):
						SaORF[i+15] = True
					if i+16 < len(ORF_1000_seq):
						SaORF[i+16] = True
			#AsCas12a
			if i+3 < len(ORF_1000_seq):
				#AsCas12 Forward PAM
				if ORF_1000_seq[i]=='T' and ORF_1000_seq[i+1]=='T' and ORF_1000_seq[i+2]=='T' and (ORF_1000_seq[i+3]=='A' or ORF_1000_seq[i+3]=='G' or ORF_1000_seq[i+3]=='C'):
					As_PAM = As_PAM + 1
					if i < len(ORF_1000_seq):
						AsORF[i] = True
					if i+1 < len(ORF_1000_seq):
						AsORF[i+1] = True
					if i+2 < len(ORF_1000_seq):
						AsORF[i+2] = True
					if i+4 < len(ORF_1000_seq):
						AsORF[i+4] = True
					if i+5 < len(ORF_1000_seq):
						AsORF[i+5] = True
					if i+6 < len(ORF_1000_seq):
						AsORF[i+6] = True
					if i+7 < len(ORF_1000_seq):
						AsORF[i+7] = True
					if i+8 < len(ORF_1000_seq):
						AsORF[i+8] = True
					if i+9 < len(ORF_1000_seq):
						AsORF[i+9] = True
					if i+10 < len(ORF_1000_seq):
						AsORF[i+10] = True
					if i+11 < len(ORF_1000_seq):
						AsORF[i+11] = True
					if i+12 < len(ORF_1000_seq):
						AsORF[i+12] = True
					if i+13 < len(ORF_1000_seq):
						AsORF[i+13] = True
					if i+14 < len(ORF_1000_seq):
						AsORF[i+14] = True
					if i+15 < len(ORF_1000_seq):
						AsORF[i+15] = True
					if i+16 < len(ORF_1000_seq):
						AsORF[i+16] = True
					if i+17 < len(ORF_1000_seq):
						AsORF[i+17] = True
					if i+18 < len(ORF_1000_seq):
						AsORF[i+18] = True
					if i+19 < len(ORF_1000_seq):
						AsORF[i+19] = True
					if i+20 < len(ORF_1000_seq):
						AsORF[i+20] = True

				#AsCas12a Revwerse PAM
				if (ORF_1000_seq[i]=='A' or ORF_1000_seq[i]=='G' or ORF_1000_seq[i]=='C') and ORF_1000_seq[i+1]=='A' and ORF_1000_seq[i+2]=='A' and ORF_1000_seq[i+3]=='A':
					As_PAM = As_PAM + 1
					if i-17 >= 0:
						AsORF[i-17] = True
					if i-16 >= 0:
						AsORF[i-16] = True
					if i-15 >= 0:
						AsORF[i-15] = True
					if i-14 >= 0:
						AsORF[i-14] = True
					if i-13 >= 0:
						AsORF[i-13] = True
					if i-12 >= 0:
						AsORF[i-12] = True
					if i-11 >= 0:
						AsORF[i-11] = True
					if i-10 >= 0:
						AsORF[i-10] = True
					if i-9 >= 0:
						AsORF[i-9] = True
					if i-8 >= 0:
						AsORF[i-8] = True
					if i-7 >= 0:
						AsORF[i-7] = True
					if i-6 >= 0:
						AsORF[i-6] = True
					if i-5 >= 0:
						AsORF[i-5] = True
					if i-4 >= 0:
						AsORF[i-4] = True
					if i-3 >= 0:
						AsORF[i-3] = True
					if i-2 >= 0:
						AsORF[i-2] = True
					if i-1 >= 0:
						AsORF[i-1] = True
					if i+1 < len(ORF_1000_seq):
						AsORF[i+1] = True
					if i+2 < len(ORF_1000_seq):
						AsORF[i+2] = True
					if i+3 < len(ORF_1000_seq):
						AsORF[i+3] = True

		SpNterm = (SpORF[1000] or SpORF[1001] or SpORF[1002])
		SaNterm = (SaORF[1000] or SaORF[1001] or SaORF[1002])
		AsNterm = (AsORF[1000] or AsORF[1001] or AsORF[1002])

		ORF_data.append(ORF_name)
		ORF_data.append(SpNterm)
		ORF_data.append(SaNterm)
		ORF_data.append(AsNterm)

		line = ''
		for item in ORF_data:
			line = line + '\t' +str(item)
		print(line)
