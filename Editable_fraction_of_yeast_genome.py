# -*- coding: utf-8 -*-
import sys
import re
import subprocess
from Bio import SeqIO
import numpy as np

argvs = sys.argv
argc = len(argvs)

if argc < 2:
	print("Usage: python3 {} [reference_genome_seq.fasta]".format(argvs[0]))
else:
	reference_genome_seq_name = argvs[1]

	chrI_len = 0
	chrII_len = 0
	chrIII_len = 0
	chrIV_len = 0
	chrV_len = 0
	chrVI_len = 0
	chrVII_len = 0
	chrVIII_len = 0
	chrIX_len = 0
	chrX_len = 0
	chrXI_len = 0
	chrXII_len = 0
	chrXIII_len = 0
	chrXIV_len = 0
	chrXV_len = 0
	chrXVI_len = 0
	chrmt_len = 0

	for record in SeqIO.parse(reference_genome_seq_name, 'fasta'):
		if record.id == "chrI":
			chrI_len = len(record.seq)
		elif record.id == "chrII":
			chrII_len = len(record.seq)
		elif record.id == "chrIII":
			chrIII_len = len(record.seq)
		elif record.id == "chrIV":
			chrIV_len = len(record.seq)
		elif record.id == "chrV":
			chrV_len = len(record.seq)
		elif record.id == "chrVI":
			chrVI_len = len(record.seq)
		elif record.id == "chrVII":
			chrVII_len = len(record.seq)
		elif record.id == "chrVIII":
			chrVIII_len = len(record.seq)
		elif record.id == "chrIX":
			chrIX_len = len(record.seq)
		elif record.id == "chrX":
			chrX_len = len(record.seq)
		elif record.id == "chrXI":
			chrXI_len = len(record.seq)
		elif record.id == "chrXII":
			chrXII_len = len(record.seq)
		elif record.id == "chrXIII":
			chrXIII_len = len(record.seq)
		elif record.id == "chrXIV":
			chrXIV_len = len(record.seq)
		elif record.id == "chrXV":
			chrXV_len = len(record.seq)
		elif record.id == "chrXVI":
			chrXVI_len = len(record.seq)
		elif record.id == "chrmt":
			chrmt_len = len(record.seq)
	print("Reference genome sequence loaded.")

	total_genome_size = 0
	total_genome_size += chrI_len
	total_genome_size += chrII_len
	total_genome_size += chrIII_len
	total_genome_size += chrIV_len
	total_genome_size += chrV_len
	total_genome_size += chrVI_len
	total_genome_size += chrVII_len
	total_genome_size += chrVIII_len
	total_genome_size += chrIX_len
	total_genome_size += chrX_len
	total_genome_size += chrXI_len
	total_genome_size += chrXII_len
	total_genome_size += chrXIII_len
	total_genome_size += chrXIV_len
	total_genome_size += chrXV_len
	total_genome_size += chrXVI_len

	array_Sp_chrI = np.zeros(chrI_len)
	array_Sp_chrII = np.zeros(chrII_len)
	array_Sp_chrIII = np.zeros(chrIII_len)
	array_Sp_chrIV = np.zeros(chrIV_len)
	array_Sp_chrV = np.zeros(chrV_len)
	array_Sp_chrVI = np.zeros(chrVI_len)
	array_Sp_chrVII = np.zeros(chrVII_len)
	array_Sp_chrVIII = np.zeros(chrVIII_len)
	array_Sp_chrIX = np.zeros(chrIX_len)
	array_Sp_chrX = np.zeros(chrX_len)
	array_Sp_chrXI = np.zeros(chrXI_len)
	array_Sp_chrXII = np.zeros(chrXII_len)
	array_Sp_chrXIII = np.zeros(chrXIII_len)
	array_Sp_chrXIV = np.zeros(chrXIV_len)
	array_Sp_chrXV = np.zeros(chrXV_len)
	array_Sp_chrXVI = np.zeros(chrXVI_len)
	array_Sp_chrmt = np.zeros(chrmt_len)

	array_Sa_chrI = np.zeros(chrI_len)
	array_Sa_chrII = np.zeros(chrII_len)
	array_Sa_chrIII = np.zeros(chrIII_len)
	array_Sa_chrIV = np.zeros(chrIV_len)
	array_Sa_chrV = np.zeros(chrV_len)
	array_Sa_chrVI = np.zeros(chrVI_len)
	array_Sa_chrVII = np.zeros(chrVII_len)
	array_Sa_chrVIII = np.zeros(chrVIII_len)
	array_Sa_chrIX = np.zeros(chrIX_len)
	array_Sa_chrX = np.zeros(chrX_len)
	array_Sa_chrXI = np.zeros(chrXI_len)
	array_Sa_chrXII = np.zeros(chrXII_len)
	array_Sa_chrXIII = np.zeros(chrXIII_len)
	array_Sa_chrXIV = np.zeros(chrXIV_len)
	array_Sa_chrXV = np.zeros(chrXV_len)
	array_Sa_chrXVI = np.zeros(chrXVI_len)
	array_Sa_chrmt = np.zeros(chrmt_len)

	array_As_chrI = np.zeros(chrI_len)
	array_As_chrII = np.zeros(chrII_len)
	array_As_chrIII = np.zeros(chrIII_len)
	array_As_chrIV = np.zeros(chrIV_len)
	array_As_chrV = np.zeros(chrV_len)
	array_As_chrVI = np.zeros(chrVI_len)
	array_As_chrVII = np.zeros(chrVII_len)
	array_As_chrVIII = np.zeros(chrVIII_len)
	array_As_chrIX = np.zeros(chrIX_len)
	array_As_chrX = np.zeros(chrX_len)
	array_As_chrXI = np.zeros(chrXI_len)
	array_As_chrXII = np.zeros(chrXII_len)
	array_As_chrXIII = np.zeros(chrXIII_len)
	array_As_chrXIV = np.zeros(chrXIV_len)
	array_As_chrXV = np.zeros(chrXV_len)
	array_As_chrXVI = np.zeros(chrXVI_len)
	array_As_chrmt = np.zeros(chrmt_len)

	array_SpSa_chrI = np.zeros(chrI_len)
	array_SpSa_chrII = np.zeros(chrII_len)
	array_SpSa_chrIII = np.zeros(chrIII_len)
	array_SpSa_chrIV = np.zeros(chrIV_len)
	array_SpSa_chrV = np.zeros(chrV_len)
	array_SpSa_chrVI = np.zeros(chrVI_len)
	array_SpSa_chrVII = np.zeros(chrVII_len)
	array_SpSa_chrVIII = np.zeros(chrVIII_len)
	array_SpSa_chrIX = np.zeros(chrIX_len)
	array_SpSa_chrX = np.zeros(chrX_len)
	array_SpSa_chrXI = np.zeros(chrXI_len)
	array_SpSa_chrXII = np.zeros(chrXII_len)
	array_SpSa_chrXIII = np.zeros(chrXIII_len)
	array_SpSa_chrXIV = np.zeros(chrXIV_len)
	array_SpSa_chrXV = np.zeros(chrXV_len)
	array_SpSa_chrXVI = np.zeros(chrXVI_len)
	array_SpSa_chrmt = np.zeros(chrmt_len)

	array_SpAs_chrI = np.zeros(chrI_len)
	array_SpAs_chrII = np.zeros(chrII_len)
	array_SpAs_chrIII = np.zeros(chrIII_len)
	array_SpAs_chrIV = np.zeros(chrIV_len)
	array_SpAs_chrV = np.zeros(chrV_len)
	array_SpAs_chrVI = np.zeros(chrVI_len)
	array_SpAs_chrVII = np.zeros(chrVII_len)
	array_SpAs_chrVIII = np.zeros(chrVIII_len)
	array_SpAs_chrIX = np.zeros(chrIX_len)
	array_SpAs_chrX = np.zeros(chrX_len)
	array_SpAs_chrXI = np.zeros(chrXI_len)
	array_SpAs_chrXII = np.zeros(chrXII_len)
	array_SpAs_chrXIII = np.zeros(chrXIII_len)
	array_SpAs_chrXIV = np.zeros(chrXIV_len)
	array_SpAs_chrXV = np.zeros(chrXV_len)
	array_SpAs_chrXVI = np.zeros(chrXVI_len)
	array_SpAs_chrmt = np.zeros(chrmt_len)

	array_SaAs_chrI = np.zeros(chrI_len)
	array_SaAs_chrII = np.zeros(chrII_len)
	array_SaAs_chrIII = np.zeros(chrIII_len)
	array_SaAs_chrIV = np.zeros(chrIV_len)
	array_SaAs_chrV = np.zeros(chrV_len)
	array_SaAs_chrVI = np.zeros(chrVI_len)
	array_SaAs_chrVII = np.zeros(chrVII_len)
	array_SaAs_chrVIII = np.zeros(chrVIII_len)
	array_SaAs_chrIX = np.zeros(chrIX_len)
	array_SaAs_chrX = np.zeros(chrX_len)
	array_SaAs_chrXI = np.zeros(chrXI_len)
	array_SaAs_chrXII = np.zeros(chrXII_len)
	array_SaAs_chrXIII = np.zeros(chrXIII_len)
	array_SaAs_chrXIV = np.zeros(chrXIV_len)
	array_SaAs_chrXV = np.zeros(chrXV_len)
	array_SaAs_chrXVI = np.zeros(chrXVI_len)
	array_SaAs_chrmt = np.zeros(chrmt_len)

	array_SpSaAs_chrI = np.zeros(chrI_len)
	array_SpSaAs_chrII = np.zeros(chrII_len)
	array_SpSaAs_chrIII = np.zeros(chrIII_len)
	array_SpSaAs_chrIV = np.zeros(chrIV_len)
	array_SpSaAs_chrV = np.zeros(chrV_len)
	array_SpSaAs_chrVI = np.zeros(chrVI_len)
	array_SpSaAs_chrVII = np.zeros(chrVII_len)
	array_SpSaAs_chrVIII = np.zeros(chrVIII_len)
	array_SpSaAs_chrIX = np.zeros(chrIX_len)
	array_SpSaAs_chrX = np.zeros(chrX_len)
	array_SpSaAs_chrXI = np.zeros(chrXI_len)
	array_SpSaAs_chrXII = np.zeros(chrXII_len)
	array_SpSaAs_chrXIII = np.zeros(chrXIII_len)
	array_SpSaAs_chrXIV = np.zeros(chrXIV_len)
	array_SpSaAs_chrXV = np.zeros(chrXV_len)
	array_SpSaAs_chrXVI = np.zeros(chrXVI_len)
	array_SpSaAs_chrmt = np.zeros(chrmt_len)

	for record in SeqIO.parse(reference_genome_seq_name, 'fasta'):
		if record.id == "chrI":
			chr_seq = record.seq
			Sp_PAM_chrI = 0
			Sa_PAM_chrI = 0
			As_PAM_chrI = 0
			for i in range(len(chr_seq)):
				#SpCas9
				if i+1 < len(chr_seq):
					#SpCas9 Forward PAM
					if chr_seq[i]=='G' and chr_seq[i+1]=='G':
						Sp_PAM_chrI = Sp_PAM_chrI + 1
						if i-12 >=0:
							array_Sp_chrI[i-12] = True
						if i-11 >=0:
							array_Sp_chrI[i-11] = True
						if i-10 >=0:
							array_Sp_chrI[i-10] = True
						if i-9 >= 0:
							array_Sp_chrI[i-9] = True
						if i-8 >= 0:
							array_Sp_chrI[i-8] = True
						if i-7 >= 0:
							array_Sp_chrI[i-7] = True
						if i-6 >= 0:
							array_Sp_chrI[i-6] = True
						if i-5 >= 0:
							array_Sp_chrI[i-5] = True
						if i-4 >= 0:
							array_Sp_chrI[i-4] = True
						if i-3 >= 0:
							array_Sp_chrI[i-3] = True
						if i-2 >= 0:
							array_Sp_chrI[i-2] = True
						#if i-1 >= 0:
						#	array_Sp_chrI[i-1] = True
						if i >= 0:
							array_Sp_chrI[i] = True
						if i+1 < len(chr_seq):
							array_Sp_chrI[i+1] = True

					#SpCas9 Reverse PAM
					if chr_seq[i]=='C' and chr_seq[i+1]=='C':
						Sp_PAM_chrI = Sp_PAM_chrI + 1
						if i < len(chr_seq):
							array_Sp_chrI[i] = True
						if i+1 < len(chr_seq):
							array_Sp_chrI[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sp_chrI[i+2] = True
						if i+3 < len(chr_seq):
							array_Sp_chrI[i+3] = True
						if i+4 < len(chr_seq):
							array_Sp_chrI[i+4] = True
						if i+5 < len(chr_seq):
							array_Sp_chrI[i+5] = True
						if i+6 < len(chr_seq):
							array_Sp_chrI[i+6] = True
						if i+7 < len(chr_seq):
							array_Sp_chrI[i+7] = True
						if i+8 < len(chr_seq):
							array_Sp_chrI[i+8] = True
						if i+9 < len(chr_seq):
							array_Sp_chrI[i+9] = True
						if i+10 < len(chr_seq):
							array_Sp_chrI[i+10] = True
						if i+11 < len(chr_seq):
							array_Sp_chrI[i+11] = True
						if i+12 < len(chr_seq):
							array_Sp_chrI[i+12] = True
						if i+13 < len(chr_seq):
							array_Sp_chrI[i+13] = True

				#SaCas9
				if i+3 < len(chr_seq):
					#SaCas9 Forward PAM
					if chr_seq[i]=='G' and (chr_seq[i+1]=='A' or chr_seq[i+1]=='G') and (chr_seq[i+2]=='A' or chr_seq[i+2]=='G') and chr_seq[i+3]=='T':
						Sa_PAM_chrI = Sa_PAM_chrI + 1
						if i-13 >= 0:
							array_Sa_chrI[i-13] = True
						if i-12 >= 0:
							array_Sa_chrI[i-12] = True
						if i-11 >= 0:
							array_Sa_chrI[i-11] = True
						if i-10 >= 0:
							array_Sa_chrI[i-10] = True
						if i-9 >= 0:
							array_Sa_chrI[i-9] = True
						if i-8 >= 0:
							array_Sa_chrI[i-8] = True
						if i-7 >= 0:
							array_Sa_chrI[i-7] = True
						if i-6 >= 0:
							array_Sa_chrI[i-6] = True
						if i-5 >= 0:
							array_Sa_chrI[i-5] = True
						if i-4 >= 0:
							array_Sa_chrI[i-4] = True
						if i-3 >= 0:
							array_Sa_chrI[i-3] = True
						#if i-2 >= 0:
						#	array_Sa_chrI[i-2] = True
						#if i-1 >= 0:
						#	array_Sa_chrI[i-1] = True
						if i >= 0:
							array_Sa_chrI[i] = True
						#if i+1 < len(chr_seq):
						#	array_Sa_chrI[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sa_chrI[i+2] = True
						if i+3 < len(chr_seq):
							array_Sa_chrI[i+3] = True

					#SaCas9 Reverse PAM
					if chr_seq[i]=='A' and (chr_seq[i+1]=='C' or chr_seq[i+1]=='T') and (chr_seq[i+2]=='C' or chr_seq[i+2]=='T') and chr_seq[i+3]=='C':
						Sa_PAM_chrI = Sa_PAM_chrI + 1
						if i < len(chr_seq):
							array_Sa_chrI[i] = True
						#if i+1 < len(chr_seq):
						#	array_Sa_chrI[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sa_chrI[i+2] = True
						if i+3 < len(chr_seq):
							array_Sa_chrI[i+3] = True
						#if i+4 < len(chr_seq):
						#	array_Sa_chrI[i+4] = True
						#if i+5 < len(chr_seq):
						#	array_Sa_chrI[i+5] = True
						if i+6 < len(chr_seq):
							array_Sa_chrI[i+6] = True
						if i+7 < len(chr_seq):
							array_Sa_chrI[i+7] = True
						if i+8 < len(chr_seq):
							array_Sa_chrI[i+8] = True
						if i+9 < len(chr_seq):
							array_Sa_chrI[i+9] = True
						if i+10 < len(chr_seq):
							array_Sa_chrI[i+10] = True
						if i+11 < len(chr_seq):
							array_Sa_chrI[i+11] = True
						if i+12 < len(chr_seq):
							array_Sa_chrI[i+12] = True
						if i+13 < len(chr_seq):
							array_Sa_chrI[i+13] = True
						if i+14 < len(chr_seq):
							array_Sa_chrI[i+14] = True
						if i+15 < len(chr_seq):
							array_Sa_chrI[i+15] = True
						if i+16 < len(chr_seq):
							array_Sa_chrI[i+16] = True

				#AsCas12a
				if i+3 < len(chr_seq):
					#AsCas12 Forward PAM
					if chr_seq[i]=='T' and chr_seq[i+1]=='T' and chr_seq[i+2]=='T' and (chr_seq[i+3]=='A' or chr_seq[i+3]=='G' or chr_seq[i+3]=='C'):
						As_PAM_chrI = As_PAM_chrI + 1
						if i < len(chr_seq):
							array_As_chrI[i] = True
						if i+1 < len(chr_seq):
							array_As_chrI[i+1] = True
						if i+2 < len(chr_seq):
							array_As_chrI[i+2] = True
						if i+4 < len(chr_seq):
							array_As_chrI[i+4] = True
						if i+5 < len(chr_seq):
							array_As_chrI[i+5] = True
						if i+6 < len(chr_seq):
							array_As_chrI[i+6] = True
						if i+7 < len(chr_seq):
							array_As_chrI[i+7] = True
						if i+8 < len(chr_seq):
							array_As_chrI[i+8] = True
						if i+9 < len(chr_seq):
							array_As_chrI[i+9] = True
						if i+10 < len(chr_seq):
							array_As_chrI[i+10] = True
						if i+11 < len(chr_seq):
							array_As_chrI[i+11] = True
						if i+12 < len(chr_seq):
							array_As_chrI[i+12] = True
						if i+13 < len(chr_seq):
							array_As_chrI[i+13] = True
						if i+14 < len(chr_seq):
							array_As_chrI[i+14] = True
						if i+15 < len(chr_seq):
							array_As_chrI[i+15] = True
						if i+16 < len(chr_seq):
							array_As_chrI[i+16] = True
						if i+17 < len(chr_seq):
							array_As_chrI[i+17] = True
						if i+18 < len(chr_seq):
							array_As_chrI[i+18] = True
						if i+19 < len(chr_seq):
							array_As_chrI[i+19] = True
						if i+20 < len(chr_seq):
							array_As_chrI[i+20] = True

					#AsCas12a Revwerse PAM
					if (chr_seq[i]=='A' or chr_seq[i]=='G' or chr_seq[i]=='C') and chr_seq[i+1]=='A' and chr_seq[i+2]=='A' and chr_seq[i+3]=='A':
						As_PAM_chrI = As_PAM_chrI + 1
						if i-17 >= 0:
							array_As_chrI[i-17] = True
						if i-16 >= 0:
							array_As_chrI[i-16] = True
						if i-15 >= 0:
							array_As_chrI[i-15] = True
						if i-14 >= 0:
							array_As_chrI[i-14] = True
						if i-13 >= 0:
							array_As_chrI[i-13] = True
						if i-12 >= 0:
							array_As_chrI[i-12] = True
						if i-11 >= 0:
							array_As_chrI[i-11] = True
						if i-10 >= 0:
							array_As_chrI[i-10] = True
						if i-9 >= 0:
							array_As_chrI[i-9] = True
						if i-8 >= 0:
							array_As_chrI[i-8] = True
						if i-7 >= 0:
							array_As_chrI[i-7] = True
						if i-6 >= 0:
							array_As_chrI[i-6] = True
						if i-5 >= 0:
							array_As_chrI[i-5] = True
						if i-4 >= 0:
							array_As_chrI[i-4] = True
						if i-3 >= 0:
							array_As_chrI[i-3] = True
						if i-2 >= 0:
							array_As_chrI[i-2] = True
						if i-1 >= 0:
							array_As_chrI[i-1] = True
						if i+1 < len(chr_seq):
							array_As_chrI[i+1] = True
						if i+2 < len(chr_seq):
							array_As_chrI[i+2] = True
						if i+3 < len(chr_seq):
							array_As_chrI[i+3] = True

			for i in range(len(chr_seq)):
				array_SpSa_chrI[i] = (array_Sp_chrI[i] or array_Sa_chrI[i])
				array_SpAs_chrI[i] = (array_Sp_chrI[i] or array_As_chrI[i])
				array_SaAs_chrI[i] = (array_Sa_chrI[i] or array_As_chrI[i])
				array_SpSaAs_chrI[i] = (array_Sp_chrI[i] or array_Sa_chrI[i] or array_As_chrI[i])

			Sp_target_num_chrI = 0
			Sa_target_num_chrI = 0
			As_target_num_chrI = 0
			SpSa_target_num_chrI = 0
			SpAs_target_num_chrI = 0
			SaAs_target_num_chrI = 0
			SpSaAs_target_num_chrI = 0

			for i in range(len(chr_seq)):
				Sp_target_num_chrI = Sp_target_num_chrI + array_Sp_chrI[i]
				Sa_target_num_chrI = Sa_target_num_chrI + array_Sa_chrI[i]
				As_target_num_chrI = As_target_num_chrI + array_As_chrI[i]
				SpSa_target_num_chrI = SpSa_target_num_chrI + array_SpSa_chrI[i]
				SpAs_target_num_chrI = SpAs_target_num_chrI + array_SpAs_chrI[i]
				SaAs_target_num_chrI = SaAs_target_num_chrI + array_SaAs_chrI[i]
				SpSaAs_target_num_chrI = SpSaAs_target_num_chrI + array_SpSaAs_chrI[i]

			print('chrI length: ' + str(len(chr_seq)))
			print('chrI Sp: ' + str(int(Sp_PAM_chrI)) + ' ' + str(int(Sp_target_num_chrI)) + ' ' + str(Sp_target_num_chrI/len(chr_seq)))
			print('chrI Sa: ' + str(int(Sa_PAM_chrI)) + ' ' + str(int(Sa_target_num_chrI)) + ' ' + str(Sa_target_num_chrI/len(chr_seq)))
			print('chrI As: ' + str(int(As_PAM_chrI)) + ' ' + str(int(As_target_num_chrI)) + ' ' + str(As_target_num_chrI/len(chr_seq)))
			print('chrI SpSa: ' + str(int(SpSa_target_num_chrI)) + ' ' + str(SpSa_target_num_chrI/len(chr_seq)))
			print('chrI SpAs: ' + str(int(SpAs_target_num_chrI)) + ' ' + str(SpAs_target_num_chrI/len(chr_seq)))
			print('chrI SaAs: ' + str(int(SaAs_target_num_chrI)) + ' ' + str(SaAs_target_num_chrI/len(chr_seq)))
			print('chrI SpSaAs: ' + str(int(SpSaAs_target_num_chrI)) + ' ' + str(SpSaAs_target_num_chrI/len(chr_seq)))

		elif record.id == "chrII":
			chr_seq = record.seq
			Sp_PAM_chrII = 0
			Sa_PAM_chrII = 0
			As_PAM_chrII = 0
			for i in range(len(chr_seq)):
				#SpCas9
				if i+1 < len(chr_seq):
					#SpCas9 Forward PAM
					if chr_seq[i]=='G' and chr_seq[i+1]=='G':
						Sp_PAM_chrII = Sp_PAM_chrII + 1
						if i-12 >=0:
							array_Sp_chrII[i-12] = True
						if i-11 >=0:
							array_Sp_chrII[i-11] = True
						if i-10 >=0:
							array_Sp_chrII[i-10] = True
						if i-9 >= 0:
							array_Sp_chrII[i-9] = True
						if i-8 >= 0:
							array_Sp_chrII[i-8] = True
						if i-7 >= 0:
							array_Sp_chrII[i-7] = True
						if i-6 >= 0:
							array_Sp_chrII[i-6] = True
						if i-5 >= 0:
							array_Sp_chrII[i-5] = True
						if i-4 >= 0:
							array_Sp_chrII[i-4] = True
						if i-3 >= 0:
							array_Sp_chrII[i-3] = True
						if i-2 >= 0:
							array_Sp_chrII[i-2] = True
						#if i-1 >= 0:
						#	array_Sp_chrII[i-1] = True
						if i >= 0:
							array_Sp_chrII[i] = True
						if i+1 < len(chr_seq):
							array_Sp_chrII[i+1] = True

					#SpCas9 Reverse PAM
					if chr_seq[i]=='C' and chr_seq[i+1]=='C':
						Sp_PAM_chrII = Sp_PAM_chrII + 1
						if i < len(chr_seq):
							array_Sp_chrII[i] = True
						if i+1 < len(chr_seq):
							array_Sp_chrII[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sp_chrII[i+2] = True
						if i+3 < len(chr_seq):
							array_Sp_chrII[i+3] = True
						if i+4 < len(chr_seq):
							array_Sp_chrII[i+4] = True
						if i+5 < len(chr_seq):
							array_Sp_chrII[i+5] = True
						if i+6 < len(chr_seq):
							array_Sp_chrII[i+6] = True
						if i+7 < len(chr_seq):
							array_Sp_chrII[i+7] = True
						if i+8 < len(chr_seq):
							array_Sp_chrII[i+8] = True
						if i+9 < len(chr_seq):
							array_Sp_chrII[i+9] = True
						if i+10 < len(chr_seq):
							array_Sp_chrII[i+10] = True
						if i+11 < len(chr_seq):
							array_Sp_chrII[i+11] = True
						if i+12 < len(chr_seq):
							array_Sp_chrII[i+12] = True
						if i+13 < len(chr_seq):
							array_Sp_chrII[i+13] = True

				#SaCas9
				if i+3 < len(chr_seq):
					#SaCas9 Forward PAM
					if chr_seq[i]=='G' and (chr_seq[i+1]=='A' or chr_seq[i+1]=='G') and (chr_seq[i+2]=='A' or chr_seq[i+2]=='G') and chr_seq[i+3]=='T':
						Sa_PAM_chrII = Sa_PAM_chrII + 1
						if i-13 >= 0:
							array_Sa_chrII[i-13] = True
						if i-12 >= 0:
							array_Sa_chrII[i-12] = True
						if i-11 >= 0:
							array_Sa_chrII[i-11] = True
						if i-10 >= 0:
							array_Sa_chrII[i-10] = True
						if i-9 >= 0:
							array_Sa_chrII[i-9] = True
						if i-8 >= 0:
							array_Sa_chrII[i-8] = True
						if i-7 >= 0:
							array_Sa_chrII[i-7] = True
						if i-6 >= 0:
							array_Sa_chrII[i-6] = True
						if i-5 >= 0:
							array_Sa_chrII[i-5] = True
						if i-4 >= 0:
							array_Sa_chrII[i-4] = True
						if i-3 >= 0:
							array_Sa_chrII[i-3] = True
						#if i-2 >= 0:
						#	array_Sa_chrII[i-2] = True
						#if i-1 >= 0:
						#	array_Sa_chrII[i-1] = True
						if i >= 0:
							array_Sa_chrII[i] = True
						#if i+1 < len(chr_seq):
						#	array_Sa_chrII[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sa_chrII[i+2] = True
						if i+3 < len(chr_seq):
							array_Sa_chrII[i+3] = True

					#SaCas9 Reverse PAM
					if chr_seq[i]=='A' and (chr_seq[i+1]=='C' or chr_seq[i+1]=='T') and (chr_seq[i+2]=='C' or chr_seq[i+2]=='T') and chr_seq[i+3]=='C':
						Sa_PAM_chrII = Sa_PAM_chrII + 1
						if i < len(chr_seq):
							array_Sa_chrII[i] = True
						#if i+1 < len(chr_seq):
						#	array_Sa_chrII[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sa_chrII[i+2] = True
						if i+3 < len(chr_seq):
							array_Sa_chrII[i+3] = True
						#if i+4 < len(chr_seq):
						#	array_Sa_chrII[i+4] = True
						#if i+5 < len(chr_seq):
						#	array_Sa_chrII[i+5] = True
						if i+6 < len(chr_seq):
							array_Sa_chrII[i+6] = True
						if i+7 < len(chr_seq):
							array_Sa_chrII[i+7] = True
						if i+8 < len(chr_seq):
							array_Sa_chrII[i+8] = True
						if i+9 < len(chr_seq):
							array_Sa_chrII[i+9] = True
						if i+10 < len(chr_seq):
							array_Sa_chrII[i+10] = True
						if i+11 < len(chr_seq):
							array_Sa_chrII[i+11] = True
						if i+12 < len(chr_seq):
							array_Sa_chrII[i+12] = True
						if i+13 < len(chr_seq):
							array_Sa_chrII[i+13] = True
						if i+14 < len(chr_seq):
							array_Sa_chrII[i+14] = True
						if i+15 < len(chr_seq):
							array_Sa_chrII[i+15] = True
						if i+16 < len(chr_seq):
							array_Sa_chrII[i+16] = True

				#AsCas12a
				if i+3 < len(chr_seq):
					#AsCas12 Forward PAM
					if chr_seq[i]=='T' and chr_seq[i+1]=='T' and chr_seq[i+2]=='T' and (chr_seq[i+3]=='A' or chr_seq[i+3]=='G' or chr_seq[i+3]=='C'):
						As_PAM_chrII = As_PAM_chrII + 1
						if i < len(chr_seq):
							array_As_chrII[i] = True
						if i+1 < len(chr_seq):
							array_As_chrII[i+1] = True
						if i+2 < len(chr_seq):
							array_As_chrII[i+2] = True
						if i+4 < len(chr_seq):
							array_As_chrII[i+4] = True
						if i+5 < len(chr_seq):
							array_As_chrII[i+5] = True
						if i+6 < len(chr_seq):
							array_As_chrII[i+6] = True
						if i+7 < len(chr_seq):
							array_As_chrII[i+7] = True
						if i+8 < len(chr_seq):
							array_As_chrII[i+8] = True
						if i+9 < len(chr_seq):
							array_As_chrII[i+9] = True
						if i+10 < len(chr_seq):
							array_As_chrII[i+10] = True
						if i+11 < len(chr_seq):
							array_As_chrII[i+11] = True
						if i+12 < len(chr_seq):
							array_As_chrII[i+12] = True
						if i+13 < len(chr_seq):
							array_As_chrII[i+13] = True
						if i+14 < len(chr_seq):
							array_As_chrII[i+14] = True
						if i+15 < len(chr_seq):
							array_As_chrII[i+15] = True
						if i+16 < len(chr_seq):
							array_As_chrII[i+16] = True
						if i+17 < len(chr_seq):
							array_As_chrII[i+17] = True
						if i+18 < len(chr_seq):
							array_As_chrII[i+18] = True
						if i+19 < len(chr_seq):
							array_As_chrII[i+19] = True
						if i+20 < len(chr_seq):
							array_As_chrII[i+20] = True

					#AsCas12a Revwerse PAM
					if (chr_seq[i]=='A' or chr_seq[i]=='G' or chr_seq[i]=='C') and chr_seq[i+1]=='A' and chr_seq[i+2]=='A' and chr_seq[i+3]=='A':
						As_PAM_chrII = As_PAM_chrII + 1
						if i-17 >= 0:
							array_As_chrII[i-17] = True
						if i-16 >= 0:
							array_As_chrII[i-16] = True
						if i-15 >= 0:
							array_As_chrII[i-15] = True
						if i-14 >= 0:
							array_As_chrII[i-14] = True
						if i-13 >= 0:
							array_As_chrII[i-13] = True
						if i-12 >= 0:
							array_As_chrII[i-12] = True
						if i-11 >= 0:
							array_As_chrII[i-11] = True
						if i-10 >= 0:
							array_As_chrII[i-10] = True
						if i-9 >= 0:
							array_As_chrII[i-9] = True
						if i-8 >= 0:
							array_As_chrII[i-8] = True
						if i-7 >= 0:
							array_As_chrII[i-7] = True
						if i-6 >= 0:
							array_As_chrII[i-6] = True
						if i-5 >= 0:
							array_As_chrII[i-5] = True
						if i-4 >= 0:
							array_As_chrII[i-4] = True
						if i-3 >= 0:
							array_As_chrII[i-3] = True
						if i-2 >= 0:
							array_As_chrII[i-2] = True
						if i-1 >= 0:
							array_As_chrII[i-1] = True
						if i+1 < len(chr_seq):
							array_As_chrII[i+1] = True
						if i+2 < len(chr_seq):
							array_As_chrII[i+2] = True
						if i+3 < len(chr_seq):
							array_As_chrII[i+3] = True

			#複数のCasについてtargetabilityの和を計算する
			for i in range(len(chr_seq)):
				array_SpSa_chrII[i] = (array_Sp_chrII[i] or array_Sa_chrII[i])
				array_SpAs_chrII[i] = (array_Sp_chrII[i] or array_As_chrII[i])
				array_SaAs_chrII[i] = (array_Sa_chrII[i] or array_As_chrII[i])
				array_SpSaAs_chrII[i] = (array_Sp_chrII[i] or array_Sa_chrII[i] or array_As_chrII[i])

			Sp_target_num_chrII = 0
			Sa_target_num_chrII = 0
			As_target_num_chrII = 0
			SpSa_target_num_chrII = 0
			SpAs_target_num_chrII = 0
			SaAs_target_num_chrII = 0
			SpSaAs_target_num_chrII = 0

			for i in range(len(chr_seq)):
				Sp_target_num_chrII = Sp_target_num_chrII + array_Sp_chrII[i]
				Sa_target_num_chrII = Sa_target_num_chrII + array_Sa_chrII[i]
				As_target_num_chrII = As_target_num_chrII + array_As_chrII[i]
				SpSa_target_num_chrII = SpSa_target_num_chrII + array_SpSa_chrII[i]
				SpAs_target_num_chrII = SpAs_target_num_chrII + array_SpAs_chrII[i]
				SaAs_target_num_chrII = SaAs_target_num_chrII + array_SaAs_chrII[i]
				SpSaAs_target_num_chrII = SpSaAs_target_num_chrII + array_SpSaAs_chrII[i]

			print('chrII length: ' + str(len(chr_seq)))
			print('chrII Sp: ' + str(int(Sp_PAM_chrII)) + ' ' + str(int(Sp_target_num_chrII)) + ' ' + str(Sp_target_num_chrII/len(chr_seq)))
			print('chrII Sa: ' + str(int(Sa_PAM_chrII)) + ' ' + str(int(Sa_target_num_chrII)) + ' ' + str(Sa_target_num_chrII/len(chr_seq)))
			print('chrII As: ' + str(int(As_PAM_chrII)) + ' ' + str(int(As_target_num_chrII)) + ' ' + str(As_target_num_chrII/len(chr_seq)))
			print('chrII SpSa: ' + str(int(SpSa_target_num_chrII)) + ' ' + str(SpSa_target_num_chrII/len(chr_seq)))
			print('chrII SpAs: ' + str(int(SpAs_target_num_chrII)) + ' ' + str(SpAs_target_num_chrII/len(chr_seq)))
			print('chrII SaAs: ' + str(int(SaAs_target_num_chrII)) + ' ' + str(SaAs_target_num_chrII/len(chr_seq)))
			print('chrII SpSaAs: ' + str(int(SpSaAs_target_num_chrII)) + ' ' + str(SpSaAs_target_num_chrII/len(chr_seq)))

		elif record.id == "chrIII":
			chr_seq = record.seq
			Sp_PAM_chrIII = 0
			Sa_PAM_chrIII = 0
			As_PAM_chrIII = 0
			for i in range(len(chr_seq)):
								#SpCas9
				if i+1 < len(chr_seq):
					#SpCas9 Forward PAM
					if chr_seq[i]=='G' and chr_seq[i+1]=='G':
						Sp_PAM_chrIII = Sp_PAM_chrIII + 1
						if i-12 >=0:
							array_Sp_chrIII[i-12] = True
						if i-11 >=0:
							array_Sp_chrIII[i-11] = True
						if i-10 >=0:
							array_Sp_chrIII[i-10] = True
						if i-9 >= 0:
							array_Sp_chrIII[i-9] = True
						if i-8 >= 0:
							array_Sp_chrIII[i-8] = True
						if i-7 >= 0:
							array_Sp_chrIII[i-7] = True
						if i-6 >= 0:
							array_Sp_chrIII[i-6] = True
						if i-5 >= 0:
							array_Sp_chrIII[i-5] = True
						if i-4 >= 0:
							array_Sp_chrIII[i-4] = True
						if i-3 >= 0:
							array_Sp_chrIII[i-3] = True
						if i-2 >= 0:
							array_Sp_chrIII[i-2] = True
						#if i-1 >= 0:
						#	array_Sp_chrIII[i-1] = True
						if i >= 0:
							array_Sp_chrIII[i] = True
						if i+1 < len(chr_seq):
							array_Sp_chrIII[i+1] = True

					#SpCas9 Reverse PAM
					if chr_seq[i]=='C' and chr_seq[i+1]=='C':
						Sp_PAM_chrIII = Sp_PAM_chrIII + 1
						if i < len(chr_seq):
							array_Sp_chrIII[i] = True
						if i+1 < len(chr_seq):
							array_Sp_chrIII[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sp_chrIII[i+2] = True
						if i+3 < len(chr_seq):
							array_Sp_chrIII[i+3] = True
						if i+4 < len(chr_seq):
							array_Sp_chrIII[i+4] = True
						if i+5 < len(chr_seq):
							array_Sp_chrIII[i+5] = True
						if i+6 < len(chr_seq):
							array_Sp_chrIII[i+6] = True
						if i+7 < len(chr_seq):
							array_Sp_chrIII[i+7] = True
						if i+8 < len(chr_seq):
							array_Sp_chrIII[i+8] = True
						if i+9 < len(chr_seq):
							array_Sp_chrIII[i+9] = True
						if i+10 < len(chr_seq):
							array_Sp_chrIII[i+10] = True
						if i+11 < len(chr_seq):
							array_Sp_chrIII[i+11] = True
						if i+12 < len(chr_seq):
							array_Sp_chrIII[i+12] = True
						if i+13 < len(chr_seq):
							array_Sp_chrIII[i+13] = True

				#SaCas9
				if i+3 < len(chr_seq):
					#SaCas9 Forward PAM
					if chr_seq[i]=='G' and (chr_seq[i+1]=='A' or chr_seq[i+1]=='G') and (chr_seq[i+2]=='A' or chr_seq[i+2]=='G') and chr_seq[i+3]=='T':
						Sa_PAM_chrIII = Sa_PAM_chrIII + 1
						if i-13 >= 0:
							array_Sa_chrIII[i-13] = True
						if i-12 >= 0:
							array_Sa_chrIII[i-12] = True
						if i-11 >= 0:
							array_Sa_chrIII[i-11] = True
						if i-10 >= 0:
							array_Sa_chrIII[i-10] = True
						if i-9 >= 0:
							array_Sa_chrIII[i-9] = True
						if i-8 >= 0:
							array_Sa_chrIII[i-8] = True
						if i-7 >= 0:
							array_Sa_chrIII[i-7] = True
						if i-6 >= 0:
							array_Sa_chrIII[i-6] = True
						if i-5 >= 0:
							array_Sa_chrIII[i-5] = True
						if i-4 >= 0:
							array_Sa_chrIII[i-4] = True
						if i-3 >= 0:
							array_Sa_chrIII[i-3] = True
						#if i-2 >= 0:
						#	array_Sa_chrIII[i-2] = True
						#if i-1 >= 0:
						#	array_Sa_chrIII[i-1] = True
						if i >= 0:
							array_Sa_chrIII[i] = True
						#if i+1 < len(chr_seq):
						#	array_Sa_chrIII[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sa_chrIII[i+2] = True
						if i+3 < len(chr_seq):
							array_Sa_chrIII[i+3] = True

					#SaCas9 Reverse PAM
					if chr_seq[i]=='A' and (chr_seq[i+1]=='C' or chr_seq[i+1]=='T') and (chr_seq[i+2]=='C' or chr_seq[i+2]=='T') and chr_seq[i+3]=='C':
						Sa_PAM_chrIII = Sa_PAM_chrIII + 1
						if i < len(chr_seq):
							array_Sa_chrIII[i] = True
						#if i+1 < len(chr_seq):
						#	array_Sa_chrIII[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sa_chrIII[i+2] = True
						if i+3 < len(chr_seq):
							array_Sa_chrIII[i+3] = True
						#if i+4 < len(chr_seq):
						#	array_Sa_chrIII[i+4] = True
						#if i+5 < len(chr_seq):
						#	array_Sa_chrIII[i+5] = True
						if i+6 < len(chr_seq):
							array_Sa_chrIII[i+6] = True
						if i+7 < len(chr_seq):
							array_Sa_chrIII[i+7] = True
						if i+8 < len(chr_seq):
							array_Sa_chrIII[i+8] = True
						if i+9 < len(chr_seq):
							array_Sa_chrIII[i+9] = True
						if i+10 < len(chr_seq):
							array_Sa_chrIII[i+10] = True
						if i+11 < len(chr_seq):
							array_Sa_chrIII[i+11] = True
						if i+12 < len(chr_seq):
							array_Sa_chrIII[i+12] = True
						if i+13 < len(chr_seq):
							array_Sa_chrIII[i+13] = True
						if i+14 < len(chr_seq):
							array_Sa_chrIII[i+14] = True
						if i+15 < len(chr_seq):
							array_Sa_chrIII[i+15] = True
						if i+16 < len(chr_seq):
							array_Sa_chrIII[i+16] = True

				#AsCas12a
				if i+3 < len(chr_seq):
					#AsCas12 Forward PAM
					if chr_seq[i]=='T' and chr_seq[i+1]=='T' and chr_seq[i+2]=='T' and (chr_seq[i+3]=='A' or chr_seq[i+3]=='G' or chr_seq[i+3]=='C'):
						As_PAM_chrIII = As_PAM_chrIII + 1
						if i < len(chr_seq):
							array_As_chrIII[i] = True
						if i+1 < len(chr_seq):
							array_As_chrIII[i+1] = True
						if i+2 < len(chr_seq):
							array_As_chrIII[i+2] = True
						if i+4 < len(chr_seq):
							array_As_chrIII[i+4] = True
						if i+5 < len(chr_seq):
							array_As_chrIII[i+5] = True
						if i+6 < len(chr_seq):
							array_As_chrIII[i+6] = True
						if i+7 < len(chr_seq):
							array_As_chrIII[i+7] = True
						if i+8 < len(chr_seq):
							array_As_chrIII[i+8] = True
						if i+9 < len(chr_seq):
							array_As_chrIII[i+9] = True
						if i+10 < len(chr_seq):
							array_As_chrIII[i+10] = True
						if i+11 < len(chr_seq):
							array_As_chrIII[i+11] = True
						if i+12 < len(chr_seq):
							array_As_chrIII[i+12] = True
						if i+13 < len(chr_seq):
							array_As_chrIII[i+13] = True
						if i+14 < len(chr_seq):
							array_As_chrIII[i+14] = True
						if i+15 < len(chr_seq):
							array_As_chrIII[i+15] = True
						if i+16 < len(chr_seq):
							array_As_chrIII[i+16] = True
						if i+17 < len(chr_seq):
							array_As_chrIII[i+17] = True
						if i+18 < len(chr_seq):
							array_As_chrIII[i+18] = True
						if i+19 < len(chr_seq):
							array_As_chrIII[i+19] = True
						if i+20 < len(chr_seq):
							array_As_chrIII[i+20] = True

					#AsCas12a Revwerse PAM
					if (chr_seq[i]=='A' or chr_seq[i]=='G' or chr_seq[i]=='C') and chr_seq[i+1]=='A' and chr_seq[i+2]=='A' and chr_seq[i+3]=='A':
						As_PAM_chrIII = As_PAM_chrIII + 1
						if i-17 >= 0:
							array_As_chrIII[i-17] = True
						if i-16 >= 0:
							array_As_chrIII[i-16] = True
						if i-15 >= 0:
							array_As_chrIII[i-15] = True
						if i-14 >= 0:
							array_As_chrIII[i-14] = True
						if i-13 >= 0:
							array_As_chrIII[i-13] = True
						if i-12 >= 0:
							array_As_chrIII[i-12] = True
						if i-11 >= 0:
							array_As_chrIII[i-11] = True
						if i-10 >= 0:
							array_As_chrIII[i-10] = True
						if i-9 >= 0:
							array_As_chrIII[i-9] = True
						if i-8 >= 0:
							array_As_chrIII[i-8] = True
						if i-7 >= 0:
							array_As_chrIII[i-7] = True
						if i-6 >= 0:
							array_As_chrIII[i-6] = True
						if i-5 >= 0:
							array_As_chrIII[i-5] = True
						if i-4 >= 0:
							array_As_chrIII[i-4] = True
						if i-3 >= 0:
							array_As_chrIII[i-3] = True
						if i-2 >= 0:
							array_As_chrIII[i-2] = True
						if i-1 >= 0:
							array_As_chrIII[i-1] = True
						if i+1 < len(chr_seq):
							array_As_chrIII[i+1] = True
						if i+2 < len(chr_seq):
							array_As_chrIII[i+2] = True
						if i+3 < len(chr_seq):
							array_As_chrIII[i+3] = True

			for i in range(len(chr_seq)):
				array_SpSa_chrIII[i] = (array_Sp_chrIII[i] or array_Sa_chrIII[i])
				array_SpAs_chrIII[i] = (array_Sp_chrIII[i] or array_As_chrIII[i])
				array_SaAs_chrIII[i] = (array_Sa_chrIII[i] or array_As_chrIII[i])
				array_SpSaAs_chrIII[i] = (array_Sp_chrIII[i] or array_Sa_chrIII[i] or array_As_chrIII[i])

			Sp_target_num_chrIII = 0
			Sa_target_num_chrIII = 0
			As_target_num_chrIII = 0
			SpSa_target_num_chrIII = 0
			SpAs_target_num_chrIII = 0
			SaAs_target_num_chrIII = 0
			SpSaAs_target_num_chrIII = 0

			for i in range(len(chr_seq)):
				Sp_target_num_chrIII = Sp_target_num_chrIII + array_Sp_chrIII[i]
				Sa_target_num_chrIII = Sa_target_num_chrIII + array_Sa_chrIII[i]
				As_target_num_chrIII = As_target_num_chrIII + array_As_chrIII[i]
				SpSa_target_num_chrIII = SpSa_target_num_chrIII + array_SpSa_chrIII[i]
				SpAs_target_num_chrIII = SpAs_target_num_chrIII + array_SpAs_chrIII[i]
				SaAs_target_num_chrIII = SaAs_target_num_chrIII + array_SaAs_chrIII[i]
				SpSaAs_target_num_chrIII = SpSaAs_target_num_chrIII + array_SpSaAs_chrIII[i]

			print('chrIII length: ' + str(len(chr_seq)))
			print('chrIII Sp: ' + str(int(Sp_PAM_chrIII)) + ' ' + str(int(Sp_target_num_chrIII)) + ' ' + str(Sp_target_num_chrIII/len(chr_seq)))
			print('chrIII Sa: ' + str(int(Sa_PAM_chrIII)) + ' ' + str(int(Sa_target_num_chrIII)) + ' ' + str(Sa_target_num_chrIII/len(chr_seq)))
			print('chrIII As: ' + str(int(As_PAM_chrIII)) + ' ' + str(int(As_target_num_chrIII)) + ' ' + str(As_target_num_chrIII/len(chr_seq)))
			print('chrIII SpSa: ' + str(int(SpSa_target_num_chrIII)) + ' ' + str(SpSa_target_num_chrIII/len(chr_seq)))
			print('chrIII SpAs: ' + str(int(SpAs_target_num_chrIII)) + ' ' + str(SpAs_target_num_chrIII/len(chr_seq)))
			print('chrIII SaAs: ' + str(int(SaAs_target_num_chrIII)) + ' ' + str(SaAs_target_num_chrIII/len(chr_seq)))
			print('chrIII SpSaAs: ' + str(int(SpSaAs_target_num_chrIII)) + ' ' + str(SpSaAs_target_num_chrIII/len(chr_seq)))

		elif record.id == "chrIV":
			chr_seq = record.seq
			Sp_PAM_chrIV = 0
			Sa_PAM_chrIV = 0
			As_PAM_chrIV = 0
			for i in range(len(chr_seq)):
				#SpCas9
				if i+1 < len(chr_seq):
					#SpCas9 Forward PAM
					if chr_seq[i]=='G' and chr_seq[i+1]=='G':
						Sp_PAM_chrIV = Sp_PAM_chrIV + 1
						if i-12 >=0:
							array_Sp_chrIV[i-12] = True
						if i-11 >=0:
							array_Sp_chrIV[i-11] = True
						if i-10 >=0:
							array_Sp_chrIV[i-10] = True
						if i-9 >= 0:
							array_Sp_chrIV[i-9] = True
						if i-8 >= 0:
							array_Sp_chrIV[i-8] = True
						if i-7 >= 0:
							array_Sp_chrIV[i-7] = True
						if i-6 >= 0:
							array_Sp_chrIV[i-6] = True
						if i-5 >= 0:
							array_Sp_chrIV[i-5] = True
						if i-4 >= 0:
							array_Sp_chrIV[i-4] = True
						if i-3 >= 0:
							array_Sp_chrIV[i-3] = True
						if i-2 >= 0:
							array_Sp_chrIV[i-2] = True
						#if i-1 >= 0:
						#	array_Sp_chrIV[i-1] = True
						if i >= 0:
							array_Sp_chrIV[i] = True
						if i+1 < len(chr_seq):
							array_Sp_chrIV[i+1] = True

					#SpCas9 Reverse PAM
					if chr_seq[i]=='C' and chr_seq[i+1]=='C':
						Sp_PAM_chrIV = Sp_PAM_chrIV + 1
						if i < len(chr_seq):
							array_Sp_chrIV[i] = True
						if i+1 < len(chr_seq):
							array_Sp_chrIV[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sp_chrIV[i+2] = True
						if i+3 < len(chr_seq):
							array_Sp_chrIV[i+3] = True
						if i+4 < len(chr_seq):
							array_Sp_chrIV[i+4] = True
						if i+5 < len(chr_seq):
							array_Sp_chrIV[i+5] = True
						if i+6 < len(chr_seq):
							array_Sp_chrIV[i+6] = True
						if i+7 < len(chr_seq):
							array_Sp_chrIV[i+7] = True
						if i+8 < len(chr_seq):
							array_Sp_chrIV[i+8] = True
						if i+9 < len(chr_seq):
							array_Sp_chrIV[i+9] = True
						if i+10 < len(chr_seq):
							array_Sp_chrIV[i+10] = True
						if i+11 < len(chr_seq):
							array_Sp_chrIV[i+11] = True
						if i+12 < len(chr_seq):
							array_Sp_chrIV[i+12] = True
						if i+13 < len(chr_seq):
							array_Sp_chrIV[i+13] = True

				#SaCas9
				if i+3 < len(chr_seq):
					#SaCas9 Forward PAM
					if chr_seq[i]=='G' and (chr_seq[i+1]=='A' or chr_seq[i+1]=='G') and (chr_seq[i+2]=='A' or chr_seq[i+2]=='G') and chr_seq[i+3]=='T':
						Sa_PAM_chrIV = Sa_PAM_chrIV + 1
						if i-13 >= 0:
							array_Sa_chrIV[i-13] = True
						if i-12 >= 0:
							array_Sa_chrIV[i-12] = True
						if i-11 >= 0:
							array_Sa_chrIV[i-11] = True
						if i-10 >= 0:
							array_Sa_chrIV[i-10] = True
						if i-9 >= 0:
							array_Sa_chrIV[i-9] = True
						if i-8 >= 0:
							array_Sa_chrIV[i-8] = True
						if i-7 >= 0:
							array_Sa_chrIV[i-7] = True
						if i-6 >= 0:
							array_Sa_chrIV[i-6] = True
						if i-5 >= 0:
							array_Sa_chrIV[i-5] = True
						if i-4 >= 0:
							array_Sa_chrIV[i-4] = True
						if i-3 >= 0:
							array_Sa_chrIV[i-3] = True
						#if i-2 >= 0:
						#	array_Sa_chrIV[i-2] = True
						#if i-1 >= 0:
						#	array_Sa_chrIV[i-1] = True
						if i >= 0:
							array_Sa_chrIV[i] = True
						#if i+1 < len(chr_seq):
						#	array_Sa_chrIV[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sa_chrIV[i+2] = True
						if i+3 < len(chr_seq):
							array_Sa_chrIV[i+3] = True

					#SaCas9 Reverse PAM
					if chr_seq[i]=='A' and (chr_seq[i+1]=='C' or chr_seq[i+1]=='T') and (chr_seq[i+2]=='C' or chr_seq[i+2]=='T') and chr_seq[i+3]=='C':
						Sa_PAM_chrIV = Sa_PAM_chrIV + 1
						if i < len(chr_seq):
							array_Sa_chrIV[i] = True
						#if i+1 < len(chr_seq):
						#	array_Sa_chrIV[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sa_chrIV[i+2] = True
						if i+3 < len(chr_seq):
							array_Sa_chrIV[i+3] = True
						#if i+4 < len(chr_seq):
						#	array_Sa_chrIV[i+4] = True
						#if i+5 < len(chr_seq):
						#	array_Sa_chrIV[i+5] = True
						if i+6 < len(chr_seq):
							array_Sa_chrIV[i+6] = True
						if i+7 < len(chr_seq):
							array_Sa_chrIV[i+7] = True
						if i+8 < len(chr_seq):
							array_Sa_chrIV[i+8] = True
						if i+9 < len(chr_seq):
							array_Sa_chrIV[i+9] = True
						if i+10 < len(chr_seq):
							array_Sa_chrIV[i+10] = True
						if i+11 < len(chr_seq):
							array_Sa_chrIV[i+11] = True
						if i+12 < len(chr_seq):
							array_Sa_chrIV[i+12] = True
						if i+13 < len(chr_seq):
							array_Sa_chrIV[i+13] = True
						if i+14 < len(chr_seq):
							array_Sa_chrIV[i+14] = True
						if i+15 < len(chr_seq):
							array_Sa_chrIV[i+15] = True
						if i+16 < len(chr_seq):
							array_Sa_chrIV[i+16] = True

				#AsCas12a
				if i+3 < len(chr_seq):
					#AsCas12 Forward PAM
					if chr_seq[i]=='T' and chr_seq[i+1]=='T' and chr_seq[i+2]=='T' and (chr_seq[i+3]=='A' or chr_seq[i+3]=='G' or chr_seq[i+3]=='C'):
						As_PAM_chrIV = As_PAM_chrIV + 1
						if i < len(chr_seq):
							array_As_chrIV[i] = True
						if i+1 < len(chr_seq):
							array_As_chrIV[i+1] = True
						if i+2 < len(chr_seq):
							array_As_chrIV[i+2] = True
						if i+4 < len(chr_seq):
							array_As_chrIV[i+4] = True
						if i+5 < len(chr_seq):
							array_As_chrIV[i+5] = True
						if i+6 < len(chr_seq):
							array_As_chrIV[i+6] = True
						if i+7 < len(chr_seq):
							array_As_chrIV[i+7] = True
						if i+8 < len(chr_seq):
							array_As_chrIV[i+8] = True
						if i+9 < len(chr_seq):
							array_As_chrIV[i+9] = True
						if i+10 < len(chr_seq):
							array_As_chrIV[i+10] = True
						if i+11 < len(chr_seq):
							array_As_chrIV[i+11] = True
						if i+12 < len(chr_seq):
							array_As_chrIV[i+12] = True
						if i+13 < len(chr_seq):
							array_As_chrIV[i+13] = True
						if i+14 < len(chr_seq):
							array_As_chrIV[i+14] = True
						if i+15 < len(chr_seq):
							array_As_chrIV[i+15] = True
						if i+16 < len(chr_seq):
							array_As_chrIV[i+16] = True
						if i+17 < len(chr_seq):
							array_As_chrIV[i+17] = True
						if i+18 < len(chr_seq):
							array_As_chrIV[i+18] = True
						if i+19 < len(chr_seq):
							array_As_chrIV[i+19] = True
						if i+20 < len(chr_seq):
							array_As_chrIV[i+20] = True

					#AsCas12a Revwerse PAM
					if (chr_seq[i]=='A' or chr_seq[i]=='G' or chr_seq[i]=='C') and chr_seq[i+1]=='A' and chr_seq[i+2]=='A' and chr_seq[i+3]=='A':
						As_PAM_chrIV = As_PAM_chrIV + 1
						if i-17 >= 0:
							array_As_chrIV[i-17] = True
						if i-16 >= 0:
							array_As_chrIV[i-16] = True
						if i-15 >= 0:
							array_As_chrIV[i-15] = True
						if i-14 >= 0:
							array_As_chrIV[i-14] = True
						if i-13 >= 0:
							array_As_chrIV[i-13] = True
						if i-12 >= 0:
							array_As_chrIV[i-12] = True
						if i-11 >= 0:
							array_As_chrIV[i-11] = True
						if i-10 >= 0:
							array_As_chrIV[i-10] = True
						if i-9 >= 0:
							array_As_chrIV[i-9] = True
						if i-8 >= 0:
							array_As_chrIV[i-8] = True
						if i-7 >= 0:
							array_As_chrIV[i-7] = True
						if i-6 >= 0:
							array_As_chrIV[i-6] = True
						if i-5 >= 0:
							array_As_chrIV[i-5] = True
						if i-4 >= 0:
							array_As_chrIV[i-4] = True
						if i-3 >= 0:
							array_As_chrIV[i-3] = True
						if i-2 >= 0:
							array_As_chrIV[i-2] = True
						if i-1 >= 0:
							array_As_chrIV[i-1] = True
						if i+1 < len(chr_seq):
							array_As_chrIV[i+1] = True
						if i+2 < len(chr_seq):
							array_As_chrIV[i+2] = True
						if i+3 < len(chr_seq):
							array_As_chrIV[i+3] = True

			for i in range(len(chr_seq)):
				array_SpSa_chrIV[i] = (array_Sp_chrIV[i] or array_Sa_chrIV[i])
				array_SpAs_chrIV[i] = (array_Sp_chrIV[i] or array_As_chrIV[i])
				array_SaAs_chrIV[i] = (array_Sa_chrIV[i] or array_As_chrIV[i])
				array_SpSaAs_chrIV[i] = (array_Sp_chrIV[i] or array_Sa_chrIV[i] or array_As_chrIV[i])

			Sp_target_num_chrIV = 0
			Sa_target_num_chrIV = 0
			As_target_num_chrIV = 0
			SpSa_target_num_chrIV = 0
			SpAs_target_num_chrIV = 0
			SaAs_target_num_chrIV = 0
			SpSaAs_target_num_chrIV = 0

			for i in range(len(chr_seq)):
				Sp_target_num_chrIV = Sp_target_num_chrIV + array_Sp_chrIV[i]
				Sa_target_num_chrIV = Sa_target_num_chrIV + array_Sa_chrIV[i]
				As_target_num_chrIV = As_target_num_chrIV + array_As_chrIV[i]
				SpSa_target_num_chrIV = SpSa_target_num_chrIV + array_SpSa_chrIV[i]
				SpAs_target_num_chrIV = SpAs_target_num_chrIV + array_SpAs_chrIV[i]
				SaAs_target_num_chrIV = SaAs_target_num_chrIV + array_SaAs_chrIV[i]
				SpSaAs_target_num_chrIV = SpSaAs_target_num_chrIV + array_SpSaAs_chrIV[i]

			print('chrIV length: ' + str(len(chr_seq)))
			print('chrIV Sp: ' + str(int(Sp_PAM_chrIV)) + ' ' + str(int(Sp_target_num_chrIV)) + ' ' + str(Sp_target_num_chrIV/len(chr_seq)))
			print('chrIV Sa: ' + str(int(Sa_PAM_chrIV)) + ' ' + str(int(Sa_target_num_chrIV)) + ' ' + str(Sa_target_num_chrIV/len(chr_seq)))
			print('chrIV As: ' + str(int(As_PAM_chrIV)) + ' ' + str(int(As_target_num_chrIV)) + ' ' + str(As_target_num_chrIV/len(chr_seq)))
			print('chrIV SpSa: ' + str(int(SpSa_target_num_chrIV)) + ' ' + str(SpSa_target_num_chrIV/len(chr_seq)))
			print('chrIV SpAs: ' + str(int(SpAs_target_num_chrIV)) + ' ' + str(SpAs_target_num_chrIV/len(chr_seq)))
			print('chrIV SaAs: ' + str(int(SaAs_target_num_chrIV)) + ' ' + str(SaAs_target_num_chrIV/len(chr_seq)))
			print('chrIV SpSaAs: ' + str(int(SpSaAs_target_num_chrIV)) + ' ' + str(SpSaAs_target_num_chrIV/len(chr_seq)))

		elif record.id == "chrV":
			chr_seq = record.seq
			Sp_PAM_chrV = 0
			Sa_PAM_chrV = 0
			As_PAM_chrV = 0
			for i in range(len(chr_seq)):
				#SpCas9
				if i+1 < len(chr_seq):
					#SpCas9 Forward PAM
					if chr_seq[i]=='G' and chr_seq[i+1]=='G':
						Sp_PAM_chrV = Sp_PAM_chrV + 1
						if i-12 >=0:
							array_Sp_chrV[i-12] = True
						if i-11 >=0:
							array_Sp_chrV[i-11] = True
						if i-10 >=0:
							array_Sp_chrV[i-10] = True
						if i-9 >= 0:
							array_Sp_chrV[i-9] = True
						if i-8 >= 0:
							array_Sp_chrV[i-8] = True
						if i-7 >= 0:
							array_Sp_chrV[i-7] = True
						if i-6 >= 0:
							array_Sp_chrV[i-6] = True
						if i-5 >= 0:
							array_Sp_chrV[i-5] = True
						if i-4 >= 0:
							array_Sp_chrV[i-4] = True
						if i-3 >= 0:
							array_Sp_chrV[i-3] = True
						if i-2 >= 0:
							array_Sp_chrV[i-2] = True
						#if i-1 >= 0:
						#	array_Sp_chrV[i-1] = True
						if i >= 0:
							array_Sp_chrV[i] = True
						if i+1 < len(chr_seq):
							array_Sp_chrV[i+1] = True

					#SpCas9 Reverse PAM
					if chr_seq[i]=='C' and chr_seq[i+1]=='C':
						Sp_PAM_chrV = Sp_PAM_chrV + 1
						if i < len(chr_seq):
							array_Sp_chrV[i] = True
						if i+1 < len(chr_seq):
							array_Sp_chrV[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sp_chrV[i+2] = True
						if i+3 < len(chr_seq):
							array_Sp_chrV[i+3] = True
						if i+4 < len(chr_seq):
							array_Sp_chrV[i+4] = True
						if i+5 < len(chr_seq):
							array_Sp_chrV[i+5] = True
						if i+6 < len(chr_seq):
							array_Sp_chrV[i+6] = True
						if i+7 < len(chr_seq):
							array_Sp_chrV[i+7] = True
						if i+8 < len(chr_seq):
							array_Sp_chrV[i+8] = True
						if i+9 < len(chr_seq):
							array_Sp_chrV[i+9] = True
						if i+10 < len(chr_seq):
							array_Sp_chrV[i+10] = True
						if i+11 < len(chr_seq):
							array_Sp_chrV[i+11] = True
						if i+12 < len(chr_seq):
							array_Sp_chrV[i+12] = True
						if i+13 < len(chr_seq):
							array_Sp_chrV[i+13] = True

				#SaCas9
				if i+3 < len(chr_seq):
					#SaCas9 Forward PAM
					if chr_seq[i]=='G' and (chr_seq[i+1]=='A' or chr_seq[i+1]=='G') and (chr_seq[i+2]=='A' or chr_seq[i+2]=='G') and chr_seq[i+3]=='T':
						Sa_PAM_chrV = Sa_PAM_chrV + 1
						if i-13 >= 0:
							array_Sa_chrV[i-13] = True
						if i-12 >= 0:
							array_Sa_chrV[i-12] = True
						if i-11 >= 0:
							array_Sa_chrV[i-11] = True
						if i-10 >= 0:
							array_Sa_chrV[i-10] = True
						if i-9 >= 0:
							array_Sa_chrV[i-9] = True
						if i-8 >= 0:
							array_Sa_chrV[i-8] = True
						if i-7 >= 0:
							array_Sa_chrV[i-7] = True
						if i-6 >= 0:
							array_Sa_chrV[i-6] = True
						if i-5 >= 0:
							array_Sa_chrV[i-5] = True
						if i-4 >= 0:
							array_Sa_chrV[i-4] = True
						if i-3 >= 0:
							array_Sa_chrV[i-3] = True
						#if i-2 >= 0:
						#	array_Sa_chrV[i-2] = True
						#if i-1 >= 0:
						#	array_Sa_chrV[i-1] = True
						if i >= 0:
							array_Sa_chrV[i] = True
						#if i+1 < len(chr_seq):
						#	array_Sa_chrV[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sa_chrV[i+2] = True
						if i+3 < len(chr_seq):
							array_Sa_chrV[i+3] = True

					#SaCas9 Reverse PAM
					if chr_seq[i]=='A' and (chr_seq[i+1]=='C' or chr_seq[i+1]=='T') and (chr_seq[i+2]=='C' or chr_seq[i+2]=='T') and chr_seq[i+3]=='C':
						Sa_PAM_chrV = Sa_PAM_chrV + 1
						if i < len(chr_seq):
							array_Sa_chrV[i] = True
						#if i+1 < len(chr_seq):
						#	array_Sa_chrV[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sa_chrV[i+2] = True
						if i+3 < len(chr_seq):
							array_Sa_chrV[i+3] = True
						#if i+4 < len(chr_seq):
						#	array_Sa_chrV[i+4] = True
						#if i+5 < len(chr_seq):
						#	array_Sa_chrV[i+5] = True
						if i+6 < len(chr_seq):
							array_Sa_chrV[i+6] = True
						if i+7 < len(chr_seq):
							array_Sa_chrV[i+7] = True
						if i+8 < len(chr_seq):
							array_Sa_chrV[i+8] = True
						if i+9 < len(chr_seq):
							array_Sa_chrV[i+9] = True
						if i+10 < len(chr_seq):
							array_Sa_chrV[i+10] = True
						if i+11 < len(chr_seq):
							array_Sa_chrV[i+11] = True
						if i+12 < len(chr_seq):
							array_Sa_chrV[i+12] = True
						if i+13 < len(chr_seq):
							array_Sa_chrV[i+13] = True
						if i+14 < len(chr_seq):
							array_Sa_chrV[i+14] = True
						if i+15 < len(chr_seq):
							array_Sa_chrV[i+15] = True
						if i+16 < len(chr_seq):
							array_Sa_chrV[i+16] = True

				#AsCas12a
				if i+3 < len(chr_seq):
					#AsCas12 Forward PAM
					if chr_seq[i]=='T' and chr_seq[i+1]=='T' and chr_seq[i+2]=='T' and (chr_seq[i+3]=='A' or chr_seq[i+3]=='G' or chr_seq[i+3]=='C'):
						As_PAM_chrV = As_PAM_chrV + 1
						if i < len(chr_seq):
							array_As_chrV[i] = True
						if i+1 < len(chr_seq):
							array_As_chrV[i+1] = True
						if i+2 < len(chr_seq):
							array_As_chrV[i+2] = True
						if i+4 < len(chr_seq):
							array_As_chrV[i+4] = True
						if i+5 < len(chr_seq):
							array_As_chrV[i+5] = True
						if i+6 < len(chr_seq):
							array_As_chrV[i+6] = True
						if i+7 < len(chr_seq):
							array_As_chrV[i+7] = True
						if i+8 < len(chr_seq):
							array_As_chrV[i+8] = True
						if i+9 < len(chr_seq):
							array_As_chrV[i+9] = True
						if i+10 < len(chr_seq):
							array_As_chrV[i+10] = True
						if i+11 < len(chr_seq):
							array_As_chrV[i+11] = True
						if i+12 < len(chr_seq):
							array_As_chrV[i+12] = True
						if i+13 < len(chr_seq):
							array_As_chrV[i+13] = True
						if i+14 < len(chr_seq):
							array_As_chrV[i+14] = True
						if i+15 < len(chr_seq):
							array_As_chrV[i+15] = True
						if i+16 < len(chr_seq):
							array_As_chrV[i+16] = True
						if i+17 < len(chr_seq):
							array_As_chrV[i+17] = True
						if i+18 < len(chr_seq):
							array_As_chrV[i+18] = True
						if i+19 < len(chr_seq):
							array_As_chrV[i+19] = True
						if i+20 < len(chr_seq):
							array_As_chrV[i+20] = True

					#AsCas12a Revwerse PAM
					if (chr_seq[i]=='A' or chr_seq[i]=='G' or chr_seq[i]=='C') and chr_seq[i+1]=='A' and chr_seq[i+2]=='A' and chr_seq[i+3]=='A':
						As_PAM_chrV = As_PAM_chrV + 1
						if i-17 >= 0:
							array_As_chrV[i-17] = True
						if i-16 >= 0:
							array_As_chrV[i-16] = True
						if i-15 >= 0:
							array_As_chrV[i-15] = True
						if i-14 >= 0:
							array_As_chrV[i-14] = True
						if i-13 >= 0:
							array_As_chrV[i-13] = True
						if i-12 >= 0:
							array_As_chrV[i-12] = True
						if i-11 >= 0:
							array_As_chrV[i-11] = True
						if i-10 >= 0:
							array_As_chrV[i-10] = True
						if i-9 >= 0:
							array_As_chrV[i-9] = True
						if i-8 >= 0:
							array_As_chrV[i-8] = True
						if i-7 >= 0:
							array_As_chrV[i-7] = True
						if i-6 >= 0:
							array_As_chrV[i-6] = True
						if i-5 >= 0:
							array_As_chrV[i-5] = True
						if i-4 >= 0:
							array_As_chrV[i-4] = True
						if i-3 >= 0:
							array_As_chrV[i-3] = True
						if i-2 >= 0:
							array_As_chrV[i-2] = True
						if i-1 >= 0:
							array_As_chrV[i-1] = True
						if i+1 < len(chr_seq):
							array_As_chrV[i+1] = True
						if i+2 < len(chr_seq):
							array_As_chrV[i+2] = True
						if i+3 < len(chr_seq):
							array_As_chrV[i+3] = True

			for i in range(len(chr_seq)):
				array_SpSa_chrV[i] = (array_Sp_chrV[i] or array_Sa_chrV[i])
				array_SpAs_chrV[i] = (array_Sp_chrV[i] or array_As_chrV[i])
				array_SaAs_chrV[i] = (array_Sa_chrV[i] or array_As_chrV[i])
				array_SpSaAs_chrV[i] = (array_Sp_chrV[i] or array_Sa_chrV[i] or array_As_chrV[i])

			Sp_target_num_chrV = 0
			Sa_target_num_chrV = 0
			As_target_num_chrV = 0
			SpSa_target_num_chrV = 0
			SpAs_target_num_chrV = 0
			SaAs_target_num_chrV = 0
			SpSaAs_target_num_chrV = 0

			for i in range(len(chr_seq)):
				Sp_target_num_chrV = Sp_target_num_chrV + array_Sp_chrV[i]
				Sa_target_num_chrV = Sa_target_num_chrV + array_Sa_chrV[i]
				As_target_num_chrV = As_target_num_chrV + array_As_chrV[i]
				SpSa_target_num_chrV = SpSa_target_num_chrV + array_SpSa_chrV[i]
				SpAs_target_num_chrV = SpAs_target_num_chrV + array_SpAs_chrV[i]
				SaAs_target_num_chrV = SaAs_target_num_chrV + array_SaAs_chrV[i]
				SpSaAs_target_num_chrV = SpSaAs_target_num_chrV + array_SpSaAs_chrV[i]

			print('chrV length: ' + str(len(chr_seq)))
			print('chrV Sp: ' + str(int(Sp_PAM_chrV)) + ' ' + str(int(Sp_target_num_chrV)) + ' ' + str(Sp_target_num_chrV/len(chr_seq)))
			print('chrV Sa: ' + str(int(Sa_PAM_chrV)) + ' ' + str(int(Sa_target_num_chrV)) + ' ' + str(Sa_target_num_chrV/len(chr_seq)))
			print('chrV As: ' + str(int(As_PAM_chrV)) + ' ' + str(int(As_target_num_chrV)) + ' ' + str(As_target_num_chrV/len(chr_seq)))
			print('chrV SpSa: ' + str(int(SpSa_target_num_chrV)) + ' ' + str(SpSa_target_num_chrV/len(chr_seq)))
			print('chrV SpAs: ' + str(int(SpAs_target_num_chrV)) + ' ' + str(SpAs_target_num_chrV/len(chr_seq)))
			print('chrV SaAs: ' + str(int(SaAs_target_num_chrV)) + ' ' + str(SaAs_target_num_chrV/len(chr_seq)))
			print('chrV SpSaAs: ' + str(int(SpSaAs_target_num_chrV)) + ' ' + str(SpSaAs_target_num_chrV/len(chr_seq)))

		elif record.id == "chrVI":
			chr_seq = record.seq
			Sp_PAM_chrVI = 0
			Sa_PAM_chrVI = 0
			As_PAM_chrVI = 0
			for i in range(len(chr_seq)):
				#SpCas9
				if i+1 < len(chr_seq):
					#SpCas9 Forward PAM
					if chr_seq[i]=='G' and chr_seq[i+1]=='G':
						Sp_PAM_chrVI = Sp_PAM_chrVI + 1
						if i-12 >=0:
							array_Sp_chrVI[i-12] = True
						if i-11 >=0:
							array_Sp_chrVI[i-11] = True
						if i-10 >=0:
							array_Sp_chrVI[i-10] = True
						if i-9 >= 0:
							array_Sp_chrVI[i-9] = True
						if i-8 >= 0:
							array_Sp_chrVI[i-8] = True
						if i-7 >= 0:
							array_Sp_chrVI[i-7] = True
						if i-6 >= 0:
							array_Sp_chrVI[i-6] = True
						if i-5 >= 0:
							array_Sp_chrVI[i-5] = True
						if i-4 >= 0:
							array_Sp_chrVI[i-4] = True
						if i-3 >= 0:
							array_Sp_chrVI[i-3] = True
						if i-2 >= 0:
							array_Sp_chrVI[i-2] = True
						#if i-1 >= 0:
						#	array_Sp_chrVI[i-1] = True
						if i >= 0:
							array_Sp_chrVI[i] = True
						if i+1 < len(chr_seq):
							array_Sp_chrVI[i+1] = True

					#SpCas9 Reverse PAM
					if chr_seq[i]=='C' and chr_seq[i+1]=='C':
						Sp_PAM_chrVI = Sp_PAM_chrVI + 1
						if i < len(chr_seq):
							array_Sp_chrVI[i] = True
						if i+1 < len(chr_seq):
							array_Sp_chrVI[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sp_chrVI[i+2] = True
						if i+3 < len(chr_seq):
							array_Sp_chrVI[i+3] = True
						if i+4 < len(chr_seq):
							array_Sp_chrVI[i+4] = True
						if i+5 < len(chr_seq):
							array_Sp_chrVI[i+5] = True
						if i+6 < len(chr_seq):
							array_Sp_chrVI[i+6] = True
						if i+7 < len(chr_seq):
							array_Sp_chrVI[i+7] = True
						if i+8 < len(chr_seq):
							array_Sp_chrVI[i+8] = True
						if i+9 < len(chr_seq):
							array_Sp_chrVI[i+9] = True
						if i+10 < len(chr_seq):
							array_Sp_chrVI[i+10] = True
						if i+11 < len(chr_seq):
							array_Sp_chrVI[i+11] = True
						if i+12 < len(chr_seq):
							array_Sp_chrVI[i+12] = True
						if i+13 < len(chr_seq):
							array_Sp_chrVI[i+13] = True

				#SaCas9
				if i+3 < len(chr_seq):
					#SaCas9 Forward PAM
					if chr_seq[i]=='G' and (chr_seq[i+1]=='A' or chr_seq[i+1]=='G') and (chr_seq[i+2]=='A' or chr_seq[i+2]=='G') and chr_seq[i+3]=='T':
						Sa_PAM_chrVI = Sa_PAM_chrVI + 1
						if i-13 >= 0:
							array_Sa_chrVI[i-13] = True
						if i-12 >= 0:
							array_Sa_chrVI[i-12] = True
						if i-11 >= 0:
							array_Sa_chrVI[i-11] = True
						if i-10 >= 0:
							array_Sa_chrVI[i-10] = True
						if i-9 >= 0:
							array_Sa_chrVI[i-9] = True
						if i-8 >= 0:
							array_Sa_chrVI[i-8] = True
						if i-7 >= 0:
							array_Sa_chrVI[i-7] = True
						if i-6 >= 0:
							array_Sa_chrVI[i-6] = True
						if i-5 >= 0:
							array_Sa_chrVI[i-5] = True
						if i-4 >= 0:
							array_Sa_chrVI[i-4] = True
						if i-3 >= 0:
							array_Sa_chrVI[i-3] = True
						#if i-2 >= 0:
						#	array_Sa_chrVI[i-2] = True
						#if i-1 >= 0:
						#	array_Sa_chrVI[i-1] = True
						if i >= 0:
							array_Sa_chrVI[i] = True
						#if i+1 < len(chr_seq):
						#	array_Sa_chrVI[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sa_chrVI[i+2] = True
						if i+3 < len(chr_seq):
							array_Sa_chrVI[i+3] = True

					#SaCas9 Reverse PAM
					if chr_seq[i]=='A' and (chr_seq[i+1]=='C' or chr_seq[i+1]=='T') and (chr_seq[i+2]=='C' or chr_seq[i+2]=='T') and chr_seq[i+3]=='C':
						Sa_PAM_chrVI = Sa_PAM_chrVI + 1
						if i < len(chr_seq):
							array_Sa_chrVI[i] = True
						#if i+1 < len(chr_seq):
						#	array_Sa_chrVI[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sa_chrVI[i+2] = True
						if i+3 < len(chr_seq):
							array_Sa_chrVI[i+3] = True
						#if i+4 < len(chr_seq):
						#	array_Sa_chrVI[i+4] = True
						#if i+5 < len(chr_seq):
						#	array_Sa_chrVI[i+5] = True
						if i+6 < len(chr_seq):
							array_Sa_chrVI[i+6] = True
						if i+7 < len(chr_seq):
							array_Sa_chrVI[i+7] = True
						if i+8 < len(chr_seq):
							array_Sa_chrVI[i+8] = True
						if i+9 < len(chr_seq):
							array_Sa_chrVI[i+9] = True
						if i+10 < len(chr_seq):
							array_Sa_chrVI[i+10] = True
						if i+11 < len(chr_seq):
							array_Sa_chrVI[i+11] = True
						if i+12 < len(chr_seq):
							array_Sa_chrVI[i+12] = True
						if i+13 < len(chr_seq):
							array_Sa_chrVI[i+13] = True
						if i+14 < len(chr_seq):
							array_Sa_chrVI[i+14] = True
						if i+15 < len(chr_seq):
							array_Sa_chrVI[i+15] = True
						if i+16 < len(chr_seq):
							array_Sa_chrVI[i+16] = True

				#AsCas12a
				if i+3 < len(chr_seq):
					#AsCas12 Forward PAM
					if chr_seq[i]=='T' and chr_seq[i+1]=='T' and chr_seq[i+2]=='T' and (chr_seq[i+3]=='A' or chr_seq[i+3]=='G' or chr_seq[i+3]=='C'):
						As_PAM_chrVI = As_PAM_chrVI + 1
						if i < len(chr_seq):
							array_As_chrVI[i] = True
						if i+1 < len(chr_seq):
							array_As_chrVI[i+1] = True
						if i+2 < len(chr_seq):
							array_As_chrVI[i+2] = True
						if i+4 < len(chr_seq):
							array_As_chrVI[i+4] = True
						if i+5 < len(chr_seq):
							array_As_chrVI[i+5] = True
						if i+6 < len(chr_seq):
							array_As_chrVI[i+6] = True
						if i+7 < len(chr_seq):
							array_As_chrVI[i+7] = True
						if i+8 < len(chr_seq):
							array_As_chrVI[i+8] = True
						if i+9 < len(chr_seq):
							array_As_chrVI[i+9] = True
						if i+10 < len(chr_seq):
							array_As_chrVI[i+10] = True
						if i+11 < len(chr_seq):
							array_As_chrVI[i+11] = True
						if i+12 < len(chr_seq):
							array_As_chrVI[i+12] = True
						if i+13 < len(chr_seq):
							array_As_chrVI[i+13] = True
						if i+14 < len(chr_seq):
							array_As_chrVI[i+14] = True
						if i+15 < len(chr_seq):
							array_As_chrVI[i+15] = True
						if i+16 < len(chr_seq):
							array_As_chrVI[i+16] = True
						if i+17 < len(chr_seq):
							array_As_chrVI[i+17] = True
						if i+18 < len(chr_seq):
							array_As_chrVI[i+18] = True
						if i+19 < len(chr_seq):
							array_As_chrVI[i+19] = True
						if i+20 < len(chr_seq):
							array_As_chrVI[i+20] = True

					#AsCas12a Revwerse PAM
					if (chr_seq[i]=='A' or chr_seq[i]=='G' or chr_seq[i]=='C') and chr_seq[i+1]=='A' and chr_seq[i+2]=='A' and chr_seq[i+3]=='A':
						As_PAM_chrVI = As_PAM_chrVI + 1
						if i-17 >= 0:
							array_As_chrVI[i-17] = True
						if i-16 >= 0:
							array_As_chrVI[i-16] = True
						if i-15 >= 0:
							array_As_chrVI[i-15] = True
						if i-14 >= 0:
							array_As_chrVI[i-14] = True
						if i-13 >= 0:
							array_As_chrVI[i-13] = True
						if i-12 >= 0:
							array_As_chrVI[i-12] = True
						if i-11 >= 0:
							array_As_chrVI[i-11] = True
						if i-10 >= 0:
							array_As_chrVI[i-10] = True
						if i-9 >= 0:
							array_As_chrVI[i-9] = True
						if i-8 >= 0:
							array_As_chrVI[i-8] = True
						if i-7 >= 0:
							array_As_chrVI[i-7] = True
						if i-6 >= 0:
							array_As_chrVI[i-6] = True
						if i-5 >= 0:
							array_As_chrVI[i-5] = True
						if i-4 >= 0:
							array_As_chrVI[i-4] = True
						if i-3 >= 0:
							array_As_chrVI[i-3] = True
						if i-2 >= 0:
							array_As_chrVI[i-2] = True
						if i-1 >= 0:
							array_As_chrVI[i-1] = True
						if i+1 < len(chr_seq):
							array_As_chrVI[i+1] = True
						if i+2 < len(chr_seq):
							array_As_chrVI[i+2] = True
						if i+3 < len(chr_seq):
							array_As_chrVI[i+3] = True

			for i in range(len(chr_seq)):
				array_SpSa_chrVI[i] = (array_Sp_chrVI[i] or array_Sa_chrVI[i])
				array_SpAs_chrVI[i] = (array_Sp_chrVI[i] or array_As_chrVI[i])
				array_SaAs_chrVI[i] = (array_Sa_chrVI[i] or array_As_chrVI[i])
				array_SpSaAs_chrVI[i] = (array_Sp_chrVI[i] or array_Sa_chrVI[i] or array_As_chrVI[i])

			Sp_target_num_chrVI = 0
			Sa_target_num_chrVI = 0
			As_target_num_chrVI = 0
			SpSa_target_num_chrVI = 0
			SpAs_target_num_chrVI = 0
			SaAs_target_num_chrVI = 0
			SpSaAs_target_num_chrVI = 0

			for i in range(len(chr_seq)):
				Sp_target_num_chrVI = Sp_target_num_chrVI + array_Sp_chrVI[i]
				Sa_target_num_chrVI = Sa_target_num_chrVI + array_Sa_chrVI[i]
				As_target_num_chrVI = As_target_num_chrVI + array_As_chrVI[i]
				SpSa_target_num_chrVI = SpSa_target_num_chrVI + array_SpSa_chrVI[i]
				SpAs_target_num_chrVI = SpAs_target_num_chrVI + array_SpAs_chrVI[i]
				SaAs_target_num_chrVI = SaAs_target_num_chrVI + array_SaAs_chrVI[i]
				SpSaAs_target_num_chrVI = SpSaAs_target_num_chrVI + array_SpSaAs_chrVI[i]

			print('chrVI length: ' + str(len(chr_seq)))
			print('chrVI Sp: ' + str(int(Sp_PAM_chrVI)) + ' ' + str(int(Sp_target_num_chrVI)) + ' ' + str(Sp_target_num_chrVI/len(chr_seq)))
			print('chrVI Sa: ' + str(int(Sa_PAM_chrVI)) + ' ' + str(int(Sa_target_num_chrVI)) + ' ' + str(Sa_target_num_chrVI/len(chr_seq)))
			print('chrVI As: ' + str(int(As_PAM_chrVI)) + ' ' + str(int(As_target_num_chrVI)) + ' ' + str(As_target_num_chrVI/len(chr_seq)))
			print('chrVI SpSa: ' + str(int(SpSa_target_num_chrVI)) + ' ' + str(SpSa_target_num_chrVI/len(chr_seq)))
			print('chrVI SpAs: ' + str(int(SpAs_target_num_chrVI)) + ' ' + str(SpAs_target_num_chrVI/len(chr_seq)))
			print('chrVI SaAs: ' + str(int(SaAs_target_num_chrVI)) + ' ' + str(SaAs_target_num_chrVI/len(chr_seq)))
			print('chrVI SpSaAs: ' + str(int(SpSaAs_target_num_chrVI)) + ' ' + str(SpSaAs_target_num_chrVI/len(chr_seq)))

		elif record.id == "chrVII":
			chr_seq = record.seq
			Sp_PAM_chrVII = 0
			Sa_PAM_chrVII = 0
			As_PAM_chrVII = 0
			for i in range(len(chr_seq)):
								#SpCas9
				if i+1 < len(chr_seq):
					#SpCas9 Forward PAM
					if chr_seq[i]=='G' and chr_seq[i+1]=='G':
						Sp_PAM_chrVII = Sp_PAM_chrVII + 1
						if i-12 >=0:
							array_Sp_chrVII[i-12] = True
						if i-11 >=0:
							array_Sp_chrVII[i-11] = True
						if i-10 >=0:
							array_Sp_chrVII[i-10] = True
						if i-9 >= 0:
							array_Sp_chrVII[i-9] = True
						if i-8 >= 0:
							array_Sp_chrVII[i-8] = True
						if i-7 >= 0:
							array_Sp_chrVII[i-7] = True
						if i-6 >= 0:
							array_Sp_chrVII[i-6] = True
						if i-5 >= 0:
							array_Sp_chrVII[i-5] = True
						if i-4 >= 0:
							array_Sp_chrVII[i-4] = True
						if i-3 >= 0:
							array_Sp_chrVII[i-3] = True
						if i-2 >= 0:
							array_Sp_chrVII[i-2] = True
						#if i-1 >= 0:
						#	array_Sp_chrVII[i-1] = True
						if i >= 0:
							array_Sp_chrVII[i] = True
						if i+1 < len(chr_seq):
							array_Sp_chrVII[i+1] = True

					#SpCas9 Reverse PAM
					if chr_seq[i]=='C' and chr_seq[i+1]=='C':
						Sp_PAM_chrVII = Sp_PAM_chrVII + 1
						if i < len(chr_seq):
							array_Sp_chrVII[i] = True
						if i+1 < len(chr_seq):
							array_Sp_chrVII[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sp_chrVII[i+2] = True
						if i+3 < len(chr_seq):
							array_Sp_chrVII[i+3] = True
						if i+4 < len(chr_seq):
							array_Sp_chrVII[i+4] = True
						if i+5 < len(chr_seq):
							array_Sp_chrVII[i+5] = True
						if i+6 < len(chr_seq):
							array_Sp_chrVII[i+6] = True
						if i+7 < len(chr_seq):
							array_Sp_chrVII[i+7] = True
						if i+8 < len(chr_seq):
							array_Sp_chrVII[i+8] = True
						if i+9 < len(chr_seq):
							array_Sp_chrVII[i+9] = True
						if i+10 < len(chr_seq):
							array_Sp_chrVII[i+10] = True
						if i+11 < len(chr_seq):
							array_Sp_chrVII[i+11] = True
						if i+12 < len(chr_seq):
							array_Sp_chrVII[i+12] = True
						if i+13 < len(chr_seq):
							array_Sp_chrVII[i+13] = True

				#SaCas9
				if i+3 < len(chr_seq):
					#SaCas9 Forward PAM
					if chr_seq[i]=='G' and (chr_seq[i+1]=='A' or chr_seq[i+1]=='G') and (chr_seq[i+2]=='A' or chr_seq[i+2]=='G') and chr_seq[i+3]=='T':
						Sa_PAM_chrVII = Sa_PAM_chrVII + 1
						if i-13 >= 0:
							array_Sa_chrVII[i-13] = True
						if i-12 >= 0:
							array_Sa_chrVII[i-12] = True
						if i-11 >= 0:
							array_Sa_chrVII[i-11] = True
						if i-10 >= 0:
							array_Sa_chrVII[i-10] = True
						if i-9 >= 0:
							array_Sa_chrVII[i-9] = True
						if i-8 >= 0:
							array_Sa_chrVII[i-8] = True
						if i-7 >= 0:
							array_Sa_chrVII[i-7] = True
						if i-6 >= 0:
							array_Sa_chrVII[i-6] = True
						if i-5 >= 0:
							array_Sa_chrVII[i-5] = True
						if i-4 >= 0:
							array_Sa_chrVII[i-4] = True
						if i-3 >= 0:
							array_Sa_chrVII[i-3] = True
						#if i-2 >= 0:
						#	array_Sa_chrVII[i-2] = True
						#if i-1 >= 0:
						#	array_Sa_chrVII[i-1] = True
						if i >= 0:
							array_Sa_chrVII[i] = True
						#if i+1 < len(chr_seq):
						#	array_Sa_chrVII[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sa_chrVII[i+2] = True
						if i+3 < len(chr_seq):
							array_Sa_chrVII[i+3] = True

					#SaCas9 Reverse PAM
					if chr_seq[i]=='A' and (chr_seq[i+1]=='C' or chr_seq[i+1]=='T') and (chr_seq[i+2]=='C' or chr_seq[i+2]=='T') and chr_seq[i+3]=='C':
						Sa_PAM_chrVII = Sa_PAM_chrVII + 1
						if i < len(chr_seq):
							array_Sa_chrVII[i] = True
						#if i+1 < len(chr_seq):
						#	array_Sa_chrVII[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sa_chrVII[i+2] = True
						if i+3 < len(chr_seq):
							array_Sa_chrVII[i+3] = True
						#if i+4 < len(chr_seq):
						#	array_Sa_chrVII[i+4] = True
						#if i+5 < len(chr_seq):
						#	array_Sa_chrVII[i+5] = True
						if i+6 < len(chr_seq):
							array_Sa_chrVII[i+6] = True
						if i+7 < len(chr_seq):
							array_Sa_chrVII[i+7] = True
						if i+8 < len(chr_seq):
							array_Sa_chrVII[i+8] = True
						if i+9 < len(chr_seq):
							array_Sa_chrVII[i+9] = True
						if i+10 < len(chr_seq):
							array_Sa_chrVII[i+10] = True
						if i+11 < len(chr_seq):
							array_Sa_chrVII[i+11] = True
						if i+12 < len(chr_seq):
							array_Sa_chrVII[i+12] = True
						if i+13 < len(chr_seq):
							array_Sa_chrVII[i+13] = True
						if i+14 < len(chr_seq):
							array_Sa_chrVII[i+14] = True
						if i+15 < len(chr_seq):
							array_Sa_chrVII[i+15] = True
						if i+16 < len(chr_seq):
							array_Sa_chrVII[i+16] = True

				#AsCas12a
				if i+3 < len(chr_seq):
					#AsCas12 Forward PAM
					if chr_seq[i]=='T' and chr_seq[i+1]=='T' and chr_seq[i+2]=='T' and (chr_seq[i+3]=='A' or chr_seq[i+3]=='G' or chr_seq[i+3]=='C'):
						As_PAM_chrVII = As_PAM_chrVII + 1
						if i < len(chr_seq):
							array_As_chrVII[i] = True
						if i+1 < len(chr_seq):
							array_As_chrVII[i+1] = True
						if i+2 < len(chr_seq):
							array_As_chrVII[i+2] = True
						if i+4 < len(chr_seq):
							array_As_chrVII[i+4] = True
						if i+5 < len(chr_seq):
							array_As_chrVII[i+5] = True
						if i+6 < len(chr_seq):
							array_As_chrVII[i+6] = True
						if i+7 < len(chr_seq):
							array_As_chrVII[i+7] = True
						if i+8 < len(chr_seq):
							array_As_chrVII[i+8] = True
						if i+9 < len(chr_seq):
							array_As_chrVII[i+9] = True
						if i+10 < len(chr_seq):
							array_As_chrVII[i+10] = True
						if i+11 < len(chr_seq):
							array_As_chrVII[i+11] = True
						if i+12 < len(chr_seq):
							array_As_chrVII[i+12] = True
						if i+13 < len(chr_seq):
							array_As_chrVII[i+13] = True
						if i+14 < len(chr_seq):
							array_As_chrVII[i+14] = True
						if i+15 < len(chr_seq):
							array_As_chrVII[i+15] = True
						if i+16 < len(chr_seq):
							array_As_chrVII[i+16] = True
						if i+17 < len(chr_seq):
							array_As_chrVII[i+17] = True
						if i+18 < len(chr_seq):
							array_As_chrVII[i+18] = True
						if i+19 < len(chr_seq):
							array_As_chrVII[i+19] = True
						if i+20 < len(chr_seq):
							array_As_chrVII[i+20] = True

					#AsCas12a Revwerse PAM
					if (chr_seq[i]=='A' or chr_seq[i]=='G' or chr_seq[i]=='C') and chr_seq[i+1]=='A' and chr_seq[i+2]=='A' and chr_seq[i+3]=='A':
						As_PAM_chrVII = As_PAM_chrVII + 1
						if i-17 >= 0:
							array_As_chrVII[i-17] = True
						if i-16 >= 0:
							array_As_chrVII[i-16] = True
						if i-15 >= 0:
							array_As_chrVII[i-15] = True
						if i-14 >= 0:
							array_As_chrVII[i-14] = True
						if i-13 >= 0:
							array_As_chrVII[i-13] = True
						if i-12 >= 0:
							array_As_chrVII[i-12] = True
						if i-11 >= 0:
							array_As_chrVII[i-11] = True
						if i-10 >= 0:
							array_As_chrVII[i-10] = True
						if i-9 >= 0:
							array_As_chrVII[i-9] = True
						if i-8 >= 0:
							array_As_chrVII[i-8] = True
						if i-7 >= 0:
							array_As_chrVII[i-7] = True
						if i-6 >= 0:
							array_As_chrVII[i-6] = True
						if i-5 >= 0:
							array_As_chrVII[i-5] = True
						if i-4 >= 0:
							array_As_chrVII[i-4] = True
						if i-3 >= 0:
							array_As_chrVII[i-3] = True
						if i-2 >= 0:
							array_As_chrVII[i-2] = True
						if i-1 >= 0:
							array_As_chrVII[i-1] = True
						if i+1 < len(chr_seq):
							array_As_chrVII[i+1] = True
						if i+2 < len(chr_seq):
							array_As_chrVII[i+2] = True
						if i+3 < len(chr_seq):
							array_As_chrVII[i+3] = True

			for i in range(len(chr_seq)):
				array_SpSa_chrVII[i] = (array_Sp_chrVII[i] or array_Sa_chrVII[i])
				array_SpAs_chrVII[i] = (array_Sp_chrVII[i] or array_As_chrVII[i])
				array_SaAs_chrVII[i] = (array_Sa_chrVII[i] or array_As_chrVII[i])
				array_SpSaAs_chrVII[i] = (array_Sp_chrVII[i] or array_Sa_chrVII[i] or array_As_chrVII[i])

			Sp_target_num_chrVII = 0
			Sa_target_num_chrVII = 0
			As_target_num_chrVII = 0
			SpSa_target_num_chrVII = 0
			SpAs_target_num_chrVII = 0
			SaAs_target_num_chrVII = 0
			SpSaAs_target_num_chrVII = 0

			for i in range(len(chr_seq)):
				Sp_target_num_chrVII = Sp_target_num_chrVII + array_Sp_chrVII[i]
				Sa_target_num_chrVII = Sa_target_num_chrVII + array_Sa_chrVII[i]
				As_target_num_chrVII = As_target_num_chrVII + array_As_chrVII[i]
				SpSa_target_num_chrVII = SpSa_target_num_chrVII + array_SpSa_chrVII[i]
				SpAs_target_num_chrVII = SpAs_target_num_chrVII + array_SpAs_chrVII[i]
				SaAs_target_num_chrVII = SaAs_target_num_chrVII + array_SaAs_chrVII[i]
				SpSaAs_target_num_chrVII = SpSaAs_target_num_chrVII + array_SpSaAs_chrVII[i]

			print('chrVII length: ' + str(len(chr_seq)))
			print('chrVII Sp: ' + str(int(Sp_PAM_chrVII)) + ' ' + str(int(Sp_target_num_chrVII)) + ' ' + str(Sp_target_num_chrVII/len(chr_seq)))
			print('chrVII Sa: ' + str(int(Sa_PAM_chrVII)) + ' ' + str(int(Sa_target_num_chrVII)) + ' ' + str(Sa_target_num_chrVII/len(chr_seq)))
			print('chrVII As: ' + str(int(As_PAM_chrVII)) + ' ' + str(int(As_target_num_chrVII)) + ' ' + str(As_target_num_chrVII/len(chr_seq)))
			print('chrVII SpSa: ' + str(int(SpSa_target_num_chrVII)) + ' ' + str(SpSa_target_num_chrVII/len(chr_seq)))
			print('chrVII SpAs: ' + str(int(SpAs_target_num_chrVII)) + ' ' + str(SpAs_target_num_chrVII/len(chr_seq)))
			print('chrVII SaAs: ' + str(int(SaAs_target_num_chrVII)) + ' ' + str(SaAs_target_num_chrVII/len(chr_seq)))
			print('chrVII SpSaAs: ' + str(int(SpSaAs_target_num_chrVII)) + ' ' + str(SpSaAs_target_num_chrVII/len(chr_seq)))

		elif record.id == "chrVIII":
			chr_seq = record.seq
			Sp_PAM_chrVIII = 0
			Sa_PAM_chrVIII = 0
			As_PAM_chrVIII = 0
			for i in range(len(chr_seq)):
								#SpCas9
				if i+1 < len(chr_seq):
					#SpCas9 Forward PAM
					if chr_seq[i]=='G' and chr_seq[i+1]=='G':
						Sp_PAM_chrVIII = Sp_PAM_chrVIII + 1
						if i-12 >=0:
							array_Sp_chrVIII[i-12] = True
						if i-11 >=0:
							array_Sp_chrVIII[i-11] = True
						if i-10 >=0:
							array_Sp_chrVIII[i-10] = True
						if i-9 >= 0:
							array_Sp_chrVIII[i-9] = True
						if i-8 >= 0:
							array_Sp_chrVIII[i-8] = True
						if i-7 >= 0:
							array_Sp_chrVIII[i-7] = True
						if i-6 >= 0:
							array_Sp_chrVIII[i-6] = True
						if i-5 >= 0:
							array_Sp_chrVIII[i-5] = True
						if i-4 >= 0:
							array_Sp_chrVIII[i-4] = True
						if i-3 >= 0:
							array_Sp_chrVIII[i-3] = True
						if i-2 >= 0:
							array_Sp_chrVIII[i-2] = True
						#if i-1 >= 0:
						#	array_Sp_chrVIII[i-1] = True
						if i >= 0:
							array_Sp_chrVIII[i] = True
						if i+1 < len(chr_seq):
							array_Sp_chrVIII[i+1] = True

					#SpCas9 Reverse PAM
					if chr_seq[i]=='C' and chr_seq[i+1]=='C':
						Sp_PAM_chrVIII = Sp_PAM_chrVIII + 1
						if i < len(chr_seq):
							array_Sp_chrVIII[i] = True
						if i+1 < len(chr_seq):
							array_Sp_chrVIII[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sp_chrVIII[i+2] = True
						if i+3 < len(chr_seq):
							array_Sp_chrVIII[i+3] = True
						if i+4 < len(chr_seq):
							array_Sp_chrVIII[i+4] = True
						if i+5 < len(chr_seq):
							array_Sp_chrVIII[i+5] = True
						if i+6 < len(chr_seq):
							array_Sp_chrVIII[i+6] = True
						if i+7 < len(chr_seq):
							array_Sp_chrVIII[i+7] = True
						if i+8 < len(chr_seq):
							array_Sp_chrVIII[i+8] = True
						if i+9 < len(chr_seq):
							array_Sp_chrVIII[i+9] = True
						if i+10 < len(chr_seq):
							array_Sp_chrVIII[i+10] = True
						if i+11 < len(chr_seq):
							array_Sp_chrVIII[i+11] = True
						if i+12 < len(chr_seq):
							array_Sp_chrVIII[i+12] = True
						if i+13 < len(chr_seq):
							array_Sp_chrVIII[i+13] = True

				#SaCas9
				if i+3 < len(chr_seq):
					#SaCas9 Forward PAM
					if chr_seq[i]=='G' and (chr_seq[i+1]=='A' or chr_seq[i+1]=='G') and (chr_seq[i+2]=='A' or chr_seq[i+2]=='G') and chr_seq[i+3]=='T':
						Sa_PAM_chrVIII = Sa_PAM_chrVIII + 1
						if i-13 >= 0:
							array_Sa_chrVIII[i-13] = True
						if i-12 >= 0:
							array_Sa_chrVIII[i-12] = True
						if i-11 >= 0:
							array_Sa_chrVIII[i-11] = True
						if i-10 >= 0:
							array_Sa_chrVIII[i-10] = True
						if i-9 >= 0:
							array_Sa_chrVIII[i-9] = True
						if i-8 >= 0:
							array_Sa_chrVIII[i-8] = True
						if i-7 >= 0:
							array_Sa_chrVIII[i-7] = True
						if i-6 >= 0:
							array_Sa_chrVIII[i-6] = True
						if i-5 >= 0:
							array_Sa_chrVIII[i-5] = True
						if i-4 >= 0:
							array_Sa_chrVIII[i-4] = True
						if i-3 >= 0:
							array_Sa_chrVIII[i-3] = True
						#if i-2 >= 0:
						#	array_Sa_chrVIII[i-2] = True
						#if i-1 >= 0:
						#	array_Sa_chrVIII[i-1] = True
						if i >= 0:
							array_Sa_chrVIII[i] = True
						#if i+1 < len(chr_seq):
						#	array_Sa_chrVIII[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sa_chrVIII[i+2] = True
						if i+3 < len(chr_seq):
							array_Sa_chrVIII[i+3] = True

					#SaCas9 Reverse PAM
					if chr_seq[i]=='A' and (chr_seq[i+1]=='C' or chr_seq[i+1]=='T') and (chr_seq[i+2]=='C' or chr_seq[i+2]=='T') and chr_seq[i+3]=='C':
						Sa_PAM_chrVIII = Sa_PAM_chrVIII + 1
						if i < len(chr_seq):
							array_Sa_chrVIII[i] = True
						#if i+1 < len(chr_seq):
						#	array_Sa_chrVIII[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sa_chrVIII[i+2] = True
						if i+3 < len(chr_seq):
							array_Sa_chrVIII[i+3] = True
						#if i+4 < len(chr_seq):
						#	array_Sa_chrVIII[i+4] = True
						#if i+5 < len(chr_seq):
						#	array_Sa_chrVIII[i+5] = True
						if i+6 < len(chr_seq):
							array_Sa_chrVIII[i+6] = True
						if i+7 < len(chr_seq):
							array_Sa_chrVIII[i+7] = True
						if i+8 < len(chr_seq):
							array_Sa_chrVIII[i+8] = True
						if i+9 < len(chr_seq):
							array_Sa_chrVIII[i+9] = True
						if i+10 < len(chr_seq):
							array_Sa_chrVIII[i+10] = True
						if i+11 < len(chr_seq):
							array_Sa_chrVIII[i+11] = True
						if i+12 < len(chr_seq):
							array_Sa_chrVIII[i+12] = True
						if i+13 < len(chr_seq):
							array_Sa_chrVIII[i+13] = True
						if i+14 < len(chr_seq):
							array_Sa_chrVIII[i+14] = True
						if i+15 < len(chr_seq):
							array_Sa_chrVIII[i+15] = True
						if i+16 < len(chr_seq):
							array_Sa_chrVIII[i+16] = True

				#AsCas12a
				if i+3 < len(chr_seq):
					#AsCas12 Forward PAM
					if chr_seq[i]=='T' and chr_seq[i+1]=='T' and chr_seq[i+2]=='T' and (chr_seq[i+3]=='A' or chr_seq[i+3]=='G' or chr_seq[i+3]=='C'):
						As_PAM_chrVIII = As_PAM_chrVIII + 1
						if i < len(chr_seq):
							array_As_chrVIII[i] = True
						if i+1 < len(chr_seq):
							array_As_chrVIII[i+1] = True
						if i+2 < len(chr_seq):
							array_As_chrVIII[i+2] = True
						if i+4 < len(chr_seq):
							array_As_chrVIII[i+4] = True
						if i+5 < len(chr_seq):
							array_As_chrVIII[i+5] = True
						if i+6 < len(chr_seq):
							array_As_chrVIII[i+6] = True
						if i+7 < len(chr_seq):
							array_As_chrVIII[i+7] = True
						if i+8 < len(chr_seq):
							array_As_chrVIII[i+8] = True
						if i+9 < len(chr_seq):
							array_As_chrVIII[i+9] = True
						if i+10 < len(chr_seq):
							array_As_chrVIII[i+10] = True
						if i+11 < len(chr_seq):
							array_As_chrVIII[i+11] = True
						if i+12 < len(chr_seq):
							array_As_chrVIII[i+12] = True
						if i+13 < len(chr_seq):
							array_As_chrVIII[i+13] = True
						if i+14 < len(chr_seq):
							array_As_chrVIII[i+14] = True
						if i+15 < len(chr_seq):
							array_As_chrVIII[i+15] = True
						if i+16 < len(chr_seq):
							array_As_chrVIII[i+16] = True
						if i+17 < len(chr_seq):
							array_As_chrVIII[i+17] = True
						if i+18 < len(chr_seq):
							array_As_chrVIII[i+18] = True
						if i+19 < len(chr_seq):
							array_As_chrVIII[i+19] = True
						if i+20 < len(chr_seq):
							array_As_chrVIII[i+20] = True

					#AsCas12a Revwerse PAM
					if (chr_seq[i]=='A' or chr_seq[i]=='G' or chr_seq[i]=='C') and chr_seq[i+1]=='A' and chr_seq[i+2]=='A' and chr_seq[i+3]=='A':
						As_PAM_chrVIII = As_PAM_chrVIII + 1
						if i-17 >= 0:
							array_As_chrVIII[i-17] = True
						if i-16 >= 0:
							array_As_chrVIII[i-16] = True
						if i-15 >= 0:
							array_As_chrVIII[i-15] = True
						if i-14 >= 0:
							array_As_chrVIII[i-14] = True
						if i-13 >= 0:
							array_As_chrVIII[i-13] = True
						if i-12 >= 0:
							array_As_chrVIII[i-12] = True
						if i-11 >= 0:
							array_As_chrVIII[i-11] = True
						if i-10 >= 0:
							array_As_chrVIII[i-10] = True
						if i-9 >= 0:
							array_As_chrVIII[i-9] = True
						if i-8 >= 0:
							array_As_chrVIII[i-8] = True
						if i-7 >= 0:
							array_As_chrVIII[i-7] = True
						if i-6 >= 0:
							array_As_chrVIII[i-6] = True
						if i-5 >= 0:
							array_As_chrVIII[i-5] = True
						if i-4 >= 0:
							array_As_chrVIII[i-4] = True
						if i-3 >= 0:
							array_As_chrVIII[i-3] = True
						if i-2 >= 0:
							array_As_chrVIII[i-2] = True
						if i-1 >= 0:
							array_As_chrVIII[i-1] = True
						if i+1 < len(chr_seq):
							array_As_chrVIII[i+1] = True
						if i+2 < len(chr_seq):
							array_As_chrVIII[i+2] = True
						if i+3 < len(chr_seq):
							array_As_chrVIII[i+3] = True

			for i in range(len(chr_seq)):
				array_SpSa_chrVIII[i] = (array_Sp_chrVIII[i] or array_Sa_chrVIII[i])
				array_SpAs_chrVIII[i] = (array_Sp_chrVIII[i] or array_As_chrVIII[i])
				array_SaAs_chrVIII[i] = (array_Sa_chrVIII[i] or array_As_chrVIII[i])
				array_SpSaAs_chrVIII[i] = (array_Sp_chrVIII[i] or array_Sa_chrVIII[i] or array_As_chrVIII[i])

			Sp_target_num_chrVIII = 0
			Sa_target_num_chrVIII = 0
			As_target_num_chrVIII = 0
			SpSa_target_num_chrVIII = 0
			SpAs_target_num_chrVIII = 0
			SaAs_target_num_chrVIII = 0
			SpSaAs_target_num_chrVIII = 0

			for i in range(len(chr_seq)):
				Sp_target_num_chrVIII = Sp_target_num_chrVIII + array_Sp_chrVIII[i]
				Sa_target_num_chrVIII = Sa_target_num_chrVIII + array_Sa_chrVIII[i]
				As_target_num_chrVIII = As_target_num_chrVIII + array_As_chrVIII[i]
				SpSa_target_num_chrVIII = SpSa_target_num_chrVIII + array_SpSa_chrVIII[i]
				SpAs_target_num_chrVIII = SpAs_target_num_chrVIII + array_SpAs_chrVIII[i]
				SaAs_target_num_chrVIII = SaAs_target_num_chrVIII + array_SaAs_chrVIII[i]
				SpSaAs_target_num_chrVIII = SpSaAs_target_num_chrVIII + array_SpSaAs_chrVIII[i]

			print('chrVIII length: ' + str(len(chr_seq)))
			print('chrVIII Sp: ' + str(int(Sp_PAM_chrVIII)) + ' ' + str(int(Sp_target_num_chrVIII)) + ' ' + str(Sp_target_num_chrVIII/len(chr_seq)))
			print('chrVIII Sa: ' + str(int(Sa_PAM_chrVIII)) + ' ' + str(int(Sa_target_num_chrVIII)) + ' ' + str(Sa_target_num_chrVIII/len(chr_seq)))
			print('chrVIII As: ' + str(int(As_PAM_chrVIII)) + ' ' + str(int(As_target_num_chrVIII)) + ' ' + str(As_target_num_chrVIII/len(chr_seq)))
			print('chrVIII SpSa: ' + str(int(SpSa_target_num_chrVIII)) + ' ' + str(SpSa_target_num_chrVIII/len(chr_seq)))
			print('chrVIII SpAs: ' + str(int(SpAs_target_num_chrVIII)) + ' ' + str(SpAs_target_num_chrVIII/len(chr_seq)))
			print('chrVIII SaAs: ' + str(int(SaAs_target_num_chrVIII)) + ' ' + str(SaAs_target_num_chrVIII/len(chr_seq)))
			print('chrVIII SpSaAs: ' + str(int(SpSaAs_target_num_chrVIII)) + ' ' + str(SpSaAs_target_num_chrVIII/len(chr_seq)))

		elif record.id == "chrIX":
			chr_seq = record.seq
			Sp_PAM_chrIX = 0
			Sa_PAM_chrIX = 0
			As_PAM_chrIX = 0
			for i in range(len(chr_seq)):
								#SpCas9
				if i+1 < len(chr_seq):
					#SpCas9 Forward PAM
					if chr_seq[i]=='G' and chr_seq[i+1]=='G':
						Sp_PAM_chrIX = Sp_PAM_chrIX + 1
						if i-12 >=0:
							array_Sp_chrIX[i-12] = True
						if i-11 >=0:
							array_Sp_chrIX[i-11] = True
						if i-10 >=0:
							array_Sp_chrIX[i-10] = True
						if i-9 >= 0:
							array_Sp_chrIX[i-9] = True
						if i-8 >= 0:
							array_Sp_chrIX[i-8] = True
						if i-7 >= 0:
							array_Sp_chrIX[i-7] = True
						if i-6 >= 0:
							array_Sp_chrIX[i-6] = True
						if i-5 >= 0:
							array_Sp_chrIX[i-5] = True
						if i-4 >= 0:
							array_Sp_chrIX[i-4] = True
						if i-3 >= 0:
							array_Sp_chrIX[i-3] = True
						if i-2 >= 0:
							array_Sp_chrIX[i-2] = True
						#if i-1 >= 0:
						#	array_Sp_chrIX[i-1] = True
						if i >= 0:
							array_Sp_chrIX[i] = True
						if i+1 < len(chr_seq):
							array_Sp_chrIX[i+1] = True

					#SpCas9 Reverse PAM
					if chr_seq[i]=='C' and chr_seq[i+1]=='C':
						Sp_PAM_chrIX = Sp_PAM_chrIX + 1
						if i < len(chr_seq):
							array_Sp_chrIX[i] = True
						if i+1 < len(chr_seq):
							array_Sp_chrIX[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sp_chrIX[i+2] = True
						if i+3 < len(chr_seq):
							array_Sp_chrIX[i+3] = True
						if i+4 < len(chr_seq):
							array_Sp_chrIX[i+4] = True
						if i+5 < len(chr_seq):
							array_Sp_chrIX[i+5] = True
						if i+6 < len(chr_seq):
							array_Sp_chrIX[i+6] = True
						if i+7 < len(chr_seq):
							array_Sp_chrIX[i+7] = True
						if i+8 < len(chr_seq):
							array_Sp_chrIX[i+8] = True
						if i+9 < len(chr_seq):
							array_Sp_chrIX[i+9] = True
						if i+10 < len(chr_seq):
							array_Sp_chrIX[i+10] = True
						if i+11 < len(chr_seq):
							array_Sp_chrIX[i+11] = True
						if i+12 < len(chr_seq):
							array_Sp_chrIX[i+12] = True
						if i+13 < len(chr_seq):
							array_Sp_chrIX[i+13] = True

				#SaCas9
				if i+3 < len(chr_seq):
					#SaCas9 Forward PAM
					if chr_seq[i]=='G' and (chr_seq[i+1]=='A' or chr_seq[i+1]=='G') and (chr_seq[i+2]=='A' or chr_seq[i+2]=='G') and chr_seq[i+3]=='T':
						Sa_PAM_chrIX = Sa_PAM_chrIX + 1
						if i-13 >= 0:
							array_Sa_chrIX[i-13] = True
						if i-12 >= 0:
							array_Sa_chrIX[i-12] = True
						if i-11 >= 0:
							array_Sa_chrIX[i-11] = True
						if i-10 >= 0:
							array_Sa_chrIX[i-10] = True
						if i-9 >= 0:
							array_Sa_chrIX[i-9] = True
						if i-8 >= 0:
							array_Sa_chrIX[i-8] = True
						if i-7 >= 0:
							array_Sa_chrIX[i-7] = True
						if i-6 >= 0:
							array_Sa_chrIX[i-6] = True
						if i-5 >= 0:
							array_Sa_chrIX[i-5] = True
						if i-4 >= 0:
							array_Sa_chrIX[i-4] = True
						if i-3 >= 0:
							array_Sa_chrIX[i-3] = True
						#if i-2 >= 0:
						#	array_Sa_chrIX[i-2] = True
						#if i-1 >= 0:
						#	array_Sa_chrIX[i-1] = True
						if i >= 0:
							array_Sa_chrIX[i] = True
						#if i+1 < len(chr_seq):
						#	array_Sa_chrIX[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sa_chrIX[i+2] = True
						if i+3 < len(chr_seq):
							array_Sa_chrIX[i+3] = True

					#SaCas9 Reverse PAM
					if chr_seq[i]=='A' and (chr_seq[i+1]=='C' or chr_seq[i+1]=='T') and (chr_seq[i+2]=='C' or chr_seq[i+2]=='T') and chr_seq[i+3]=='C':
						Sa_PAM_chrIX = Sa_PAM_chrIX + 1
						if i < len(chr_seq):
							array_Sa_chrIX[i] = True
						#if i+1 < len(chr_seq):
						#	array_Sa_chrIX[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sa_chrIX[i+2] = True
						if i+3 < len(chr_seq):
							array_Sa_chrIX[i+3] = True
						#if i+4 < len(chr_seq):
						#	array_Sa_chrIX[i+4] = True
						#if i+5 < len(chr_seq):
						#	array_Sa_chrIX[i+5] = True
						if i+6 < len(chr_seq):
							array_Sa_chrIX[i+6] = True
						if i+7 < len(chr_seq):
							array_Sa_chrIX[i+7] = True
						if i+8 < len(chr_seq):
							array_Sa_chrIX[i+8] = True
						if i+9 < len(chr_seq):
							array_Sa_chrIX[i+9] = True
						if i+10 < len(chr_seq):
							array_Sa_chrIX[i+10] = True
						if i+11 < len(chr_seq):
							array_Sa_chrIX[i+11] = True
						if i+12 < len(chr_seq):
							array_Sa_chrIX[i+12] = True
						if i+13 < len(chr_seq):
							array_Sa_chrIX[i+13] = True
						if i+14 < len(chr_seq):
							array_Sa_chrIX[i+14] = True
						if i+15 < len(chr_seq):
							array_Sa_chrIX[i+15] = True
						if i+16 < len(chr_seq):
							array_Sa_chrIX[i+16] = True

				#AsCas12a
				if i+3 < len(chr_seq):
					#AsCas12 Forward PAM
					if chr_seq[i]=='T' and chr_seq[i+1]=='T' and chr_seq[i+2]=='T' and (chr_seq[i+3]=='A' or chr_seq[i+3]=='G' or chr_seq[i+3]=='C'):
						As_PAM_chrIX = As_PAM_chrIX + 1
						if i < len(chr_seq):
							array_As_chrIX[i] = True
						if i+1 < len(chr_seq):
							array_As_chrIX[i+1] = True
						if i+2 < len(chr_seq):
							array_As_chrIX[i+2] = True
						if i+4 < len(chr_seq):
							array_As_chrIX[i+4] = True
						if i+5 < len(chr_seq):
							array_As_chrIX[i+5] = True
						if i+6 < len(chr_seq):
							array_As_chrIX[i+6] = True
						if i+7 < len(chr_seq):
							array_As_chrIX[i+7] = True
						if i+8 < len(chr_seq):
							array_As_chrIX[i+8] = True
						if i+9 < len(chr_seq):
							array_As_chrIX[i+9] = True
						if i+10 < len(chr_seq):
							array_As_chrIX[i+10] = True
						if i+11 < len(chr_seq):
							array_As_chrIX[i+11] = True
						if i+12 < len(chr_seq):
							array_As_chrIX[i+12] = True
						if i+13 < len(chr_seq):
							array_As_chrIX[i+13] = True
						if i+14 < len(chr_seq):
							array_As_chrIX[i+14] = True
						if i+15 < len(chr_seq):
							array_As_chrIX[i+15] = True
						if i+16 < len(chr_seq):
							array_As_chrIX[i+16] = True
						if i+17 < len(chr_seq):
							array_As_chrIX[i+17] = True
						if i+18 < len(chr_seq):
							array_As_chrIX[i+18] = True
						if i+19 < len(chr_seq):
							array_As_chrIX[i+19] = True
						if i+20 < len(chr_seq):
							array_As_chrIX[i+20] = True

					#AsCas12a Revwerse PAM
					if (chr_seq[i]=='A' or chr_seq[i]=='G' or chr_seq[i]=='C') and chr_seq[i+1]=='A' and chr_seq[i+2]=='A' and chr_seq[i+3]=='A':
						As_PAM_chrIX = As_PAM_chrIX + 1
						if i-17 >= 0:
							array_As_chrIX[i-17] = True
						if i-16 >= 0:
							array_As_chrIX[i-16] = True
						if i-15 >= 0:
							array_As_chrIX[i-15] = True
						if i-14 >= 0:
							array_As_chrIX[i-14] = True
						if i-13 >= 0:
							array_As_chrIX[i-13] = True
						if i-12 >= 0:
							array_As_chrIX[i-12] = True
						if i-11 >= 0:
							array_As_chrIX[i-11] = True
						if i-10 >= 0:
							array_As_chrIX[i-10] = True
						if i-9 >= 0:
							array_As_chrIX[i-9] = True
						if i-8 >= 0:
							array_As_chrIX[i-8] = True
						if i-7 >= 0:
							array_As_chrIX[i-7] = True
						if i-6 >= 0:
							array_As_chrIX[i-6] = True
						if i-5 >= 0:
							array_As_chrIX[i-5] = True
						if i-4 >= 0:
							array_As_chrIX[i-4] = True
						if i-3 >= 0:
							array_As_chrIX[i-3] = True
						if i-2 >= 0:
							array_As_chrIX[i-2] = True
						if i-1 >= 0:
							array_As_chrIX[i-1] = True
						if i+1 < len(chr_seq):
							array_As_chrIX[i+1] = True
						if i+2 < len(chr_seq):
							array_As_chrIX[i+2] = True
						if i+3 < len(chr_seq):
							array_As_chrIX[i+3] = True

			for i in range(len(chr_seq)):
				array_SpSa_chrIX[i] = (array_Sp_chrIX[i] or array_Sa_chrIX[i])
				array_SpAs_chrIX[i] = (array_Sp_chrIX[i] or array_As_chrIX[i])
				array_SaAs_chrIX[i] = (array_Sa_chrIX[i] or array_As_chrIX[i])
				array_SpSaAs_chrIX[i] = (array_Sp_chrIX[i] or array_Sa_chrIX[i] or array_As_chrIX[i])

			Sp_target_num_chrIX = 0
			Sa_target_num_chrIX = 0
			As_target_num_chrIX = 0
			SpSa_target_num_chrIX = 0
			SpAs_target_num_chrIX = 0
			SaAs_target_num_chrIX = 0
			SpSaAs_target_num_chrIX = 0

			for i in range(len(chr_seq)):
				Sp_target_num_chrIX = Sp_target_num_chrIX + array_Sp_chrIX[i]
				Sa_target_num_chrIX = Sa_target_num_chrIX + array_Sa_chrIX[i]
				As_target_num_chrIX = As_target_num_chrIX + array_As_chrIX[i]
				SpSa_target_num_chrIX = SpSa_target_num_chrIX + array_SpSa_chrIX[i]
				SpAs_target_num_chrIX = SpAs_target_num_chrIX + array_SpAs_chrIX[i]
				SaAs_target_num_chrIX = SaAs_target_num_chrIX + array_SaAs_chrIX[i]
				SpSaAs_target_num_chrIX = SpSaAs_target_num_chrIX + array_SpSaAs_chrIX[i]

			print('chrIX length: ' + str(len(chr_seq)))
			print('chrIX Sp: ' + str(int(Sp_PAM_chrIX)) + ' ' + str(int(Sp_target_num_chrIX)) + ' ' + str(Sp_target_num_chrIX/len(chr_seq)))
			print('chrIX Sa: ' + str(int(Sa_PAM_chrIX)) + ' ' + str(int(Sa_target_num_chrIX)) + ' ' + str(Sa_target_num_chrIX/len(chr_seq)))
			print('chrIX As: ' + str(int(As_PAM_chrIX)) + ' ' + str(int(As_target_num_chrIX)) + ' ' + str(As_target_num_chrIX/len(chr_seq)))
			print('chrIX SpSa: ' + str(int(SpSa_target_num_chrIX)) + ' ' + str(SpSa_target_num_chrIX/len(chr_seq)))
			print('chrIX SpAs: ' + str(int(SpAs_target_num_chrIX)) + ' ' + str(SpAs_target_num_chrIX/len(chr_seq)))
			print('chrIX SaAs: ' + str(int(SaAs_target_num_chrIX)) + ' ' + str(SaAs_target_num_chrIX/len(chr_seq)))
			print('chrIX SpSaAs: ' + str(int(SpSaAs_target_num_chrIX)) + ' ' + str(SpSaAs_target_num_chrIX/len(chr_seq)))

		elif record.id == "chrX":
			chr_seq = record.seq
			Sp_PAM_chrX = 0
			Sa_PAM_chrX = 0
			As_PAM_chrX = 0
			for i in range(len(chr_seq)):
								#SpCas9
				if i+1 < len(chr_seq):
					#SpCas9 Forward PAM
					if chr_seq[i]=='G' and chr_seq[i+1]=='G':
						Sp_PAM_chrX = Sp_PAM_chrX + 1
						if i-12 >=0:
							array_Sp_chrX[i-12] = True
						if i-11 >=0:
							array_Sp_chrX[i-11] = True
						if i-10 >=0:
							array_Sp_chrX[i-10] = True
						if i-9 >= 0:
							array_Sp_chrX[i-9] = True
						if i-8 >= 0:
							array_Sp_chrX[i-8] = True
						if i-7 >= 0:
							array_Sp_chrX[i-7] = True
						if i-6 >= 0:
							array_Sp_chrX[i-6] = True
						if i-5 >= 0:
							array_Sp_chrX[i-5] = True
						if i-4 >= 0:
							array_Sp_chrX[i-4] = True
						if i-3 >= 0:
							array_Sp_chrX[i-3] = True
						if i-2 >= 0:
							array_Sp_chrX[i-2] = True
						#if i-1 >= 0:
						#	array_Sp_chrX[i-1] = True
						if i >= 0:
							array_Sp_chrX[i] = True
						if i+1 < len(chr_seq):
							array_Sp_chrX[i+1] = True

					#SpCas9 Reverse PAM
					if chr_seq[i]=='C' and chr_seq[i+1]=='C':
						Sp_PAM_chrX = Sp_PAM_chrX + 1
						if i < len(chr_seq):
							array_Sp_chrX[i] = True
						if i+1 < len(chr_seq):
							array_Sp_chrX[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sp_chrX[i+2] = True
						if i+3 < len(chr_seq):
							array_Sp_chrX[i+3] = True
						if i+4 < len(chr_seq):
							array_Sp_chrX[i+4] = True
						if i+5 < len(chr_seq):
							array_Sp_chrX[i+5] = True
						if i+6 < len(chr_seq):
							array_Sp_chrX[i+6] = True
						if i+7 < len(chr_seq):
							array_Sp_chrX[i+7] = True
						if i+8 < len(chr_seq):
							array_Sp_chrX[i+8] = True
						if i+9 < len(chr_seq):
							array_Sp_chrX[i+9] = True
						if i+10 < len(chr_seq):
							array_Sp_chrX[i+10] = True
						if i+11 < len(chr_seq):
							array_Sp_chrX[i+11] = True
						if i+12 < len(chr_seq):
							array_Sp_chrX[i+12] = True
						if i+13 < len(chr_seq):
							array_Sp_chrX[i+13] = True

				#SaCas9
				if i+3 < len(chr_seq):
					#SaCas9 Forward PAM
					if chr_seq[i]=='G' and (chr_seq[i+1]=='A' or chr_seq[i+1]=='G') and (chr_seq[i+2]=='A' or chr_seq[i+2]=='G') and chr_seq[i+3]=='T':
						Sa_PAM_chrX = Sa_PAM_chrX + 1
						if i-13 >= 0:
							array_Sa_chrX[i-13] = True
						if i-12 >= 0:
							array_Sa_chrX[i-12] = True
						if i-11 >= 0:
							array_Sa_chrX[i-11] = True
						if i-10 >= 0:
							array_Sa_chrX[i-10] = True
						if i-9 >= 0:
							array_Sa_chrX[i-9] = True
						if i-8 >= 0:
							array_Sa_chrX[i-8] = True
						if i-7 >= 0:
							array_Sa_chrX[i-7] = True
						if i-6 >= 0:
							array_Sa_chrX[i-6] = True
						if i-5 >= 0:
							array_Sa_chrX[i-5] = True
						if i-4 >= 0:
							array_Sa_chrX[i-4] = True
						if i-3 >= 0:
							array_Sa_chrX[i-3] = True
						#if i-2 >= 0:
						#	array_Sa_chrX[i-2] = True
						#if i-1 >= 0:
						#	array_Sa_chrX[i-1] = True
						if i >= 0:
							array_Sa_chrX[i] = True
						#if i+1 < len(chr_seq):
						#	array_Sa_chrX[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sa_chrX[i+2] = True
						if i+3 < len(chr_seq):
							array_Sa_chrX[i+3] = True

					#SaCas9 Reverse PAM
					if chr_seq[i]=='A' and (chr_seq[i+1]=='C' or chr_seq[i+1]=='T') and (chr_seq[i+2]=='C' or chr_seq[i+2]=='T') and chr_seq[i+3]=='C':
						Sa_PAM_chrX = Sa_PAM_chrX + 1
						if i < len(chr_seq):
							array_Sa_chrX[i] = True
						#if i+1 < len(chr_seq):
						#	array_Sa_chrX[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sa_chrX[i+2] = True
						if i+3 < len(chr_seq):
							array_Sa_chrX[i+3] = True
						#if i+4 < len(chr_seq):
						#	array_Sa_chrX[i+4] = True
						#if i+5 < len(chr_seq):
						#	array_Sa_chrX[i+5] = True
						if i+6 < len(chr_seq):
							array_Sa_chrX[i+6] = True
						if i+7 < len(chr_seq):
							array_Sa_chrX[i+7] = True
						if i+8 < len(chr_seq):
							array_Sa_chrX[i+8] = True
						if i+9 < len(chr_seq):
							array_Sa_chrX[i+9] = True
						if i+10 < len(chr_seq):
							array_Sa_chrX[i+10] = True
						if i+11 < len(chr_seq):
							array_Sa_chrX[i+11] = True
						if i+12 < len(chr_seq):
							array_Sa_chrX[i+12] = True
						if i+13 < len(chr_seq):
							array_Sa_chrX[i+13] = True
						if i+14 < len(chr_seq):
							array_Sa_chrX[i+14] = True
						if i+15 < len(chr_seq):
							array_Sa_chrX[i+15] = True
						if i+16 < len(chr_seq):
							array_Sa_chrX[i+16] = True

				#AsCas12a
				if i+3 < len(chr_seq):
					#AsCas12 Forward PAM
					if chr_seq[i]=='T' and chr_seq[i+1]=='T' and chr_seq[i+2]=='T' and (chr_seq[i+3]=='A' or chr_seq[i+3]=='G' or chr_seq[i+3]=='C'):
						As_PAM_chrX = As_PAM_chrX + 1
						if i < len(chr_seq):
							array_As_chrX[i] = True
						if i+1 < len(chr_seq):
							array_As_chrX[i+1] = True
						if i+2 < len(chr_seq):
							array_As_chrX[i+2] = True
						if i+4 < len(chr_seq):
							array_As_chrX[i+4] = True
						if i+5 < len(chr_seq):
							array_As_chrX[i+5] = True
						if i+6 < len(chr_seq):
							array_As_chrX[i+6] = True
						if i+7 < len(chr_seq):
							array_As_chrX[i+7] = True
						if i+8 < len(chr_seq):
							array_As_chrX[i+8] = True
						if i+9 < len(chr_seq):
							array_As_chrX[i+9] = True
						if i+10 < len(chr_seq):
							array_As_chrX[i+10] = True
						if i+11 < len(chr_seq):
							array_As_chrX[i+11] = True
						if i+12 < len(chr_seq):
							array_As_chrX[i+12] = True
						if i+13 < len(chr_seq):
							array_As_chrX[i+13] = True
						if i+14 < len(chr_seq):
							array_As_chrX[i+14] = True
						if i+15 < len(chr_seq):
							array_As_chrX[i+15] = True
						if i+16 < len(chr_seq):
							array_As_chrX[i+16] = True
						if i+17 < len(chr_seq):
							array_As_chrX[i+17] = True
						if i+18 < len(chr_seq):
							array_As_chrX[i+18] = True
						if i+19 < len(chr_seq):
							array_As_chrX[i+19] = True
						if i+20 < len(chr_seq):
							array_As_chrX[i+20] = True

					#AsCas12a Revwerse PAM
					if (chr_seq[i]=='A' or chr_seq[i]=='G' or chr_seq[i]=='C') and chr_seq[i+1]=='A' and chr_seq[i+2]=='A' and chr_seq[i+3]=='A':
						As_PAM_chrX = As_PAM_chrX + 1
						if i-17 >= 0:
							array_As_chrX[i-17] = True
						if i-16 >= 0:
							array_As_chrX[i-16] = True
						if i-15 >= 0:
							array_As_chrX[i-15] = True
						if i-14 >= 0:
							array_As_chrX[i-14] = True
						if i-13 >= 0:
							array_As_chrX[i-13] = True
						if i-12 >= 0:
							array_As_chrX[i-12] = True
						if i-11 >= 0:
							array_As_chrX[i-11] = True
						if i-10 >= 0:
							array_As_chrX[i-10] = True
						if i-9 >= 0:
							array_As_chrX[i-9] = True
						if i-8 >= 0:
							array_As_chrX[i-8] = True
						if i-7 >= 0:
							array_As_chrX[i-7] = True
						if i-6 >= 0:
							array_As_chrX[i-6] = True
						if i-5 >= 0:
							array_As_chrX[i-5] = True
						if i-4 >= 0:
							array_As_chrX[i-4] = True
						if i-3 >= 0:
							array_As_chrX[i-3] = True
						if i-2 >= 0:
							array_As_chrX[i-2] = True
						if i-1 >= 0:
							array_As_chrX[i-1] = True
						if i+1 < len(chr_seq):
							array_As_chrX[i+1] = True
						if i+2 < len(chr_seq):
							array_As_chrX[i+2] = True
						if i+3 < len(chr_seq):
							array_As_chrX[i+3] = True

			for i in range(len(chr_seq)):
				array_SpSa_chrX[i] = (array_Sp_chrX[i] or array_Sa_chrX[i])
				array_SpAs_chrX[i] = (array_Sp_chrX[i] or array_As_chrX[i])
				array_SaAs_chrX[i] = (array_Sa_chrX[i] or array_As_chrX[i])
				array_SpSaAs_chrX[i] = (array_Sp_chrX[i] or array_Sa_chrX[i] or array_As_chrX[i])

			Sp_target_num_chrX = 0
			Sa_target_num_chrX = 0
			As_target_num_chrX = 0
			SpSa_target_num_chrX = 0
			SpAs_target_num_chrX = 0
			SaAs_target_num_chrX = 0
			SpSaAs_target_num_chrX = 0

			for i in range(len(chr_seq)):
				Sp_target_num_chrX = Sp_target_num_chrX + array_Sp_chrX[i]
				Sa_target_num_chrX = Sa_target_num_chrX + array_Sa_chrX[i]
				As_target_num_chrX = As_target_num_chrX + array_As_chrX[i]
				SpSa_target_num_chrX = SpSa_target_num_chrX + array_SpSa_chrX[i]
				SpAs_target_num_chrX = SpAs_target_num_chrX + array_SpAs_chrX[i]
				SaAs_target_num_chrX = SaAs_target_num_chrX + array_SaAs_chrX[i]
				SpSaAs_target_num_chrX = SpSaAs_target_num_chrX + array_SpSaAs_chrX[i]

			print('chrX length: ' + str(len(chr_seq)))
			print('chrX Sp: ' + str(int(Sp_PAM_chrX)) + ' ' + str(int(Sp_target_num_chrX)) + ' ' + str(Sp_target_num_chrX/len(chr_seq)))
			print('chrX Sa: ' + str(int(Sa_PAM_chrX)) + ' ' + str(int(Sa_target_num_chrX)) + ' ' + str(Sa_target_num_chrX/len(chr_seq)))
			print('chrX As: ' + str(int(As_PAM_chrX)) + ' ' + str(int(As_target_num_chrX)) + ' ' + str(As_target_num_chrX/len(chr_seq)))
			print('chrX SpSa: ' + str(int(SpSa_target_num_chrX)) + ' ' + str(SpSa_target_num_chrX/len(chr_seq)))
			print('chrX SpAs: ' + str(int(SpAs_target_num_chrX)) + ' ' + str(SpAs_target_num_chrX/len(chr_seq)))
			print('chrX SaAs: ' + str(int(SaAs_target_num_chrX)) + ' ' + str(SaAs_target_num_chrX/len(chr_seq)))
			print('chrX SpSaAs: ' + str(int(SpSaAs_target_num_chrX)) + ' ' + str(SpSaAs_target_num_chrX/len(chr_seq)))

		elif record.id == "chrXI":
			chr_seq = record.seq
			Sp_PAM_chrXI = 0
			Sa_PAM_chrXI = 0
			As_PAM_chrXI = 0
			for i in range(len(chr_seq)):
								#SpCas9
				if i+1 < len(chr_seq):
					#SpCas9 Forward PAM
					if chr_seq[i]=='G' and chr_seq[i+1]=='G':
						Sp_PAM_chrXI = Sp_PAM_chrXI + 1
						if i-12 >=0:
							array_Sp_chrXI[i-12] = True
						if i-11 >=0:
							array_Sp_chrXI[i-11] = True
						if i-10 >=0:
							array_Sp_chrXI[i-10] = True
						if i-9 >= 0:
							array_Sp_chrXI[i-9] = True
						if i-8 >= 0:
							array_Sp_chrXI[i-8] = True
						if i-7 >= 0:
							array_Sp_chrXI[i-7] = True
						if i-6 >= 0:
							array_Sp_chrXI[i-6] = True
						if i-5 >= 0:
							array_Sp_chrXI[i-5] = True
						if i-4 >= 0:
							array_Sp_chrXI[i-4] = True
						if i-3 >= 0:
							array_Sp_chrXI[i-3] = True
						if i-2 >= 0:
							array_Sp_chrXI[i-2] = True
						#if i-1 >= 0:
						#	array_Sp_chrXI[i-1] = True
						if i >= 0:
							array_Sp_chrXI[i] = True
						if i+1 < len(chr_seq):
							array_Sp_chrXI[i+1] = True

					#SpCas9 Reverse PAM
					if chr_seq[i]=='C' and chr_seq[i+1]=='C':
						Sp_PAM_chrXI = Sp_PAM_chrXI + 1
						if i < len(chr_seq):
							array_Sp_chrXI[i] = True
						if i+1 < len(chr_seq):
							array_Sp_chrXI[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sp_chrXI[i+2] = True
						if i+3 < len(chr_seq):
							array_Sp_chrXI[i+3] = True
						if i+4 < len(chr_seq):
							array_Sp_chrXI[i+4] = True
						if i+5 < len(chr_seq):
							array_Sp_chrXI[i+5] = True
						if i+6 < len(chr_seq):
							array_Sp_chrXI[i+6] = True
						if i+7 < len(chr_seq):
							array_Sp_chrXI[i+7] = True
						if i+8 < len(chr_seq):
							array_Sp_chrXI[i+8] = True
						if i+9 < len(chr_seq):
							array_Sp_chrXI[i+9] = True
						if i+10 < len(chr_seq):
							array_Sp_chrXI[i+10] = True
						if i+11 < len(chr_seq):
							array_Sp_chrXI[i+11] = True
						if i+12 < len(chr_seq):
							array_Sp_chrXI[i+12] = True
						if i+13 < len(chr_seq):
							array_Sp_chrXI[i+13] = True

				#SaCas9
				if i+3 < len(chr_seq):
					#SaCas9 Forward PAM
					if chr_seq[i]=='G' and (chr_seq[i+1]=='A' or chr_seq[i+1]=='G') and (chr_seq[i+2]=='A' or chr_seq[i+2]=='G') and chr_seq[i+3]=='T':
						Sa_PAM_chrXI = Sa_PAM_chrXI + 1
						if i-13 >= 0:
							array_Sa_chrXI[i-13] = True
						if i-12 >= 0:
							array_Sa_chrXI[i-12] = True
						if i-11 >= 0:
							array_Sa_chrXI[i-11] = True
						if i-10 >= 0:
							array_Sa_chrXI[i-10] = True
						if i-9 >= 0:
							array_Sa_chrXI[i-9] = True
						if i-8 >= 0:
							array_Sa_chrXI[i-8] = True
						if i-7 >= 0:
							array_Sa_chrXI[i-7] = True
						if i-6 >= 0:
							array_Sa_chrXI[i-6] = True
						if i-5 >= 0:
							array_Sa_chrXI[i-5] = True
						if i-4 >= 0:
							array_Sa_chrXI[i-4] = True
						if i-3 >= 0:
							array_Sa_chrXI[i-3] = True
						#if i-2 >= 0:
						#	array_Sa_chrXI[i-2] = True
						#if i-1 >= 0:
						#	array_Sa_chrXI[i-1] = True
						if i >= 0:
							array_Sa_chrXI[i] = True
						#if i+1 < len(chr_seq):
						#	array_Sa_chrXI[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sa_chrXI[i+2] = True
						if i+3 < len(chr_seq):
							array_Sa_chrXI[i+3] = True

					#SaCas9 Reverse PAM
					if chr_seq[i]=='A' and (chr_seq[i+1]=='C' or chr_seq[i+1]=='T') and (chr_seq[i+2]=='C' or chr_seq[i+2]=='T') and chr_seq[i+3]=='C':
						Sa_PAM_chrXI = Sa_PAM_chrXI + 1
						if i < len(chr_seq):
							array_Sa_chrXI[i] = True
						#if i+1 < len(chr_seq):
						#	array_Sa_chrXI[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sa_chrXI[i+2] = True
						if i+3 < len(chr_seq):
							array_Sa_chrXI[i+3] = True
						#if i+4 < len(chr_seq):
						#	array_Sa_chrXI[i+4] = True
						#if i+5 < len(chr_seq):
						#	array_Sa_chrXI[i+5] = True
						if i+6 < len(chr_seq):
							array_Sa_chrXI[i+6] = True
						if i+7 < len(chr_seq):
							array_Sa_chrXI[i+7] = True
						if i+8 < len(chr_seq):
							array_Sa_chrXI[i+8] = True
						if i+9 < len(chr_seq):
							array_Sa_chrXI[i+9] = True
						if i+10 < len(chr_seq):
							array_Sa_chrXI[i+10] = True
						if i+11 < len(chr_seq):
							array_Sa_chrXI[i+11] = True
						if i+12 < len(chr_seq):
							array_Sa_chrXI[i+12] = True
						if i+13 < len(chr_seq):
							array_Sa_chrXI[i+13] = True
						if i+14 < len(chr_seq):
							array_Sa_chrXI[i+14] = True
						if i+15 < len(chr_seq):
							array_Sa_chrXI[i+15] = True
						if i+16 < len(chr_seq):
							array_Sa_chrXI[i+16] = True

				#AsCas12a
				if i+3 < len(chr_seq):
					#AsCas12 Forward PAM
					if chr_seq[i]=='T' and chr_seq[i+1]=='T' and chr_seq[i+2]=='T' and (chr_seq[i+3]=='A' or chr_seq[i+3]=='G' or chr_seq[i+3]=='C'):
						As_PAM_chrXI = As_PAM_chrXI + 1
						if i < len(chr_seq):
							array_As_chrXI[i] = True
						if i+1 < len(chr_seq):
							array_As_chrXI[i+1] = True
						if i+2 < len(chr_seq):
							array_As_chrXI[i+2] = True
						if i+4 < len(chr_seq):
							array_As_chrXI[i+4] = True
						if i+5 < len(chr_seq):
							array_As_chrXI[i+5] = True
						if i+6 < len(chr_seq):
							array_As_chrXI[i+6] = True
						if i+7 < len(chr_seq):
							array_As_chrXI[i+7] = True
						if i+8 < len(chr_seq):
							array_As_chrXI[i+8] = True
						if i+9 < len(chr_seq):
							array_As_chrXI[i+9] = True
						if i+10 < len(chr_seq):
							array_As_chrXI[i+10] = True
						if i+11 < len(chr_seq):
							array_As_chrXI[i+11] = True
						if i+12 < len(chr_seq):
							array_As_chrXI[i+12] = True
						if i+13 < len(chr_seq):
							array_As_chrXI[i+13] = True
						if i+14 < len(chr_seq):
							array_As_chrXI[i+14] = True
						if i+15 < len(chr_seq):
							array_As_chrXI[i+15] = True
						if i+16 < len(chr_seq):
							array_As_chrXI[i+16] = True
						if i+17 < len(chr_seq):
							array_As_chrXI[i+17] = True
						if i+18 < len(chr_seq):
							array_As_chrXI[i+18] = True
						if i+19 < len(chr_seq):
							array_As_chrXI[i+19] = True
						if i+20 < len(chr_seq):
							array_As_chrXI[i+20] = True

					#AsCas12a Revwerse PAM
					if (chr_seq[i]=='A' or chr_seq[i]=='G' or chr_seq[i]=='C') and chr_seq[i+1]=='A' and chr_seq[i+2]=='A' and chr_seq[i+3]=='A':
						As_PAM_chrXI = As_PAM_chrXI + 1
						if i-17 >= 0:
							array_As_chrXI[i-17] = True
						if i-16 >= 0:
							array_As_chrXI[i-16] = True
						if i-15 >= 0:
							array_As_chrXI[i-15] = True
						if i-14 >= 0:
							array_As_chrXI[i-14] = True
						if i-13 >= 0:
							array_As_chrXI[i-13] = True
						if i-12 >= 0:
							array_As_chrXI[i-12] = True
						if i-11 >= 0:
							array_As_chrXI[i-11] = True
						if i-10 >= 0:
							array_As_chrXI[i-10] = True
						if i-9 >= 0:
							array_As_chrXI[i-9] = True
						if i-8 >= 0:
							array_As_chrXI[i-8] = True
						if i-7 >= 0:
							array_As_chrXI[i-7] = True
						if i-6 >= 0:
							array_As_chrXI[i-6] = True
						if i-5 >= 0:
							array_As_chrXI[i-5] = True
						if i-4 >= 0:
							array_As_chrXI[i-4] = True
						if i-3 >= 0:
							array_As_chrXI[i-3] = True
						if i-2 >= 0:
							array_As_chrXI[i-2] = True
						if i-1 >= 0:
							array_As_chrXI[i-1] = True
						if i+1 < len(chr_seq):
							array_As_chrXI[i+1] = True
						if i+2 < len(chr_seq):
							array_As_chrXI[i+2] = True
						if i+3 < len(chr_seq):
							array_As_chrXI[i+3] = True

			for i in range(len(chr_seq)):
				array_SpSa_chrXI[i] = (array_Sp_chrXI[i] or array_Sa_chrXI[i])
				array_SpAs_chrXI[i] = (array_Sp_chrXI[i] or array_As_chrXI[i])
				array_SaAs_chrXI[i] = (array_Sa_chrXI[i] or array_As_chrXI[i])
				array_SpSaAs_chrXI[i] = (array_Sp_chrXI[i] or array_Sa_chrXI[i] or array_As_chrXI[i])

			Sp_target_num_chrXI = 0
			Sa_target_num_chrXI = 0
			As_target_num_chrXI = 0
			SpSa_target_num_chrXI = 0
			SpAs_target_num_chrXI = 0
			SaAs_target_num_chrXI = 0
			SpSaAs_target_num_chrXI = 0

			for i in range(len(chr_seq)):
				Sp_target_num_chrXI = Sp_target_num_chrXI + array_Sp_chrXI[i]
				Sa_target_num_chrXI = Sa_target_num_chrXI + array_Sa_chrXI[i]
				As_target_num_chrXI = As_target_num_chrXI + array_As_chrXI[i]
				SpSa_target_num_chrXI = SpSa_target_num_chrXI + array_SpSa_chrXI[i]
				SpAs_target_num_chrXI = SpAs_target_num_chrXI + array_SpAs_chrXI[i]
				SaAs_target_num_chrXI = SaAs_target_num_chrXI + array_SaAs_chrXI[i]
				SpSaAs_target_num_chrXI = SpSaAs_target_num_chrXI + array_SpSaAs_chrXI[i]

			print('chrXI length: ' + str(len(chr_seq)))
			print('chrXI Sp: ' + str(int(Sp_PAM_chrXI)) + ' ' + str(int(Sp_target_num_chrXI)) + ' ' + str(Sp_target_num_chrXI/len(chr_seq)))
			print('chrXI Sa: ' + str(int(Sa_PAM_chrXI)) + ' ' + str(int(Sa_target_num_chrXI)) + ' ' + str(Sa_target_num_chrXI/len(chr_seq)))
			print('chrXI As: ' + str(int(As_PAM_chrXI)) + ' ' + str(int(As_target_num_chrXI)) + ' ' + str(As_target_num_chrXI/len(chr_seq)))
			print('chrXI SpSa: ' + str(int(SpSa_target_num_chrXI)) + ' ' + str(SpSa_target_num_chrXI/len(chr_seq)))
			print('chrXI SpAs: ' + str(int(SpAs_target_num_chrXI)) + ' ' + str(SpAs_target_num_chrXI/len(chr_seq)))
			print('chrXI SaAs: ' + str(int(SaAs_target_num_chrXI)) + ' ' + str(SaAs_target_num_chrXI/len(chr_seq)))
			print('chrXI SpSaAs: ' + str(int(SpSaAs_target_num_chrXI)) + ' ' + str(SpSaAs_target_num_chrXI/len(chr_seq)))

		elif record.id == "chrXII":
			chr_seq = record.seq
			Sp_PAM_chrXII = 0
			Sa_PAM_chrXII = 0
			As_PAM_chrXII = 0
			for i in range(len(chr_seq)):
								#SpCas9
				if i+1 < len(chr_seq):
					#SpCas9 Forward PAM
					if chr_seq[i]=='G' and chr_seq[i+1]=='G':
						Sp_PAM_chrXII = Sp_PAM_chrXII + 1
						if i-12 >=0:
							array_Sp_chrXII[i-12] = True
						if i-11 >=0:
							array_Sp_chrXII[i-11] = True
						if i-10 >=0:
							array_Sp_chrXII[i-10] = True
						if i-9 >= 0:
							array_Sp_chrXII[i-9] = True
						if i-8 >= 0:
							array_Sp_chrXII[i-8] = True
						if i-7 >= 0:
							array_Sp_chrXII[i-7] = True
						if i-6 >= 0:
							array_Sp_chrXII[i-6] = True
						if i-5 >= 0:
							array_Sp_chrXII[i-5] = True
						if i-4 >= 0:
							array_Sp_chrXII[i-4] = True
						if i-3 >= 0:
							array_Sp_chrXII[i-3] = True
						if i-2 >= 0:
							array_Sp_chrXII[i-2] = True
						#if i-1 >= 0:
						#	array_Sp_chrXII[i-1] = True
						if i >= 0:
							array_Sp_chrXII[i] = True
						if i+1 < len(chr_seq):
							array_Sp_chrXII[i+1] = True

					#SpCas9 Reverse PAM
					if chr_seq[i]=='C' and chr_seq[i+1]=='C':
						Sp_PAM_chrXII = Sp_PAM_chrXII + 1
						if i < len(chr_seq):
							array_Sp_chrXII[i] = True
						if i+1 < len(chr_seq):
							array_Sp_chrXII[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sp_chrXII[i+2] = True
						if i+3 < len(chr_seq):
							array_Sp_chrXII[i+3] = True
						if i+4 < len(chr_seq):
							array_Sp_chrXII[i+4] = True
						if i+5 < len(chr_seq):
							array_Sp_chrXII[i+5] = True
						if i+6 < len(chr_seq):
							array_Sp_chrXII[i+6] = True
						if i+7 < len(chr_seq):
							array_Sp_chrXII[i+7] = True
						if i+8 < len(chr_seq):
							array_Sp_chrXII[i+8] = True
						if i+9 < len(chr_seq):
							array_Sp_chrXII[i+9] = True
						if i+10 < len(chr_seq):
							array_Sp_chrXII[i+10] = True
						if i+11 < len(chr_seq):
							array_Sp_chrXII[i+11] = True
						if i+12 < len(chr_seq):
							array_Sp_chrXII[i+12] = True
						if i+13 < len(chr_seq):
							array_Sp_chrXII[i+13] = True

				#SaCas9
				if i+3 < len(chr_seq):
					#SaCas9 Forward PAM
					if chr_seq[i]=='G' and (chr_seq[i+1]=='A' or chr_seq[i+1]=='G') and (chr_seq[i+2]=='A' or chr_seq[i+2]=='G') and chr_seq[i+3]=='T':
						Sa_PAM_chrXII = Sa_PAM_chrXII + 1
						if i-13 >= 0:
							array_Sa_chrXII[i-13] = True
						if i-12 >= 0:
							array_Sa_chrXII[i-12] = True
						if i-11 >= 0:
							array_Sa_chrXII[i-11] = True
						if i-10 >= 0:
							array_Sa_chrXII[i-10] = True
						if i-9 >= 0:
							array_Sa_chrXII[i-9] = True
						if i-8 >= 0:
							array_Sa_chrXII[i-8] = True
						if i-7 >= 0:
							array_Sa_chrXII[i-7] = True
						if i-6 >= 0:
							array_Sa_chrXII[i-6] = True
						if i-5 >= 0:
							array_Sa_chrXII[i-5] = True
						if i-4 >= 0:
							array_Sa_chrXII[i-4] = True
						if i-3 >= 0:
							array_Sa_chrXII[i-3] = True
						#if i-2 >= 0:
						#	array_Sa_chrXII[i-2] = True
						#if i-1 >= 0:
						#	array_Sa_chrXII[i-1] = True
						if i >= 0:
							array_Sa_chrXII[i] = True
						#if i+1 < len(chr_seq):
						#	array_Sa_chrXII[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sa_chrXII[i+2] = True
						if i+3 < len(chr_seq):
							array_Sa_chrXII[i+3] = True

					#SaCas9 Reverse PAM
					if chr_seq[i]=='A' and (chr_seq[i+1]=='C' or chr_seq[i+1]=='T') and (chr_seq[i+2]=='C' or chr_seq[i+2]=='T') and chr_seq[i+3]=='C':
						Sa_PAM_chrXII = Sa_PAM_chrXII + 1
						if i < len(chr_seq):
							array_Sa_chrXII[i] = True
						#if i+1 < len(chr_seq):
						#	array_Sa_chrXII[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sa_chrXII[i+2] = True
						if i+3 < len(chr_seq):
							array_Sa_chrXII[i+3] = True
						#if i+4 < len(chr_seq):
						#	array_Sa_chrXII[i+4] = True
						#if i+5 < len(chr_seq):
						#	array_Sa_chrXII[i+5] = True
						if i+6 < len(chr_seq):
							array_Sa_chrXII[i+6] = True
						if i+7 < len(chr_seq):
							array_Sa_chrXII[i+7] = True
						if i+8 < len(chr_seq):
							array_Sa_chrXII[i+8] = True
						if i+9 < len(chr_seq):
							array_Sa_chrXII[i+9] = True
						if i+10 < len(chr_seq):
							array_Sa_chrXII[i+10] = True
						if i+11 < len(chr_seq):
							array_Sa_chrXII[i+11] = True
						if i+12 < len(chr_seq):
							array_Sa_chrXII[i+12] = True
						if i+13 < len(chr_seq):
							array_Sa_chrXII[i+13] = True
						if i+14 < len(chr_seq):
							array_Sa_chrXII[i+14] = True
						if i+15 < len(chr_seq):
							array_Sa_chrXII[i+15] = True
						if i+16 < len(chr_seq):
							array_Sa_chrXII[i+16] = True

				#AsCas12a
				if i+3 < len(chr_seq):
					#AsCas12 Forward PAM
					if chr_seq[i]=='T' and chr_seq[i+1]=='T' and chr_seq[i+2]=='T' and (chr_seq[i+3]=='A' or chr_seq[i+3]=='G' or chr_seq[i+3]=='C'):
						As_PAM_chrXII = As_PAM_chrXII + 1
						if i < len(chr_seq):
							array_As_chrXII[i] = True
						if i+1 < len(chr_seq):
							array_As_chrXII[i+1] = True
						if i+2 < len(chr_seq):
							array_As_chrXII[i+2] = True
						if i+4 < len(chr_seq):
							array_As_chrXII[i+4] = True
						if i+5 < len(chr_seq):
							array_As_chrXII[i+5] = True
						if i+6 < len(chr_seq):
							array_As_chrXII[i+6] = True
						if i+7 < len(chr_seq):
							array_As_chrXII[i+7] = True
						if i+8 < len(chr_seq):
							array_As_chrXII[i+8] = True
						if i+9 < len(chr_seq):
							array_As_chrXII[i+9] = True
						if i+10 < len(chr_seq):
							array_As_chrXII[i+10] = True
						if i+11 < len(chr_seq):
							array_As_chrXII[i+11] = True
						if i+12 < len(chr_seq):
							array_As_chrXII[i+12] = True
						if i+13 < len(chr_seq):
							array_As_chrXII[i+13] = True
						if i+14 < len(chr_seq):
							array_As_chrXII[i+14] = True
						if i+15 < len(chr_seq):
							array_As_chrXII[i+15] = True
						if i+16 < len(chr_seq):
							array_As_chrXII[i+16] = True
						if i+17 < len(chr_seq):
							array_As_chrXII[i+17] = True
						if i+18 < len(chr_seq):
							array_As_chrXII[i+18] = True
						if i+19 < len(chr_seq):
							array_As_chrXII[i+19] = True
						if i+20 < len(chr_seq):
							array_As_chrXII[i+20] = True

					#AsCas12a Revwerse PAM
					if (chr_seq[i]=='A' or chr_seq[i]=='G' or chr_seq[i]=='C') and chr_seq[i+1]=='A' and chr_seq[i+2]=='A' and chr_seq[i+3]=='A':
						As_PAM_chrXII = As_PAM_chrXII + 1
						if i-17 >= 0:
							array_As_chrXII[i-17] = True
						if i-16 >= 0:
							array_As_chrXII[i-16] = True
						if i-15 >= 0:
							array_As_chrXII[i-15] = True
						if i-14 >= 0:
							array_As_chrXII[i-14] = True
						if i-13 >= 0:
							array_As_chrXII[i-13] = True
						if i-12 >= 0:
							array_As_chrXII[i-12] = True
						if i-11 >= 0:
							array_As_chrXII[i-11] = True
						if i-10 >= 0:
							array_As_chrXII[i-10] = True
						if i-9 >= 0:
							array_As_chrXII[i-9] = True
						if i-8 >= 0:
							array_As_chrXII[i-8] = True
						if i-7 >= 0:
							array_As_chrXII[i-7] = True
						if i-6 >= 0:
							array_As_chrXII[i-6] = True
						if i-5 >= 0:
							array_As_chrXII[i-5] = True
						if i-4 >= 0:
							array_As_chrXII[i-4] = True
						if i-3 >= 0:
							array_As_chrXII[i-3] = True
						if i-2 >= 0:
							array_As_chrXII[i-2] = True
						if i-1 >= 0:
							array_As_chrXII[i-1] = True
						if i+1 < len(chr_seq):
							array_As_chrXII[i+1] = True
						if i+2 < len(chr_seq):
							array_As_chrXII[i+2] = True
						if i+3 < len(chr_seq):
							array_As_chrXII[i+3] = True

			for i in range(len(chr_seq)):
				array_SpSa_chrXII[i] = (array_Sp_chrXII[i] or array_Sa_chrXII[i])
				array_SpAs_chrXII[i] = (array_Sp_chrXII[i] or array_As_chrXII[i])
				array_SaAs_chrXII[i] = (array_Sa_chrXII[i] or array_As_chrXII[i])
				array_SpSaAs_chrXII[i] = (array_Sp_chrXII[i] or array_Sa_chrXII[i] or array_As_chrXII[i])

			Sp_target_num_chrXII = 0
			Sa_target_num_chrXII = 0
			As_target_num_chrXII = 0
			SpSa_target_num_chrXII = 0
			SpAs_target_num_chrXII = 0
			SaAs_target_num_chrXII = 0
			SpSaAs_target_num_chrXII = 0

			for i in range(len(chr_seq)):
				Sp_target_num_chrXII = Sp_target_num_chrXII + array_Sp_chrXII[i]
				Sa_target_num_chrXII = Sa_target_num_chrXII + array_Sa_chrXII[i]
				As_target_num_chrXII = As_target_num_chrXII + array_As_chrXII[i]
				SpSa_target_num_chrXII = SpSa_target_num_chrXII + array_SpSa_chrXII[i]
				SpAs_target_num_chrXII = SpAs_target_num_chrXII + array_SpAs_chrXII[i]
				SaAs_target_num_chrXII = SaAs_target_num_chrXII + array_SaAs_chrXII[i]
				SpSaAs_target_num_chrXII = SpSaAs_target_num_chrXII + array_SpSaAs_chrXII[i]

			print('chrXII length: ' + str(len(chr_seq)))
			print('chrXII Sp: ' + str(int(Sp_PAM_chrXII)) + ' ' + str(int(Sp_target_num_chrXII)) + ' ' + str(Sp_target_num_chrXII/len(chr_seq)))
			print('chrXII Sa: ' + str(int(Sa_PAM_chrXII)) + ' ' + str(int(Sa_target_num_chrXII)) + ' ' + str(Sa_target_num_chrXII/len(chr_seq)))
			print('chrXII As: ' + str(int(As_PAM_chrXII)) + ' ' + str(int(As_target_num_chrXII)) + ' ' + str(As_target_num_chrXII/len(chr_seq)))
			print('chrXII SpSa: ' + str(int(SpSa_target_num_chrXII)) + ' ' + str(SpSa_target_num_chrXII/len(chr_seq)))
			print('chrXII SpAs: ' + str(int(SpAs_target_num_chrXII)) + ' ' + str(SpAs_target_num_chrXII/len(chr_seq)))
			print('chrXII SaAs: ' + str(int(SaAs_target_num_chrXII)) + ' ' + str(SaAs_target_num_chrXII/len(chr_seq)))
			print('chrXII SpSaAs: ' + str(int(SpSaAs_target_num_chrXII)) + ' ' + str(SpSaAs_target_num_chrXII/len(chr_seq)))

		elif record.id == "chrXIII":
			chr_seq = record.seq
			Sp_PAM_chrXIII = 0
			Sa_PAM_chrXIII = 0
			As_PAM_chrXIII = 0
			for i in range(len(chr_seq)):
								#SpCas9
				if i+1 < len(chr_seq):
					#SpCas9 Forward PAM
					if chr_seq[i]=='G' and chr_seq[i+1]=='G':
						Sp_PAM_chrXIII = Sp_PAM_chrXIII + 1
						if i-12 >=0:
							array_Sp_chrXIII[i-12] = True
						if i-11 >=0:
							array_Sp_chrXIII[i-11] = True
						if i-10 >=0:
							array_Sp_chrXIII[i-10] = True
						if i-9 >= 0:
							array_Sp_chrXIII[i-9] = True
						if i-8 >= 0:
							array_Sp_chrXIII[i-8] = True
						if i-7 >= 0:
							array_Sp_chrXIII[i-7] = True
						if i-6 >= 0:
							array_Sp_chrXIII[i-6] = True
						if i-5 >= 0:
							array_Sp_chrXIII[i-5] = True
						if i-4 >= 0:
							array_Sp_chrXIII[i-4] = True
						if i-3 >= 0:
							array_Sp_chrXIII[i-3] = True
						if i-2 >= 0:
							array_Sp_chrXIII[i-2] = True
						#if i-1 >= 0:
						#	array_Sp_chrXIII[i-1] = True
						if i >= 0:
							array_Sp_chrXIII[i] = True
						if i+1 < len(chr_seq):
							array_Sp_chrXIII[i+1] = True

					#SpCas9 Reverse PAM
					if chr_seq[i]=='C' and chr_seq[i+1]=='C':
						Sp_PAM_chrXIII = Sp_PAM_chrXIII + 1
						if i < len(chr_seq):
							array_Sp_chrXIII[i] = True
						if i+1 < len(chr_seq):
							array_Sp_chrXIII[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sp_chrXIII[i+2] = True
						if i+3 < len(chr_seq):
							array_Sp_chrXIII[i+3] = True
						if i+4 < len(chr_seq):
							array_Sp_chrXIII[i+4] = True
						if i+5 < len(chr_seq):
							array_Sp_chrXIII[i+5] = True
						if i+6 < len(chr_seq):
							array_Sp_chrXIII[i+6] = True
						if i+7 < len(chr_seq):
							array_Sp_chrXIII[i+7] = True
						if i+8 < len(chr_seq):
							array_Sp_chrXIII[i+8] = True
						if i+9 < len(chr_seq):
							array_Sp_chrXIII[i+9] = True
						if i+10 < len(chr_seq):
							array_Sp_chrXIII[i+10] = True
						if i+11 < len(chr_seq):
							array_Sp_chrXIII[i+11] = True
						if i+12 < len(chr_seq):
							array_Sp_chrXIII[i+12] = True
						if i+13 < len(chr_seq):
							array_Sp_chrXIII[i+13] = True

				#SaCas9
				if i+3 < len(chr_seq):
					#SaCas9 Forward PAM
					if chr_seq[i]=='G' and (chr_seq[i+1]=='A' or chr_seq[i+1]=='G') and (chr_seq[i+2]=='A' or chr_seq[i+2]=='G') and chr_seq[i+3]=='T':
						Sa_PAM_chrXIII = Sa_PAM_chrXIII + 1
						if i-13 >= 0:
							array_Sa_chrXIII[i-13] = True
						if i-12 >= 0:
							array_Sa_chrXIII[i-12] = True
						if i-11 >= 0:
							array_Sa_chrXIII[i-11] = True
						if i-10 >= 0:
							array_Sa_chrXIII[i-10] = True
						if i-9 >= 0:
							array_Sa_chrXIII[i-9] = True
						if i-8 >= 0:
							array_Sa_chrXIII[i-8] = True
						if i-7 >= 0:
							array_Sa_chrXIII[i-7] = True
						if i-6 >= 0:
							array_Sa_chrXIII[i-6] = True
						if i-5 >= 0:
							array_Sa_chrXIII[i-5] = True
						if i-4 >= 0:
							array_Sa_chrXIII[i-4] = True
						if i-3 >= 0:
							array_Sa_chrXIII[i-3] = True
						#if i-2 >= 0:
						#	array_Sa_chrXIII[i-2] = True
						#if i-1 >= 0:
						#	array_Sa_chrXIII[i-1] = True
						if i >= 0:
							array_Sa_chrXIII[i] = True
						#if i+1 < len(chr_seq):
						#	array_Sa_chrXIII[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sa_chrXIII[i+2] = True
						if i+3 < len(chr_seq):
							array_Sa_chrXIII[i+3] = True

					#SaCas9 Reverse PAM
					if chr_seq[i]=='A' and (chr_seq[i+1]=='C' or chr_seq[i+1]=='T') and (chr_seq[i+2]=='C' or chr_seq[i+2]=='T') and chr_seq[i+3]=='C':
						Sa_PAM_chrXIII = Sa_PAM_chrXIII + 1
						if i < len(chr_seq):
							array_Sa_chrXIII[i] = True
						#if i+1 < len(chr_seq):
						#	array_Sa_chrXIII[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sa_chrXIII[i+2] = True
						if i+3 < len(chr_seq):
							array_Sa_chrXIII[i+3] = True
						#if i+4 < len(chr_seq):
						#	array_Sa_chrXIII[i+4] = True
						#if i+5 < len(chr_seq):
						#	array_Sa_chrXIII[i+5] = True
						if i+6 < len(chr_seq):
							array_Sa_chrXIII[i+6] = True
						if i+7 < len(chr_seq):
							array_Sa_chrXIII[i+7] = True
						if i+8 < len(chr_seq):
							array_Sa_chrXIII[i+8] = True
						if i+9 < len(chr_seq):
							array_Sa_chrXIII[i+9] = True
						if i+10 < len(chr_seq):
							array_Sa_chrXIII[i+10] = True
						if i+11 < len(chr_seq):
							array_Sa_chrXIII[i+11] = True
						if i+12 < len(chr_seq):
							array_Sa_chrXIII[i+12] = True
						if i+13 < len(chr_seq):
							array_Sa_chrXIII[i+13] = True
						if i+14 < len(chr_seq):
							array_Sa_chrXIII[i+14] = True
						if i+15 < len(chr_seq):
							array_Sa_chrXIII[i+15] = True
						if i+16 < len(chr_seq):
							array_Sa_chrXIII[i+16] = True

				#AsCas12a
				if i+3 < len(chr_seq):
					#AsCas12 Forward PAM
					if chr_seq[i]=='T' and chr_seq[i+1]=='T' and chr_seq[i+2]=='T' and (chr_seq[i+3]=='A' or chr_seq[i+3]=='G' or chr_seq[i+3]=='C'):
						As_PAM_chrXIII = As_PAM_chrXIII + 1
						if i < len(chr_seq):
							array_As_chrXIII[i] = True
						if i+1 < len(chr_seq):
							array_As_chrXIII[i+1] = True
						if i+2 < len(chr_seq):
							array_As_chrXIII[i+2] = True
						if i+4 < len(chr_seq):
							array_As_chrXIII[i+4] = True
						if i+5 < len(chr_seq):
							array_As_chrXIII[i+5] = True
						if i+6 < len(chr_seq):
							array_As_chrXIII[i+6] = True
						if i+7 < len(chr_seq):
							array_As_chrXIII[i+7] = True
						if i+8 < len(chr_seq):
							array_As_chrXIII[i+8] = True
						if i+9 < len(chr_seq):
							array_As_chrXIII[i+9] = True
						if i+10 < len(chr_seq):
							array_As_chrXIII[i+10] = True
						if i+11 < len(chr_seq):
							array_As_chrXIII[i+11] = True
						if i+12 < len(chr_seq):
							array_As_chrXIII[i+12] = True
						if i+13 < len(chr_seq):
							array_As_chrXIII[i+13] = True
						if i+14 < len(chr_seq):
							array_As_chrXIII[i+14] = True
						if i+15 < len(chr_seq):
							array_As_chrXIII[i+15] = True
						if i+16 < len(chr_seq):
							array_As_chrXIII[i+16] = True
						if i+17 < len(chr_seq):
							array_As_chrXIII[i+17] = True
						if i+18 < len(chr_seq):
							array_As_chrXIII[i+18] = True
						if i+19 < len(chr_seq):
							array_As_chrXIII[i+19] = True
						if i+20 < len(chr_seq):
							array_As_chrXIII[i+20] = True

					#AsCas12a Revwerse PAM
					if (chr_seq[i]=='A' or chr_seq[i]=='G' or chr_seq[i]=='C') and chr_seq[i+1]=='A' and chr_seq[i+2]=='A' and chr_seq[i+3]=='A':
						As_PAM_chrXIII = As_PAM_chrXIII + 1
						if i-17 >= 0:
							array_As_chrXIII[i-17] = True
						if i-16 >= 0:
							array_As_chrXIII[i-16] = True
						if i-15 >= 0:
							array_As_chrXIII[i-15] = True
						if i-14 >= 0:
							array_As_chrXIII[i-14] = True
						if i-13 >= 0:
							array_As_chrXIII[i-13] = True
						if i-12 >= 0:
							array_As_chrXIII[i-12] = True
						if i-11 >= 0:
							array_As_chrXIII[i-11] = True
						if i-10 >= 0:
							array_As_chrXIII[i-10] = True
						if i-9 >= 0:
							array_As_chrXIII[i-9] = True
						if i-8 >= 0:
							array_As_chrXIII[i-8] = True
						if i-7 >= 0:
							array_As_chrXIII[i-7] = True
						if i-6 >= 0:
							array_As_chrXIII[i-6] = True
						if i-5 >= 0:
							array_As_chrXIII[i-5] = True
						if i-4 >= 0:
							array_As_chrXIII[i-4] = True
						if i-3 >= 0:
							array_As_chrXIII[i-3] = True
						if i-2 >= 0:
							array_As_chrXIII[i-2] = True
						if i-1 >= 0:
							array_As_chrXIII[i-1] = True
						if i+1 < len(chr_seq):
							array_As_chrXIII[i+1] = True
						if i+2 < len(chr_seq):
							array_As_chrXIII[i+2] = True
						if i+3 < len(chr_seq):
							array_As_chrXIII[i+3] = True

			for i in range(len(chr_seq)):
				array_SpSa_chrXIII[i] = (array_Sp_chrXIII[i] or array_Sa_chrXIII[i])
				array_SpAs_chrXIII[i] = (array_Sp_chrXIII[i] or array_As_chrXIII[i])
				array_SaAs_chrXIII[i] = (array_Sa_chrXIII[i] or array_As_chrXIII[i])
				array_SpSaAs_chrXIII[i] = (array_Sp_chrXIII[i] or array_Sa_chrXIII[i] or array_As_chrXIII[i])

			Sp_target_num_chrXIII = 0
			Sa_target_num_chrXIII = 0
			As_target_num_chrXIII = 0
			SpSa_target_num_chrXIII = 0
			SpAs_target_num_chrXIII = 0
			SaAs_target_num_chrXIII = 0
			SpSaAs_target_num_chrXIII = 0

			for i in range(len(chr_seq)):
				Sp_target_num_chrXIII = Sp_target_num_chrXIII + array_Sp_chrXIII[i]
				Sa_target_num_chrXIII = Sa_target_num_chrXIII + array_Sa_chrXIII[i]
				As_target_num_chrXIII = As_target_num_chrXIII + array_As_chrXIII[i]
				SpSa_target_num_chrXIII = SpSa_target_num_chrXIII + array_SpSa_chrXIII[i]
				SpAs_target_num_chrXIII = SpAs_target_num_chrXIII + array_SpAs_chrXIII[i]
				SaAs_target_num_chrXIII = SaAs_target_num_chrXIII + array_SaAs_chrXIII[i]
				SpSaAs_target_num_chrXIII = SpSaAs_target_num_chrXIII + array_SpSaAs_chrXIII[i]

			print('chrXIII length: ' + str(len(chr_seq)))
			print('chrXIII Sp: ' + str(int(Sp_PAM_chrXIII)) + ' ' + str(int(Sp_target_num_chrXIII)) + ' ' + str(Sp_target_num_chrXIII/len(chr_seq)))
			print('chrXIII Sa: ' + str(int(Sa_PAM_chrXIII)) + ' ' + str(int(Sa_target_num_chrXIII)) + ' ' + str(Sa_target_num_chrXIII/len(chr_seq)))
			print('chrXIII As: ' + str(int(As_PAM_chrXIII)) + ' ' + str(int(As_target_num_chrXIII)) + ' ' + str(As_target_num_chrXIII/len(chr_seq)))
			print('chrXIII SpSa: ' + str(int(SpSa_target_num_chrXIII)) + ' ' + str(SpSa_target_num_chrXIII/len(chr_seq)))
			print('chrXIII SpAs: ' + str(int(SpAs_target_num_chrXIII)) + ' ' + str(SpAs_target_num_chrXIII/len(chr_seq)))
			print('chrXIII SaAs: ' + str(int(SaAs_target_num_chrXIII)) + ' ' + str(SaAs_target_num_chrXIII/len(chr_seq)))
			print('chrXIII SpSaAs: ' + str(int(SpSaAs_target_num_chrXIII)) + ' ' + str(SpSaAs_target_num_chrXIII/len(chr_seq)))

		elif record.id == "chrXIV":
			chr_seq = record.seq
			Sp_PAM_chrXIV = 0
			Sa_PAM_chrXIV = 0
			As_PAM_chrXIV = 0
			for i in range(len(chr_seq)):
								#SpCas9
				if i+1 < len(chr_seq):
					#SpCas9 Forward PAM
					if chr_seq[i]=='G' and chr_seq[i+1]=='G':
						Sp_PAM_chrXIV = Sp_PAM_chrXIV + 1
						if i-12 >=0:
							array_Sp_chrXIV[i-12] = True
						if i-11 >=0:
							array_Sp_chrXIV[i-11] = True
						if i-10 >=0:
							array_Sp_chrXIV[i-10] = True
						if i-9 >= 0:
							array_Sp_chrXIV[i-9] = True
						if i-8 >= 0:
							array_Sp_chrXIV[i-8] = True
						if i-7 >= 0:
							array_Sp_chrXIV[i-7] = True
						if i-6 >= 0:
							array_Sp_chrXIV[i-6] = True
						if i-5 >= 0:
							array_Sp_chrXIV[i-5] = True
						if i-4 >= 0:
							array_Sp_chrXIV[i-4] = True
						if i-3 >= 0:
							array_Sp_chrXIV[i-3] = True
						if i-2 >= 0:
							array_Sp_chrXIV[i-2] = True
						#if i-1 >= 0:
						#	array_Sp_chrXIV[i-1] = True
						if i >= 0:
							array_Sp_chrXIV[i] = True
						if i+1 < len(chr_seq):
							array_Sp_chrXIV[i+1] = True

					#SpCas9 Reverse PAM
					if chr_seq[i]=='C' and chr_seq[i+1]=='C':
						Sp_PAM_chrXIV = Sp_PAM_chrXIV + 1
						if i < len(chr_seq):
							array_Sp_chrXIV[i] = True
						if i+1 < len(chr_seq):
							array_Sp_chrXIV[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sp_chrXIV[i+2] = True
						if i+3 < len(chr_seq):
							array_Sp_chrXIV[i+3] = True
						if i+4 < len(chr_seq):
							array_Sp_chrXIV[i+4] = True
						if i+5 < len(chr_seq):
							array_Sp_chrXIV[i+5] = True
						if i+6 < len(chr_seq):
							array_Sp_chrXIV[i+6] = True
						if i+7 < len(chr_seq):
							array_Sp_chrXIV[i+7] = True
						if i+8 < len(chr_seq):
							array_Sp_chrXIV[i+8] = True
						if i+9 < len(chr_seq):
							array_Sp_chrXIV[i+9] = True
						if i+10 < len(chr_seq):
							array_Sp_chrXIV[i+10] = True
						if i+11 < len(chr_seq):
							array_Sp_chrXIV[i+11] = True
						if i+12 < len(chr_seq):
							array_Sp_chrXIV[i+12] = True
						if i+13 < len(chr_seq):
							array_Sp_chrXIV[i+13] = True

				#SaCas9
				if i+3 < len(chr_seq):
					#SaCas9 Forward PAM
					if chr_seq[i]=='G' and (chr_seq[i+1]=='A' or chr_seq[i+1]=='G') and (chr_seq[i+2]=='A' or chr_seq[i+2]=='G') and chr_seq[i+3]=='T':
						Sa_PAM_chrXIV = Sa_PAM_chrXIV + 1
						if i-13 >= 0:
							array_Sa_chrXIV[i-13] = True
						if i-12 >= 0:
							array_Sa_chrXIV[i-12] = True
						if i-11 >= 0:
							array_Sa_chrXIV[i-11] = True
						if i-10 >= 0:
							array_Sa_chrXIV[i-10] = True
						if i-9 >= 0:
							array_Sa_chrXIV[i-9] = True
						if i-8 >= 0:
							array_Sa_chrXIV[i-8] = True
						if i-7 >= 0:
							array_Sa_chrXIV[i-7] = True
						if i-6 >= 0:
							array_Sa_chrXIV[i-6] = True
						if i-5 >= 0:
							array_Sa_chrXIV[i-5] = True
						if i-4 >= 0:
							array_Sa_chrXIV[i-4] = True
						if i-3 >= 0:
							array_Sa_chrXIV[i-3] = True
						#if i-2 >= 0:
						#	array_Sa_chrXIV[i-2] = True
						#if i-1 >= 0:
						#	array_Sa_chrXIV[i-1] = True
						if i >= 0:
							array_Sa_chrXIV[i] = True
						#if i+1 < len(chr_seq):
						#	array_Sa_chrXIV[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sa_chrXIV[i+2] = True
						if i+3 < len(chr_seq):
							array_Sa_chrXIV[i+3] = True

					#SaCas9 Reverse PAM
					if chr_seq[i]=='A' and (chr_seq[i+1]=='C' or chr_seq[i+1]=='T') and (chr_seq[i+2]=='C' or chr_seq[i+2]=='T') and chr_seq[i+3]=='C':
						Sa_PAM_chrXIV = Sa_PAM_chrXIV + 1
						if i < len(chr_seq):
							array_Sa_chrXIV[i] = True
						#if i+1 < len(chr_seq):
						#	array_Sa_chrXIV[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sa_chrXIV[i+2] = True
						if i+3 < len(chr_seq):
							array_Sa_chrXIV[i+3] = True
						#if i+4 < len(chr_seq):
						#	array_Sa_chrXIV[i+4] = True
						#if i+5 < len(chr_seq):
						#	array_Sa_chrXIV[i+5] = True
						if i+6 < len(chr_seq):
							array_Sa_chrXIV[i+6] = True
						if i+7 < len(chr_seq):
							array_Sa_chrXIV[i+7] = True
						if i+8 < len(chr_seq):
							array_Sa_chrXIV[i+8] = True
						if i+9 < len(chr_seq):
							array_Sa_chrXIV[i+9] = True
						if i+10 < len(chr_seq):
							array_Sa_chrXIV[i+10] = True
						if i+11 < len(chr_seq):
							array_Sa_chrXIV[i+11] = True
						if i+12 < len(chr_seq):
							array_Sa_chrXIV[i+12] = True
						if i+13 < len(chr_seq):
							array_Sa_chrXIV[i+13] = True
						if i+14 < len(chr_seq):
							array_Sa_chrXIV[i+14] = True
						if i+15 < len(chr_seq):
							array_Sa_chrXIV[i+15] = True
						if i+16 < len(chr_seq):
							array_Sa_chrXIV[i+16] = True

				#AsCas12a
				if i+3 < len(chr_seq):
					#AsCas12 Forward PAM
					if chr_seq[i]=='T' and chr_seq[i+1]=='T' and chr_seq[i+2]=='T' and (chr_seq[i+3]=='A' or chr_seq[i+3]=='G' or chr_seq[i+3]=='C'):
						As_PAM_chrXIV = As_PAM_chrXIV + 1
						if i < len(chr_seq):
							array_As_chrXIV[i] = True
						if i+1 < len(chr_seq):
							array_As_chrXIV[i+1] = True
						if i+2 < len(chr_seq):
							array_As_chrXIV[i+2] = True
						if i+4 < len(chr_seq):
							array_As_chrXIV[i+4] = True
						if i+5 < len(chr_seq):
							array_As_chrXIV[i+5] = True
						if i+6 < len(chr_seq):
							array_As_chrXIV[i+6] = True
						if i+7 < len(chr_seq):
							array_As_chrXIV[i+7] = True
						if i+8 < len(chr_seq):
							array_As_chrXIV[i+8] = True
						if i+9 < len(chr_seq):
							array_As_chrXIV[i+9] = True
						if i+10 < len(chr_seq):
							array_As_chrXIV[i+10] = True
						if i+11 < len(chr_seq):
							array_As_chrXIV[i+11] = True
						if i+12 < len(chr_seq):
							array_As_chrXIV[i+12] = True
						if i+13 < len(chr_seq):
							array_As_chrXIV[i+13] = True
						if i+14 < len(chr_seq):
							array_As_chrXIV[i+14] = True
						if i+15 < len(chr_seq):
							array_As_chrXIV[i+15] = True
						if i+16 < len(chr_seq):
							array_As_chrXIV[i+16] = True
						if i+17 < len(chr_seq):
							array_As_chrXIV[i+17] = True
						if i+18 < len(chr_seq):
							array_As_chrXIV[i+18] = True
						if i+19 < len(chr_seq):
							array_As_chrXIV[i+19] = True
						if i+20 < len(chr_seq):
							array_As_chrXIV[i+20] = True

					#AsCas12a Revwerse PAM
					if (chr_seq[i]=='A' or chr_seq[i]=='G' or chr_seq[i]=='C') and chr_seq[i+1]=='A' and chr_seq[i+2]=='A' and chr_seq[i+3]=='A':
						As_PAM_chrXIV = As_PAM_chrXIV + 1
						if i-17 >= 0:
							array_As_chrXIV[i-17] = True
						if i-16 >= 0:
							array_As_chrXIV[i-16] = True
						if i-15 >= 0:
							array_As_chrXIV[i-15] = True
						if i-14 >= 0:
							array_As_chrXIV[i-14] = True
						if i-13 >= 0:
							array_As_chrXIV[i-13] = True
						if i-12 >= 0:
							array_As_chrXIV[i-12] = True
						if i-11 >= 0:
							array_As_chrXIV[i-11] = True
						if i-10 >= 0:
							array_As_chrXIV[i-10] = True
						if i-9 >= 0:
							array_As_chrXIV[i-9] = True
						if i-8 >= 0:
							array_As_chrXIV[i-8] = True
						if i-7 >= 0:
							array_As_chrXIV[i-7] = True
						if i-6 >= 0:
							array_As_chrXIV[i-6] = True
						if i-5 >= 0:
							array_As_chrXIV[i-5] = True
						if i-4 >= 0:
							array_As_chrXIV[i-4] = True
						if i-3 >= 0:
							array_As_chrXIV[i-3] = True
						if i-2 >= 0:
							array_As_chrXIV[i-2] = True
						if i-1 >= 0:
							array_As_chrXIV[i-1] = True
						if i+1 < len(chr_seq):
							array_As_chrXIV[i+1] = True
						if i+2 < len(chr_seq):
							array_As_chrXIV[i+2] = True
						if i+3 < len(chr_seq):
							array_As_chrXIV[i+3] = True

			for i in range(len(chr_seq)):
				array_SpSa_chrXIV[i] = (array_Sp_chrXIV[i] or array_Sa_chrXIV[i])
				array_SpAs_chrXIV[i] = (array_Sp_chrXIV[i] or array_As_chrXIV[i])
				array_SaAs_chrXIV[i] = (array_Sa_chrXIV[i] or array_As_chrXIV[i])
				array_SpSaAs_chrXIV[i] = (array_Sp_chrXIV[i] or array_Sa_chrXIV[i] or array_As_chrXIV[i])

			Sp_target_num_chrXIV = 0
			Sa_target_num_chrXIV = 0
			As_target_num_chrXIV = 0
			SpSa_target_num_chrXIV = 0
			SpAs_target_num_chrXIV = 0
			SaAs_target_num_chrXIV = 0
			SpSaAs_target_num_chrXIV = 0

			for i in range(len(chr_seq)):
				Sp_target_num_chrXIV = Sp_target_num_chrXIV + array_Sp_chrXIV[i]
				Sa_target_num_chrXIV = Sa_target_num_chrXIV + array_Sa_chrXIV[i]
				As_target_num_chrXIV = As_target_num_chrXIV + array_As_chrXIV[i]
				SpSa_target_num_chrXIV = SpSa_target_num_chrXIV + array_SpSa_chrXIV[i]
				SpAs_target_num_chrXIV = SpAs_target_num_chrXIV + array_SpAs_chrXIV[i]
				SaAs_target_num_chrXIV = SaAs_target_num_chrXIV + array_SaAs_chrXIV[i]
				SpSaAs_target_num_chrXIV = SpSaAs_target_num_chrXIV + array_SpSaAs_chrXIV[i]

			print('chrXIV length: ' + str(len(chr_seq)))
			print('chrXIV Sp: ' + str(int(Sp_PAM_chrXIV)) + ' ' + str(int(Sp_target_num_chrXIV)) + ' ' + str(Sp_target_num_chrXIV/len(chr_seq)))
			print('chrXIV Sa: ' + str(int(Sa_PAM_chrXIV)) + ' ' + str(int(Sa_target_num_chrXIV)) + ' ' + str(Sa_target_num_chrXIV/len(chr_seq)))
			print('chrXIV As: ' + str(int(As_PAM_chrXIV)) + ' ' + str(int(As_target_num_chrXIV)) + ' ' + str(As_target_num_chrXIV/len(chr_seq)))
			print('chrXIV SpSa: ' + str(int(SpSa_target_num_chrXIV)) + ' ' + str(SpSa_target_num_chrXIV/len(chr_seq)))
			print('chrXIV SpAs: ' + str(int(SpAs_target_num_chrXIV)) + ' ' + str(SpAs_target_num_chrXIV/len(chr_seq)))
			print('chrXIV SaAs: ' + str(int(SaAs_target_num_chrXIV)) + ' ' + str(SaAs_target_num_chrXIV/len(chr_seq)))
			print('chrXIV SpSaAs: ' + str(int(SpSaAs_target_num_chrXIV)) + ' ' + str(SpSaAs_target_num_chrXIV/len(chr_seq)))

		elif record.id == "chrXV":
			chr_seq = record.seq
			Sp_PAM_chrXV = 0
			Sa_PAM_chrXV = 0
			As_PAM_chrXV = 0
			for i in range(len(chr_seq)):
								#SpCas9
				if i+1 < len(chr_seq):
					#SpCas9 Forward PAM
					if chr_seq[i]=='G' and chr_seq[i+1]=='G':
						Sp_PAM_chrXV = Sp_PAM_chrXV + 1
						if i-12 >=0:
							array_Sp_chrXV[i-12] = True
						if i-11 >=0:
							array_Sp_chrXV[i-11] = True
						if i-10 >=0:
							array_Sp_chrXV[i-10] = True
						if i-9 >= 0:
							array_Sp_chrXV[i-9] = True
						if i-8 >= 0:
							array_Sp_chrXV[i-8] = True
						if i-7 >= 0:
							array_Sp_chrXV[i-7] = True
						if i-6 >= 0:
							array_Sp_chrXV[i-6] = True
						if i-5 >= 0:
							array_Sp_chrXV[i-5] = True
						if i-4 >= 0:
							array_Sp_chrXV[i-4] = True
						if i-3 >= 0:
							array_Sp_chrXV[i-3] = True
						if i-2 >= 0:
							array_Sp_chrXV[i-2] = True
						#if i-1 >= 0:
						#	array_Sp_chrXV[i-1] = True
						if i >= 0:
							array_Sp_chrXV[i] = True
						if i+1 < len(chr_seq):
							array_Sp_chrXV[i+1] = True

					#SpCas9 Reverse PAM
					if chr_seq[i]=='C' and chr_seq[i+1]=='C':
						Sp_PAM_chrXV = Sp_PAM_chrXV + 1
						if i < len(chr_seq):
							array_Sp_chrXV[i] = True
						if i+1 < len(chr_seq):
							array_Sp_chrXV[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sp_chrXV[i+2] = True
						if i+3 < len(chr_seq):
							array_Sp_chrXV[i+3] = True
						if i+4 < len(chr_seq):
							array_Sp_chrXV[i+4] = True
						if i+5 < len(chr_seq):
							array_Sp_chrXV[i+5] = True
						if i+6 < len(chr_seq):
							array_Sp_chrXV[i+6] = True
						if i+7 < len(chr_seq):
							array_Sp_chrXV[i+7] = True
						if i+8 < len(chr_seq):
							array_Sp_chrXV[i+8] = True
						if i+9 < len(chr_seq):
							array_Sp_chrXV[i+9] = True
						if i+10 < len(chr_seq):
							array_Sp_chrXV[i+10] = True
						if i+11 < len(chr_seq):
							array_Sp_chrXV[i+11] = True
						if i+12 < len(chr_seq):
							array_Sp_chrXV[i+12] = True
						if i+13 < len(chr_seq):
							array_Sp_chrXV[i+13] = True

				#SaCas9
				if i+3 < len(chr_seq):
					#SaCas9 Forward PAM
					if chr_seq[i]=='G' and (chr_seq[i+1]=='A' or chr_seq[i+1]=='G') and (chr_seq[i+2]=='A' or chr_seq[i+2]=='G') and chr_seq[i+3]=='T':
						Sa_PAM_chrXV = Sa_PAM_chrXV + 1
						if i-13 >= 0:
							array_Sa_chrXV[i-13] = True
						if i-12 >= 0:
							array_Sa_chrXV[i-12] = True
						if i-11 >= 0:
							array_Sa_chrXV[i-11] = True
						if i-10 >= 0:
							array_Sa_chrXV[i-10] = True
						if i-9 >= 0:
							array_Sa_chrXV[i-9] = True
						if i-8 >= 0:
							array_Sa_chrXV[i-8] = True
						if i-7 >= 0:
							array_Sa_chrXV[i-7] = True
						if i-6 >= 0:
							array_Sa_chrXV[i-6] = True
						if i-5 >= 0:
							array_Sa_chrXV[i-5] = True
						if i-4 >= 0:
							array_Sa_chrXV[i-4] = True
						if i-3 >= 0:
							array_Sa_chrXV[i-3] = True
						#if i-2 >= 0:
						#	array_Sa_chrXV[i-2] = True
						#if i-1 >= 0:
						#	array_Sa_chrXV[i-1] = True
						if i >= 0:
							array_Sa_chrXV[i] = True
						#if i+1 < len(chr_seq):
						#	array_Sa_chrXV[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sa_chrXV[i+2] = True
						if i+3 < len(chr_seq):
							array_Sa_chrXV[i+3] = True

					#SaCas9 Reverse PAM
					if chr_seq[i]=='A' and (chr_seq[i+1]=='C' or chr_seq[i+1]=='T') and (chr_seq[i+2]=='C' or chr_seq[i+2]=='T') and chr_seq[i+3]=='C':
						Sa_PAM_chrXV = Sa_PAM_chrXV + 1
						if i < len(chr_seq):
							array_Sa_chrXV[i] = True
						#if i+1 < len(chr_seq):
						#	array_Sa_chrXV[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sa_chrXV[i+2] = True
						if i+3 < len(chr_seq):
							array_Sa_chrXV[i+3] = True
						#if i+4 < len(chr_seq):
						#	array_Sa_chrXV[i+4] = True
						#if i+5 < len(chr_seq):
						#	array_Sa_chrXV[i+5] = True
						if i+6 < len(chr_seq):
							array_Sa_chrXV[i+6] = True
						if i+7 < len(chr_seq):
							array_Sa_chrXV[i+7] = True
						if i+8 < len(chr_seq):
							array_Sa_chrXV[i+8] = True
						if i+9 < len(chr_seq):
							array_Sa_chrXV[i+9] = True
						if i+10 < len(chr_seq):
							array_Sa_chrXV[i+10] = True
						if i+11 < len(chr_seq):
							array_Sa_chrXV[i+11] = True
						if i+12 < len(chr_seq):
							array_Sa_chrXV[i+12] = True
						if i+13 < len(chr_seq):
							array_Sa_chrXV[i+13] = True
						if i+14 < len(chr_seq):
							array_Sa_chrXV[i+14] = True
						if i+15 < len(chr_seq):
							array_Sa_chrXV[i+15] = True
						if i+16 < len(chr_seq):
							array_Sa_chrXV[i+16] = True

				#AsCas12a
				if i+3 < len(chr_seq):
					#AsCas12 Forward PAM
					if chr_seq[i]=='T' and chr_seq[i+1]=='T' and chr_seq[i+2]=='T' and (chr_seq[i+3]=='A' or chr_seq[i+3]=='G' or chr_seq[i+3]=='C'):
						As_PAM_chrXV = As_PAM_chrXV + 1
						if i < len(chr_seq):
							array_As_chrXV[i] = True
						if i+1 < len(chr_seq):
							array_As_chrXV[i+1] = True
						if i+2 < len(chr_seq):
							array_As_chrXV[i+2] = True
						if i+4 < len(chr_seq):
							array_As_chrXV[i+4] = True
						if i+5 < len(chr_seq):
							array_As_chrXV[i+5] = True
						if i+6 < len(chr_seq):
							array_As_chrXV[i+6] = True
						if i+7 < len(chr_seq):
							array_As_chrXV[i+7] = True
						if i+8 < len(chr_seq):
							array_As_chrXV[i+8] = True
						if i+9 < len(chr_seq):
							array_As_chrXV[i+9] = True
						if i+10 < len(chr_seq):
							array_As_chrXV[i+10] = True
						if i+11 < len(chr_seq):
							array_As_chrXV[i+11] = True
						if i+12 < len(chr_seq):
							array_As_chrXV[i+12] = True
						if i+13 < len(chr_seq):
							array_As_chrXV[i+13] = True
						if i+14 < len(chr_seq):
							array_As_chrXV[i+14] = True
						if i+15 < len(chr_seq):
							array_As_chrXV[i+15] = True
						if i+16 < len(chr_seq):
							array_As_chrXV[i+16] = True
						if i+17 < len(chr_seq):
							array_As_chrXV[i+17] = True
						if i+18 < len(chr_seq):
							array_As_chrXV[i+18] = True
						if i+19 < len(chr_seq):
							array_As_chrXV[i+19] = True
						if i+20 < len(chr_seq):
							array_As_chrXV[i+20] = True

					#AsCas12a Revwerse PAM
					if (chr_seq[i]=='A' or chr_seq[i]=='G' or chr_seq[i]=='C') and chr_seq[i+1]=='A' and chr_seq[i+2]=='A' and chr_seq[i+3]=='A':
						As_PAM_chrXV = As_PAM_chrXV + 1
						if i-17 >= 0:
							array_As_chrXV[i-17] = True
						if i-16 >= 0:
							array_As_chrXV[i-16] = True
						if i-15 >= 0:
							array_As_chrXV[i-15] = True
						if i-14 >= 0:
							array_As_chrXV[i-14] = True
						if i-13 >= 0:
							array_As_chrXV[i-13] = True
						if i-12 >= 0:
							array_As_chrXV[i-12] = True
						if i-11 >= 0:
							array_As_chrXV[i-11] = True
						if i-10 >= 0:
							array_As_chrXV[i-10] = True
						if i-9 >= 0:
							array_As_chrXV[i-9] = True
						if i-8 >= 0:
							array_As_chrXV[i-8] = True
						if i-7 >= 0:
							array_As_chrXV[i-7] = True
						if i-6 >= 0:
							array_As_chrXV[i-6] = True
						if i-5 >= 0:
							array_As_chrXV[i-5] = True
						if i-4 >= 0:
							array_As_chrXV[i-4] = True
						if i-3 >= 0:
							array_As_chrXV[i-3] = True
						if i-2 >= 0:
							array_As_chrXV[i-2] = True
						if i-1 >= 0:
							array_As_chrXV[i-1] = True
						if i+1 < len(chr_seq):
							array_As_chrXV[i+1] = True
						if i+2 < len(chr_seq):
							array_As_chrXV[i+2] = True
						if i+3 < len(chr_seq):
							array_As_chrXV[i+3] = True

			for i in range(len(chr_seq)):
				array_SpSa_chrXV[i] = (array_Sp_chrXV[i] or array_Sa_chrXV[i])
				array_SpAs_chrXV[i] = (array_Sp_chrXV[i] or array_As_chrXV[i])
				array_SaAs_chrXV[i] = (array_Sa_chrXV[i] or array_As_chrXV[i])
				array_SpSaAs_chrXV[i] = (array_Sp_chrXV[i] or array_Sa_chrXV[i] or array_As_chrXV[i])

			Sp_target_num_chrXV = 0
			Sa_target_num_chrXV = 0
			As_target_num_chrXV = 0
			SpSa_target_num_chrXV = 0
			SpAs_target_num_chrXV = 0
			SaAs_target_num_chrXV = 0
			SpSaAs_target_num_chrXV = 0

			for i in range(len(chr_seq)):
				Sp_target_num_chrXV = Sp_target_num_chrXV + array_Sp_chrXV[i]
				Sa_target_num_chrXV = Sa_target_num_chrXV + array_Sa_chrXV[i]
				As_target_num_chrXV = As_target_num_chrXV + array_As_chrXV[i]
				SpSa_target_num_chrXV = SpSa_target_num_chrXV + array_SpSa_chrXV[i]
				SpAs_target_num_chrXV = SpAs_target_num_chrXV + array_SpAs_chrXV[i]
				SaAs_target_num_chrXV = SaAs_target_num_chrXV + array_SaAs_chrXV[i]
				SpSaAs_target_num_chrXV = SpSaAs_target_num_chrXV + array_SpSaAs_chrXV[i]

			print('chrXV length: ' + str(len(chr_seq)))
			print('chrXV Sp: ' + str(int(Sp_PAM_chrXV)) + ' ' + str(int(Sp_target_num_chrXV)) + ' ' + str(Sp_target_num_chrXV/len(chr_seq)))
			print('chrXV Sa: ' + str(int(Sa_PAM_chrXV)) + ' ' + str(int(Sa_target_num_chrXV)) + ' ' + str(Sa_target_num_chrXV/len(chr_seq)))
			print('chrXV As: ' + str(int(As_PAM_chrXV)) + ' ' + str(int(As_target_num_chrXV)) + ' ' + str(As_target_num_chrXV/len(chr_seq)))
			print('chrXV SpSa: ' + str(int(SpSa_target_num_chrXV)) + ' ' + str(SpSa_target_num_chrXV/len(chr_seq)))
			print('chrXV SpAs: ' + str(int(SpAs_target_num_chrXV)) + ' ' + str(SpAs_target_num_chrXV/len(chr_seq)))
			print('chrXV SaAs: ' + str(int(SaAs_target_num_chrXV)) + ' ' + str(SaAs_target_num_chrXV/len(chr_seq)))
			print('chrXV SpSaAs: ' + str(int(SpSaAs_target_num_chrXV)) + ' ' + str(SpSaAs_target_num_chrXV/len(chr_seq)))

		elif record.id == "chrXVI":
			chr_seq = record.seq
			Sp_PAM_chrXVI = 0
			Sa_PAM_chrXVI = 0
			As_PAM_chrXVI = 0
			for i in range(len(chr_seq)):
								#SpCas9
				if i+1 < len(chr_seq):
					#SpCas9 Forward PAM
					if chr_seq[i]=='G' and chr_seq[i+1]=='G':
						Sp_PAM_chrXVI = Sp_PAM_chrXVI + 1
						if i-12 >=0:
							array_Sp_chrXVI[i-12] = True
						if i-11 >=0:
							array_Sp_chrXVI[i-11] = True
						if i-10 >=0:
							array_Sp_chrXVI[i-10] = True
						if i-9 >= 0:
							array_Sp_chrXVI[i-9] = True
						if i-8 >= 0:
							array_Sp_chrXVI[i-8] = True
						if i-7 >= 0:
							array_Sp_chrXVI[i-7] = True
						if i-6 >= 0:
							array_Sp_chrXVI[i-6] = True
						if i-5 >= 0:
							array_Sp_chrXVI[i-5] = True
						if i-4 >= 0:
							array_Sp_chrXVI[i-4] = True
						if i-3 >= 0:
							array_Sp_chrXVI[i-3] = True
						if i-2 >= 0:
							array_Sp_chrXVI[i-2] = True
						#if i-1 >= 0:
						#	array_Sp_chrXVI[i-1] = True
						if i >= 0:
							array_Sp_chrXVI[i] = True
						if i+1 < len(chr_seq):
							array_Sp_chrXVI[i+1] = True

					#SpCas9 Reverse PAM
					if chr_seq[i]=='C' and chr_seq[i+1]=='C':
						Sp_PAM_chrXVI = Sp_PAM_chrXVI + 1
						if i < len(chr_seq):
							array_Sp_chrXVI[i] = True
						if i+1 < len(chr_seq):
							array_Sp_chrXVI[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sp_chrXVI[i+2] = True
						if i+3 < len(chr_seq):
							array_Sp_chrXVI[i+3] = True
						if i+4 < len(chr_seq):
							array_Sp_chrXVI[i+4] = True
						if i+5 < len(chr_seq):
							array_Sp_chrXVI[i+5] = True
						if i+6 < len(chr_seq):
							array_Sp_chrXVI[i+6] = True
						if i+7 < len(chr_seq):
							array_Sp_chrXVI[i+7] = True
						if i+8 < len(chr_seq):
							array_Sp_chrXVI[i+8] = True
						if i+9 < len(chr_seq):
							array_Sp_chrXVI[i+9] = True
						if i+10 < len(chr_seq):
							array_Sp_chrXVI[i+10] = True
						if i+11 < len(chr_seq):
							array_Sp_chrXVI[i+11] = True
						if i+12 < len(chr_seq):
							array_Sp_chrXVI[i+12] = True
						if i+13 < len(chr_seq):
							array_Sp_chrXVI[i+13] = True

				#SaCas9
				if i+3 < len(chr_seq):
					#SaCas9 Forward PAM
					if chr_seq[i]=='G' and (chr_seq[i+1]=='A' or chr_seq[i+1]=='G') and (chr_seq[i+2]=='A' or chr_seq[i+2]=='G') and chr_seq[i+3]=='T':
						Sa_PAM_chrXVI = Sa_PAM_chrXVI + 1
						if i-13 >= 0:
							array_Sa_chrXVI[i-13] = True
						if i-12 >= 0:
							array_Sa_chrXVI[i-12] = True
						if i-11 >= 0:
							array_Sa_chrXVI[i-11] = True
						if i-10 >= 0:
							array_Sa_chrXVI[i-10] = True
						if i-9 >= 0:
							array_Sa_chrXVI[i-9] = True
						if i-8 >= 0:
							array_Sa_chrXVI[i-8] = True
						if i-7 >= 0:
							array_Sa_chrXVI[i-7] = True
						if i-6 >= 0:
							array_Sa_chrXVI[i-6] = True
						if i-5 >= 0:
							array_Sa_chrXVI[i-5] = True
						if i-4 >= 0:
							array_Sa_chrXVI[i-4] = True
						if i-3 >= 0:
							array_Sa_chrXVI[i-3] = True
						#if i-2 >= 0:
						#	array_Sa_chrXVI[i-2] = True
						#if i-1 >= 0:
						#	array_Sa_chrXVI[i-1] = True
						if i >= 0:
							array_Sa_chrXVI[i] = True
						#if i+1 < len(chr_seq):
						#	array_Sa_chrXVI[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sa_chrXVI[i+2] = True
						if i+3 < len(chr_seq):
							array_Sa_chrXVI[i+3] = True

					#SaCas9 Reverse PAM
					if chr_seq[i]=='A' and (chr_seq[i+1]=='C' or chr_seq[i+1]=='T') and (chr_seq[i+2]=='C' or chr_seq[i+2]=='T') and chr_seq[i+3]=='C':
						Sa_PAM_chrXVI = Sa_PAM_chrXVI + 1
						if i < len(chr_seq):
							array_Sa_chrXVI[i] = True
						#if i+1 < len(chr_seq):
						#	array_Sa_chrXVI[i+1] = True
						#if i+2 < len(chr_seq):
						#	array_Sa_chrXVI[i+2] = True
						if i+3 < len(chr_seq):
							array_Sa_chrXVI[i+3] = True
						#if i+4 < len(chr_seq):
						#	array_Sa_chrXVI[i+4] = True
						#if i+5 < len(chr_seq):
						#	array_Sa_chrXVI[i+5] = True
						if i+6 < len(chr_seq):
							array_Sa_chrXVI[i+6] = True
						if i+7 < len(chr_seq):
							array_Sa_chrXVI[i+7] = True
						if i+8 < len(chr_seq):
							array_Sa_chrXVI[i+8] = True
						if i+9 < len(chr_seq):
							array_Sa_chrXVI[i+9] = True
						if i+10 < len(chr_seq):
							array_Sa_chrXVI[i+10] = True
						if i+11 < len(chr_seq):
							array_Sa_chrXVI[i+11] = True
						if i+12 < len(chr_seq):
							array_Sa_chrXVI[i+12] = True
						if i+13 < len(chr_seq):
							array_Sa_chrXVI[i+13] = True
						if i+14 < len(chr_seq):
							array_Sa_chrXVI[i+14] = True
						if i+15 < len(chr_seq):
							array_Sa_chrXVI[i+15] = True
						if i+16 < len(chr_seq):
							array_Sa_chrXVI[i+16] = True

				#AsCas12a
				if i+3 < len(chr_seq):
					#AsCas12 Forward PAM
					if chr_seq[i]=='T' and chr_seq[i+1]=='T' and chr_seq[i+2]=='T' and (chr_seq[i+3]=='A' or chr_seq[i+3]=='G' or chr_seq[i+3]=='C'):
						As_PAM_chrXVI = As_PAM_chrXVI + 1
						if i < len(chr_seq):
							array_As_chrXVI[i] = True
						if i+1 < len(chr_seq):
							array_As_chrXVI[i+1] = True
						if i+2 < len(chr_seq):
							array_As_chrXVI[i+2] = True
						if i+4 < len(chr_seq):
							array_As_chrXVI[i+4] = True
						if i+5 < len(chr_seq):
							array_As_chrXVI[i+5] = True
						if i+6 < len(chr_seq):
							array_As_chrXVI[i+6] = True
						if i+7 < len(chr_seq):
							array_As_chrXVI[i+7] = True
						if i+8 < len(chr_seq):
							array_As_chrXVI[i+8] = True
						if i+9 < len(chr_seq):
							array_As_chrXVI[i+9] = True
						if i+10 < len(chr_seq):
							array_As_chrXVI[i+10] = True
						if i+11 < len(chr_seq):
							array_As_chrXVI[i+11] = True
						if i+12 < len(chr_seq):
							array_As_chrXVI[i+12] = True
						if i+13 < len(chr_seq):
							array_As_chrXVI[i+13] = True
						if i+14 < len(chr_seq):
							array_As_chrXVI[i+14] = True
						if i+15 < len(chr_seq):
							array_As_chrXVI[i+15] = True
						if i+16 < len(chr_seq):
							array_As_chrXVI[i+16] = True
						if i+17 < len(chr_seq):
							array_As_chrXVI[i+17] = True
						if i+18 < len(chr_seq):
							array_As_chrXVI[i+18] = True
						if i+19 < len(chr_seq):
							array_As_chrXVI[i+19] = True
						if i+20 < len(chr_seq):
							array_As_chrXVI[i+20] = True

					#AsCas12a Revwerse PAM
					if (chr_seq[i]=='A' or chr_seq[i]=='G' or chr_seq[i]=='C') and chr_seq[i+1]=='A' and chr_seq[i+2]=='A' and chr_seq[i+3]=='A':
						As_PAM_chrXVI = As_PAM_chrXVI + 1
						if i-17 >= 0:
							array_As_chrXVI[i-17] = True
						if i-16 >= 0:
							array_As_chrXVI[i-16] = True
						if i-15 >= 0:
							array_As_chrXVI[i-15] = True
						if i-14 >= 0:
							array_As_chrXVI[i-14] = True
						if i-13 >= 0:
							array_As_chrXVI[i-13] = True
						if i-12 >= 0:
							array_As_chrXVI[i-12] = True
						if i-11 >= 0:
							array_As_chrXVI[i-11] = True
						if i-10 >= 0:
							array_As_chrXVI[i-10] = True
						if i-9 >= 0:
							array_As_chrXVI[i-9] = True
						if i-8 >= 0:
							array_As_chrXVI[i-8] = True
						if i-7 >= 0:
							array_As_chrXVI[i-7] = True
						if i-6 >= 0:
							array_As_chrXVI[i-6] = True
						if i-5 >= 0:
							array_As_chrXVI[i-5] = True
						if i-4 >= 0:
							array_As_chrXVI[i-4] = True
						if i-3 >= 0:
							array_As_chrXVI[i-3] = True
						if i-2 >= 0:
							array_As_chrXVI[i-2] = True
						if i-1 >= 0:
							array_As_chrXVI[i-1] = True
						if i+1 < len(chr_seq):
							array_As_chrXVI[i+1] = True
						if i+2 < len(chr_seq):
							array_As_chrXVI[i+2] = True
						if i+3 < len(chr_seq):
							array_As_chrXVI[i+3] = True

			for i in range(len(chr_seq)):
				array_SpSa_chrXVI[i] = (array_Sp_chrXVI[i] or array_Sa_chrXVI[i])
				array_SpAs_chrXVI[i] = (array_Sp_chrXVI[i] or array_As_chrXVI[i])
				array_SaAs_chrXVI[i] = (array_Sa_chrXVI[i] or array_As_chrXVI[i])
				array_SpSaAs_chrXVI[i] = (array_Sp_chrXVI[i] or array_Sa_chrXVI[i] or array_As_chrXVI[i])

			Sp_target_num_chrXVI = 0
			Sa_target_num_chrXVI = 0
			As_target_num_chrXVI = 0
			SpSa_target_num_chrXVI = 0
			SpAs_target_num_chrXVI = 0
			SaAs_target_num_chrXVI = 0
			SpSaAs_target_num_chrXVI = 0

			for i in range(len(chr_seq)):
				Sp_target_num_chrXVI = Sp_target_num_chrXVI + array_Sp_chrXVI[i]
				Sa_target_num_chrXVI = Sa_target_num_chrXVI + array_Sa_chrXVI[i]
				As_target_num_chrXVI = As_target_num_chrXVI + array_As_chrXVI[i]
				SpSa_target_num_chrXVI = SpSa_target_num_chrXVI + array_SpSa_chrXVI[i]
				SpAs_target_num_chrXVI = SpAs_target_num_chrXVI + array_SpAs_chrXVI[i]
				SaAs_target_num_chrXVI = SaAs_target_num_chrXVI + array_SaAs_chrXVI[i]
				SpSaAs_target_num_chrXVI = SpSaAs_target_num_chrXVI + array_SpSaAs_chrXVI[i]

			print('chrXVI length: ' + str(len(chr_seq)))
			print('chrXVI Sp: ' + str(int(Sp_PAM_chrXVI)) + ' ' + str(int(Sp_target_num_chrXVI)) + ' ' + str(Sp_target_num_chrXVI/len(chr_seq)))
			print('chrXVI Sa: ' + str(int(Sa_PAM_chrXVI)) + ' ' + str(int(Sa_target_num_chrXVI)) + ' ' + str(Sa_target_num_chrXVI/len(chr_seq)))
			print('chrXVI As: ' + str(int(As_PAM_chrXVI)) + ' ' + str(int(As_target_num_chrXVI)) + ' ' + str(As_target_num_chrXVI/len(chr_seq)))
			print('chrXVI SpSa: ' + str(int(SpSa_target_num_chrXVI)) + ' ' + str(SpSa_target_num_chrXVI/len(chr_seq)))
			print('chrXVI SpAs: ' + str(int(SpAs_target_num_chrXVI)) + ' ' + str(SpAs_target_num_chrXVI/len(chr_seq)))
			print('chrXVI SaAs: ' + str(int(SaAs_target_num_chrXVI)) + ' ' + str(SaAs_target_num_chrXVI/len(chr_seq)))
			print('chrXVI SpSaAs: ' + str(int(SpSaAs_target_num_chrXVI)) + ' ' + str(SpSaAs_target_num_chrXVI/len(chr_seq)))
