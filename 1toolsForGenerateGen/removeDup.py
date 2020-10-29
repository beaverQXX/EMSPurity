#!/usr/bin/env python
# -*- coding: utf-8 -*-
#python removeDup.py nor.sim s1.sim diff_nor_s1.txt
#python removeDup.py s100.snpincnv.txt train11.txt train12.txt
#python removeDup.py snpcnv.txt snp.txt  cnv.txt
import sys
input_file1 = sys.argv[1]
input_file2 = sys.argv[2]
output_file = sys.argv[3]
snp_file1 = []
snp_file2 = []
snp_dump = []
#input_file1 = "nor.sim"
#input_file2 = "s1.sim"
#output_file = "diff_nor_s1.txt"
def load_snps(file):
	snps = []
	with open(file, 'r') as fin:
		for line in fin.readlines():
			snps.append(line)
	return snps
if __name__ == '__main__':
	fout = open(output_file,'w')
	snp_file1 = load_snps(input_file1)
	snp_file2 = load_snps(input_file2)
	snp_dump = snp_file2
	for i in snp_file1:
		if i in snp_file2:
			snp_dump.remove(i)
	for i in range(len(snp_dump)):
		fout.write(snp_dump[i])
	fout.close()
