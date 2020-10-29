## -*- coding: utf-8 -*-
#python getCNVInf.py  s100.snp.sv.vcf
#生成all_snp.txt
import vcf
#染色体 起始位置 长度 终止位置

import sys
assert len(sys.argv) == 2
file1= sys.argv[1]

vcf_reader = vcf.Reader(open(file1,"r"))
str1 = "CHROM"+'\t'+'POS'+'\t'+'LENGTH'+'\t'+"END"+'\n'
for record in vcf_reader:
	str1 += str(record.CHROM)+'\t'+str(record.POS)+'\t'+str(record.INFO['SVLEN'][0])+'\t'+str(record.INFO['END'])+'\n' 


txtName = "all_CNV.txt"
f = open(txtName,'w')
f.write(str1)
f.close()
