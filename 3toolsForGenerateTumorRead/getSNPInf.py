## -*- coding: utf-8 -*-
#python getSNPInf.py  s100.snp.vcf
#生成all_snp.txt
import vcf

import sys
assert len(sys.argv) == 2
file1= sys.argv[1]

vcf_reader = vcf.Reader(open(file1,"r"))
str1 = "CHROM"+'\t'+'POS'+'\t'+'REF'+'\t'+"ALT"+'\t'+"AB"+'\t'+"AF"+'\t'+"DP"+'\t'+"RO"+'\t'+"AO"+"\n"
for record in vcf_reader:
	str1 += str(record.CHROM)+'\t'+str(record.POS)+'\t'+str(record.REF)+'\t'+str(record.ALT[0])+'\t'+str(record.INFO['AB'][0])+'\t'+str(record.INFO['AF'][0])+'\t'+str(record.INFO['DP'])+'\t'+str(record.INFO['RO'])+'\t'+str(record.INFO['AO'][0])+'\n'
txtName = "all_snp.txt"
f = open(txtName,'w')
f.write(str1)
f.close()

