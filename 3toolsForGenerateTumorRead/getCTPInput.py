## -*- coding: utf-8 -*-
#python getCTPInput.py  all_snp.txt
#生成CTPsingle的输入，CTPInput.txt

import sys
assert len(sys.argv) == 2
file1= sys.argv[1]

with open(file1) as fb:
	lines = fb.readlines()
str1 = "Chromosome"+"\t"+"Position"+"\t"+"Mutant"+"\t"+"Reference"+"\t"+"Mcount"+"\t"+"Rcount"+"\t"+"Multiplier"+"\t"+"Gender"+"\n"
for line in lines[1:]:
	line = line.strip().split('\t')
	Chromosome = line[0][3]
	Position = line[1]
	Mutant = line[3]
	Reference = line[2]
	Mcount = line[8]
	Rcount = line[7]
	Multiplier = 2
	Gender = "unknown"
	str1 += Chromosome+"\t"+str(Position)+"\t"+Mutant+"\t"+Reference+"\t"+str(Mcount)+"\t"+str(Rcount)+"\t"+str(Multiplier)+"\t"+Gender+"\n"
txtName = "CTPInput.txt"
f = open(txtName,'w')
f.write(str1)
f.close()



