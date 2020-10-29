## -*- coding: utf-8 -*-
#python getCNV_SNP.py    all_snp.txt    all_CNV.txt


import sys
assert len(sys.argv) == 3
file1,file2 = sys.argv[1],sys.argv[2]

def getSNP(file2,cnv_interval):
    cross_CNV_SNP=[]
    not_cross_CNV_SNP = []
    with open(file2) as fb:
        lines = fb.readlines()
    for line in lines[1:]:
        line = line.strip().split('\t')
        pos = int(line[1])
        #print(pos)
        find = False
        for j in range(len(cnv_interval)):
            if(not find and pos < int(cnv_interval[j][0])):
                break
            elif(not find and int(cnv_interval[j][0]) <= pos < int(cnv_interval[j][1])):
                cross_CNV_SNP.append(line)
                find = True
                break
            elif(not find and pos >= int(cnv_interval[j][1])):
                continue
        if(not find):
        	not_cross_CNV_SNP.append(line)
    return cross_CNV_SNP,not_cross_CNV_SNP


def getCNV(file2):
    cnv_interval = []
    with open(file2) as fb:
        lines = fb.readlines()
    for line in lines[1:]:
        line = line.strip().split('\t')
        cnv_interval.append([line[1],line[3]])#前闭后开
    return cnv_interval
def getCrossCNVSNP(cross_CNV_SNP):
	str1 = "CHROM"+"\t"+"POS"+"\t"+"REF"+"\t"+"ALT"+"\t"+"AB"+"\t"+"AF"+"\t"+"DP"+"\t"+"RO"+"\t"+"AO"+"\n"
	for i in range(len(cross_CNV_SNP)):
		str1 += cross_CNV_SNP[i][0]+'\t'+cross_CNV_SNP[i][1]+'\t'+cross_CNV_SNP[i][2]+'\t'+cross_CNV_SNP[i][3]+'\t'+cross_CNV_SNP[i][4]+'\t'+cross_CNV_SNP[i][5]+'\t'+cross_CNV_SNP[i][6]+'\t'+cross_CNV_SNP[i][7]+'\t'+cross_CNV_SNP[i][8]+'\n'
		#print(cross_CNV_SNP[i])
	txtName = "cross_CNV_SNP.txt"
	f = open(txtName,'w')
	f.write(str1)
	f.close()
def getNotCrossCNVSnp(not_cross_CNV_SNP):
	str1 = "CHROM"+"\t"+"POS"+"\t"+"REF"+"\t"+"ALT"+"\t"+"AB"+"\t"+"AF"+"\t"+"DP"+"\t"+"RO"+"\t"+"AO"+"\n"
	for i in range(len(not_cross_CNV_SNP)):
		str1 += not_cross_CNV_SNP[i][0]+'\t'+not_cross_CNV_SNP[i][1]+'\t'+not_cross_CNV_SNP[i][2]+'\t'+not_cross_CNV_SNP[i][3]+'\t'+not_cross_CNV_SNP[i][4]+'\t'+not_cross_CNV_SNP[i][5]+'\t'+not_cross_CNV_SNP[i][6]+'\t'+not_cross_CNV_SNP[i][7]+'\t'+not_cross_CNV_SNP[i][8]+'\n'
		#print(not_cross_CNV_SNP[i])
	txtName = "netural_SNP.txt"
	f = open(txtName,'w')
	f.write(str1)
	f.close()

cnv_interval = getCNV(file2)
cross_CNV_SNP,not_cross_CNV_SNP = getSNP(file1,cnv_interval)
getCrossCNVSNP(cross_CNV_SNP)
getNotCrossCNVSnp(not_cross_CNV_SNP)
