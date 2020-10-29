# encoding:utf-8
'''
1.读入 diff_s2_s3.txt
2.读入inspos1.txt和inspos2.txt
3.生成cnvpos3.txt 
4.生成ref
5。截取一段写入inspos3.txt
'''
#python genCNV3.py diff_s2_s3.txt inspos1.txt inspos2.txt s33.sim chr1.fa
import sys
import random
assert len(sys.argv) == 6
file1,file2,file3,simFile,refFile = sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5]

def getDiff_s2_s3(file1):
	data = []
	with open(file1) as fb:
		lines = fb.readlines()
	for line in lines:
		line = line.strip().replace(" ","").split('\t')
		data.append([int(line[0]),int(line[1]),int(line[2])])
	return data
def getInsertPos(file2,file3):
	pos = []
	with open(file2) as fb:
		lines = fb.readlines()
	for line in lines:
		line = line.strip().split(" ")
		pos.append(int(line[1]))
	with open(file3) as fb:
		lines = fb.readlines()
	for line in lines:
		line = line.strip().split(" ")
		pos.append(int(line[1]))
	
	pos.sort()
	return pos
def getCnvPos(snpPos,CNVPos):
	candidateSnp = []
	for snp in snpPos:
		cross = False
		for cnv in CNVPos:
			if(cnv-9999<=snp[0]<=cnv or cnv-9999<=snp[0]+9999 < cnv):
				cross = True
				break
		if(cross == False):
			candidateSnp.append(snp)
	return candidateSnp
def getRef1(refFile):
	with open(refFile) as fb:
		lines = fb.readlines()
		ref = ""
	for line in lines[1:]:
		line = line.strip()
		ref += line
	return ref


def getRef(refFile,simFile):
	dict1={80:'A',81:'C',82:'G',83:'T'}
	ref = getRef1(refFile)
	ref=list(ref)
	with open(simFile) as fb:
		lines = fb.readlines()
	pos_mutation=[]
	for line in lines[1:]:
		line = line.strip().replace(" ","").split('\t')
		pos_mutation.append([int(line[0]),int(line[1])])
		if(int(line[1])>4 and int(line[1])<90):
			ref[int(line[0])-1] = dict1[int(line[1])]
	return ref,pos_mutation
def getInsposTxt(ref,candidateSnp):
	selectNum= 5
	random_list = []
	str1 = ""
	for i in range(selectNum):
	    random_list.append(random.randint(0, len(candidateSnp)))
	for i in range(len(random_list)):
		pos = candidateSnp[random_list[i]][0]
		copyNum = random.randint(3, 7)
		copyPart = "".join(ref[pos-1:pos+9999]) *copyNum
		str1 += str(i+1)+' '+str(pos+9999)+" "+str(copyNum*10000)+" 1 1 0 "+copyPart+'\n'
	f = open("inspos3.txt",'w')
	f.write(str1)



snpPos=getDiff_s2_s3(file1)
CNVPos=getInsertPos(file2,file3)
candidateSnp=getCnvPos(snpPos,CNVPos)
ref,pos_mutation = getRef(refFile,simFile)
getInsposTxt(ref,candidateSnp)

