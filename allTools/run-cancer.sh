#生成fastq文件，合并，建索引
#生成read 2:!
./readgen -l 250 -r 250 -d 500 -s 20 -e 0  -c 200 -I nor.sim.idx chr1.fa nor.sim n_left100.fq n_right100.fq &&
./readgen -l 250 -r 250 -d 500 -s 20 -e 0  -c 100 -I s1.sim.idx chr1.fa s1.sim s1_left100.fq s1_right100.fq &&
./readgen -l 250 -r 250 -d 500 -s 20 -e 0  -c 100 -I s2.sim.idx chr1.fa s2.sim s2_left100.fq s2_right100.fq &&
./readgen -l 250 -r 250 -d 500 -s 20 -e 0  -c 100 -I s3.sim.idx chr1.fa s3.sim s3_left100.fq s3_right100.fq &&

#合并read tumor
./mergeread n_left100.fq s1_left100.fq s2_left100.fq s3_left100.fq s_left100.fq &&
./mergeread n_right100.fq s1_right100.fq s2_right100.fq s3_right100.fq  s_right100.fq &&

#index ref
samtools faidx chr1.fa &&

#speedseq处理
speedseq align -o s100.snp -M 3 -R "@RG\tID:id\tSM:sample\tLB:lib" chr1.fa s_left100.fq s_right100.fq &&
speedseq sv -o s100.snp -B s100.snp.bam -S s100.snp.splitters.bam -D s100.snp.discordants.bam -R chr1.fa -d &&
speedseq var -o s100.snp chr1.fa s100.snp.bam &&

gunzip s100.snp.sv.vcf.gz &&
gunzip s100.snp.vcf.gz &&
#生成all_snp.txt
python getSNPInf.py  s100.snp.vcf &&
#生成all_CNV.txt
python getCNVInf.py  s100.snp.sv.vcf &&
#找到横跨CNV的SNP，生成cross_CNV_SNP.txt 和netural_SNP.txt
python getCNV_SNP.py    all_snp.txt    all_CNV.txt 
#中性的位点中去除异常点，得到real_netural_SNP.txt
#python removeOutlierFromNetural.py    netural_SNP.txt

