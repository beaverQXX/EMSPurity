#g++ MergeRead.cpp -o mergeread &&
#g++ NorSim_2.0.cpp -o norsim &&
#g++ TumSim_2.0.cpp -o tumsim &&
#g++ ReadGen_2.0.cpp -o readgen &&
#生成正常样本sim
./norsim -r 0.001 -X 0 -D 0 -B 0.333330 -I inspos.txt -A nor_AB.idx chr1.fa nor.sim &&


#生成肿瘤样本sim s1初代克隆
./tumsim -r 0.00007 -X 0 -D 0 -B 0.333330 -b 0 -A 0 -p 0 -n 0 -l 0 -I inspos.txt chr1.fa nor.sim nor_AB.idx s11.sim &&
#generate diff_nor_S1
python removeDup.py nor.sim s11.sim diff_nor_s1.txt &&
#生成拷贝数变异克隆s1:inspos1.txt cnvpos_s1.txt
#g++ GeneCNVs1.cpp -o gencnvs1 &&
./gencnvs1 diff_nor_s1.txt chr1.fa &&
#生成肿瘤样本sim s1初代克隆 cnv 
./tumsim -r 0 -X 0 -D 0 -B 0 -b 0 -A 0 -p 0 -n 0 -l 0 -I inspos1.txt chr1.fa s11.sim s11_AB.idx s1.sim &&



#生成肿瘤样本sim s2为s1的子克隆
./tumsim -r 0.00002 -X 0 -D 0 -B 0.333330 -b 0 -A 0 -p 0 -n 0 -l 0 -I inspos.txt chr1.fa s1.sim s1_AB.idx s22.sim &&
#generate diff_S1_S2
python removeDup.py s1.sim s22.sim diff_s1_s2.txt &&
#生成拷贝数变异克隆s1子克隆s2:inspos2.txt cnvpos_s2.txt
#g++ GeneCNVs2.cpp -o gencnvs2 &&
./gencnvs2 diff_s1_s2.txt cnvpos_s1.txt chr1.fa &&
#生成肿瘤样本sim s2子克隆 cnv 
./tumsim -r 0 -X 0 -D 0 -B 0 -b 0 -A 0 -p 0 -n 0 -l 0 -I inspos2.txt chr1.fa s22.sim s22_AB.idx s2.sim &&


#生成肿瘤样本sim s3为s2的子克隆
./tumsim -r 0.00001  -X 0 -D 0 -B 0.333330 -b 0 -A 0 -p 0 -n 0 -l 0 -I inspos.txt chr1.fa s2.sim s2_AB.idx s33.sim &&
#generate diff_S2_S3
python removeDup.py s2.sim s33.sim diff_s2_s3.txt &&
#生成拷贝数变异克隆s1子克隆s2:inspos3.txt
python genCNV3.py diff_s2_s3.txt inspos1.txt inspos2.txt s33.sim chr1.fa &&
#生成肿瘤样本sim s2子克隆 cnv 
./tumsim -r 0 -X 0 -D 0 -B 0 -b 0 -A 0 -p 0 -n 0 -l 0 -I inspos3.txt chr1.fa s33.sim s33_AB.idx s3.sim &&



#生成肿瘤样本sim s4为s3的子克隆
./tumsim -r 0.00001 -X 0 -D 0 -B 0.333330 -b 0 -A 0 -p 0 -n 0 -l 0 -I inspos.txt chr1.fa s3.sim s3_AB.idx s44.sim &&
#generate diff_S3_S4
python removeDup.py s3.sim s44.sim diff_s3_s4.txt &&
#生成拷贝数变异克隆s4子克隆s4:inspos4.txt
python genCNV4.py diff_s3_s4.txt inspos1.txt inspos2.txt inspos3.txt s44.sim chr1.fa &&
#生成肿瘤样本sim s4子克隆 cnv
./tumsim -r 0 -X 0 -D 0 -B 0 -b 0 -A 0 -p 0 -n 0 -l 0 -I inspos4.txt chr1.fa s44.sim s44_AB.idx s4.sim 

