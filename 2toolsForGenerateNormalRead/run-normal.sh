./readgen -l 250 -r 250 -d 500 -s 20 -e 0  -c 1500 -I nor.sim.idx chr1.fa nor.sim n_left100.fq n_right100.fq &&

samtools faidx chr1.fa &&

speedseq align -o s100.snp -M 3 -R "@RG\tID:id\tSM:sample\tLB:lib" chr1.fa n_left100.fq n_right100.fq &&
speedseq var -o s100.snp chr1.fa s100.snp.bam &&

samtools sort s100.snp.bam normal.sorted
