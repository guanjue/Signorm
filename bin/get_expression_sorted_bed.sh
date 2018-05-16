### get cmp_ery reads count table
tail -n+2 rnaHtseqCountsall.txt | awk -F '\t' -v OFS='\t' '{print $1, ($6+$7)/2, ($8+$9)/2}' > rnaHtseqCountsall.cmp_ery.txt

### get protein coding gene only
cat gencode.vM4.annotation.bed | awk -F '\t' -v OFS='\t' '{if ($5=="protein_coding") print $0}' > gencode.vM4.annotation.pc.bed

### sort cmp_ery reads count table by matched bed id
time python ~/group/software/signorm/bin/vlookup.py -t rnaHtseqCountsall.cmp_ery.txt -m 1 -s gencode.vM4.annotation.pc.bed -n 4 -o rnaHtseqCountsall.cmp_ery.bedmatched.tmp1.txt

### get expression signal only
cut -f2,3 rnaHtseqCountsall.cmp_ery.bedmatched.tmp1.txt > rnaHtseqCountsall.cmp_ery.bedmatched.tmp2.txt

### paste bed info & sort by fold-change between cmp & ery_fl
#paste gencode.vM4.annotation.pc.bed rnaHtseqCountsall.cmp_ery.bedmatched.tmp2.txt | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8, ($7+100)/($8+100)}' | sort -k9,9n > rnaHtseqCountsall.cmp_ery.bedmatched.txt
#paste gencode.vM4.annotation.pc.bed rnaHtseqCountsall.cmp_ery.bedmatched.tmp2.txt | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6, -(($7)/($3-$2)+100)/(($8)/($3-$2)/21653542*51105650+100)}' | sort -k7,7n > rnaHtseqCountsall.cmp_ery.bedmatched.txt
#paste gencode.vM4.annotation.pc.bed rnaHtseqCountsall.cmp_ery.bedmatched.tmp2.txt | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6, -($8)/($3-$2)}' > rnaHtseqCountsall.cmp_ery.bedmatched.txt

paste gencode.vM4.annotation.pc.bed rnaHtseqCountsall.cmp_ery.bedmatched.tmp2.txt > rnaHtseqCountsall.cmp_ery.bedmatched.txt

Rscript sort_exp_ery.R



rm rnaHtseqCountsall.cmp_ery.bedmatched.tmp1.txt
rm rnaHtseqCountsall.cmp_ery.bedmatched.tmp2.txt
rm rnaHtseqCountsall.cmp_ery.txt
#rm rnaHtseqCountsall.cmp_ery.bedmatched.txt


