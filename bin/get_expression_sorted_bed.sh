### get cmp_ery reads count table
tail -n+2 rnaHtseqCountsall.txt | awk -F '\t' -v OFS='\t' '{print $1, ($6+$7)/2, ($8+$9)/2}' > rnaHtseqCountsall.cmp_ery.txt

### get protein coding gene only
cat gencode.vM4.annotation.bed | awk -F '\t' -v OFS='\t' '{if ($5=="protein_coding") print $0}' > gencode.vM4.annotation.pc.bed

### sort cmp_ery reads count table by matched bed id
time python ~/group/software/signorm/bin/vlookup.py -t rnaHtseqCountsall.cmp_ery.txt -m 1 -s gencode.vM4.annotation.pc.bed -n 4 -o rnaHtseqCountsall.cmp_ery.bedmatched.tmp1.txt

### get expression signal only
cut -f2,3 rnaHtseqCountsall.cmp_ery.bedmatched.tmp1.txt > rnaHtseqCountsall.cmp_ery.bedmatched.tmp2.txt

### paste bed info & sort by fold-change between cmp & ery_fl
paste gencode.vM4.annotation.pc.bed rnaHtseqCountsall.cmp_ery.bedmatched.tmp2.txt | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6, ($7+10)/($8+10)}' | sort -k7,7n > rnaHtseqCountsall.cmp_ery.bedmatched.txt

### get sorted bed file
cut -f1,2,3,4,5,6 rnaHtseqCountsall.cmp_ery.bedmatched.txt > gencode.vM4.annotation.pc.cmp_ery_fc.bed

rm rnaHtseqCountsall.cmp_ery.bedmatched.tmp1.txt
rm rnaHtseqCountsall.cmp_ery.bedmatched.tmp2.txt
rm rnaHtseqCountsall.cmp_ery.txt
rm rnaHtseqCountsall.cmp_ery.bedmatched.txt


