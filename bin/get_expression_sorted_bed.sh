
tail -n+2 rnaHtseqCountsall.txt | awk -F '\t' -v OFS='\t' '{print $1, ($6+$7)/2, ($8+$9)/2}' > rnaHtseqCountsall.cmp_ery.txt

time python ~/group/software/signorm/bin/vlookup.py -t rnaHtseqCountsall.cmp_ery.txt -m 1 -s gencode.vM4.annotation.bed -n 4 -o rnaHtseqCountsall.cmp_ery.bedmatched.tmp1.txt

cut -f2,3 rnaHtseqCountsall.cmp_ery.bedmatched.tmp1.txt > rnaHtseqCountsall.cmp_ery.bedmatched.tmp2.txt

paste gencode.vM4.annotation.bed rnaHtseqCountsall.cmp_ery.bedmatched.tmp2.txt | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6, ($7+10)/($8+10)}' | sort -k7,7n > rnaHtseqCountsall.cmp_ery.bedmatched.txt

cut -f1,2,3,4,5,6 rnaHtseqCountsall.cmp_ery.bedmatched.txt > gencode.vM4.annotation.cmp_ery_fc.bed

rm rnaHtseqCountsall.cmp_ery.bedmatched.tmp1.txt
rm rnaHtseqCountsall.cmp_ery.bedmatched.tmp2.txt
rm rnaHtseqCountsall.cmp_ery.txt
rm rnaHtseqCountsall.cmp_ery.bedmatched.txt


