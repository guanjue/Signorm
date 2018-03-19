
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed
chmod 777 bigBedToBed

for i in {1..20}
do
        echo $i
        #wget 'http://bx.psu.edu/~yuzhang/tmp/pknorm_16lim_ref1mo_new.'$i'.bb' ./
done

mkdir ideas_bb
mv *bb ideas_bb/

for i in {1..20}
do
        echo $i
        ./bigBedToBed 'ideas_bb/pknorm_16lim_ref1mo_new.'$i'.bb' 'ideas_bb/pknorm_16lim_ref1mo_new.'$i'.bed'
done


