##### script for assembling metaTs from OXR Yep 1 with IDBA

# set cast/bottle for files
CB=$1

if [[ -z $CB ]];
then echo "Missing input to run command. Please provide an ID to identify samples"; exit;
else echo "All variables provided. Proceeding to IDBA."; fi

# create sampleset directory and do into
mkdir ${CB}_idba && cd ${CB}_idba

# convert to single fasta format required by IDBA-UD
for j in {1..3}; \
do /data/app/idba_mt/idba-1.1.3/bin/fq2fa --merge --filter ../${CB}F${j}_QC_R1.fq ../${CB}F${j}_QC_R2.fq ${CB}F${j}_merged.fa; done
cat ${CB}F*_merged.fa > ${CB}_reads4idba-ud.fa 

# run IDBA-UD for initial assembly
/data/app/idba_mt/idba_ud -r ${CB}_reads4idba-ud.fa -o ${CB}_idba-ud --num_threads 12
/data/app/seqkit stats ./${CB}_idba-ud/contig.fa >${CB}_idba-ud.stats
echo "IDBA-UD complete"

# run IDBA-MT on contigs from IDBA-MT

# convert individual R1/R2 fastqs to fastas
for i in {1..3}; do /data/app/idba_mt/idba-1.1.3/bin/fq2fa ../${CB}F${i}_QC_R1.fq ${CB}F${i}_R1.fa; done
for i in {1..3}; do /data/app/idba_mt/idba-1.1.3/bin/fq2fa ../${CB}F${i}_QC_R2.fq ${CB}F${i}_R2.fa; done
cat ${CB}F*R1.fa > ${CB}_R1.fa
cat ${CB}F*R2.fa > ${CB}_R2.fa

/data/app/idba_mt/idba-mt -t ${CB}_R1.fa -f ${CB}_R2.fa -O ${CB}_idba-mt.fa -c ./${CB}_idba-ud/contig.fa -r 88
/data/app/seqkit stats ${CB}_idba-mt.fa >${CB}_idba-mt.stats

echo $CB " complete"
