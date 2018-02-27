inputfile=$1
perl /gpfs/users/yanghao/project/somatic_motif_from_bam/SomaticExtract.pl $inputfile > $inputfile.somatic
echo "Sample chr pos ref alt" | sed -e 's/ /\t/g'  > $inputfile.somatic.forMotif
cat $inputfile.somatic |grep -v "#"| cut -f 1,2,4,5 |awk -v input=$inputfile '{OFS="\t";print input, "chr"$1, $2, $3, $4}' >> $inputfile.somatic.forMotif
