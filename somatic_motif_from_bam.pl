#!/apps/perl/perl5/perls/perl-5.10.1/bin/perl
use strict;
use warnings;

use lib '/home/yanghao/perl5/lib/perl5/';
use Term::ANSIColor;
use Getopt::Long qw/GetOptions/;
use Cwd qw/getcwd abs_path/;
use File::Basename;
use Parallel::ForkManager;
use feature 'say';

my $base = dirname abs_path $0;
my ($sample,$bam);
my $eof = <<EOF;

----------------------------------------------------------
        {体细胞突变检测流程}
. 从任意一个样本的bam文件以及NA12878，得到somatic motif突变谱信息

使用方法：

perl $0 -s <sample ID> -b <sample BAM>

---------------------------------------------------------

EOF

GetOptions(
    's=s' => \$sample,
    'b=s' => \$bam,
);

unless($sample && $bam && -e $bam){
    say colored($eof, "bright_yellow");
    exit;
}

my $pm = new Parallel::ForkManager (24);

for my $chr (1..22,"X","Y"){
	my $pid = $pm -> start and next;
	calling($chr);
	$pm -> finish;
}
$pm -> wait_all_children;

`/gpfs/users/yanghao/software/anaconda2/bin/bcftools concat $sample.NA12878.*.somatic.raw.vcf.gz > $sample.vcf`;
`/gpfs/users/yanghao/software/annovar/table_annovar.pl $sample.vcf /gpfs/users/yanghao/database/annotation/humandb/ -buildver hg19 -protocol gnomad_genome,gnomad_exome -operation f,f -vcfinput -nastring . --thread 6 --maxgenethread 6 -remove`;
`grep -v ^# $sample.vcf.hg19_multianno.vcf|grep SOMATIC|perl -lane'(\$a)=\$F[7]=~m|gnomAD_genome_ALL=([^;]+)|;\$a=\$a eq "."?0:\$a;(\$b)=\$F[7]=~m|gnomAD_genome_ALL=([^;]+)|;\$a=\$b eq "."?0:\$b;print if \$a <= 0.001 && \$b<=0.001' > $sample.filter.vcf`;
`/gpfs/users/yanghao/software/anaconda2/bin/Rscript $base/SomaticMotif.R $sample.filter.vcf $sample`;


sub calling{
	my $chr = shift @_;
	`bash -c \"java -jar /gpfs/users/yanghao/software/VarScan.v2.3.9.jar somatic <\(/gpfs/bin/samtools-1.3.1/samtools mpileup -L 10000 -d 10000 -q 1 -f /gpfs/genomedb/b37/human_g1k_v37.fasta -l /gpfs/users/yanghao/bed/xgen.idt.bed.anno.sorted.bed -r $chr /gpfs/users/yanghao/database/giab/NA12878/NIST-hg001-7001-ready.bam\) <\(/gpfs/bin/samtools-1.3.1/samtools mpileup -L 10000 -d 10000 -q 1 -f /gpfs/genomedb/b37/human_g1k_v37_decoy.fasta -l /gpfs/users/yanghao/bed/xgen.idt.bed.anno.sorted.bed -r $chr $bam\) $sample.NA12878.$chr.somatic --min-var-freq 0.009 --output-vcf 1 --strand-filter 1\"`;
	`/gpfs/users/yanghao/software/anaconda2/bin/bgzip -c $sample.NA12878.$chr.somatic.snp.vcf > $sample.NA12878.$chr.somatic.snp.vcf.gz`;
	`/gpfs/users/yanghao/software/anaconda2/bin/tabix -p vcf -f $sample.NA12878.$chr.somatic.snp.vcf.gz`;
	`/gpfs/users/yanghao/software/anaconda2/bin/bgzip -c $sample.NA12878.$chr.somatic.indel.vcf > $sample.NA12878.$chr.somatic.indel.vcf.gz`;
	`/gpfs/users/yanghao/software/anaconda2/bin/tabix -p vcf -f $sample.NA12878.$chr.somatic.indel.vcf.gz`;
	`/gpfs/users/yanghao/software/anaconda2/bin/bcftools concat -a $sample.NA12878.$chr.somatic.snp.vcf.gz $sample.NA12878.$chr.somatic.indel.vcf.gz | /gpfs/users/yanghao/software/anaconda2/bin/bgzip -c >  $sample.NA12878.$chr.somatic.raw.vcf.gz`;
	`/gpfs/users/yanghao/software/anaconda2/bin/tabix -p vcf -f $sample.NA12878.$chr.somatic.raw.vcf.gz`;
}
