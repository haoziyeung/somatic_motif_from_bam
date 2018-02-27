
open (INPUT, "<$ARGV[0]");

while(<INPUT>){
	chomp;
	@array = split "\t", $_;
	@gtdetails = split ":", $array[10];
	my $vaf = $gtdetails[5];
	$vaf =~ s/\%//;
	print "$_\n" if ($vaf >= 1 && $vaf <= 3 ) || $_=~/^#/;
}
