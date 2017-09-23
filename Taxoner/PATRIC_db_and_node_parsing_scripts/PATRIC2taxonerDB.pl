use strict; use warnings;
use Cwd;
use Cwd 'abs_path';

my $verbose = 0;
my $version = "1.2";

my %p2t = ();	#stores PATRIC ID to taxon id (hash) 
my $maxdbid = 0; #highest db index
my @genomeindex; #store genome index
my @contig; #store contig indexes

my @fasta; #store fasta files

my $argu = $ARGV[0];
my $gindex = $ARGV[2];

my $dir = getcwd;
my $abs_path = abs_path();

my $GB = 4;
my $GBchars = 1048576 * 1000 * $GB;
my $printedChars = 0;
my $mb = 0;

my $out = $ARGV[1];
my $currDB = 0;
my $currfile;

if(!defined $ARGV[0] || !defined $ARGV[1] || !defined $ARGV[2]) {
	PrintUsage();
	exit;
}

if(defined $ARGV[2]) {
	if($ARGV[2] =~ /(^\d|\d+$)/) {
		$GB = $ARGV[2];
	}

	else {
		print "$ARGV[2] is not an integer\n";
	}
}


ImportIDs();
GetFastaFiles();
ImportFastaFiles();

unless(-d $out) {
	mkdir($out);
}

$currfile = $out."/".$currDB.".fasta";
open OUT, ">$currfile" or die;


for(my $f = 0; $f <= $#fasta; $f++) {
	my @temp = split /\//, $fasta[$f];

	if($temp[$#temp] =~ m/(.+)(\.fna$)/) {
		printf("Parsing file [ %d / %d ]\n", $f + 1, $#fasta + 1);
		my $pgenome = $1;

		if($printedChars >= $GBchars) {
			$mb = $printedChars / (1000 * 1000);
			print "Printed $mb MB to $currfile\n";

			close(OUT);
			$currDB++;
			$currfile = $out."/".$currDB.".fasta";
			open OUT, ">$currfile" or die;
			$printedChars = 0;
		}

		if(defined $p2t{$pgenome}) {
			open FH, $fasta[$f] or die;
			while(<FH>) {
				chomp($_);
				if($_ =~ m/(^\>)(.+)/) {
					my @toks = split " ", $2;

					if(!defined $p2t{$toks[0]}) {
						push @contig, "$toks[0]\t$p2t{$pgenome}\t$maxdbid\t$pgenome";
						$p2t{$toks[0]} = $p2t{$pgenome};
					}

					else {
						if($verbose == 1) {print "Error, found duplicate with ID: $pgenome\n";}
					}
					
					my $title = sprintf("%s", ">$maxdbid;$p2t{$pgenome}");
					print OUT "$title\n";
					$printedChars = $printedChars + length($title) + 2;

					$maxdbid++;
				}

				else {
					if(length($_) > 0) {
						print OUT "$_\n";
						$printedChars = $printedChars + length($_) + 1;
					}
				}
			}
		}
	}
}

close(OUT);

$mb = $printedChars / (1000 * 1000);
print "Printed $mb MB to $currfile\n";

open OUT, ">PATRIC_final_genome_index.txt" or die;

foreach(@genomeindex) {
	print OUT "$_\n";
}

foreach(@contig) {
	print OUT "$_\n";
}

close(OUT);


############ SUBS #############

sub ImportFastaFiles {
	open FH, "fasta_files.txt" or die;

	while(<FH>) {
		chomp($_);
		if($_ =~ m/(^\.\/)(.+)/) {
			push @fasta, $argu."/".$2;
		}
	}

	close(FH);

	
}

sub ImportIDs {
	open FH, $gindex or die;
	@genomeindex = <FH>;
	close(FH);
	chomp(@genomeindex);

	for(my $i = 1; $i <= $#genomeindex; $i++) {
		my @line = split "\t", $genomeindex[$i];
		if(defined $line[0] && defined $line[1] && defined $line[2]) {
			if($line[2] > $maxdbid) {
				$maxdbid = $line[2];
			}

			$p2t{$line[0]} = $line[1];
		}
	}

	$maxdbid++;
}

sub GetFastaFiles {
	my $folder = $abs_path."/".$argu;
	print "Getting fasta file paths\n";
	chdir($folder);
	system("find -follow | grep \".fna\" > $abs_path/fasta_files.txt");
	print "Finished seaching for fasta files\n";
	chdir($abs_path);
}

sub PrintUsage {
	print "Create taxoner database from PATRIC datasets\n";
	print "Program: $0\n";
	print "Version: v$version\n\n";

	print "Usage: $0 <PATRIC input> <OUTPUT folder> <Index file> <optional: genom size>\n";
	print "===============================================================\n\n";

	print "PATRIC input\tfolder that contains fasta files from PATRIC database\n";
	print "\t\tFasta files can be in any subdir inside specified folder\n";
	print "OUTPUT folder\tFolder name, where results should be printed\n";
	print "Index file\tIndex file created from PATRIC (something like: PATRIC_genomeIndex.txt)\n";
	print "Genome size\tOptional, set the size of taxoner DB size (in Gb). Default: $GB\n";
}
