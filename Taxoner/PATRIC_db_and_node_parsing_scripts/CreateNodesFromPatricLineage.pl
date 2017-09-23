use strict; use warnings;

my %id = ();
my @taxon;

open FH, $ARGV[0] or die;
my @data = <FH>;
close(FH);
chomp(@data);

for(my $i = 1; $i <= $#data; $i++) {
	my @line = split "\t", $data[$i];
	
	if(defined $line[10]) {
		my @toks = split /\;/, $line[10];
		for(my $j = $#toks; $j > 0; $j--) {
			if(!defined $id{$toks[$j]}) {
				push @taxon, $toks[$j];
				$id{$toks[$j]} = $toks[$j-1];
			}
		}
	}
}

open OUT, ">Patric_nodes.dmp" or die;

print OUT "131567\t1";
foreach(@taxon) {
	if(!defined $id{$_}) {
		print "Error: missing parent for $_?\n";
		print "Exiting\n";
		exit;
	}
	
	print OUT "\n$_\t$id{$_}";
}