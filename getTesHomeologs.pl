#!/share/apps/bin/perl
my $VERSION = '0.1';
my $lastModif = '2017-6-19';
use strict;
use warnings;
use Getopt::Long;
use Statistics::Descriptive;
use File::Basename;
use Data::Dumper;
use Text::CSV;

my $homeologFile;
#my $outBed;
&GetOptions("h|help" => \&help,
            "f=s"    => \$homeologFile);
@ARGV or &help;
&main(\@ARGV);

sub main {
	my $self = {};
	bless $self ;
	$self->{inputFiles} = shift;
	$self->setOptions();
	$self->readEnvTeFile();
	$self->readHomeologs();
	exit(0);
}

sub setOptions {
	my $self = shift;
	-e $homeologFile or die "Cannot find file: " . $homeologFile . "\n";
	$self->{param}->{homeologFile} = $homeologFile;
#	$self->{param}->{outBed} = $outBed;
	return 1;
}

sub readEnvTeFile {
	my $self = shift;
	foreach my $file (@{$self->{inputFiles}}){
		my $csv = Text::CSV->new({sep_char => "\t", 
								  quote_char => ''});
		open(my $data, "<".$file) or die "Cannot open file: " . $file . "\n";
		my $k = 0;
		while (my $col = $csv->getline($data)){
			++$k==1 and next;
			# pick the left or right TE depending on the strand of the gene to consider the closest TE in 5'
			my $colIndex = ($col->[4] eq '+' ? 6 : 11); 
			my ($fam,$sf) = ($col->[$colIndex],'NA');
			if($fam=~/^([A-Z]{3})_.*$/){
				$sf = $1;
				$fam=~s/\.\d+$//;
			}
			my $promoter = $fam eq 'NA' ? ['NA','NA'] : getPromoterCoords($col);
			$self->{data}->{closestTe}->{$col->[1]} = {'fam' => $fam, 
			                                           'sf' => $sf,
			                                           'chr' => $col->[0],
			                                           'promoterStart' => $promoter->[0],
			                                           'promoterStop' => $promoter->[1],
			                                           'promoterStrand' => $col->[4]};
		}
		close $data;
		print STDERR $file . "->ok\n";
	}
	print STDERR '->Nb of genes analyzed: ', scalar(keys %{$self->{data}->{closestTe}}), "\n";
	return 1;
}


sub getPromoterCoords {
	my $col = shift;
	my ($start, $stop);
	if($col->[4] eq '+'){
		$start = $col->[8]-10000;
		$stop = $col->[3];
	}
	else{
		$start = $col->[2];
		$stop = $col->[13]+10000;
	}
	$start<1 and $start=1;
	return [$start, $stop];
}


sub readHomeologs {
	my $self = shift;
	$self->printHeader();
	my $csv = Text::CSV->new({sep_char => ",", 
	                          quote_char => '"'});
	open(my $file, "<".$self->{param}->{homeologFile}) or die "Cannot open CSV file: " . $self->{param}->{homeologFile} . "\n";
	my $k = 0;
	while (my $col = $csv->getline($file)){
		if($col->[5] eq '1:1:1'){
			my ($h,$score) = ({},{});
			foreach('fam','sf'){
				$score->{$_} = 0;
				$h->{$_} = [$self->{data}->{closestTe}->{$col->[6]}->{$_},
			                $self->{data}->{closestTe}->{$col->[7]}->{$_},
			                $self->{data}->{closestTe}->{$col->[8]}->{$_}];
			    if($h->{$_}->[0] eq $h->{$_}->[1]){ $score->{$_}+=1; }
			    if($h->{$_}->[0] eq $h->{$_}->[2]){ $score->{$_}+=1; }
			    if($h->{$_}->[1] eq $h->{$_}->[2]){ $score->{$_}+=1; }
			}
			
			print join("\t", @{$col}[6,7,8], $h->{sf}->[0], $h->{fam}->[0], 
			                                 $h->{sf}->[1], $h->{fam}->[1], 
			                                 $h->{sf}->[2], $h->{fam}->[2],
			                                 $score->{sf},
			                                 $score->{fam},
			                                 $self->{data}->{closestTe}->{$col->[6]}->{chr},
			                                 $self->{data}->{closestTe}->{$col->[6]}->{promoterStart},
			                                 $self->{data}->{closestTe}->{$col->[6]}->{promoterStop},
			                                 $self->{data}->{closestTe}->{$col->[6]}->{promoterStrand},
			                                 $self->{data}->{closestTe}->{$col->[7]}->{chr},
			                                 $self->{data}->{closestTe}->{$col->[7]}->{promoterStart},
			                                 $self->{data}->{closestTe}->{$col->[7]}->{promoterStop},
			                                 $self->{data}->{closestTe}->{$col->[7]}->{promoterStrand},
			                                 $self->{data}->{closestTe}->{$col->[8]}->{chr},
			                                 $self->{data}->{closestTe}->{$col->[8]}->{promoterStart},
			                                 $self->{data}->{closestTe}->{$col->[8]}->{promoterStop},
			                                 $self->{data}->{closestTe}->{$col->[8]}->{promoterStrand}
			                 ),"\n";
		}
	}
	close $file;
	print STDERR "readHomeologs ->done\n";
	return 1;
}


sub printHeader {
	my $self = shift;
	my @header = qw(idA idB idD 
	                teSfamA teFamA 
	                teSfamB teFamB 
	                teSfamD teFamD 
	                sfamScore famScore 
	                chrA promoterStartA promoterStopA promoterStrandA
	                chrB promoterStartB promoterStopB promoterStrandB 
	                chrD promoterStartD promoterStopD promoterStrandD);
	print join("\t", @header), "\n";
	return 1;
}

sub help() {
	my $prog = basename($0) ;
	print STDERR <<EOF ;
### $prog $VERSION ###
# LAST MODIF:  $lastModif
# PURPOSE:     This script is used to get TEs in the promoter of homeogous triplets

USAGE: 
       $prog  -f <homeologFile>  <list of file *closestTe.csv>

       ### OPTIONS ###
       -h, --help                        Print this help
       -f <file>                         File with homeologs ("wheat.homeolog_groups.release.nonTE_segmental.csv")


EOF
	exit(1) ;
}
