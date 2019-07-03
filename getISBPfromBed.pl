#!/share/apps/bin/perl

use strict;
use warnings;
use Getopt::Long ;
use File::Basename ;

my $VERSION = '2.0';
my $lastmodif = '2019/07/03';
my $time = localtime;
my $overlap = 99;
my $baseFile;

&GetOptions ( "h|help"      => \&help,
	          "o|overlap=i" => \$overlap) ;

$#ARGV==0 or &help;
&main($ARGV[0]);

# Main routine

sub main{
	my $self={};
	bless $self;
	$self->{inputFile}=shift@_;
	$self-> setOption();							 # Setting the options (overlap, output file format)
	unless(-e $self->{inputFile}){
		print STDERR "[Error] Could not find".$self->{inputFile}."\n";
		exit;
	}
	$self-> getAllIsbp(); 							 # getting ISBP flanking postions (+75/-75bp)
	if($overlap<150){								 # ISBP with allowed overlap (0 <= overlap < 150bp)
		$self-> getSelectedIsbpAccordingToOverlap(); # Saved ISBP (when overlap is under the threshold)
			$self-> printResultsWithOverlap();       # Print the results (ie ISBP posistions, bed file format)
	}
	else{											 # All ISBP ( overlap == 150 )
			$self-> printResultsWithoutOverlap();    # Print the results (ie ISBP posistions, bed file format)
	}		
 	exit(0);
}



# Functions

sub setOption{
	my $self=shift;
	$overlap>=0 or die "Overlap must be >= 0 (option -o)\n";
	$overlap<151 or die "Overlap max value : 150 bp (option -o)\n";
	$self->{param}->{overlap} = $overlap;
	$self->{param}->{isbpSize} = $overlap;
	return 1;
}


sub getAllIsbp{
	my$self=shift;
	open (F,"<".$self->{inputFile}) or die ("Can't open the reference bed input file\n" . 
	                                        "For more information print the help (option -h|--help)\n" . 
	                                        $self->{inputFile} . "\n"); # Input file (bed file format)

	#Assessing ISBP postions from TEs insertion sites
	while(<F>){
		chomp;
		my ($seqId,$start,$stop,$teId) = split/\t/;
		foreach (['start',$start], ['stop',$stop]) {
			my $type = $_->[0];
			my $pos = $_->[1];
			$pos<75 and next;
			my $isbpId = ($type eq 'start' ? $teId.'_5prime;' : $teId.'_3prime;');
			$self->{isbp}->{$seqId}->{all_ISBP}->{$isbpId}->{start}=$pos-75;
			$self->{isbp}->{$seqId}->{all_ISBP}->{$isbpId}->{stop}=$pos+75;
			$self->{isbp}->{$seqId}->{all_ISBP}->{$isbpId}->{teId}=$teId;
		}
	}
	print STDERR "getAllIsbp -> done\n" ;
	foreach my $seqId (sort keys %{$self->{isbp}}){
		print STDERR $seqId, ":", scalar keys %{$self->{isbp}->{$seqId}->{all_ISBP}}, " ISBPs analyzed\n";
	}
	close(F);
	return 1;
}

sub getSelectedIsbpAccordingToOverlap{
	my$self=shift;
	foreach my $seqId ( sort keys (%{$self->{isbp}}) ){ # For each chromosome
		my ($k, $kSelected) = (0, 0);
		foreach my $isbpId (sort { $self->{isbp}->{$seqId}->{all_ISBP}->{$a}->{start} <=> $self->{isbp}->{$seqId}->{all_ISBP}->{$b}->{start} } keys (%{ $self->{isbp}->{$seqId}->{all_ISBP} })){ # For each ISBP start position sorted (ISBP are sorted by start positions)
			$k++;
			if($k==1){
				$self->{isbp}->{$seqId}->{selected}->{++$kSelected}->{selectedIsbpId} = $isbpId;
				next;		
			}
			my $previousIsbpId = $self->{isbp}->{$seqId}->{selected}->{$kSelected}->{selectedIsbpId};
			my $startCurrent = $self->{isbp}->{$seqId}->{all_ISBP}->{$isbpId}->{start};
			my $stopPrevious = $self->{isbp}->{$seqId}->{all_ISBP}->{$previousIsbpId}->{stop};
			# in case of overlap
			if ($startCurrent >= $stopPrevious - $self->{param}->{overlap} + 1 ){ 
				my$resStop=$stopPrevious - $self->{param}->{overlap} + 1;
				$self->{isbp}->{$seqId}->{selected}->{++$kSelected}->{selectedIsbpId} = $isbpId;
			}
			# no overlap
			else {
				my$resStopnotok=$stopPrevious - $self->{param}->{overlap} + 1;
				push(@{$self->{isbp}->{$seqId}->{selected}->{$kSelected}->{rejectedIsbpId}}, $isbpId);
			}
		}
	}	
	return 1;	
}


sub printResultsWithOverlap {
	my $self = shift;
	foreach my $seqId (sort keys(%{$self->{isbp}})){
		foreach my $k (sort { $a <=> $b } keys %{$self->{isbp}->{$seqId}->{selected}}){
			my $isbpId = $self->{isbp}->{$seqId}->{selected}->{$k}->{selectedIsbpId};
			my $rejectedIds="";
			if ( exists  $self->{isbp}->{$seqId}->{selected}->{$k}->{rejectedIsbpId}  ) {
				$rejectedIds = join("", @{$self->{isbp}->{$seqId}->{selected}->{$k}->{rejectedIsbpId}});
			}		
			my @out = ($seqId, 
			           $self->{isbp}->{$seqId}->{all_ISBP}->{$isbpId}->{start}, 
			           $self->{isbp}->{$seqId}->{all_ISBP}->{$isbpId}->{stop},
			           $isbpId . $rejectedIds);		
			print join("\t", @out), "\n";
		}
	}
	return 1;
}




sub printResultsWithoutOverlap {
	my $self = shift;
	foreach my $seqId (sort keys(%{$self->{isbp}})){
		foreach my $isbpId (sort { $self->{isbp}->{$seqId}->{all_ISBP}->{$a}->{start} <=> $self->{isbp}->{$seqId}->{all_ISBP}->{$b}->{start} } keys %{ $self->{isbp}->{$seqId}->{all_ISBP} } ){
			print $seqId."\t".$self->{isbp}->{$seqId}->{all_ISBP}->{$isbpId}->{start},"\t",
			      $self->{isbp}->{$seqId}->{all_ISBP}->{$isbpId}->{stop},"\t",
			      $isbpId,"\n";
		}
	}
	return 1;
}


sub help {
my $prog = basename($0);
print STDERR <<EOF ;
# # # $prog # # #
#
# CREATED:    2016/12/28
# LAST MODIF: $lastmodif
# AUTHOR:     Romain DE OLIVEIRA (INRA Clermont-Ferrand)
# VERSION:    $VERSION
# PURPOSE:    This script is used to design 150 bp ISBPs (=junctions between TE/TE or TE/LowCopyDNA)
#             with 75 nt spanning each side of the junction.
#             Input: bed file of the TE positions. Must be ordered by chromosome and start position 
#             (e.g., sort -k1,1 -k2,2n in.bed > in.sorted.bed for BED files)

USAGE:
             $prog  [options] bed_file

OPTIONS:
          -h               Print this help
          -o <integer>     Overlap allowed between ISBP (in bp) [Default: 99 bp]
EOF
        exit(1) ;
}

__END__

