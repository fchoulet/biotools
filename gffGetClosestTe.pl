#!/share/apps/bin/perl
my $VERSION = '0.1';
my $lastModif = '2017-2-17';
use strict;
use warnings;
use Getopt::Long;
use Statistics::Descriptive;
use File::Basename;
#use Data::Dumper;

my $gffVersion = 3;
my $optMode = 'closest';
my $optTag = 'ID';
my $range = 5000;
#my $bedSeqLength;
&GetOptions ( "h|help" => \&help,
              "m|mode=s" => \$optMode,
              "r|range=i" => \$range,
              "gv|gffversion=i" => \$gffVersion);
$#ARGV==0 or &help;
&main($ARGV[0]);

sub main {
	my $self = {};
	bless $self ;
	$self->{inputFile} = shift;
	$self->setOptions();
	$self->colNames();
	$self->readGff();
	$self->indexFeatures();
	$self->getTeAroundGenes();
	print STDERR $self->{inputFile} . " ->done\n";
	exit(0);
}


sub setOptions {
	my $self = shift;
#	$self->readBed();
	$self->{param}->{gffVersion} = $gffVersion;
	$self->{param}->{mode} = ($optMode =~/^closest$|^env$/ ? $optMode : 'closest') ;
	unless($range >= 200 and $range <= 20000){
		die "Range must be between 200 and 20000 bp\n";
	}
	$self->{param}->{range} = $range;
	return 1;
}


sub colNames {
	my $self = shift;
	if($self->{param}->{mode} eq 'closest'){ 
		@{$self->{colNames}} = qw(seqId geneId geneStart geneStop geneStrand 
		                          leftTeId leftTeCompo leftTeStatus leftTeStop leftDistanceToTe 
		                          rightTeId rightTeCompo rightTeStatus rightTeStart rightDistanceToTe);
	}
	else{
		@{$self->{colNames}} = qw(seqId geneId geneStart geneStop geneStrand 
		                          teId teCompo teStatus tePosition distanceToTe);
	}
	return 1;
}

# deprecated -> used to get the length of the sequence -> now replaced by "range" option
=head readBed
sub readBed {
	my $self = shift;
	-e $bedSeqLength or die "*** Enter BED file with sequence length with option: -bed <file> ***\n";
	open(BED, "<".$bedSeqLength) or die "Cannot open file: " . $bedSeqLength . "\n";
	while (<BED>){
		chomp;
		/^\s*$/ and next;
		/^#/ and next;
		my @col = split(/\t/,$_);
		$self->{seqLength}->{$col[0]} = $col[2]+1;
	}
	close BED;
	return 1;
}
=cut

sub readGff {
	my $self = shift;
	open(GFF, "<".$self->{inputFile}) or die "Cannot open file: " . $self->{inputFile} . "\n";
	my $kFt = 0;
	my ($kIndexStart, $kIndexStop) = 0;
	while (<GFF>){
		chomp;
		/^\s*$/ and next;
		/^#/ and next;
		$self->{_currentLine} = $_;
		my @col = split(/\t/,$_);
		scalar @col == 9 or die "Not a 9-column format at line " . $_ . "\n";
		$kFt++;

		my ($seqid,$source,$pritag,$start,$stop,$score,$strand,$frame,$note) = @col;
		$self->{feat}->{$seqid}->{$kFt}->{source} = $source;
		$self->{feat}->{$seqid}->{$kFt}->{pritag} = $pritag;
		$self->{feat}->{$seqid}->{$kFt}->{start} = $start;
		$self->{feat}->{$seqid}->{$kFt}->{stop} = $stop;
		$self->{feat}->{$seqid}->{$kFt}->{score} = $score;
		$self->{feat}->{$seqid}->{$kFt}->{strand} = $strand;
		$self->{feat}->{$seqid}->{$kFt}->{note} = $note;

		$self->{featBySeqidAndPritag}->{$seqid}->{$pritag}->{$kFt} = $start;
		$self->setParent2Child($seqid, $kFt);
	}
	close GFF;
	print STDERR "readGff ->done\n";
	return 1;
}


sub setParent2Child {
	my $self = shift;
	my ($seqid, $kFt) = @_;
	foreach (split(";", $self->{feat}->{$seqid}->{$kFt}->{note})){
		my ($tag,$value) = $self->getTagValue($_);
		if(uc($tag) eq 'ID'){
			$self->{feat}->{$seqid}->{$kFt}->{ID} and die "Feature with several IDs at line:\n" . $self->{_currentLine} . "\n";
			$self->{feat}->{$seqid}->{$kFt}->{ID} = $value;
#			$self->{id2k}->{$value} = $kFt;
		}
		elsif(uc($tag) eq 'PARENT'){
			$self->{feat}->{$seqid}->{$kFt}->{Parent} and die "Feature with several Parents at line:\n" . $self->{_currentLine} . "\n";
			$self->{feat}->{$seqid}->{$kFt}->{Parent} = $value;
		}
		elsif(uc($tag) eq 'COMPO'){
			$self->{feat}->{$seqid}->{$kFt}->{$tag} = (split/\s+/, $value)[0];
		}
		else{
			$self->{feat}->{$seqid}->{$kFt}->{$tag} = $value;
		}
	}
	$self->{feat}->{$seqid}->{$kFt}->{ID} or die "No ID found for feature:\n" . $self->{_currentLine} . "\n";
	return 1;
}


sub getTagValue {
	my $self = shift;
	my $field = shift;
	$field =~s/^\s+|\s+$//g;
	my ($tag,$value);
	if($self->{param}->{gffVersion}==3){
		($tag,$value) = split(/=/,$field,2);
	}
	else{
		($tag,$value) = $field =~/^(\S+) \"*(.*)\"*$/;
	}
	($tag and $value) or die "Could not understand: " . $field . "\n";
	return ($tag, $value);
}

 
sub indexFeatures {
	my $self = shift;
	foreach my $seqid (sort {$a cmp $b} keys %{$self->{feat}}){
		my $index = 0;
		foreach my $kFt (sort {$self->{feat}->{$seqid}->{$a}->{start} <=> $self->{feat}->{$seqid}->{$b}->{start}}
		                 keys %{$self->{feat}->{$seqid}}){
			$self->{index}->{start}->{$seqid}->{++$index} = $kFt;
			$self->{feat}->{$seqid}->{$kFt}->{index}->{start} = $index;
		}
		$index = 0;
		foreach my $kFt (sort {$self->{feat}->{$seqid}->{$a}->{stop} <=> $self->{feat}->{$seqid}->{$b}->{stop}}
		                 keys %{$self->{feat}->{$seqid}}){
			$self->{index}->{stop}->{$seqid}->{++$index} = $kFt;
			$self->{feat}->{$seqid}->{$kFt}->{index}->{stop} = $index;
		}
	}
	return 1;
}


sub getTeAroundGenes {
	my $self = shift;
	print join("\t", @{$self->{colNames}}), "\n";
	foreach my $seqid (sort {$a cmp $b} keys %{$self->{featBySeqidAndPritag}}){
#		$self->{seqLength}->{$seqid} or die "no feature \"region\" for seqid=" . $seqid . " =>exit\n";
		foreach my $kFt (sort {$self->{featBySeqidAndPritag}->{$seqid}->{gene}->{$a} <=> $self->{featBySeqidAndPritag}->{$seqid}->{gene}->{$b}}
		                 keys %{$self->{featBySeqidAndPritag}->{$seqid}->{gene}}){

			# set all values to NA.
			foreach (@{$self->{colNames}}) { $self->{print}->{$_} = 'NA'; }
			my $id         = $self->{feat}->{$seqid}->{$kFt}->{ID};
			my $geneStart  = $self->{feat}->{$seqid}->{$kFt}->{start};
			my $geneStop   = $self->{feat}->{$seqid}->{$kFt}->{stop};
			my $geneStrand = $self->{feat}->{$seqid}->{$kFt}->{strand};
			my $indexStart = $self->{feat}->{$seqid}->{$kFt}->{index}->{start};
			my $indexStop  = $self->{feat}->{$seqid}->{$kFt}->{index}->{stop};

			$self->{print}->{seqId} = $seqid;
			$self->{print}->{geneId} = $id;
			$self->{print}->{geneStart} = $geneStart;
			$self->{print}->{geneStop} = $geneStop;
			$self->{print}->{geneStrand} = $geneStrand;
			
			if($self->{param}->{mode} eq 'closest'){
				$self->searchClosestTe($seqid, $geneStart, $geneStop, $indexStart, $indexStop);
			}
			else{
				if($geneStrand eq '+'){ $self->getTeLeft($seqid, $geneStart, $indexStop); }
				else                  { $self->getTeRight($seqid, $geneStop, $indexStart); }
			}
		}
	}
	return 1;
}


sub searchClosestTe {
	my $self = shift;
	my ($seqid, $geneStart, $geneStop, $indexStart, $indexStop) = @_;
	my $kLeftTe = $self->searchByIndex($seqid, $geneStart, $indexStop, 'left');
	my $kRightTe = $self->searchByIndex($seqid, $geneStop, $indexStart, 'right');
	foreach(['left',$kLeftTe], ['right',$kRightTe]){
		$_->[1] eq 'NA' and next;
		$self->{print}->{$_->[0].'TeId'} = $self->{feat}->{$seqid}->{$_->[1]}->{ID};
		$self->{print}->{$_->[0].'TeCompo'} = $self->{feat}->{$seqid}->{$_->[1]}->{compo};
		$self->{print}->{$_->[0].'TeStatus'} = $self->{feat}->{$seqid}->{$_->[1]}->{status};

		if ($_->[0] eq 'left'){
			$self->{print}->{$_->[0].'TeStop'} = $self->{feat}->{$seqid}->{$_->[1]}->{stop};
			$self->{print}->{$_->[0].'DistanceToTe'} = $geneStart - $self->{feat}->{$seqid}->{$_->[1]}->{stop};
		}
		else{
			$self->{print}->{$_->[0].'TeStart'} = $self->{feat}->{$seqid}->{$_->[1]}->{start};
			$self->{print}->{$_->[0].'DistanceToTe'} = $self->{feat}->{$seqid}->{$_->[1]}->{start} - $geneStop;
		}
	}

	### PRINT ###
	$self->printData();

	return 1;
}


sub searchByIndex {
	my $self = shift;
	my ($seqid, $genePosition, $index, $leftOrRight) = @_;
	
	my $candidate = 'NA';
	my $i = $index;
	my $startOrStop;
	while($i>$index-1000 and $i<$index+1000){
		if($leftOrRight eq 'left'){
			$i--;
			$startOrStop = 'stop';
		}
		else{
			$i++;
			$startOrStop = 'start';
		}
		$self->{index}->{$startOrStop}->{$seqid}->{$i} or last;
		my $kFt = $self->{index}->{$startOrStop}->{$seqid}->{$i};
		my $pritag = $self->{feat}->{$seqid}->{$kFt}->{pritag};
		next unless ($pritag eq 'repeat_region' || $pritag eq 'nested_repeat');

		if($leftOrRight eq 'left'){
			next if($self->{feat}->{$seqid}->{$kFt}->{$startOrStop} > $genePosition);
		}
		else{
			next if($self->{feat}->{$seqid}->{$kFt}->{$startOrStop} < $genePosition);
		}
		$candidate = $kFt;
		last;
	}
	return $candidate;
}


sub getTeLeft {
	my $self = shift;
	my ($seqid, $genePosition, $index) = @_;
	my $candidates = 0;
	my $i = $index;
	while($i>$index-1000){
		$i--;
		$self->{index}->{stop}->{$seqid}->{$i} or last;
		my $kFt = $self->{index}->{stop}->{$seqid}->{$i};
		my $pritag = $self->{feat}->{$seqid}->{$kFt}->{pritag};
		next unless ($pritag eq 'repeat_region' || $pritag eq 'nested_repeat');
		next if($self->{feat}->{$seqid}->{$kFt}->{stop} > $genePosition);
		last if($self->{feat}->{$seqid}->{$kFt}->{stop} < $genePosition-$self->{param}->{range});
		
		$candidates++;
		$self->{print}->{teId} = $self->{feat}->{$seqid}->{$kFt}->{ID};
		$self->{print}->{teCompo} = $self->{feat}->{$seqid}->{$kFt}->{compo};
		$self->{print}->{teStatus} = $self->{feat}->{$seqid}->{$kFt}->{status};
		$self->{print}->{tePosition} = $self->{feat}->{$seqid}->{$kFt}->{stop};
		$self->{print}->{distanceToTe} = $genePosition-$self->{feat}->{$seqid}->{$kFt}->{stop};
		$self->printData();
	}
	$candidates==0 and $self->printData();
	return 1;
}


sub getTeRight {
	my $self = shift;
	my ($seqid, $genePosition, $index) = @_;
	my $candidates = 0;
	my $i = $index;
	while($i<$index+1000){
		$i++;
		$self->{index}->{start}->{$seqid}->{$i} or last;
		my $kFt = $self->{index}->{start}->{$seqid}->{$i};
		my $pritag = $self->{feat}->{$seqid}->{$kFt}->{pritag};
		next unless ($pritag eq 'repeat_region' || $pritag eq 'nested_repeat');
		next if($self->{feat}->{$seqid}->{$kFt}->{start} < $genePosition);
		last if($self->{feat}->{$seqid}->{$kFt}->{start} > $genePosition+$self->{param}->{range});
		
		$candidates++;
		$self->{print}->{teId} = $self->{feat}->{$seqid}->{$kFt}->{ID};
		$self->{print}->{teCompo} = $self->{feat}->{$seqid}->{$kFt}->{compo};
		$self->{print}->{teStatus} = $self->{feat}->{$seqid}->{$kFt}->{status};
		$self->{print}->{tePosition} = $self->{feat}->{$seqid}->{$kFt}->{start};
		$self->{print}->{distanceToTe} = $self->{feat}->{$seqid}->{$kFt}->{start}-$genePosition;
		$self->printData();
	}
	$candidates==0 and $self->printData();
	return 1;
}


sub printData {
	my $self = shift;
	my @out; 
	foreach (@{$self->{colNames}}) {
		push(@out, $self->{print}->{$_});
	}
	print join("\t", @out), "\n";
	return 1;
}


sub help() {
	my $prog = basename($0) ;
	print STDERR <<EOF ;
### $prog $VERSION ###
# LAST MODIF:  $lastModif
# PURPOSE:     This script is used to get the clostest 5' and 3' TEs in the vicinity of genes
USAGE: 
       $prog  <gffFile>

       ### OPTIONS ###
       -m, --mode <closest|env>         Analysis mode: 'closest' or 'env' [Default: closest]
                                        closest: find the closest TE for each gene (5' and 3')
                                        env:     find all TEs present in the 5' region of each gene (custom range with -r option)
       -r, --range <integer>            Range for search mode "env" [Default: 5000 bp]
                                        Must be between 200 and 20000 bps
       -gv, --gffversion <integer>      GFF version [Default: 3]
       -h, --help                       Print this help


EOF
	exit(1) ;
}
