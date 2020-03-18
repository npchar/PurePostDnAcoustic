#!/usr/bin/perl

# Intend to transform DNA into a list of character readable by PureData - from sequence as fasta and annotation store as bed format

use warnings;
use strict;
use Getopt::Long;
use POSIX;


my $FASTAfile ;
my $BEDfile;
my $help ;
my $outputName ;
my $composition = 0 ;
my $SlidingWindowSize = 100 ;

GetOptions(
           'f=s'          => \$FASTAfile,
           'b=s'          => \$BEDfile,
           'h'        => \$help,
	   'o=s'	=>\$outputName,
	   'c'	=>\$composition, 
          );

my $usage = "\nUsage: $0 -f inputFastaFile -b inputBedFile-o outputFileName\n";
$usage .= "Description : Convert Genome (Fasta, annotations, ... ) into a Pure data readable file (a list)\n";
$usage .= "Option :    \n";
$usage .= "            -c Base composition (\%GC) \n";
$usage .= "            -s Size of the Sliding Window [def.=$SlidingWindowSize]\n";
$usage .= "            -h   Displays this help and exit\n";

if($help){
    print $usage;
    exit 0;
}
unless(defined($FASTAfile) && -f $FASTAfile ){
    print STDERR "File(s) [$FASTAfile] not accessible : $!\n";
    print STDERR $usage;
    exit 1;
}
unless(defined($BEDfile) && -f $BEDfile ){
    print STDERR "File(s) [$BEDfile] not accessible : $!\n";
    print STDERR $usage;
    exit 1;
}


#Open Fasta file and output file
open(IN, "<$FASTAfile" ) or die "Unable to open $FASTAfile !" ;
open(OUT, ">$outputName" ) or die "Unable to open $outputName !" ;
my $sequence ;
#Reading and Storing
while (my $line = <IN>){
	chomp $line ;
	if($line =~ /^>/ ){next ;} 
	$sequence .= $line ;
}
close IN;

#Open Bed file and store informations
open(IN, "<$BEDfile" ) or die "Unable to open $FASTAfile !" ;
my %Bed ; 
my $FeatureCounter = 0 ;
#Reading and Storing
while (my $line = <IN>){
	chomp $line ;
	my @words = split "\t", $line ;
	$Bed{"$FeatureCounter"}{"start"} = $words[1] ;
	$Bed{"$FeatureCounter"}{"stop"} = $words[2] ;
	$Bed{"$FeatureCounter"}{"strand"} = $words[5] ;
	$Bed{"$FeatureCounter"}{"type"} = $words[3] ;
	
$FeatureCounter+=1;
}
close IN;

#Wrapping Sequence
my $lengthSeq = length($sequence) ;
my @words = split '', $sequence ;
my $counter = 0 ;
foreach my $base ( @words ){
	my $GCcontent ;
	print $counter ;
	print "\n" ;
	
	### Dealing with composition (if -c Flag)
	if($composition){
		my $SlidingWindow ;
		
		if($counter < ($SlidingWindowSize/2) ){
			my $limitSupsup = $counter + ($SlidingWindowSize/2) ;
			my $limitInfinf = $lengthSeq - ($SlidingWindowSize/2) + $counter ;
			my $limitInfsup = ($lengthSeq-1) ;
			$SlidingWindow = join('', @words[0..$limitSupsup]) ;
			$SlidingWindow .= join('', @words[$limitInfinf..$limitInfsup]) ;
#			print $counter." ".$limitInfinf." ".$limitInfsup." ".$limitSupsup." " ;
#			print "Length:".length($SlidingWindow)." ".join('', $words[0..$limitSupsup])." ".$SlidingWindow."\n";
		
		}elsif($counter>= ( $lengthSeq - ($SlidingWindowSize/2) ) ){
			my $limitSupsup = ($lengthSeq-1) ; 
			my $limitSupinf = $counter - ($SlidingWindowSize/2) ;
			my $limitInfinf = 0 ;
			my $limitInfsup = 0 + ($SlidingWindowSize/2) - ($lengthSeq-$counter) ;
			$SlidingWindow = join('', @words[0..$limitInfsup]) ;
			$SlidingWindow .= join('', @words[$limitSupinf..$limitSupsup]) ;
#			print $counter." ".$limitInfinf." ".$limitInfsup." ".$limitSupinf." ".$limitSupsup." " ;
#			print "Length:".length($SlidingWindow)." ".join('', $words[0..$limitInfsup])." ".$SlidingWindow."\n";
#			print $counter." End part...\n";
		}
		
		else{
			my $limitSup= $counter + ( $SlidingWindowSize / 2 ) ;
			my $limitInf= $counter - ( $SlidingWindowSize / 2 ) ;
			$SlidingWindow = join('', @words[$limitInf .. $limitSup]) ;
#			print $counter." ".$limitInf." ".$limitSup." " ;
#			print "Length:".length($SlidingWindow)." ".$SlidingWindow."\n";
		}
	# compute GC prop:
		my @Bases = split '', $SlidingWindow ;
		my %count;
		foreach my $str (@Bases) {
			$count{$str}++;
		}
		$GCcontent = ( $count{'G'} + $count{'C'} ) / ( $SlidingWindowSize + 1 ) ;
		
	}
	
	
	print OUT $base ;
	if($composition){
		print OUT " ".sprintf("%.2f", $GCcontent) ;
	}
	
	# Dealing with Bed Storing structure
	my $infosFeature ='' ;
#	print "$counter " ;
	foreach my $feature (keys %Bed ){
#		print $Bed{$feature}{"strand"} ;
#		print " " ;
#		print $Bed{$feature}{"start"} ;
#		print "\n" ;
		if( $counter ge $Bed{$feature}{"start"} && $counter le $Bed{$feature}{"stop"} && $Bed{$feature}{"strand"} eq '+' ){
			$infosFeature = $Bed{$feature}{"type"} ;
		}
	}
	print OUT " ".$infosFeature ;
		
	# Preparing next round
	print OUT "\n" ;
	$counter+=1;
}



close(IN) ;
close(OUT) ;
