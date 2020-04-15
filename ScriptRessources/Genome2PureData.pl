#!/usr/bin/perl

# Intend to transform DNA into a list of character readable by PureData - from sequence as fasta and annotations as bed format


### To do
# Store Informations (Bed, structure)
# Check informations (Annotated files)
# Verbose version for debug
# Translate positions for every features
# Create list


use warnings;
use strict;
use Getopt::Long;
use POSIX;


my $FASTAfile ;
my $BEDfile;
my $debug ;
my $verbose;
my $help ;
my $outputName ;
my $composition = 0 ;
my $SlidingWindowSize = 99 ;

GetOptions(
           'f=s'          => \$FASTAfile,
           'b=s'          => \$BEDfile,
           'h'        => \$help,
	   'o=s'	=>\$outputName,
	   's=i'	=> \$SlidingWindowSize,
	   'c!'	=>\$composition, 
	   'd!'	=>\$debug,
	   'v!'	=> \$verbose
          );

my $usage = "\nUsage: $0 -f inputFastaFile -b inputBedFile -o outputFileName\n";
$usage .= "Description : Convert Genome (Fasta, annotations, ... ) into a Pure data readable file (a list)\n";
$usage .= "Option :    \n";
$usage .= "            -c Base composition (\%GC) \n";
$usage .= "            -s Size of the Sliding Window [def.=$SlidingWindowSize]\n";
$usage .= "            -v Verbose mode\n";
$usage .= "            -d Ultra-verbose mode (debug mode)\n";
$usage .= "            -h Displays this help and exit\n";

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
unless(defined($outputName) ){
    print STDERR "Output file name is missing $!\n";
    print STDERR $usage;
    exit 1;
}

#########################################################################
### Script Structure :
#########################################################################
### i. Opens Fasta, stores bases in a hash (position are keys -> begin at 1)
### ii. Opens Bed File, stores Features in a hash (Features are keys)
### iii. Checks availability of annotated files for every features
### iv. Opens annotated files, stores informations into something...
### v. Loops through DNA, get information into Feature's Hash and ....
### vi print
# Store Informations (Bed, structure)
# Check informations (Annotated files)
# Verbose version for debug
# Translate positions for every features
# Create list



### i. Dealing with Fasta file
#==========================
#Open Fasta file 
open(IN, "<$FASTAfile" ) or die "Unable to open $FASTAfile !" ;
my $sequence ;
my $GenomeFastaLength ;
my %DNAgenome ;
my $NbrKeysDNAgenomeFasta ;
#Reading and Storing
my $ReadingPosition = 1;
while (my $line = <IN>){
	chomp $line ;
	#Next if header
	if($line =~ /^>/ ){next ;}
	
	# last version - Genome stored as scalar
	$sequence .= $line ;
	
	my @bases = split('', $line) ;
	foreach my $ba (@bases){
		# Storing base
		$DNAgenome{$ReadingPosition}{'Base'} = $ba ;
		# Incrementing counter
		$ReadingPosition+=1;
	}
}
$NbrKeysDNAgenomeFasta = keys %DNAgenome ;
$GenomeFastaLength =$ReadingPosition-1 ;
close IN;


### ii. Dealing with Features structure (Bed file)
#=================================================
#Open Bed file and store informations
open(IN, "<$BEDfile" ) or die "Unable to open $BEDfile !" ;
my %Bed ; # last version
my %Features ;
my $FeatureCounter = 0 ;
#Reading and Storing
while (my $line = <IN>){
	chomp $line ;
	my @words = split "\t", $line ;
	if($words[0] =~ /chr/i){next;}
	# last version
	$Bed{"$FeatureCounter"}{"start"} = $words[1] ;
	$Bed{"$FeatureCounter"}{"stop"} = $words[2] ;
	$Bed{"$FeatureCounter"}{"strand"} = $words[5] ;
	$Bed{"$FeatureCounter"}{"type"} = $words[3] ;
	
	# dev version
	my @NameFeatureDetails = split( ' ', $words[3]) ;
	my $NameFeature = $NameFeatureDetails[0] ;
	$Features{$NameFeature}{"start"} = $words[1] ;
	$Features{$NameFeature}{"stop"} = $words[2] ;
	$Features{$NameFeature}{"strand"} = $words[5] ;
	$Features{$NameFeature}{"file"} = $words[0].'_'.${NameFeature}.'.CFss1D.tab' ; # Error prone -> my own logic
	
	# Incrementing counter
	$FeatureCounter+=1;
}
close IN;

if ($verbose){
	print "Old way\n" ;
	foreach my $feature (keys %Bed ){
		print $feature , "\t", $Bed{$feature}{"type"}, "\t" , $Bed{$feature}{"start"} , "\t", $Bed{$feature}{"stop"}, "\n" ;
	}
	print "\nNew way\n" ;
	foreach my $feature (keys %Features ){
		print $feature, "\t", $Features{$feature}{"start"}, "\t", $Features{$feature}{"stop"}, "\t", $Features{$feature}{"strand"}, "\t", $Features{$feature}{"file"}, "\n" ;
	}
	
}

### iii. Checking availibility of annotated files
#=================================================
foreach my $feature (keys %Features ){
	if ( -f $Features{$feature}{"file"}){
		$Features{$feature}{"fileAvail"} = "Found" ;
	} else {
		$Features{$feature}{"fileAvail"} = "Not Found" ;
		if($verbose){
			print $Features{$feature}{"file"}, " not found !\n" ; 
		}
	}
}

### iv. Opens annotated files, Stores inforamtions
#=================================================
foreach my $feature (keys %Features ){
	unless ( -f $Features{$feature}{"file"}){
		$Features{$feature}{"fileAvail"} = "Not Found" ;
		if($verbose or $debug){
			print $Features{$feature}{"file"}, " not found !\n" ; 
		}
	} else {
		open(IN, "<$Features{$feature}{'file'}" ) or die "Unable to open $Features{$feature}{'file'} !\n" ;
		my $PositionInFeature = 0 ;
		while(my $line = <IN>){
			chomp $line ;
			my @words = split "[ \*]+", $line ;
			if($words[1] !~ /^[0-9]/){next;}
			if ($debug) { 
				print $words[1], " AND ", $words[2], " AND ", $words[3], " AND ", $words[12],"\n" ;
			}
			# Find position on genome
			my $DNApos = $Features{$feature}{'start'} + $PositionInFeature * 3 ;
			my $DNAposSc = $DNApos + 1 ;
			my $DNAposTh = $DNApos + 2 ;
			# Get AA
			$DNAgenome{$DNApos}{'AA'} = $words[2] ;
			$DNAgenome{$DNAposSc}{'AA'} = $words[2] ;
			$DNAgenome{$DNAposTh}{'AA'} = $words[2] ;
			# Get Secondary Structure 
			$DNAgenome{$DNApos}{'SecStru'} = $words[12] ;
			$DNAgenome{$DNAposSc}{'SecStru'} = $words[12] ;
			$DNAgenome{$DNAposTh}{'SecStru'} = $words[12] ;
			# Get Feature
			$DNAgenome{$DNApos}{'Feature'} = $feature ;
			$DNAgenome{$DNAposSc}{'Feature'} = $feature ;
			$DNAgenome{$DNAposTh}{'Feature'} = $feature ;
			
			
			if ($debug) { 
				print $words[1], " AND ", $words[2], " AND ", $words[3], " AND ", $words[12], " : SecStruct of postion $DNAposSc  following $DNApos = ", $DNAgenome{$DNAposSc}{'SecStru'},"\n" ; 
			}
			
			# Incrementing counter
			$PositionInFeature+= 1 ;
		}
		close IN
	}
}

### v. Compute GC content in a sliding windows
#==============================================
foreach my $pos (keys %DNAgenome){
	my %GCcount ;
	$GCcount{'G'} = 0 ;
	$GCcount{'C'} = 0 ;
	my $GCcontent ; 
	my $DNAslideStart ;
	my $DNAslideStop ;
	my $TheoriticalStart =  $pos - int( $SlidingWindowSize / 2 ) ;
	my $TheoriticalStop = $pos + int( $SlidingWindowSize / 2 ) ;
	
	#Finds limits :
	if($TheoriticalStart <= 1){
		my $diff = $pos - 1 ;
		$DNAslideStart = 1 ;
		$DNAslideStop = $pos + $diff;
	}elsif ($TheoriticalStop >= $GenomeFastaLength ){
		my $diff = $GenomeFastaLength - $pos ;
		$DNAslideStart = $pos - $diff ;
		$DNAslideStop = $GenomeFastaLength ;
	}else{
		$DNAslideStart = $TheoriticalStart ;
		$DNAslideStop = $TheoriticalStop ;
	}
	
	#Computes GCcontent :
	my $slideSeq ;
	for(my $i=$DNAslideStart; $i <= $DNAslideStop; $i++) {
		my $baseChar = $DNAgenome{$i}{'Base'} ;
		$slideSeq .= $baseChar ;
		$GCcount{$baseChar}++;
	}
	if ($debug) {print "Position : $pos ; SlideSeq =  $slideSeq \n"; }
	$GCcontent = ( $GCcount{'G'} + $GCcount{'C'} ) / ( $DNAslideStop - $DNAslideStart + 1 ) ;
	
	# Stores results :
	$DNAgenome{$pos}{'GC'} = $GCcontent ;
}
my $NbrKeysDNAgenomeAnnotation = keys %DNAgenome ;


### vi. Formats and prints
#==============================================
open(OUT, ">$outputName" ) or die "Unable to open $outputName !" ;
for (my $pos=1 ; $pos <= $NbrKeysDNAgenomeFasta; $pos++){
#foreach my $pos (keys %DNAgenome){
	#position, base, GC content, Feature, AA, Sec Struct.
	my $line = $pos."\t" ;
	if ( $debug){
		print join(" ", keys($DNAgenome{$pos}))."\t".join(" ", values($DNAgenome{$pos}))."\n" ;
	}
	$line .= $DNAgenome{$pos}{'Base'}."\t" ;
	unless( defined $DNAgenome{$pos}{'GC'}){ print "$pos \n"; next ;}
	$line .= sprintf("%.2f", $DNAgenome{$pos}{'GC'})."\t" ;
	if ($DNAgenome{$pos}{'Feature'}){
		$line .= $DNAgenome{$pos}{'Feature'}."\t" ;
		$line .= $DNAgenome{$pos}{'AA'}."\t" ;
		$line .= $DNAgenome{$pos}{'SecStru'}."\n" ;
	}
	else{
		$line .= "\t" ;
		$line .= "\t" ;
		$line .= "\n" ;
	}

	if ($verbose){ print $line ; }
	print OUT $line ;
}
close OUT ;
if ($debug) {print "Nbr of keys: $NbrKeysDNAgenomeFasta $NbrKeysDNAgenomeAnnotation \n";}

### Here I am
exit 1;
### Here I am

