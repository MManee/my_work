#!/usr/bin/perl
use strict;
use warnings;
no warnings 'uninitialized';
use Archive::Extract;
use LWP::Simple;
use Statistics::R;


##############################
#### PhastCons processing ####
##############################

# download UCSC PC dataset
my $PC_url = 'http://hgdownload.soe.ucsc.edu/goldenPath/dm3/database/phastConsElements15way.txt.gz';
my $PC_file = 'phastConsElements15way.txt.gz';
getstore($PC_url, $PC_file);
my $ae = Archive::Extract->new( archive => 'phastConsElements15way.txt.gz' );
my $ok = $ae->extract;

# remove PC gzip file
`rm -rf phastConsElements15way.txt.gz`;

# perform overlapSelect to identify non-exonic blocks for PC
system("overlapSelect -selectCoordCols=0,1,2 -inCoordCols=1,2,3 dm3_exons.txt phastConsElements15way.txt PC_filter.txt");

my $datafile = "phastConsElements15way.txt";
my $filterfile = "PC_filter.txt";

my %compdata;
my $line;

open FILEDATA, "< " . $datafile or die "Cannot open the phastcons file: $!";
open FILEFILTER, "< " . $filterfile or die "Cannot open the filter file: $!";

# read both the files into a hash for quick checking
my @data = <FILEDATA>;
while ($line = <FILEFILTER>)
{
	$line =~ s/\s+//g;
	$compdata{$line}=1;
}
close FILEDATA;
close FILEFILTER;
# check if given line in pashtcons file exists in filter file 
foreach (@data)
{
	# if the line exists write EXO as a new column
	chomp($_);
	$line = $_;
	$line =~ s/\s+//g;
	if (exists $compdata{$line})
	{
		$_ .= "\tEXO\n";
	}
        # if the line does not exist write INT as a new column
	else
	{
		$_ .= "\tINT\n";
	}
}
#overwrite the phastcons file with the new data
open FILEDATA, "> ". $datafile;
foreach (@data)
{
    print FILEDATA $_;
}
close FILEDATA;

open FILEOUTPUT1, "> PC_nonexonic_exonic_spacer.txt";
open FILEOUTPUT2, "> PC_nonexonic_spacer.txt";
open FILEOUTPUT3, "> PC_exonic_spacer.txt";
open FILEOUTPUT4, "> PC_nonexonic_block.txt";
open FILEOUTPUT5, "> PC_exonic_block.txt";

my @currLine;
my @nextLine;
my $exonic_blocks_count=0;
my $nonexonic_blocks_count=0;
my $exonic_spacer_count=0;
my $nonexonic_spacer_count=0;
my $nonexonic_exonic_spacer_count=0;
foreach (0 .. $#data)
{
	
	@currLine = split /\s+/, $data[$_];
	@nextLine = split /\s+/, $data[$_+1];

	if (((($data[$_] =~ /INT/) && ($data[$_+1] =~ /EXO/)) || (($data[$_] =~ /EXO/) && ($data[$_+1] =~ /INT/))) && ($nextLine[1] eq $currLine[1]) && (($currLine[1] eq 'chr2L') || ($currLine[1] eq 'chr2R') || ($currLine[1] eq 'chr3L') || ($currLine[1] eq 'chr3R') || ($currLine[1] eq 'chrX')))
	{
    		# exonic/nonexonic spacers
    		print FILEOUTPUT1 $currLine[1], "\t", $currLine[3], "\t", $nextLine[2], "\n";
    		$nonexonic_exonic_spacer_count++;
	}

	if (($data[$_] =~ /INT/) && ($data[$_+1] =~ /INT/) && ($nextLine[1] eq $currLine[1]) && (($currLine[1] eq 'chr2L') || ($currLine[1] eq 'chr2R') || ($currLine[1] eq 'chr3L') || ($currLine[1] eq 'chr3R') || ($currLine[1] eq 'chrX')))
    	{ 
        	# nonexonic spacers
        	print FILEOUTPUT2 $currLine[1], "\t", $currLine[3], "\t", $nextLine[2], "\n";
		$nonexonic_spacer_count++;
    	}

    	if (($data[$_] =~ /EXO/) && ($data[$_+1] =~ /EXO/) && ($nextLine[1] eq $currLine[1]) && (($currLine[1] eq 'chr2L') || ($currLine[1] eq 'chr2R') || ($currLine[1] eq 'chr3L') || ($currLine[1] eq 'chr3R') || ($currLine[1] eq 'chrX')))
    	{
        	# exonic spacers
        	print FILEOUTPUT3 $currLine[1], "\t", $currLine[3], "\t", $nextLine[2], "\n";
		$exonic_spacer_count++;
    	}
		   
    	if (($data[$_] =~ /INT/) && (($currLine[1] eq 'chr2L') || ($currLine[1] eq 'chr2R') || ($currLine[1] eq 'chr3L') || ($currLine[1] eq 'chr3R') || ($currLine[1] eq 'chrX')))
    	{
		# nonexonic blocks
		print FILEOUTPUT4 $currLine[0], "\t", $currLine[1], "\t", $currLine[2], "\t", $currLine[3], "\t", $currLine[4], "\t", $currLine[5], "\n";
		$nonexonic_blocks_count++;
    	}

    	if (($data[$_] =~ /EXO/) && (($currLine[1] eq 'chr2L') || ($currLine[1] eq 'chr2R') || ($currLine[1] eq 'chr3L') || ($currLine[1] eq 'chr3R') || ($currLine[1] eq 'chrX')))
    	{
		# exonic blocks
		print FILEOUTPUT5 $currLine[0], "\t", $currLine[1], "\t", $currLine[2], "\t", $currLine[3], "\t", $currLine[4], "\t", $currLine[5], "\n";
		$exonic_blocks_count++;
    	}

}
close FILEOUTPUT1;
close FILEOUTPUT2;
close FILEOUTPUT3;
close FILEOUTPUT4;
close FILEOUTPUT5;

print "\nNumber of PhastCons exonic blocks: ". $exonic_blocks_count. "\n";
print "\nNumber of PhastCons nonexonic blocks: ". $nonexonic_blocks_count. "\n";  
print "\nNumber of PhastCons exonic spacers: ". $exonic_spacer_count. "\n";
print "\nNumber of PhastCons nonexonic spacers: ". $nonexonic_spacer_count. "\n";
print "\nNumber of PhastCons exonic/nonexonic spacers: ". $nonexonic_exonic_spacer_count. "\n\n\n";


##############################
###### Tanay processing ######
##############################

# download KT dataset
my $KT_url = 'http://www.wisdom.weizmann.ac.il/~effi/plosgen2013/CE.txt';
my $KT_file ='KT_CEs.txt';
getstore($KT_url, $KT_file);

# replace spaces with tabs
system("perl -p -i -e 's/ /\t/g' $KT_file");

# remove header line form KT file
open KTFILE, "< " . $KT_file or die "Cannot open the KT file: $!";
my @lines = <KTFILE>;
foreach (@lines)
{
	shift @lines if $lines[0] =~ /^chrom	start	end/;
}
# overwrite the KT file without header line
open KTFILE, "> ". $KT_file;
foreach (@lines)
{
    print KTFILE $_;
}
close KTFILE;


# perform overlapSelect to identify non-exonic blocks for KT
system("overlapSelect -selectCoordCols=0,1,2 -inCoordCols=0,1,2 dm3_exons.txt KT_CEs.txt KT_filter.txt");

my $tdatafile = "KT_CEs.txt";
my $tfilterfile = "KT_filter.txt";

my %tcompdata;
my $tline;

open FILEDATA, "< " . $tdatafile or die "Cannot open the phastcons file: $!";
open FILEFILTER, "< " . $tfilterfile or die "Cannot open the filter file: $!";

# read both the files into a hash for quick checking
my @tdata = <FILEDATA>;
while ($tline = <FILEFILTER>)
{
	$tline =~ s/\s+//g;
	$tcompdata{$tline}=1;
}
close FILEDATA;
close FILEFILTER;
# check if given line in KT file exists in filter file 
foreach (@tdata)
{
	# if the line exists write EXO as a new column
	chomp($_);
	$tline = $_;
	$tline =~ s/\s+//g;
	if (exists $tcompdata{$tline})
	{
		$_ .= "\tEXO\n";
	}
        # if the line does not exist write INT as a new column
	else
	{
		$_ .= "\tINT\n";
	}
}
#overwrite the KT file with the new data
open FILEDATA, "> ". $tdatafile;
foreach (@tdata)
{
    	print FILEDATA $_;
}
close FILEDATA;

open FILEOUTPUT1, "> KT_nonexonic_exonic_spacer.txt";
open FILEOUTPUT2, "> KT_nonexonic_spacer.txt";
open FILEOUTPUT3, "> KT_exonic_spacer.txt";
open FILEOUTPUT4, "> KT_nonexonic_block.txt";
open FILEOUTPUT5, "> KT_exonic_block.txt";

my @tcurrLine;
my @tnextLine;
my $texonic_blocks_count=0;
my $tnonexonic_blocks_count=0;
my $texonic_spacer_count=0;
my $tnonexonic_spacer_count=0;
my $tnonexonic_exonic_spacer_count=0;
foreach (0 .. $#tdata)
{

    	@tcurrLine = split /\s+/, $tdata[$_];
    	@tnextLine = split /\s+/, $tdata[$_+1];

    	if (((($tdata[$_] =~ /INT/) && ($tdata[$_+1] =~ /EXO/)) || (($tdata[$_] =~ /EXO/) && ($tdata[$_+1] =~ /INT/))) && ($tnextLine[0] eq $tcurrLine[0]))
    	{
    		# exonic/nonexonic spacers	
    		print FILEOUTPUT1 $tcurrLine[0], "\t", $tcurrLine[2], "\t", $tnextLine[1], "\n";
    		$tnonexonic_exonic_spacer_count++;
    	}

    	if (($tdata[$_] =~ /INT/) && ($tdata[$_+1] =~ /INT/) && ($tnextLine[0] eq $tcurrLine[0]))
    	{ 
        	# nonexonic spacers
        	print FILEOUTPUT2 $tcurrLine[0], "\t", $tcurrLine[2], "\t", $tnextLine[1], "\n";
		$tnonexonic_spacer_count++;
    	}

    	if (($tdata[$_] =~ /EXO/) && ($tdata[$_+1] =~ /EXO/) && ($tnextLine[0] eq $tcurrLine[0]))
    	{
        	# exonic spacers
        	print FILEOUTPUT3 $tcurrLine[0], "\t", $tcurrLine[2], "\t", $tnextLine[1], "\n";
		$texonic_spacer_count++;
    	}
		   
    	if ($tdata[$_] =~ /INT/)
    	{
		# nonexonic blocks
		print FILEOUTPUT4 $tcurrLine[0], "\t", $tcurrLine[1], "\t", $tcurrLine[2], "\n";
		$tnonexonic_blocks_count++;
    	}

    	if ($tdata[$_] =~ /EXO/)
    	{
		# exonic blocks
		print FILEOUTPUT5 $tcurrLine[0], "\t", $tcurrLine[1], "\t", $tcurrLine[2], "\n";
		$texonic_blocks_count++;
    	}

}

close FILEOUTPUT1;
close FILEOUTPUT2;
close FILEOUTPUT3;
close FILEOUTPUT4;
close FILEOUTPUT5;

print "\nNumber of Tanay exonic blocks: ". $texonic_blocks_count. "\n";
print "\nNumber of Tanay nonexonic blocks: ". $tnonexonic_blocks_count. "\n";  
print "\nNumber of Tanay exonic spacers: ". $texonic_spacer_count. "\n";
print "\nNumber of Tanay nonexonic spacers: ". $tnonexonic_spacer_count. "\n";
print "\nNumber of Tanay exonic/nonexonic spacers: ". $tnonexonic_exonic_spacer_count. "\n\n\n";



