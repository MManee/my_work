#!/usr/bin/perl

open FFILE, ">>results.txt" or die $!;
print FFILE "Run\tA\tB\tRoot\tB-A\tA-Root\tB-Root\n";

$file = 'outputname_TRUE.fas';
$i=1;
open(FILE, $file) || die("Couldn't read file $file\n");

local $/ = "\n>";        # read by FASTA record

while (my $seq = <FILE>)
{

    chomp $seq;
    $seq =~ s/^>*.+\n//;        # remove FASTA headers
    $seq =~ s/\n//g;            # remove endlines
    $seq =~ s/-//g;             # remove "-" characters
    $seq =~ s/ //g;             # remove spaces
    $seq =~ s/[^a-zA-Z0-9]//g;  # remove special characters
    
    $R = length $seq;
    push(@T,$R);
    #unshift   
 }
 
 for($k=1;$k<=@T;$k=$k+3) 
 {
 $BA = $T[$k+0] - $T[$k-1]; # B-A
 $AR = $T[$k-1] - $T[$k+1]; # A-Root
 $BR = $T[$k+0] - $T[$k+1]; # B-Root
 
 print FFILE "$i\t$T[$k-1]\t$T[$k+0]\t$T[$k+1]\t$BA\t$AR\t$BR\n";
 
 $i++; 
}

exit;
