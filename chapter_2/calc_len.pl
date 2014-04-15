#! /usr/bin/perl
use strict;
use warnings;

open BASE_COMP, "> " . "pha_blc_len.txt" or die "Cannot open the block file: $!";
print BASE_COMP "CE_A\tCE_T\tCE_AT\tCE_G\tCE_C\tCE_CG\tTotal\n";
		
		
		my $file = "BLC_SEQ_PHA.txt";
		open(FILE, $file) || die("Couldn't read file $file\n");
		my $a=0;
		my $t=0;
		my $c=0;
		my $g=0;
		my @T;
		my $T;
		my $k;
		my $cg;
		my $at;
		my $total;
		while (my $seq = <FILE>)
		{
			chomp $seq;
			$seq =~ s/ //g; 
			$seq =~ s/\n//g;
			push(@T,$seq);
		}
			
		for($k=0;$k<@T;$k++)
		{	
			$a = $T[$k] =~ tr/a//;
			$t = $T[$k] =~ tr/t//;
			$c = $T[$k] =~ tr/c//;
			$g = $T[$k] =~ tr/g//;
			$at=$a+$t;$cg=$c+$g;
			$total=$a+$t+$c+$g;			
			print BASE_COMP $a, "\t", $t, "\t", $at, "\t", $g, "\t", $c, "\t", $cg, "\t", $total, "\n";
		}
		close (FILE);		
		
close (BASE_COMP);
