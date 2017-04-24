#!/usr/bin/perl

use strict;
use warnings;

#
# Name:    merge_rpkms.pl
#
# Purpose: Merges multiple RKPM files into a single matrix, removing
#          lines where more than one of the values is a zero.
#
# INPUT:   RPKM files output from RSubread analysis
# OUTPUT:  Matrix of RPKM values written to stdout
#
# Usage:   ./merge_rpkms.pl [file1.rpkmlog] [file2.rpkmlog] ... > merged.rpkmlog
#
# Author:   David Jenkins
# Date:     Created 20141216
#
############################################################################

#read in options
die("You must submit multiple rpkm files\n") if scalar(@ARGV) < 2;

my @rpkm;
$rpkm[0][0] = 'Gene'; 
my $curr_col = 1;
foreach my $infile(@ARGV){
    open(my $fi, $infile) or die("ERROR: $infile is not a valid infile!\n");
    my $s_name = $infile;
    $s_name =~ s/\.[^\.]+$//g;
    $rpkm[0][$curr_col] = $s_name;
    my $curr_row = 1;
    while(my $line = <$fi>){
        chomp($line);
        my @l_a = split(/\t/, $line);
        if(defined $rpkm[$curr_row][0]){
            #die if it is not the same as the current file
            die("$l_a[0] does not match $rpkm[$curr_row][0]\n") if $l_a[0] ne $rpkm[$curr_row][0];
        }
        else{
            $rpkm[$curr_row][0] = $l_a[0];
        }
        $rpkm[$curr_row][$curr_col] = $l_a[1];
        $curr_row++;
    }
    $curr_col++;
}

for(my $i = 0; $i <= $#rpkm; $i++){
    my $outline;
    my $zero_count = 0;
    $zero_count = 100 if $rpkm[$i][0] eq '';
    for(my $j = 0; $j <= $#{$rpkm[0]} ; $j++){
#        $zero_count++ if ($j != 0 && $i != 0 && $rpkm[$i][$j] == 0);
        $outline .= "$rpkm[$i][$j]\t";
    }
    chop($outline);
    print "$outline\n" if $zero_count < 2;
}
