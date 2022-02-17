#!/usr/bin/perl

# appendTreesblock.pl v0.1.0
# Author: MH Seabolt
# Last updated: 8-13-2021

# SYNOPSIS:
# Appends a basic TREES block to the end of a NEXUS alignment file so that the user can edit the parameters
# NOTE:
# This script will APPEND A NEW TREES BLOCK to the end of the Nexus file given,
# PAUP will only consider the last trees block read, so it is fine to have multiple TREES blocks in the same NEXUS file
# This also applies to tree naming --> we dont care if we repeat names since PAUP wont consider the pre-existing ones

##################################################################################
# The MIT License
#
# Copyright (c) 2021 Matthew H. Seabolt
#
# Permission is hereby granted, free of charge, 
# to any person obtaining a copy of this software and 
# associated documentation files (the "Software"), to 
# deal in the Software without restriction, including 
# without limitation the rights to use, copy, modify, 
# merge, publish, distribute, sublicense, and/or sell 
# copies of the Software, and to permit persons to whom 
# the Software is furnished to do so, 
# subject to the following conditions:
#
# The above copyright notice and this permission notice 
# shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR 
# ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
##################################################################################

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

# Required input parameters
my $nexus = "--";
my $output = "--";
my $newick;
my $name;
my $unroot;

sub usage {
	my $usage = "appendTreesblock.pl\n
	PURPOSE: 	Appends a basic TREES block to the end of a NEXUS file.\n
	USAGE:	appendTreesblock.pl -i input.nex -t tree.nwk -n tree_name
	-i		input NEXUS file
	-o 		output file name
	-t		newick file containing the tree you want to add
	-u 		INT flag Denote trees as unrooted?
	-n 		name of the tree (no spaces!)\n";
	
	print $usage;
}

GetOptions(	'input|i=s' => \$nexus,
			'output|o=s' => \$output,
			'tree|t=s' => \$newick,
			'unroot|u=i' => \$unroot,
			'name|n=s' => \$name,
			
) or die usage();

# Defaults
$name = ( defined($name) )? $name : "tree1";
if    ( defined(($unroot) && $unroot == 1 ) )	{	$unroot = 1;	}
elsif ( defined(($unroot) && $unroot == 2 ) )	{	$unroot = 2;	}
elsif ( defined(($unroot) && $unroot == 0 ) )	{	$unroot = 0;	}
else											{	$unroot = 0;	}

# Get the new tree data to append
my @Newicks = split ",", $newick;

# Foreach new tree FILE in the list of tree files:
my %Trees = ();
my $k = 1;
my $m = $k - 1;
my $x;		# This is a counter variable that will save off the number of trees that we wrote
foreach my $nwk_tree ( @Newicks )	{
	chomp $nwk_tree;
	
	# Get the tree data
	open NWK, $nwk_tree or die "$!\n";
		my @newick_file = <NWK>;
		my $tree = $newick_file[0];	# The tree will only be on the first line, should also incl the ';' at the end 
	close NWK;
	
	$Trees{$k} = [ $tree, $Newicks[$m] ];
	$x = $k;
	$k++;
	$m++;
}

# Open and read the NEXUS file, look for an existing TREES block
my $fh = *STDIN;
my $succin = open(NEX, "<", "$nexus") if ( $nexus ne "--" && -e $nexus );
$fh = *NEX if ( $succin ); 
	my @data_in = <$fh>;
close $fh if ( $succin );

# Set the output filehandles
my $succout = open( OUT, ">", "$output" ) if $output ne "--";
my $fhout;
if ( $succout )		{	$fhout = *OUT;			}
else				{	$fhout = *STDOUT;		}

# Write the tree block
print $fhout "$_" foreach (@data_in);
print $fhout "\n\nBEGIN TREES;\n";
foreach my $tree ( sort keys %Trees )	{
	if ( $unroot == 1 )	{
		print $fhout "\ttree $tree = [&U] $Trees{$tree}[0]";
	}
	elsif ( $unroot == 2 )	{
		print $fhout "\ttree $tree = [&R] $Trees{$tree}[0]";
	}
	else	{
		print $fhout "\ttree $tree = $Trees{$tree}[0]";
	}
}
foreach my $tree ( sort keys %Trees )	{
	print $fhout "\t[tree $tree file = $Trees{$tree}[1]]\n";
}
print $fhout "\t[ntrees=$x]\n";		# Print out the number of trees as a comment line, so that we can easily extract this number later
print $fhout "END;\n\n";	

close $fhout if ( $succout );
exit;
