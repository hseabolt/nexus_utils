#!/usr/bin/perl

# appendPaupblock.pl v0.1.0
# Author: MH Seabolt
# Last updated: 8-13-2021

# SYNOPSIS:
# Appends a Paup block to the end of a NEXUS alignment file

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
no warnings 'uninitialized'; 
use Getopt::Long qw(GetOptions);

# Required input parameters
my $nexus = "--";
my $output = "--";
my $log;
my $method;
my $statistics;
my $nreps;
my $outgroup;
my $tree_scores;
my $quit;

sub usage {
	my $usage = "appendPaupblock.pl\n
	PURPOSE: 	Appends basic PAUP blocks to the end of a NEXUS file.\n
	
	USAGE:	append_parsimony_block.pl -i input.nex -l logfile.log [[ options ]]
	-i		input NEXUS file
	-b 		output file name
	-l		name of log file for PAUP
	-m 		PAUP method to use (parsimony, likelihood)
	-d 		INT flag; describe tree scores? ( only for simple parsimoney analysis, Default: ON )
	-s		compute statistics? (simple, bootstrap)
	-o 		outgroup taxon name
	-r 		n reps for tree searching (Default 1000)
	-q 		INT flag; print command to auto-close PAUP? ( Default: ON )
	\n";
	
	print $usage;
}

GetOptions(	'input|i=s' => \$nexus,
			'output|b=s' => \$output,
			'method|m=s' => \$method,
			'stats|s=s' => \$statistics,
			'nreps|r=i' => \$nreps,
			'describe_trees|d=i' => \$tree_scores,
			'log|l=s' => \$log,
			'outgroup|o=s' => \$outgroup,
			'quit|q=i' => \$quit,
			
) or die usage();

# Defaults
$outgroup = ( ($outgroup ne "NA") || ($outgroup ne "") )? $outgroup : undef; 
$nreps = ( $nreps > 0 )? $nreps : 1000;
$tree_scores = ( $tree_scores && $tree_scores == 0 )? 0 : 1; 
$quit = ( $quit && $quit == 0 )? 0 : 1; 

#######################################################################################
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

#### Print the basic PAUP block
print "$_" foreach (@data_in);

# -- log start
print $fhout "BEGIN PAUP;\n";
print $fhout "\tlog start replace file = $log\;\n";

# -- set basic parameters
print $fhout "[Tree searching parameters]\n";

# PARSIMONY ANALYSIS
if ( $statistics eq "simple" && $method eq "parsimony" )	{
	if ( $tree_scores == 1 )	{	print $fhout "\tset autoclose = yes criterion = parsimony root=outgroup maxtrees=100 increase=no storebrlens=yes;\n";		}
	else						{	print $fhout "\tset autoclose = yes criterion = parsimony root=outgroup increase=auto storebrlens=yes;\n";		}
	if ( $outgroup )	{
		print $fhout "\toutgroup $outgroup\;\n";
	}
	print $fhout "\thsearch addseq = random nreps = $nreps swap = tbr hold = 1;\n";
	print $fhout "\tsavetrees file = paup_simpleMP_trees.tre format = altnex brlens = yes;\n";
	print $fhout "\tpscores /TL HI CI RI;\n" if ( $tree_scores == 1 );
	print $fhout "\tcontree all \/ majrule = yes strict = no treefile = paup_simpleMP_consensus.tre;\n";
	
}
# PARSIMONY BOOTSTRAP ANALYSIS
elsif ( $statistics eq "bootstrap" && $method eq "parsimony" )	{
	print $fhout "\tset autoclose = yes criterion = parsimony root=outgroup increase=auto storebrlens=yes;\n";
	print $fhout "\tbootstrap nreps = $nreps search = heuristic\/ addseq = random nreps = 10 swap = tbr hold = 1;\n";
	print $fhout "\tsavetrees from = 1 to = 1 file = paup_MPboot_trees.tre format = altnex brlens = yes savebootp = NodeLabels MaxDecimals = 0;\n";
}
# MAXIMUM LIKELIHOOD
# This is going to use a GTR model
elsif ( $statistics eq "simple" && $method eq "likelihood" )	{
	print $fhout "\tset autoclose = yes criterion = distance root=outgroup increase=auto storebrlens=yes;\n";
	if ( $outgroup )	{
		print $fhout "\toutgroup $outgroup\;\n";
	}
	print $fhout "\tDSet distance = JC objective = ME base = equal rates = equal pinv = 0 subst = all negbrlen = setzero;\n";
    print $fhout "\tNJ showtree = no breakties = random;\n";
	print $fhout "\t    set criterion = like;\n";
    print $fhout "\tLset Base = (0.2892 0.2928 0.1309) Nst = 6 Rmat = (3.7285 46.5293 1.3888 2.3793 16.4374)\n";		# Need to edit the numbers in Base and Rmat
    print $fhout "\tRates = gamma Shape = 0.9350 Pinvar = 0.5691;\n";													# Need to edit rates, shape, and pinvar (run jmodeltest!)
	print $fhout "\thsearch addseq = random nreps = 5 swap = tbr;\n";													# Probably adjust nreps here too if you want
    print $fhout "\tsavetrees file = paup_simpleML_tree.tre format = altnex brlens = yes maxdecimals = 6;\n";
}
# MAXIMUM LIKELIHOOD BOOTSTRAP ANALYSIS
elsif ( $statistics eq "bootstrap" && $method eq "likelihood" )	{
	
	print $fhout "\tset autoclose = yes criterion = like root = outgroup increase = auto storebrlens = yes;\n";
	if ( $outgroup )	{
		print $fhout "\toutgroup $outgroup\;\n";
	}
    print $fhout "\tLset Base=(0.2892 0.2928 0.1309) Nst=6 Rmat=(3.7285 46.5293 1.3888 2.3793 16.4374)\n";
    print $fhout "\tRates = gamma Shape = 0.9350 Pinvar = 0.5691;\n";
    print $fhout "\tbootstrap nreps = $nreps search = heuristic\/ addseq = random swap = tbr hold = 1;\n";
    print $fhout "\tsavetrees from = 1 to = 1 file = paup_MLboot_tree.tre format = altnex brlens = yes savebootp = NodeLabels MaxDecimals = 0;\n";
}
# Default ELSE: defaults to simple parsimony
else	{
	print $fhout "\tset autoclose = yes criterion = parsimony root=outgroup increase=auto storebrlens=yes;\n";
	if ( $outgroup )	{
		print $fhout "\toutgroup $outgroup\;\n";
	}
	print $fhout "\thsearch addseq = random nreps = 1000 swap = tbr hold = 1;\n";
	print $fhout "\tsavetrees file = paup_simpleMP_trees.tre format = altnex brlens = yes;\n";
	print $fhout "\tpscores /TL HI CI RI;\n" if ( $tree_scores == 1 );
	print $fhout "\tcontree all \/ majrule = yes strict = no treefile = paup_simpleMP_consensus.tre;\n";
}

print $fhout "\tlog stop;\n";
print $fhout "\tquit;\n";
print $fhout "END;\n\n";

close $fhout if ( $succout );
exit;
