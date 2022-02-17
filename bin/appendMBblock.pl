#!/usr/bin/perl

# appendMBblock.pl v0.1.0
# Author: MH Seabolt
# Last updated: 8-13-2021

# SYNOPSIS:
# Appends a basic MrBayes block to the end of a NEXUS file.
# The MrBayes block printed is a generic K2P+I+G model, with common ngen and burninfrac settings, 
# which are meant to be reviewed/editted by the user prior to a MrBayes run.

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
my $tree_flag;
my $ngen;
my $outgroup;

sub usage {
	my $usage = "appendMBblock.pl\n
	PURPOSE: 	Appends a basic MrBayes block to the end of a NEXUS file.\n
	NOTE:		The MB block printed is a generic K2P+I+G model, with common ngen and burninfrac settings.
				*** You typically will still want to adjust the model parameters, etc after running this script.
				*** This is only meant to rapidly generate a MrBayes block for you,
				*** so you can focus on the parameter values rather than the block syntax.
	
	USAGE:	appendMBblock.pl -i input.nex -l logfile.log -n 1000000 -t 0 -b out.mb.nex
	-i		input NEXUS file
	-b 		output NEXUS file
	-l		name of log file for MrBayes
	-o 		name of outgroup taxon
	-n 		INT; number of generations to set for MrBayes (default: 10,000,000)
	-t 		use tree in TREES block? INT flag\n";
	print $usage;
}

GetOptions(	'input|i=s' => \$nexus,
			'output|b=s' => \$output,
			'log|l=s' => \$log,
			'use_tree|t=i' => \$tree_flag,
			'ngen|n=i' => \$ngen,
			'outgroup|o=s' => \$outgroup,
			
) or die usage();

# Parameter defaults
$ngen = ( $ngen && int($ngen) >= 1 )? $ngen : 10000000;
$tree_flag = ($tree_flag != 1)? 0 : 1;
$outgroup = ( ($outgroup ne "NA") || ($outgroup ne "") )? $outgroup : undef; 

#######################################################################################
# Open and read the NEXUS file, look for an existing TREES block
my $fh = *STDIN;
my $succin = open(NEX, "<", "$nexus") if ( $nexus ne "--" && -e $nexus );
$fh = *NEX if ( $succin ); 
	my @data_in = <$fh>;
close $fh if ( $succin );

# Get the first tree from the TREES block, if one exists 
#(if there are multiple and you want a different one, you should use this and then just edit the name)
my $tree;
if ( $tree_flag == 1 )	{
	my @trees = grep(/^tree /, @data_in);
	$tree = shift @trees;
	$tree =~ m/tree (.*) =/;
	$tree = $1;
}
else 	{
	$tree = undef;
}

#######################################################################################
# Set the output filehandles
my $succout = open( NEX, ">", "$output" ) if $output ne "--";
my $fhout;
if ( $succout )		{	$fhout = *NEX;			}
else				{	$fhout = *STDOUT;		}

#######################################################################################
# Print the existing NEXUS data
print $fhout "$_" foreach (@data_in);

# Now print the basic MrBayes block
# -- log start
print $fhout "BEGIN MRBAYES;\n";
print $fhout "\tlog start replace filename = $log\;\n";
print $fhout "\tset autoclose = yes nowarn = yes;\n";

# -- set basic model parameters
print $fhout "[set model parameters using K2P+I+G model followed by set priors and independently estimated parameters]\n";
print $fhout "\tlset applyto = (all) nst=2 [HKY model, but fix stationary state frequencies to use K2P model] rates = invgamma ngammacat=4;\n";
print $fhout "\tprset applyto = (all) statefreqpr = fixed(equal) [statefreqpr fixes the stationary state frequencies, converting the HKY model to K2P];\n";
print $fhout "\tunlink revmat = (all) shape = (all) pinvar = (all) statefreq = (all) tratio = (all);\n";
print $fhout "\tshowmodel;\n";

# -- starting tree if we have one
if ( $tree )	{
	print $fhout "[user specified input tree statements]\n";
	print $fhout "\tstartvals tau = $tree V = $tree\;\n";
}

# -- MB execute commands
print $fhout "[execute commands]\n";
print $fhout "\tmcmc ngen = $ngen printfreq = 1000 samplefreq = 1000 nchains = 4 temp = 0.2 checkfreq = 50000 diagnfreq = 1000 stopval = 0.01 stoprule = no;\n";
if ( $outgroup )	{
	print $fhout "\toutgroup $outgroup\;\n";
}
print $fhout "\tsumt relburnin = yes burninfrac = 0.25 contype = halfconpat;\n";
print $fhout "\tsump relburnin = yes burninfrac = 0.25;\n";
print $fhout "\tlog stop;\n";
print $fhout "END;\n\n";

close $fhout if ( $succout );
exit;
