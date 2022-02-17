#!/usr/bin/perl

# nexus2fasta.pl v0.1.0
# Author: MH Seabolt
# Last updated: 8-13-2021

# SYNOPSIS:
# Converts a NEXUS file to FASTA format

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
use Data::Dumper;

# Declare variables
my $input = "--";
my $output = "--";
my $split;
my $fetch;
my $no_dash;
my $keep_gaps;
my $revcom;
my $substring;

my $usage = "nexus2fasta.pl\n
PURPOSE: Converts a NEXUS file to FASTA format.\n
USAGE:	nexus2fasta.pl -i sequences.nexus -o output.fasta
-i		input file in NEXUS format
-o 		output file name (including extensions!)
-s		INT flag; split sequences into individual FASTA files?
-f 		STR; string pattern to match for fetching subset of sequences (use 'all' to split all seqs into unique files)
-revcom		Reverse complement
-keep_gaps	keep gap chars in the sequences ( gaps are removed by default )
-substr		STR; print out only a substring of the original sequences, given as startPosition:endPosition (eg. 1:1000 returns the first 1000 bp)
\n";

GetOptions(	'in|i=s' => \$input,
			'out|o=s' => \$output,
			'split|s=i' => \$split,
			'fetch|f=s' => \$fetch,
			'no_dash|n=i' => \$no_dash,
			'keep_gaps' => \$keep_gaps,
			'revcom' => \$revcom,
			'substr=s' => \$substring,
) or die "$usage";


$no_dash = ( defined($no_dash) && $no_dash == 1 )? 1 : 0;
$split = ( defined($split) && $split != 0 )? 1 : 0;

my $start;
my $end;
if ( $substring )	{
	my @sub = split(":", $substring);
	
	if ( $sub[1] < $sub[0] )	{
		$start = $sub[1] - 1 ;						# -1 is included to account for Perl's 0-indexing
		$end = abs($sub[0] - $start);
		$revcom = 1;
	}
	else	{
		$start = $sub[0] - 1 ;						# -1 is included to account for Perl's 0-indexing
		$end = abs($sub[1] - $start);
	}
}
	
# Open and read the NEXUS file, look for an existing TREES block
my $fh = *STDIN;
my $succin = open(NEX, "<", "$input") if ( $input ne "--" && -e $input );
$fh = *NEX if ( $succin ); 
	my @data_in = <$fh>;
close $fh if ( $succin );

my $flag = 0;
my @lines;
my $type;	my $gap;	my $missing;	my $interleave;
foreach my $line ( @data_in )	{
	chomp $line;
	next if ( $line =~ /^$/ );		# Skip blank lines
	if ( $flag == 1 )	{
		push @lines, $line;
		$flag = 0 if ( $line =~ /END;/ );
		next;
	}
	elsif ( $flag == 0 )	{
		if ( $line =~ /MATRIX/ )	{
			$flag = 1;
			push @lines, $line;
		}
		elsif ( $line =~ /^FORMAT DATATYPE/ )	{
			$line =~ /FORMAT DATATYPE = (.*) GAP = (.*) MISSING = (.*) Interleave = (.*);/;
			( $type, $gap, $missing, $interleave ) = ( $1, $2, $3, $4 );
		}
		next;
	}
}

# Get rid of the first and last elements, which will be "MATRIX" and "END;"
my $trash = shift @lines;
$trash = pop @lines;	
$trash = pop @lines if ( $lines[-1] eq ";" );

# Store the data in a hash
my %Data = ();
foreach my $record ( @lines )	{
	my ($header, $sequence) = split( /\s{2,}/, $record);			# Split the record on at least 2 spaces, just in case we have spaces in the names (which we shouldnt!)

	# Optional fiddling with the FASTA headers
	if ( $no_dash == 1 ) { $header =~ s/\.1/-1/g;  $header =~ s/-//g; } #}		# Remove any dashes and/or periods in the header name if this option is toggled on
	
	$header =~ s/^[\s]*//g;						# Remove leading whitespace from NEXUS formatting
	$header =~ s/ /_/g;							# Convert any spaces to underscores				
	$sequence =~ s/-//g if not ( $keep_gaps );		# Remove any gap chars if the option to keep them is not invoked
	$sequence = uc $sequence;					# Convert everything to uppercase 
		
	# Store the sequences as a hash
	$Data{$header} = $sequence;
}

# If split and fetch are invoked
my @heads;
if ( $split == 1 )	{
	if ( $fetch )	{
		foreach my $header ( sort keys %Data )	{
			if ( $header =~ m/$fetch/ )	{
				push @heads, $header;
			}
			else	{
				next;
			}
		}
	}
	else	{
		@heads = sort keys %Data;
	}
	# Print each sequence to its own FASTA file
	foreach my $rec ( @heads )	{
		my $rec2 = $rec;  $rec2 =~ s/xxx_//g;
		open OUT, ">", "$rec2.fasta" or die "$!\n";
		
		# Final operations -- revcom and substr
		$Data{$rec} = revcom( $Data{$rec} ) if ( $revcom );
		$Data{$rec} = substr( $Data{$rec}, $start, $end ) if ( $substring );
		
		print OUT ">$rec2\n$Data{$rec}\n";
		close OUT;
	}
}
# Else if fetch invoked but not split (ie. we just want a subset of the original fasta file):
elsif ( $split == 0 && $fetch )	{
	foreach my $header ( sort keys %Data )	{
		if ( $header =~ m/$fetch/ )	{
			push @heads, $header;
		}
		else	{
			next;
		}
	}
	# Print each sequence to its own FASTA file
	open OUT, ">", "$output.subset.fasta" or die "$!\n";
	foreach my $rec ( @heads )	{
		my $rec2 = $rec;  $rec2 =~ s/xxx_//g;
		
		# Final operations -- revcom and substr
		$Data{$rec} = revcom( $Data{$rec} ) if ( $revcom );
		$Data{$rec} = substr( $Data{$rec}, $start, $end ) if ( $substring );
		
		print OUT ">$rec2\n$Data{$rec}\n";
	}
	close OUT;
}

# Else if neither split nor fetch is invoked:
else	{
	open OUT, ">", $output or die "$!\n";
	foreach my $key ( keys %Data ) 	{
		$key =~ s/xxx_//g;
		
		# Final operations -- revcom and substr
		$Data{$key} = revcom( $Data{$key} ) if ( $revcom );
		$Data{$key} = substr( $Data{$key}, $start, $end ) if ( $substring );
		print OUT ">$key\n$Data{$key}\n";
	}
	close OUT;
}

exit;

###################### SUBROUTINES ############################

# Reverse complements a sequence
sub revcom	{
	my ( $seq ) = @_;
	my $revcom = $seq;
	
	$revcom =~ tr/ACGTacgt/TGCAtgca/;
	print "$revcom\n";
	$revcom = reverse $revcom;
	
	return ( uc $revcom );
}

