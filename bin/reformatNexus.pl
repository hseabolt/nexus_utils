#!/usr/bin/perl

# reformatNexus.pl v0.1.0
# Author: MH Seabolt
# Last updated: 8-13-2021

# SYNOPSIS:
# A script to perform multiple operations on NEXUS files, particularly DATA blocks of sequences

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

# Declare variables
my $input = "--";
my $output = "--";
my $split;
my $fetch;
my $no_dash;
my $keep_gaps;
my $revcom;
my $substring;
my $drop;
my $no_ambig;

my $usage = "reformatNexus.pl\n
PURPOSE: Reformats a NEXUS file to fetch/drop sequences from an alignment.\n
USAGE:	reformatNexus.pl -i sequences.nexus -o output.nex
-i		input file in NEXUS format
-o 		output file name (including extensions!)
-s		INT flag; split sequences into individual FASTA files?
-f 		STR; string pattern to match for fetching subset of sequences (use 'all' to split all seqs into unique files)
-d		Drop sequences with headers/names matching a given EXPR (can be given as a CSV or as a REGEX)
-revcom		Reverse complement
-keep_gaps	keep gap chars in the sequences ( gaps are removed by default )
-no_ambig 	disallow ambiguous bases in sequences? ( Converts them during processing, Default: OFF )
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
			'drop|d=s' => \$drop,
			'no_ambig' => \$no_ambig,
) or die "$usage";


$no_dash = ( defined($no_dash) && $no_dash == 1 )? 1 : 0;
$split = ( defined($split) && $split != 0 )? 1 : 0;
$output = ($output && $output ne "--" || $output ne "stdout")? $output : "--";

my $start;
my $end;
if ( $substring )	{
	my @sub = split(":", $substring);
	$start = $sub[0] - 1 ;						# -1 is included to account for Perl's 0-indexing
	$end = abs($sub[1] - $sub[0]) - 1;			# abs() needed here in case $end is less than $start (eg implies $revcom)
	
	$revcom = 1 if ( $sub[1] < $sub[0] );
}	

# If user wants to drop sequences
my @drop;
if ( $drop )	{
	@drop = split(",", $drop);
	chomp $_ foreach @drop;
}

# Open and read the NEXUS file, look for an existing TREES block
my $fh = *STDIN;
my $succin = open(NEX, "<", "$nexus") if ( $nexus ne "--" && -e $nexus );
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
$trash = pop @lines;

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
	
	# Correct ambiguity codes -- be careful here: this assumes nucleotide data and will cause major problems if the data is amino acid!
	if ( $no_ambig )	{
		$sequence =~ s/[YKDHBN]/T/g;		# Corrects any acceptable ambiguity code to T if T is an option
		$sequence =~ s/[RWMV]/A/g;			# Corrects any acceptable ambiguity code to A if A is an option
		$sequence =~ s/[S]/G/g;				# Else, call it a G
	}
	
	# Store the sequences as a hash
	$Data{$header} = $sequence;
}

# If we are dropping sequences, do it here
foreach my $drop ( @drop )	{
	foreach my $key ( sort keys %Data )	{
		if ( $key =~ /$drop/ )	{
			delete( $Data{$key} );
		}
	}
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
		
	# Print each sequence to its own NEXUS file
	foreach my $rec ( @heads )	{
		my $rec2 = $rec;  $rec2 =~ s/xxx_//g;
		my $succout = open ( OUT, ">", "$rec2.nex" ) if ( $output ne "--" );
		chomp(my $date = `date`); 
		my $head = sprintf("#NEXUS\n[written $date by reformatNexus.pl]\n\n");
		
		# Final operations -- revcom and substr
		$Data{$rec} = revcom( $Data{$rec} ) if ( $revcom );
		$Data{$rec} = substr( $Data{$rec}, $start, (2+$end) ) if ( $substring );
		
		# Set some basic info ---- in practice, a NEXUS file with a single sequence in it makes absolutely no sense, but who knows what we might need at the time...
		my $ntax = 1;
		my $nchar = length( $Data{$rec} );
		my $longest_name = length( $rec );
		
		if ( $succout )		{ 
			print OUT "$head";	
			print OUT "BEGIN DATA;\n";
			print OUT "DIMENSIONS NTAX=$ntax NCHAR=$nchar;\n";
			print OUT "FORMAT DATATYPE = $type GAP = $gap MISSING = $missing Interleave = no;\n";
			print OUT "\tMATRIX\n";
			print OUT "\t[Note: terminal stop codons removed for protein data, if detected.]\n" if ( $type eq "protein" );
		}
		else				{ 
			print STDOUT "$head";	
			print STDOUT "BEGIN DATA;\n";
			print STDOUT "DIMENSIONS NTAX=$ntax NCHAR=$nchar;\n";
			print STDOUT "FORMAT DATATYPE = $type GAP = $gap MISSING = $missing Interleave = no;\n";
			print STDOUT "\tMATRIX\n";
			print STDOUT "\t[Note: terminal stop codons removed for protein data, if detected.]\n" if ( $type eq "protein" );
		}
		
		# The sequences
		my $line = "\t";
		$line .= sprintf("$rec2" . (" " x ((10 + $longest_name) - length($rec2))) );
		my $endgaps = sprintf("-" x ($nchar - length($Data{$rec}))) if ( $nchar - length($Data{$rec}) > 0 );
		my $sequence .= $Data{$rec} . $endgaps;
		$sequence =~ s/\s+//g;
		$line .= $sequence;
		if ( $succout )		{	print OUT "$line\n";	}
		else				{	print STDOUT "$line\n";	}
		
		if ( $succout )		{	print OUT ";\nEND;\n\n";	}
		else				{	print STDOUT ";\nEND;\n\n";	}
		close OUT if ( $succout );
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
	
	# Print the subset to a NEXUS file
	my $succout = open( OUT, ">", "$output.subset.nex" ) if ( $output ne "--" );
	chomp(my $date = `date`); 
	my $head = sprintf("#NEXUS\n[written $date by reformatNexus.pl]\n\n");
	
	# Set some basic info ---- in practice, a NEXUS file with a single sequence in it makes absolutely no sense, but who knows what we might need at the time...
	my $ntax = scalar @heads;
	my $long = longest_name( @heads );
	my $nchar;
	if ( $substring )	{		$nchar = 2 + ($end - $start);				}
	else				{		$nchar = length( $Data{$long} );	}
	my $longest_name = length( $long );
	
	# Generate the DATA block with the alignment
	if ( $succout )		{ 
		print OUT "$head";	
		print OUT "BEGIN DATA;\n";
		print OUT "DIMENSIONS NTAX=$ntax NCHAR=$nchar;\n";
		print OUT "FORMAT DATATYPE = $type GAP = $gap MISSING = $missing Interleave = no;\n";
		print OUT "\tMATRIX\n";
		print OUT "\t[Note: terminal stop codons removed for protein data, if detected.]\n" if ( $type eq "protein" );
	}
	else				{ 
		print STDOUT "$head";	
		print STDOUT "BEGIN DATA;\n";
		print STDOUT "DIMENSIONS NTAX=$ntax NCHAR=$nchar;\n";
		print STDOUT "FORMAT DATATYPE = $type GAP = $gap MISSING = $missing Interleave = no;\n";
		print STDOUT "\tMATRIX\n";
		print STDOUT "\t[Note: terminal stop codons removed for protein data, if detected.]\n" if ( $type eq "protein" );
	}
	
	# Print the sequences
	foreach my $rec ( @heads )	{
		my $rec2 = $rec;  $rec2 =~ s/xxx_//g;
		
		# Final operations -- revcom and substr
		$Data{$rec} = revcom( $Data{$rec} ) if ( $revcom );
		$Data{$rec} = substr( $Data{$rec}, $start, (2+$end) ) if ( $substring );
		
		my $line = "\t";
		$line .= sprintf("$rec2" . (" " x ((10 + $longest_name) - length($rec2))) );
		my $endgaps = sprintf("-" x ($nchar - length($Data{$rec}))) if ( $nchar - length($Data{$rec}) > 0 );
		my $sequence .= $Data{$rec} . $endgaps;
		$sequence =~ s/\s+//g;
		$line .= $sequence;
		if ( $succout )		{	print OUT "$line\n";	}
		else				{	print STDOUT "$line\n";	}
	}
	if ( $succout )		{	print OUT ";\nEND;\n\n";	}
	else				{	print STDOUT ";\nEND;\n\n";	}
	close OUT if ( $succout );
}

# Else if neither split nor fetch is invoked:
else	{
	# Print a NEXUS file
	my $succout = open( OUT, ">", "$output" ) if ( $output ne "--" );
	chomp(my $date = `date`); 
	my $head = sprintf("#NEXUS\n[written $date by reformatNexus.pl]\n\n");
	
	# Get some basic info from the hash
	my $ntax = scalar ( keys %Data );
	my $nchar = length( get_longest_hash_value(\%Data) ); #- 1;
	my $longest_name = length( get_longest_hash_key(\%Data) );
	
	# Generate the DATA block with the alignment
	if ( $succout )		{ 
		print OUT "$head";	
		print OUT "BEGIN DATA;\n";
		print OUT "DIMENSIONS NTAX=$ntax NCHAR=$nchar;\n";
		print OUT "FORMAT DATATYPE = $type GAP = $gap MISSING = $missing Interleave = no;\n";
		print OUT "\tMATRIX\n";
		print OUT "\t[Note: terminal stop codons removed for protein data, if detected.]\n" if ( $type eq "protein" );
	}
	else				{ 
		print STDOUT "$head";	
		print STDOUT "BEGIN DATA;\n";
		print STDOUT "DIMENSIONS NTAX=$ntax NCHAR=$nchar;\n";
		print STDOUT "FORMAT DATATYPE = $type GAP = $gap MISSING = $missing Interleave = no;\n";
		print STDOUT "\tMATRIX\n";
		print STDOUT "\t[Note: terminal stop codons removed for protein data, if detected.]\n" if ( $type eq "protein" );
	}
	
	# Print the sequences
	foreach my $rec ( keys %Data )	{
		my $rec2 = $rec;  $rec2 =~ s/xxx_//g;
		
		# Final operations -- revcom and substr
		$Data{$rec} = revcom( $Data{$rec} ) if ( $revcom );
		$Data{$rec} = substr( $Data{$rec}, $start, (2+$end) ) if ( $substring );
		
		my $line = "\t";
		$line .= sprintf("$rec2" . (" " x ((10 + $longest_name) - length($rec2))) );
		my $endgaps = sprintf("-" x ($nchar - length($Data{$rec}))) if ( $nchar - length($Data{$rec}) > 0 );
		my $sequence .= $Data{$rec} . $endgaps;
		$sequence =~ s/\s+//g;
		$line .= $sequence;
		if ( $succout )		{	print OUT "$line\n";	}
		else				{	print STDOUT "$line\n";	}
		
		
	}
	if ( $succout )		{	print OUT ";\nEND;\n\n";	}
	else				{	print STDOUT ";\nEND;\n\n";	}
	close OUT if ( $succout );
}

exit;

###################### SUBROUTINES ############################

# Get the longest VALUE in a hash
# Accepts a hash reference as input, returns a scalar VALUE (not its length!)
sub get_longest_hash_value	{
	my ( $hash ) = @_;
	my %Hash = %{$hash};
	my $longest = (sort values %Hash)[0];
	
	foreach my $value ( sort values %Hash )	{
		if ( length($longest) > length($value) )	{
			$longest = $longest;
		}
		else		{
			$longest = $value;
		}
	}
	return($longest);	# - 1
}

# Get the longest VALUE in a hash
# Accepts a hash reference as input, returns a scalar KEY (not its length!)
sub get_longest_hash_key	{
	my ( $hash ) = @_;
	my %Hash = %{$hash};
	my $longest = (sort keys %Hash)[0];
	
	foreach my $value ( keys %Hash )	{
		$longest = length($longest) > length($value)? $longest : $value;
	}
	
	return $longest;
}

# Corrects lengthy fasta headers that are often generated by automatic programs -- ie. find_core_genome.pl // core_genome_reads.pl
# Returns a hash reference with corrected headers
sub fix_headers	{
	my $hash = shift;
	my %Hash = %{$hash};
	
	my %NewHash = ();
	foreach my $header ( sort keys %Hash )	{
		my $seq = $Hash{$header};
		my $new_header = $header;
		$new_header =~ s/\.fasta.*$//g;		
		$NewHash{$new_header} = $seq;
	}
	
	return \%NewHash;

}

# Reverse complements a sequence
sub revcom	{
	my ( $seq ) = @_;
	my $revcom = $seq;
	
	$revcom =~ tr/ACGTacgt/TGCAtgca/;
	print "$revcom\n";
	$revcom = reverse $revcom;
	
	return ( uc $revcom );
}

sub longest_name	{
	my (@names) = @_;
	my $longest = shift @names;
	
	foreach ( @names )	{
		$longest = ( length($longest) >= length($_) )? $longest : $_;
	}
	
	return $longest;
}