#!/usr/bin/perl

# nexus2phylip.pl v0.1.0
# Author: MH Seabolt
# Last updated: 8-13-2021

# SYNOPSIS:
# A converter from NEXUS to PHYLIP format.

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
my $input = "--";
my $output = "--";
my $fix_headers;
my $revcom;
my $substring;
my $split;
my $fetch;

sub usage {
	my $usage = "nexus2phylip.pl\n
	PURPOSE: 	A converter from NEXUS to PHYLIP format.\n
	USAGE:	nexus2phylip.pl -i input.nexus -o output.phylip -revcom -substr 1:1000 -keep_gaps
	-i		input alignment file in FASTA format
	-c		INT flag, correct headers? (Default OFF)
	-s		INT flag; split sequences into individual FASTA files?
	-o		output file name with .phylip extension
	-f 		STR; string pattern to match for fetching subset of sequences (use 'all' to split all seqs into unique files)
	-revcom		Reverse complement
	-substr		STR; print out only a substring of the original sequences, given as startPosition:endPosition (eg. 1:1000 returns the first 1000 bp)
	\n";
	print $usage;
}

GetOptions(	'input|i=s' => \$input,
			'out|o=s' => \$output,
			'correct|c=i' => \$fix_headers,
			'revcom' => \$revcom,
			'substr=s' => \$substring,
			'split|s=i' => \$split,
			'fetch|f=s' => \$fetch,

) or die usage();

# Parameter setups
$fix_headers = ( $fix_headers && $fix_headers == 1 )? 1 : 0;
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
my %Alignment = ();
foreach my $record ( @lines )	{
	my ($header, $sequence) = split( /\s{2,}/, $record);			# Split the record on at least 2 spaces, just in case we have spaces in the names (which we shouldnt!)
	
	$header =~ s/^[\s]*//g;						# Remove leading whitespace from NEXUS formatting
	$header =~ s/ /_/g;							# Convert any spaces to underscores				
	$sequence = uc $sequence;					# Convert everything to uppercase 
	
	$sequence = revcom( $sequence ) if ( $revcom );
	$sequence = substr( $sequence, $start, $end ) if ( $substring );
	
	# Store the sequences as a hash
	$Alignment{$header} = $sequence;
}

# Correct the weird headers if this flag is activated
my %NewHash = ();
if ( $fix_headers == 1 )	{
	my $newhash = fix_headers(\%Alignment);
	%NewHash = %{$newhash};
	%Alignment = %NewHash;
}

# If split and fetch are invoked
my @heads;
if ( $split == 1 )	{
	my %Subset = ();
	if ( $fetch )	{
		foreach my $header ( sort keys %Alignment )	{
			if ( $header =~ m/$fetch/ )	{
				$Subset{$header} = $Alignment{$header};
			}
			else	{
				next;
			}
		}
		@heads = sort keys %Subset;
	}
	else	{
		@heads = sort keys %Alignment;
	}
	# Print each sequence to its own FASTA file
	foreach my $rec ( @heads )	{
		my $ntax = 1;
		my $nchar = length( $Alignment{$rec} ); #- 1;
		my $longest_name = length( $rec );
	
		# Print each sequence to its own FASTA file
		open PHYLIP, ">", "$rec.phylip" or die;
		print PHYLIP "$ntax $nchar\n";
		my $line = sprintf("$rec" . (" " x ((10 + $longest_name) - length($rec))) );
		my $endgaps;
		if ( $nchar - length($Alignment{$rec}) > 0 )	{
			$endgaps = sprintf("-" x ($nchar - length($Alignment{$rec})));
		}
		else	{
			$endgaps = "";
		}
		my $sequence .= $Alignment{$rec} . $endgaps;
		
		$sequence =~ s/\s+//g;
		$line .= $sequence;
		print PHYLIP "$line\n";
		close PHYLIP;
	}
}
# Else if fetch invoked but not split (ie. we just want a subset of the original fasta file):
elsif ( $split == 0 && $fetch )	{
	my %Subset = ();
	foreach my $header ( sort keys %Alignment )	{
		if ( $header =~ m/$fetch/ )	{
			$Subset{$header} = $Alignment{$header};
		}
		else	{
			next;
		}
	}
	my $ntax = scalar ( keys %Subset );
	my $nchar = length( get_longest_hash_value(\%Subset) ); #- 1;
	my $longest_name = length( get_longest_hash_key(\%Subset) );
	
	# Print each sequence to its own FASTA file
	open PHYLIP, ">", "$output.subset.phylip" or die;
	print PHYLIP "$ntax $nchar\n";
	foreach my $sequence ( sort keys %Subset )	{
		my $line = sprintf("$sequence" . (" " x ((10 + $longest_name) - length($sequence))) );
		my $endgaps;
		if ( $nchar - length($Subset{$sequence}) > 0 )	{
			$endgaps = sprintf("-" x ($nchar - length($Subset{$sequence})));
		}
		else	{
			$endgaps = "";
		}
		my $sequence .= $Subset{$sequence} . $endgaps;
		
		$sequence =~ s/\s+//g;
		$line .= $sequence;
		print PHYLIP "$line\n";
	}
	close PHYLIP;
}

# Else if neither split nor fetch is invoked:
else	{
	# Get some basic info from the hash
	my $ntax = scalar ( keys %Alignment );
	my $nchar = length( get_longest_hash_value(\%Alignment) ); #- 1;
	my $longest_name = length( get_longest_hash_key(\%Alignment) );

	open PHYLIP, ">", "$output.phylip" or die;
	print PHYLIP "$ntax $nchar\n";
	foreach my $sequence ( keys %Alignment )	{
		my $line = sprintf("$sequence" . (" " x ((10 + $longest_name) - length($sequence))) );
		my $endgaps;
		if ( $nchar - length($Alignment{$sequence}) > 0 )	{
			$endgaps = sprintf("-" x ($nchar - length($Alignment{$sequence})));
		}
		else	{
			$endgaps = "";
		}
		my $sequence .= $Alignment{$sequence} . $endgaps;
		
		$sequence =~ s/\s+//g;
		$line .= $sequence;
		print PHYLIP "$line\n";
	}
	close PHYLIP;
}





exit;

########################## SUBRPHYLIPINES  ##################################

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
	my $revcom;
	
	$revcom =~ tr/ACGTacgt/TGCAtgca/;
	$revcom = reverse $revcom;
	
	return ( uc $revcom );
}

