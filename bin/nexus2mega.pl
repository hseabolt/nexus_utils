#!/usr/bin/perl

# nexus2mega.pl v0.1.0
# Author: MH Seabolt
# Last updated: 8-13-2021

# SYNOPSIS:
# Converts alignments from NEXUS format to MEGA format.

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
my $wrap;

sub usage {
	my $usage = "nexus2mega.pl\n
	PURPOSE: 	Converts alignments from NEXUS format to sequential MEGA format.\n
	USAGE:	nexus2mega.pl -i input.nex -o output.mega -revcom -substr 1:1000
	-i		input alignment file in NEXUS format
	-c		INT flag, correct headers? (Default OFF)
	-s		INT flag; split sequences into individual FASTA files?
	-o		output file name ( no extensions )
	-f 		STR; string pattern to match for fetching subset of sequences (use 'all' to split all seqs into unique files)
	-w 		INT; wrap text to a given length ( Default OFF )
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
			'wrap|w=i' => \$wrap,

) or die usage();

# Parameter setups
$fix_headers = ( $fix_headers && $fix_headers == 1 )? 1 : 0;
$split = ( defined($split) && $split != 0 )? 1 : 0;
$wrap = ( $wrap && int($wrap) > 0 )? $wrap : undef;

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
	
	if ( $wrap )	{
		$sequence = wrap_text( $sequence, $wrap );
	}
	
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
	if ( $fetch )	{
		foreach my $header ( sort keys %Alignment )	{
			if ( $header =~ m/$fetch/ )	{
				push @heads, $header;
			}
			else	{
				next;
			}
		}
	}
	else	{
		@heads = sort keys %Alignment;
	}
	# Print each sequence to its own MEGA file
	foreach my $rec ( @heads )	{
	
		# Print each sequence to its own MEGA file
		open MEGA, ">", "$rec.mega" or die;
		print MEGA "#mega\n";
		print MEGA "!Title $rec;\n";
		print MEGA "!Format DataType=DNA indel=-;\n\n";
		print MEGA "#$rec\n$Alignment{$rec}\n\n";
		close MEGA;
	}
}
# Else if fetch invoked but not split (ie. we just want a subset of the original fasta file):
elsif ( $split == 0 && $fetch )	{
	foreach my $header ( sort keys %Alignment )	{
		if ( $header =~ m/$fetch/ )	{
			push @heads, $header;
		}
		else	{
			next;
		}
	}
	
	# Print each sequence to its own FASTA file
	open MEGA, ">", "$output.subset.mega" or die;
	print MEGA "#mega\n";
	print MEGA "!Title $output\_subset--$fetch\n";
	print MEGA "!Format DataType=DNA indel=-;\n\n";
	foreach my $rec ( @heads )	{
		print MEGA "#$rec\n$Alignment{$rec}\n\n";
	}
	close MEGA;
}

# Else if neither split nor fetch is invoked:
else	{
	open MEGA, ">", "$output.mega" or die;
	print MEGA "#mega\n";
	print MEGA "!Title $output;\n";
	print MEGA "!Format DataType=DNA indel=-;\n\n";
	foreach my $rec ( keys %Alignment )	{
		print MEGA "#$rec\n$Alignment{$rec}\n\n";
	}
	close MEGA;
}

exit;

########################## SUBROUTINES  ##################################

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
	$revcom = reverse $revcom;
	
	return ( uc $revcom );
}

# Wraps text lines to a given length
sub wrap_text	{
	my ( $text, $length ) = @_;
	
	my $wrapped_text = "";
	for ( my $i=0; $i < length $text; $i += $length )	{
		my $substr = substr($text, $i, $length);
		$wrapped_text .= "$substr\n";
	}
	
	return $wrapped_text;
}


