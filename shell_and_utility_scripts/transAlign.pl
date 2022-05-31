#!/usr/bin/perl -w
#
# transAlign.pl v1.2
# Last modified June 22, 2005 9:36
# (c) Olaf R.P. Bininda-Emonds
#
# Input:
#   DNA sequence data in any of fasta, nexus, (classic or extended) phylip, or Se-Al formats
#
# Output:
#   Aligned DNA sequence data in any of fasta, nexus, (classic or extended) phylip, and/or
#	Se-Al formats. Sequences are aligned by translating them to amino acid sequences (after
#	correction for reading frame and orientation) that are then passed to Clustalw before
#	being back-translated to DNA sequences. Sequences that appear to be frame shifted are
#	either 1) deleted outright, 2) aligned using the amino acid sequence (with the associated
#	errors), or 3) subsequently profile aligned using the DNA sequences (default).
#
# Requires:
#   ClustalW that can be called up remotely. Must either be in $PATH or user must
#	specify path to program.
#
# Usage: transAlign.pl -d<filename> [-b<number> -c<number>] [-f<a|d|x>] [-g<a|f|n>] [-i<f|n|p|s>] [-m<b|g|p>] [-n<number>] [-o<n|pc|pe|s>] [-p<string>] [-r<a|c|i>] [-s<n<number>|p<number>>] [-t] [-u] [-v] [-h]
#	options: -b<number> = alpha level (from 0 to 1) to test for poorly aligning sequences using a one-tailed two-sample t-test
#                         of ClustalW pairwise alignment scores (default = 0 (off))
#            -c<number> = global genetic code (to be used unless otherwise specified)
#                         (default = standard)
#            -d<filename> = file containing raw sequence information
#            -f<a|d|x> = delete frame-shifted sequences outright (x) or align them to remaining sequences either normally as
#						 amino acid sequence (a) or subsequently as DNA (d -- default)
#            -g<a|f|n> = strip all explicit gaps (a), only those flanking the sequence (f -- default), or do not remove any gaps (n)
#            -i<f|n|p|s> = format of sequence file (fasta (f), nexus (n), phylip (p), or Se-Al (s))
#            -m<b|g|p> = protein weight matrix for ClustalW alignment (BLOSUM (b), GONNET (g -- default), or PAM (p))
#            -n<number> = maximum percentage of Ns any sequence can contain (default = 5)
#            -o<n|pc|pe|s> = output results additionally in nexus (n), classic or extended phylip (pc or pe), and/or Se-Al (s) formats (fasta is always output)
#            -p<string> = path to clustal program, if not in $PATH (default = clustalw in $PATH)
#            -r<a|c|i> = order sequences in final output alphabetically by name (a; default), according to order from ClustalW (c),
#						 or in input order from file (i)
#            -s<n<number>|p<number>> = number (n -- default = 1) or percentage (p -- suggested = 2) of stop codons (not including terminal codon)
#                                      each sequence is permitted to have
#            -t = translate input DNA sequences in all possible reading frames of all possible orientations
#                 (default = all possible reading frames of input orientation only)
#            -u = interactive user-input mode
#            -h = print this message and quit
#            -v = verbose output

use strict;
use constant PI => 3.1415926536;
use File::Basename;

# Set user-input defaults and associated parameters
	# Data set variables
		my $inputType = "";	# Options are "fasta", "nexus", "phylip", and "Se-Al"
		my $seqFile = "";
			my ($dataSource, $dir, $ext);
		my (@accNum, %nameLabel, %sequence, %geneticCode, %accPresent);
		my @AAseqs;
			my (%rawAA, %framedSeq, %alignedAA, %finalSeq);
		my $nDelCount = 0;
		my $frameDelCount = 0;
		my (@shiftedSeqs, %deletedSeq);
			my (%alignedDNA, %insertions, %insertionPoint, %deletions, %indel);
		my $seqCount;
		my %clustalOrder;

	# User input variables
		my $globalGenCode = 1;
		my $frameShift = "DNA";	# Options are "AA", "delete", and "DNA" (default)
		my $gapStrip = "flank";	# Options are "all", "flank" (default), and "none"
		my $protMatrix = "GONNET";	# Options are "BLOSUM", "GONNET" (default), and "PAM"
		my $nThreshold = 5;
		my $path = "clustalw";
		my $seqOrder = "alphabetical";	# Options are "alphabetical" (default), "clustal", and "input"
			my %clustalSpeciesCount;
		my $stopForm = "abs";	# Options are "abs" (default) and "perc"
			my $stopThreshold = 1;
		my $allOrientations = 0;
		my $badAlign = 0;	# Alpha level for rejecting badly aligning sequences in a one-tailed two-sample t-test
			my $badAlignCount = 0;

	# Output variables
		my $maxLength = 0;
		my $fastaPrint = 1;
			my $fastaOut;
		my $nexusPrint = 0;
			my $nexusOut;
		my ($phylipTradPrint, $phylipExtPrint) = (0, 0);
			my $phylipOut;
		my $sealPrint = 0;
			my $sealOut;

	# DNA translation codes
		# Ambiguity codes and complements
			my @ambigList = qw(A C G T M R W S Y K V H D B N);
			my %ambigCode = ('A' => 'A', 'C' => 'C', 'G' => 'G', 'T' => 'T',
							'AC' => 'M', 'AG' => 'R', 'AT' => 'W', 'CG' => 'S', 'CT' => 'Y', 'GT' => 'K',
							'ACG' => 'V', 'ACT' => 'H', 'AGT' => 'D', 'CGT' => 'B', 'ACGT' => 'N');
			my %complement = ('A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A',
							 'M' => 'K', 'R' => 'Y', 'W' => 'S', 'S' => 'W', 'Y' => 'R', 'K' => 'M',
							 'B' => 'V', 'D' => 'H', 'H' => 'D', 'V' => 'B', 'N' => 'N', '-' => '-', '?' => '?');
			my (%constitNTlist);
				while (my ($nt, $code) = each %ambigCode)	# Where $nt = key and $code = value
					{
					push @{$constitNTlist{$code}}, $_ foreach (split("",$nt));
					}
			my @orientationList = qw(A C R RC);

		# Genetic codes
			my %transTable = ('1' => 'standard',
							  '2' => 'vertebrate mitochondrial',
							  '3' => 'yeast mitochondrial',
							  '4' => 'mold, protozoan and colenterate mitochondrial and mycoplasam/spiroplasma',
							  '5' => 'invertebrate mitochondrial',
							  '6' => 'ciliate, dasycladacean and hexamita nuclear',
							  '9' => 'echinoderm mitochondrial',
							  '10' => 'euplotid nuclear',
							  '11' => 'bacterial and plant plastid',
							  '12' => 'alternative yeast nuclear',
							  '13' => 'ascidian mitochondrial',
							  '14' => 'alternative flatworm mitochondrial',
							  '15' => 'Blepharisma nuclear',
							  '16' => 'chlorophycean mitochondrial',
							  '21' => 'trematode mitochondrial',
							  '22' => 'Scenedesmus obliquus mitochondrial',
							  '23' => 'Thraustochytrium mitochondrial');
			my %genCodes;
				foreach my $code (qw(1 2 3 4 5 6 9 10 11 12 13 14 15 16 21 22 23))
					{
					$genCodes{$code} = 1;
					}
			my %DNAtoAA;
			my %gapAdd = ('1' => 0, '2' => 2, '3' => 1);
			
	# Miscellaneous variables
		my $verbose = 0;
		my $debug = 0;

# Read in user input
	if (not @ARGV or join(' ', @ARGV) =~ /\s-u/ or $ARGV[0] =~ /^-u/)	# Enter interactive user-input mode
		{
		print "Entering interactive user-input mode. Type \"q\" at any prompt to exit program.\n";

		# Get datafile
			until ($seqFile)
				{
				print "\tEnter name of data file [$seqFile]: ";
				$seqFile = <stdin>;
					chomp ($seqFile);
					exit(0) if ($seqFile eq "q");
				unless (-e $seqFile)
					{
					print "\t\tFile '$seqFile' does not exist\n";
					$seqFile = "";
					}
				}
		
		# Get format of datafile
			my $defaultInput = "autodetect";
				undef $inputType;
			until (defined $inputType)
				{
				print "\tEnter format of file $seqFile (fasta|nexus|phylip|Se-Al) [$defaultInput]:";
				$inputType = <stdin>;
					chomp ($inputType);
					exit(0) if ($inputType =~ /^q/i);
					if (substr($inputType, 0, 1) =~ /^a/i or $inputType eq "")
						{
						$inputType = "autodetect";
						}
					elsif (substr($inputType, 0, 1) =~ /^f/i)
						{
						$inputType = "fasta";
						}
					elsif (substr($inputType, 0, 1) =~ /^n/i)
						{
						$inputType = "nexus";
						}
					elsif (substr($inputType, 0, 1) =~ /^p/i)
						{
						$inputType = "phylip";
						}
					elsif (substr($inputType, 0, 1) =~ /^s/i)
						{
						$inputType = "Se-Al";
						}
					else
						{
						print "\t\tInvalid input ($inputType)\n";
						undef $inputType;
						}
				}
			$inputType = "" if ($inputType eq "autodetect");
		
		# Get how to handle explicit gaps
			my $defaultGap = $gapStrip;
				$gapStrip = "";
			until ($gapStrip eq "all" or $gapStrip eq "flank" or $gapStrip eq "none")
				{
				print "\tEnter whether to strip all explict gaps (all), only those flanking sequence (flank), or none (none) [$defaultGap]: ";
				$gapStrip = <stdin>;
					chomp ($gapStrip);
					exit(0) if ($gapStrip =~ /^q/i);
					if (substr($gapStrip, 0, 1) =~ /^f/i or $gapStrip eq "")
						{
						$gapStrip = "flank";
						}
					elsif (substr($gapStrip, 0, 1) =~ /^a/i)
						{
						$gapStrip = "all";
						}
					elsif (substr($gapStrip, 0, 1) =~ /^n/i)
						{
						$gapStrip = "none";
						}
					else
						{
						print "\t\tInvalid input ($gapStrip)\n";
						$gapStrip = "";
						}
				}
		
		# Get genetic code
			my $defaultCode = $globalGenCode;
				$globalGenCode = 0;
			print "\tGenetic codes available:\n";
				foreach my $code (qw(1 2 3 4 5 6 9 10 11 12 13 14 15 16 21 22 23))
					{
					print "\t\t$code: $transTable{$code}";
					print "\n";
					}
			until (defined $genCodes{$globalGenCode})
				{
				print "\tEnter global genetic code to be applied [$defaultCode]: ";
				$globalGenCode = <stdin>;
					chomp ($globalGenCode);
					exit(0) if ($globalGenCode =~ /^q/i);
					$globalGenCode = $defaultCode if ($globalGenCode eq "");
					print "\t\tInvalid input ($globalGenCode)\n" unless (defined $genCodes{$globalGenCode});
				}
		
		# Get N threshold
			my $defaultN = $nThreshold;
				$nThreshold = "zero";
			while ($nThreshold =~ /\D/)
				{
				print "\tEnter maximum percentage of Ns any sequences can contain [$defaultN]: ";
				$nThreshold = <stdin>;
					chomp ($nThreshold);
					exit(0) if ($nThreshold =~ /^q/i);
					$nThreshold = $defaultN if ($nThreshold eq "");
					print "\t\tInvalid input ($nThreshold)\n" if ($nThreshold =~ /\D/);
				}
		
		# Get stop codon threshold
			my $defaultForm = $stopForm;
				$stopForm = "";
			until ($stopForm eq "abs" or $stopForm eq "perc")
				{
				print "\tEnter whether stop codon threshold should be a number (abs) or a percentage (perc) [$defaultForm]: ";
				$stopForm = <stdin>;
					chomp ($stopForm);
					exit(0) if ($stopForm =~ /^q/i);
					if (substr($stopForm, 0, 1) =~ /^a/i or $stopForm eq "")
						{
						$stopForm = "abs";
						}
					elsif (substr($stopForm, 0, 1) =~ /^p/i)
						{
						$stopForm = "perc";
						}
					else
						{
						print "\t\tInvalid input ($stopForm)\n";
						$stopForm = "";
						}
				}

			my $defaultStop = $stopThreshold;
				$defaultStop = 2 if ($stopForm eq "perc");
				$stopThreshold = "zero";
			while ($stopThreshold =~ /\D/)
				{
				if ($stopForm eq "abs")
					{
					print "\tEnter maximum number of stop codons (excluding terminal one) input sequence can have [$defaultStop]: ";
					}
				else
					{
					print "\tEnter maximum percentage of stop codons (excluding terminal one) input sequence can have [$defaultStop]: ";
					}
				$stopThreshold = <stdin>;
					chomp ($stopThreshold);
					exit(0) if ($stopThreshold =~ /^q/i);
					$stopThreshold = $defaultStop if ($stopThreshold eq "");
					print "\t\tInvalid input ($stopThreshold)\n" if ($stopThreshold =~ /\D/);
				}
		
		# Get orientations to be checked
			my $defaultOrientations = "n";
				undef $allOrientations;
			until (defined $allOrientations)
				{
				print "\tCheck all possible orientations and reading frames for input sequences (y|n) [$defaultOrientations]: ";
				$allOrientations = <stdin>;
					chomp ($allOrientations);
					exit(0) if ($allOrientations =~ /^q/i);
					if (substr($allOrientations, 0, 1) =~ /^y/i)
						{
						$allOrientations = 1;
						}
					elsif (substr($allOrientations, 0, 1) =~ /^n/i or $allOrientations eq "")
						{
						$allOrientations = 0;
						}
					else
						{
						print "\t\tInvalid input ($allOrientations)\n";
						undef $allOrientations;
						}
				}
		
		# Get how to handle frame-shifted sequences
			my $defaultFrame = $frameShift;
				undef $frameShift;
			until (defined $frameShift)
				{
				print "\tEnter how to handle frame-shifted sequences (exclude|AA alignment|DNA alignment) [$defaultFrame]: ";
				$frameShift = <stdin>;
					chomp ($frameShift);
					exit(0) if ($frameShift =~ /^q/i);
					if (substr($frameShift, 0, 1) =~ /^e/i)
						{
						$frameShift = "delete";
						}
					elsif (substr($frameShift, 0, 1) =~ /^a/i)
						{
						$frameShift = "AA";
						}
					elsif (substr($frameShift, 0, 1) =~ /^d/i or $frameShift eq "")
						{
						$frameShift = "DNA";
						}
					else
						{
						print "\t\tInvalid input ($frameShift)\n";
						undef $frameShift;
						}
				}

		# Get protein matrix
			my $defaultProtein = $protMatrix;
				undef $protMatrix;
			until (defined $protMatrix)
				{
				print "\tEnter protein weight matrix to use for ClustalW alignment (BLOSUM|GONNET|PAM) [$defaultProtein]: ";
				$protMatrix = <stdin>;
					chomp ($protMatrix);
					exit(0) if ($protMatrix =~ /^q/i);
					if (substr($protMatrix, 0, 1) =~ /^b/i)
						{
						$protMatrix = "BLOSUM";
						}
					elsif (substr($protMatrix, 0, 1) =~ /^g/i or $protMatrix eq "")
						{
						$protMatrix = "GONNET";
						}
					elsif (substr($protMatrix, 0, 1) =~ /^d/i)
						{
						$protMatrix = "PAM";
						}
					else
						{
						print "\t\tInvalid input ($protMatrix)\n";
						undef $protMatrix;
						}
				}
		
		# Get path to ClustalW
			my $defaultPath = $path;
				undef $path;
			until ($path)
				{
				print "\tEnter full path and name to Clustal program [$defaultPath]: ";
				$path = <stdin>;
					chomp ($path);
					exit(0) if ($path =~ /^q/i);
					$path = $defaultPath if ($path eq "");
				unless (-e $path)
					{
					print "\t\tClustal program not found at $path\n";
					$path = "";
					}
				}

		# Get how to handle poorly aligning sequences
			my $defaultBad = $badAlign;
				undef $badAlign;
			until (defined $badAlign)
				{
				print "\tAlpha level with which to test for poorly aligning sequences (from 0 to 1; 0 = off) [$defaultBad]: ";
				$badAlign = <stdin>;
					chomp ($badAlign);
					exit(0) if ($badAlign =~ /^q/i);
					$badAlign = $defaultBad if ($badAlign eq "");
					if ($badAlign =~ /\d?\.?\d+/)
						{
						if ($badAlign < 0 or $badAlign > 1)
							{
							print "\t\tValue ($badAlign) must be from 0 to 1\n";
							undef $badAlign;
							}
						}
					else
						{
						print "\t\tInvalid input ($badAlign)\n";
						undef $badAlign;
						}
				}

		# Get output order of sequences
			my $defaultOrder = $seqOrder;
				undef $seqOrder;
			until (defined $seqOrder)
				{
				print "\tEnter output order for sequences (alphabetical|clustal|input file) [$defaultOrder]: ";
				$seqOrder = <stdin>;
					chomp ($seqOrder);
					exit(0) if ($seqOrder =~ /^q/i);
					if (substr($seqOrder, 0, 1) =~ /^c/i)
						{
						$seqOrder = "clustal";
						}
					elsif (substr($seqOrder, 0, 1) =~ /^i/i)
						{
						$seqOrder = "input";
						}
					elsif (substr($seqOrder, 0, 1) =~ /^a/i or $seqOrder eq "")
						{
						$seqOrder = "alphabetical";
						}
					else
						{
						print "\t\tInvalid input ($seqOrder)\n";
						undef $seqOrder;
						}
				}

		# Get output formats
			my $defaultNexus = "n";
				undef $nexusPrint;
			until (defined $nexusPrint)
				{
				print "\tOutput results in nexus format (y|n) [$defaultNexus]: ";
				$nexusPrint = <stdin>;
					chomp ($nexusPrint);
					exit(0) if ($nexusPrint =~ /^q/i);
					if (substr($nexusPrint, 0, 1) =~ /^y/i)
						{
						$nexusPrint = 1;
						}
					elsif (substr($nexusPrint, 0, 1) =~ /^n/i or $nexusPrint eq "")
						{
						$nexusPrint = 0;
						}
					else
						{
						print "\t\tInvalid input ($nexusPrint)\n";
						undef $nexusPrint;
						}
				}

			my $defaultPhylip = "n";
				undef $phylipTradPrint;
			until (defined $phylipTradPrint or $phylipExtPrint)
				{
				print "\tOutput results in traditional phylip format (y|n) [$defaultPhylip]: ";
				$phylipTradPrint = <stdin>;
					chomp ($phylipTradPrint);
					exit(0) if ($phylipTradPrint =~ /^q/i);
					if (substr($phylipTradPrint, 0, 1) =~ /^y/i)
						{
						$phylipTradPrint = 1;
						}
					elsif (substr($phylipTradPrint, 0, 1) =~ /^n/i or $phylipTradPrint eq "")
						{
						$phylipTradPrint = 0;
						}
					else
						{
						print "\t\tInvalid input ($phylipTradPrint)\n";
						undef $phylipTradPrint;
						}
				}
				
				if ($phylipTradPrint == 0)	# Check for extended format
					{
					my $defaultPhylip = "n";
						undef $phylipExtPrint;
					until (defined $phylipExtPrint or $phylipExtPrint)
						{
						print "\tOutput results in extended phylip format (y|n) [$defaultPhylip]: ";
						$phylipExtPrint = <stdin>;
							chomp ($phylipExtPrint);
							exit(0) if ($phylipExtPrint =~ /^q/i);
							if (substr($phylipExtPrint, 0, 1) =~ /^y/i)
								{
								$phylipExtPrint = 1;
								}
							elsif (substr($phylipExtPrint, 0, 1) =~ /^n/i or $phylipExtPrint eq "")
								{
								$phylipExtPrint = 0;
								}
							else
								{
								print "\t\tInvalid input ($phylipExtPrint)\n";
								undef $phylipExtPrint;
								}
						}
					}

			my $defaultSeal = "n";
				undef $sealPrint;
			until (defined $sealPrint)
				{
				print "\tOutput results in Se-Al format (y|n) [$defaultSeal]: ";
				$sealPrint = <stdin>;
					chomp ($sealPrint);
					exit(0) if ($sealPrint =~ /^q/i);
					if (substr($sealPrint, 0, 1) =~ /^y/i)
						{
						$sealPrint = 1;
						}
					elsif (substr($sealPrint, 0, 1) =~ /^n/i or $sealPrint eq "")
						{
						$sealPrint = 0;
						}
					else
						{
						print "\t\tInvalid input ($sealPrint)\n";
						undef $sealPrint;
						}
				}
		
		# Get verbose output mode
			my $defaultVerbose = "n";
				undef $verbose;
			until (defined $verbose)
				{
				print "\tOutput verbose results to screen (y|n) [$defaultVerbose]: ";
				$verbose = <stdin>;
					chomp ($verbose);
					exit(0) if ($verbose =~ /^q/i);
					if (substr($verbose, 0, 1) =~ /^y/i)
						{
						$verbose = 1;
						print "\n";
						}
					elsif (substr($verbose, 0, 1) =~ /^n/i or $verbose eq "")
						{
						$verbose = 0;
						}
					elsif (substr($verbose, 0, 1) =~ /^x/i or $verbose eq "")
						{
						$verbose = $debug = 1;
						}
					else
						{
						print "\t\tInvalid input ($verbose)\n";
						undef $verbose;
						}
				}
		}

	elsif (join(' ', @ARGV) =~ /\s-h/ or $ARGV[0] =~ /^-h/)	# Print help screen
		{
		print "Usage: transAlign.pl -d<filename> [-b<number> -c<number>] [-f<a|d|x>] [-g<a|f|n>] [-i<f|n|p|s>] [-m<b|g|p>] [-n<number>] [-o<n|pc|pe|s>] [-p<string>] [-r<a|c|i>] [-s<n<number>|p<number>>] [-t] [-u] [-v] [-h]\n";
		print "\nOptions: -b<number> = alpha level (from 0 to 1) to test for poorly aligning sequences using a\n";
		print "                      one-tailed two-sample t-test of ClustalW pairwise alignment scores (default = 0 (off))\n";
		print "         -c<number> = global genetic code (to be used unless otherwise specified)\n";
			foreach my $code (qw(1 2 3 4 5 6 9 10 11 12 13 14 15 16 21 22 23))
				{
				print "                      $code: $transTable{$code}";
				print " (default)" if ($code eq "1");
				print "\n";
				}
		print "         -d<filename> = file containing raw sequence information\n";
		print "         -f<a|d|x> = delete frame-shifted sequences outright (x) or align them to remaining sequences\n";
		print "                     either normally as amino acid sequence (a) or subsequently as DNA (d -- default)\n";
		print "         -g<a|f|n> = strip all explicit gaps (a), only those flanking the sequence (f -- default), or do not remove any gaps (n)\n";
		print "         -i<f|n|p|s> = format of sequence file (fasta (f), nexus (n), phylip (p), or Se-Al (s))\n";
		print "         -m<b|g|p> = protein weight matrix for ClustalW alignment (BLOSUM (b), GONNET (g -- default), or PAM (p))\n";
		print "         -n<number> = maximum percentage of Ns any sequence can contain (default = 5)\n";
		print "         -o<n|pc|pe|s> = output results additionally in nexus (n), classic or extended phylip (pc or pe), and/or Se-Al (s) formats (fasta is always output)\n";
		print "         -p<string> = path to clustal program, if not in \$PATH (default = clustalw)\n";
		print "         -r<a|c|i> = order sequences in final output alphabetically by name (a; default), according to order from ClustalW (c),\n";
		print "                     or in input order from file (i)\n";
		print "         -s<n<number>|p<number>> = number (n -- default = 1) or percentage (p -- suggested = 2) of stop codons (not including terminal codon)\n";
		print "                                   each sequence is permitted to have\n";
		print "         -t = translate input DNA sequences in all possible reading frames of all possible orientations\n";
		print "              (default = all possible reading frames of input orientation)\n";
		print "         -u = interactive user-input mode\n";
		print "         -h = print this message and quit\n";
		print "         -v = verbose output\n";
		print "\nCitation: Bininda-Emonds, O.R.P. 2005. transAlign: using amino acids to facilitate the multiple alignment of protein-coding DNA sequences. BMC Bioinformatics 6(156) (22Jun05).\n";
		exit(0);
		}
	else	# Process flags
		{
		for (my $i = 0; $i <= $#ARGV; $i++)
			{
			if ($ARGV[$i] =~ /-b(\d?\.?\d+)/)
				{
				$badAlign = $1;
				$badAlign = 0 if ($badAlign < 0 or $badAlign > 1);
				}
			elsif ($ARGV[$i] =~ /-c(\d+)/)
				{
				$globalGenCode = $1;
				if (not defined $genCodes{$globalGenCode})
					{
					print "Don't understand argument: $ARGV[$i]\n";
					print "Usage: transAlign.pl -d<filename> [-b<number> -c<number>] [-f<a|d|x>] [-g<a|f|n>] [-i<f|n|p|s>] [-m<b|g|p>] [-n<number>] [-o<n|pc|pe|s>] [-p<string>] [-r<a|c|i>] [-s<n<number>|p<number>>] [-t] [-u] [-v] [-h]\n";
					exit(1);
					}
				}
			elsif ($ARGV[$i] =~ /^-d(.*)/)
				{
				$seqFile = $1;
				}
			elsif ($ARGV[$i] eq "-fa")
				{
				$frameShift = "AA";
				}
			elsif ($ARGV[$i] eq "-fd")
				{
				$frameShift = "DNA";
				}
			elsif ($ARGV[$i] eq "-fx")
				{
				$frameShift = "delete";
				}
			elsif ($ARGV[$i] eq "-ga")
				{
				$gapStrip = "all";
				}
			elsif ($ARGV[$i] eq "-gf")
				{
				$gapStrip = "flank";
				}
			elsif ($ARGV[$i] eq "-gn")
				{
				$gapStrip = "none";
				}
			elsif ($ARGV[$i] eq "-if")
				{
				$inputType = "fasta";
				}
			elsif ($ARGV[$i] eq "-in")
				{
				$inputType = "nexus";
				}
			elsif ($ARGV[$i] eq "-ip")
				{
				$inputType = "phylip";
				}
			elsif ($ARGV[$i] eq "-is")
				{
				$inputType = "Se-Al";
				}
			elsif ($ARGV[$i] eq "-mb")
				{
				$protMatrix = "BLOSUM";
				}
			elsif ($ARGV[$i] eq "-mg")
				{
				$protMatrix = "GONNET";
				}
			elsif ($ARGV[$i] eq "-mp")
				{
				$protMatrix = "PAM";
				}
			elsif ($ARGV[$i] =~ /^-n(\d+)/)
				{
				$nThreshold = $1;
				}
			elsif ($ARGV[$i] eq "-on")
				{
				$nexusPrint = 1;
				}
			elsif ($ARGV[$i] eq "-opc")
				{
				$phylipTradPrint = 1;
				}
			elsif ($ARGV[$i] eq "-ope")
				{
				$phylipExtPrint = 1;
				}
			elsif ($ARGV[$i] eq "-os")
				{
				$sealPrint = 1;
				}
			elsif ($ARGV[$i] =~ /^-p(.+)/)
				{
				$path = $1;
				print "$path\n";
				}
			elsif ($ARGV[$i] eq "-ra")
				{
				$seqOrder = "alphabetical";
				}
			elsif ($ARGV[$i] eq "-rc")
				{
				$seqOrder = "clustal";
				}
			elsif ($ARGV[$i] eq "-ri")
				{
				$seqOrder = "input";
				}
			elsif ($ARGV[$i] =~ /^-sn(\d+)/)
				{
				$stopThreshold = $1;
				$stopForm = "abs";
				}
			elsif ($ARGV[$i] =~ /^-sp(\d+)/)
				{
				$stopThreshold = $1;
				$stopForm = "perc";
				}
			elsif ($ARGV[$i] eq "-t")
				{
				$allOrientations = 1;
				}
			elsif ($ARGV[$i] eq "-v")
				{
				$verbose = 1;
				}
			elsif ($ARGV[$i] eq "-x")
				{
				$debug = 1;
				$verbose = 1;
				}
			else
				{
				print "Don't understand argument: $ARGV[$i]\n";
				print "Usage: transAlign.pl -d<filename> [-b<number> -c<number>] [-f<a|d|x>] [-g<a|f|n>] [-i<f|n|p|s>] [-m<b|g|p>] [-n<number>] [-o<n|pc|pe|s>] [-p<string>] [-r<a|c|i>] [-s<n<number>|p<number>>] [-t] [-u] [-v] [-h]\n";
				exit(1); 
				}
			}
		}

# Check for input data file and process
	die "ERROR: Must supply name of file containing sequence data.\n" if (not $seqFile);
	
	($dataSource, $dir, $ext) = fileparse($seqFile, qr/\.\D.*/);	
	
	$fastaOut = $dataSource."_tAlign.fasta" if ($fastaPrint);
	$nexusOut = $dataSource."_tAlign.nex" if ($nexusPrint);
	$phylipOut = $dataSource."_tAlign.phylip" if ($phylipTradPrint or $phylipExtPrint);
	$sealOut = $dataSource."_tAlign.seal" if ($sealPrint);
	
	$phylipExtPrint = 0 if ($phylipTradPrint);
	
if ($verbose)
	{
	print "The following parameters will be used:\n";
	print "\tName of input file: $seqFile (type is $inputType)\n";
	print "\tGlobal genetic code: $globalGenCode ($transTable{$globalGenCode})\n";
	print "\tMaximum percentage of Ns allowed per sequence: $nThreshold\n";
	print "\tForm of stop codon threshold: $stopForm\n";
	print "\t\tStop codon threshold: $stopThreshold";
		print "%" if ($stopForm eq "perc");
		print "\n";
	print "\tChecking all orientations: ";
		if ($allOrientations)
			{
			print "yes\n";
			}
		else
			{
			print "no\n";
			}
	print "\tFrame-shifted sequences will be: ";
		if ($frameShift eq "delete")
			{
			print "deleted\n";
			}
		elsif ($frameShift eq "AA")
			{
			print "aligned as AAs\n";
			}
		else
			{
			print "aligned subsequently as DNA\n";
			}
	print "\tProtein matrix for AA alignment: $protMatrix\n";
	print "\tTesting for poorly aligning sequences: ";
		if ($badAlign)
			{
			print "with alpha level $badAlign\n";
			}
		else
			{
			print "off\n";
			}
	print "\tPath to ClustalW: $path\n";
	print "\tOutput order of sequences: $seqOrder\n";
	print "\tName of output file(s):\n";
		print "\t\t$fastaOut\n" if ($fastaPrint);
		print "\t\t$nexusOut\n" if ($nexusPrint);
		print "\t\t$phylipOut (classic format)\n" if ($phylipTradPrint);
		print "\t\t$phylipOut (extended format)\n" if ($phylipExtPrint);
		print "\t\t$sealOut\n" if ($sealPrint);
	print "\n";
	}

print "Establishing translation tables for all genetic codes ...\n";
	geneticCoder();

# Read in sequence data
	seqRead($seqFile);
	die "\nERROR: Could not read in sequences from file $seqFile\n" if (not @accNum);

# Process each sequence
	print "\nProcessing each sequence to determine optimal reading frame ";
		if ($allOrientations)
			{
			print "(in all possible orientations) and any frame shifts\n";
			}
		else
			{
			print "(in \"as is\" orientation only) and any frame shifts\n";
			}

	my $frameZero = time;
	foreach my $accession (@accNum)
		{
		$accPresent{$accession} = 1;
		$deletedSeq{$accession} = 0;	# Sequence defaults as included
		
		# Remove explcit gaps as desired
			unless ($gapStrip eq "none")
				{
				if ($gapStrip eq "flank")	# Remove leading and trailing gaps
					{
					$sequence{$accession} =~ s/^\-+//;
					$sequence{$accession} =~ s/\-+$//;
					}
				else	# Remove all gaps
					{
					$sequence{$accession} =~ s/\-//g;
					}
				}
		
		my $rawSeq = $sequence{$accession};
		if ($verbose)
			{
			print "\tSequence: $nameLabel{$accession}";
				print " ($accession)" unless ($accession =~ /^tAlign/);
				print "\n";
			}
		
		# Check if sequence possesses too many Ns
			if ($nThreshold < 100)
				{
				my $nCount = ($sequence{$accession} =~ tr/N//);
					$nCount = sprintf("%.2f", $nCount / length($sequence{$accession}) * 100);
				if ($nCount > $nThreshold)
					{
					print "\tWARNING: Percentage of Ns in sequence ($nCount) ";
						unless ($verbose)
							{
							print "for $nameLabel{$accession} ";
							print "($accession) " unless ($accession =~ /^tAlign/);
							}
						print "exceeds threshold; sequence will not be processed\n\n";
					$deletedSeq{$accession} = 1;
					$nDelCount++;
					next;
					}
				}

		# Determine optimal reading frame for sequence as a whole
			my (@bestAll, @bestAsIs, %DNAseq, %protSeq);
			my $minAllStops = length($rawSeq);
			
			# Translate all possible reading frames and all possible orientations for sequence and count number of stop codons to derive best reading frame
			# Order of input means that results are preferred according to simplest solutions in term of number of changes: 0N, 1N, 2N, 0C, 1C, 2C, 0R, 1R, 2R, 0RD, 1RC, 2RC (i.e., reading frame change first, followed by as is, complemented, reverse, and reverse complemented orientations)
				foreach my $orientation (@orientationList)
					{
					next if (not $allOrientations and $orientation ne "A");	# No point in checking alternative orientations
					for (my $frame = 1; $frame <= 3; $frame++)
						{
						my $frameName = $frame . $orientation;
						print "\t\tReading frame: $frame " if ($debug);
						if ($orientation eq "R")
							{
							$DNAseq{$frameName} = scalar reverse $rawSeq;
							print "(reversed)\n" if ($debug);
							}
						elsif ($orientation eq "C")
							{
							print "(complemented)\n" if ($debug);
							for (my $nt = 0; $nt < length($rawSeq); $nt++)
								{
								if (not defined $complement{substr($rawSeq, $nt, 1)})
									{
									$DNAseq{$frameName} .= "?";
									}
								else
									{
									$DNAseq{$frameName} .= $complement{substr($rawSeq, $nt, 1)};
									}
								}
							}
						elsif ($orientation eq "RC")
							{
							print "(reverse complemented)\n" if ($debug);
							$DNAseq{$frameName} = scalar reverse $DNAseq{$frame."C"};
							}
						else
							{
							print "(as is)\n" if ($debug);
							$DNAseq{$frameName} = $rawSeq;
							}

						# Adjust for reading frame
							$DNAseq{$frameName} = ("-" x $gapAdd{$frame}) . $DNAseq{$frameName};	# Create new DNA sequences for each reading frame

						# Translate sequence and count stop codons
							$protSeq{$frameName} = translate($DNAseq{$frameName}, $geneticCode{$accession});
								print "\t\t\tAA sequence: $protSeq{$frameName}\n" if ($debug);
					
							my $stopCount = ($protSeq{$frameName} =~ tr/\*//);
								$stopCount-- if (substr($protSeq{$frameName}, -1) eq "*");	# Do not penalize if last codon is stop codon
								print "\t\t\t\t$stopCount stop codons\n" if ($debug);
								
								# Find percentage of stop codons after appearance of first stop codon
									if ($stopForm eq "perc" and $stopCount > 0)
										{
										my $position = 0;
										while ($position < length ($protSeq{$frameName}))	# Get position of first stop codon
											{
											$position++;
											last if (substr($protSeq{$frameName}, $position, 1) eq "*");
											}
											
										# Truncate string and get number of informative residues
											print "\t\t\t\t\tFirst stop codon found at position $position\n" if ($debug);
											my $chopLength = 0 - (length($protSeq{$frameName}) - $position);
											my $infLength = (substr($protSeq{$frameName}, $chopLength) =~ tr/A-Z\*//);

										# Adjust stop count
											$stopCount = sprintf ("%.1f", $stopCount / $infLength * 100);
											print "\t\t\t\t\tPercentage of stop codons in remaining sequence: $stopCount\n\n" if ($debug);
										}
							
							if ($stopCount == $minAllStops) 
								{
								push @bestAll, $frameName;
								}
							if ($stopCount < $minAllStops) 
								{
								$minAllStops = $stopCount;
								undef @bestAll;
								push @bestAll, $frameName;
								}
						}
					}
			print "\t\tBest reading frame(s): ".join(" ", @bestAll)."\n" if ($allOrientations and $verbose);
			print "\t\tBest \"as is\" reading frame(s): ".join(" ", @bestAll)."\n" if (not $allOrientations and $verbose);

		# Process sequence accordingly based on whether a frame shift has apparently occurred
			my $bestFrame = shift(@bestAll);
			
			if ($minAllStops > $stopThreshold)
				{
				print "\t" if ($verbose);
				print "\tWARNING: Number of stop codons in best frame ($minAllStops) " if ($stopForm eq "abs");
				print "\tWARNING: Percentage of stop codons in best frame ($minAllStops) " if ($stopForm eq "perc");
				unless ($verbose)
					{
					print "for $nameLabel{$accession} ";
					print "($accession) " unless ($accession =~ /^tAlign/);
					}
				print "exceeds threshold ($stopThreshold); inspection of sequence is recommended\n\n";
				}
			
			if ($minAllStops <= $stopThreshold or $frameShift eq "AA")	# Process best sequence; either below threshold or frame shifted sequence to be processed regardless
				{
				push @AAseqs, $accession;
#				$frameCount{substr($bestFrame, 0, 1)}++;	# Increment reading frame count
				
				$framedSeq{$accession} = $DNAseq{$bestFrame};
					# Pad end of sequence with gaps if necessary
						my $numFullCodons = int(length($framedSeq{$accession}) / 3);
						my $numGaps = 3 - (length($framedSeq{$accession}) - ($numFullCodons * 3));
							$numGaps = 0 if ($numGaps == 3); 
						$framedSeq{$accession} .= "-" x $numGaps;
				
				$rawAA{$accession} = $protSeq{$bestFrame};

				if ($verbose)
					{
					print "\t\t\tNo frame shift correction desired by user\n" if ($minAllStops > $stopThreshold and $frameShift eq "delete");
					print "\t\t\tTranslated sequence written to file\n\n";
					}
				}
			else
				{
				push @shiftedSeqs, $accession;
				if ($frameShift eq "delete")	# Simply delete (ignore) sequence
					{
					print "\t\t\tSequence will be deleted\n\n" if ($verbose);
					$deletedSeq{$accession} = 1;
					$frameDelCount++;
					}
				else	# Store raw sequence in separate file for subsequent DNA alignment
					{
					print "\t\t\tSequence will be aligned subsequently as DNA\n\n" if ($verbose);
					}
				}
		}
	printf "\n\tTime taken to frame sequences: %s seconds\n", time - $frameZero;

# Print summary statistics and check whether ok to continue
	if ($verbose)
		{
		print "\nOf $seqCount sequences read in:\n";
		printf "\t%s will be aligned as amino acids\n", scalar(@AAseqs);
		printf "\t%s had apparent frame shifts and will be ", scalar(@shiftedSeqs);
			printf "deleted\n" if ($frameShift eq "delete");
			printf "aligned as amino acids\n" if ($frameShift eq "AA");
			printf "aligned subsequently as DNA\n" if ($frameShift eq "DNA");
		print "\t$nDelCount had more Ns than threshold and were not processed\n";
		}

	if (scalar(@AAseqs) > 1)	# Process each sequence and save to file
		{
		open (AA, ">clustalAA_fasta.txt") or die "Cannot write translated sequences to file clustalAA_fasta.txt\n";
		foreach my $accession (@AAseqs)
			{
			$rawAA{$accession} =~ s/\?/-/g;	# ClustalW will ignore ?s; code as gaps to preserve information
			$rawAA{$accession} =~ s/\*/-/g;	# ClustalW will ignore stop codons; code as gaps to preserve information

			my $fastaSeq = $rawAA{$accession};
				my $breakPoint = 79;
				until ($breakPoint > length($fastaSeq))
					{
					my $replaceString = "\n" . substr($fastaSeq, $breakPoint, 1);
					substr($fastaSeq, $breakPoint, 1) = $replaceString;
					$breakPoint += 80;
					}
			print AA ">$accession\n$fastaSeq\n";
			}
		close AA;
		}
	else
		{
		printf "\nNOTE: Insufficient sequences (%s) were found to be suitable to be aligned as amino acids. Data set should be aligned as DNA\n", scalar(@AAseqs);
		exit(0);
		}

# Call up ClustalW to align AA sequences
	print "\nPassing AA sequence data to ClustalW for alignment ...\n";
	print "\tRunning status from ClustalW saved to file clustal_AA_log.txt\n";
	my $clustalZero = time;
		my $clustalString = "-infile=clustalAA_fasta.txt -type=protein -align -matrix=$protMatrix -outfile=clustalAA_fasta.aln > clustal_AA_log.txt";
		system("$path $clustalString");

	printf "\n\tTime taken to align amino acid sequences: %s seconds\n", time - $clustalZero;
		
# Read in clustal log file to identify badly aligning sequences (if desired)
	if ($badAlign and $path =~ /clustal/i and scalar(@AAseqs) >= 20)
		{
		my $badZero = time;
		print "\nExamining for poorly aligning sequences ...\n";
		my %clustalID;
		my (%scoreSum, %scoreVar);
		
		# Read in clustal log file and process data
			open (LOG, "<clustal_AA_log.txt") or die "Cannot open clustal log file, clustal_AA_log.txt\n";
			print "\tParsing ClustalW output file ...\n" if ($debug);
				while (<LOG>)
					{
					chomp;
					next unless ($_ =~ /^Sequence/);
					
					if ($_ =~ /^Sequence\s+\d+:/)	# Relate accession numbers to clustal sequence numbers
						{
						my @logLine = split(/\s+/, $_);
							my $seqNo = $logLine[1];
								$seqNo =~ s/://;
							my $accessionID = $logLine[2];
							$clustalID{$seqNo} = $accessionID;
						print "\t\tClustalID $seqNo assigned to sequence $accessionID\n" if ($debug);
						}
					if ($_ =~ /^Sequences\s+\(/)	# Get pairwise alignment score for each pair of sequences
						{
						my @logLine = split(/\s+/, $_);
							my ($seq1, $seq2);
							if ($logLine[1] =~ /\((\d+)\:(\d+)\)/)
								{
								$seq1 = $1;
								$seq2 = $2;
								}
							my $alignScore = pop @logLine;
						print "\t\tAlignment score $alignScore assigned to sequences $clustalID{$seq1} and $clustalID{$seq2}\n" if ($debug);

						# Increment average and variance counters
							$scoreSum{$clustalID{$seq1}} += $alignScore;
							$scoreSum{$clustalID{$seq2}} += $alignScore;
							$scoreSum{all} += $alignScore;

							$scoreVar{$clustalID{$seq1}} += $alignScore**2;
							$scoreVar{$clustalID{$seq2}} += $alignScore**2;
							$scoreVar{all} += $alignScore**2;
						}
					}
			close LOG;
		
		# Perform one-tailed two-sample t-test
			# Determine relevant samples sizes
				my $numSeq = scalar(@AAseqs);
				my $numAcc = $numSeq - 1;
				my $numPop = (0.5 * $numSeq**2) - (1.5 * $numSeq) + 1;

			foreach my $accession (@AAseqs)
				{
				if ($debug)
					{
					print "\tSequence: $nameLabel{$accession}";
						print " ($accession)" unless ($accession =~ /^tAlign/);
						print "\n";
					}
				# Get mean and SS for accession
					my $avgAcc = $scoreSum{$accession} / $numAcc;
					my $ssAcc = ($scoreVar{$accession} - $avgAcc * $avgAcc * $numAcc);

				# Get mean and SS for remaining population
					$scoreSum{test} = $scoreSum{all} - $scoreSum{$accession};
						my $avgPop = $scoreSum{test} / $numPop;
					$scoreVar{test} = $scoreVar{all} - $scoreVar{$accession};
						my $ssPop =($scoreVar{test} - $avgPop * $avgPop * $numPop);
						
					my $ssPooled = ($ssAcc + $ssPop) / ($numAcc + $numPop - 2);
						$ssPooled = 0.00001 if ($ssPooled == 0);	# Safety in case ssPooled = 0
						
				if ($avgAcc >= $avgPop)	# Can never give significant result for a one-tailed test
					{
					printf "\t\tAccepted: mean score of %.3f >= mean population score of %.3f\n", $avgAcc, $avgPop if ($debug);
					next;
					}
				
				# Calculate t-statistic for avgAcc < avgPop
					my $t = ($avgAcc - $avgPop) / sqrt(($ssPooled / $numAcc) + ($ssPooled / $numPop));
					my $df = $numSeq - 2;
					
					my $lowerP = _subtprob($df, $t);
						$lowerP = 1 - $lowerP if ($t > 0);	# Actually not necessary given that stats not calculated in this instance

					if ($lowerP < $badAlign / scalar(@AAseqs))	# Sequence shows significantly worse alignment scores than remaining ones (with conservative Bonferroni correction); reject
						{
						$deletedSeq{$accession} = 1;
						$badAlignCount++;
						printf "\t\tRejected: mean score of %.3f < mean population score of %.3f (t = %.3f; p = %.3f)\n", $avgAcc, $avgPop, $t, $lowerP if ($debug);
						if ($verbose and not $debug)
							{
							printf "\tWARNING: Average alignment score (%.3f) for $nameLabel{$accession} ", $avgAcc;
							print "($accession) " unless ($accession =~ /^tAlign/);
							printf "is significantly worse than that of remaining sequences (%.3f); sequence excluded\n\n", $avgPop;
							}
						}
					else
						{
						printf "\t\tAccepted: mean score of %.3f >= mean population score of %.3f (t = %.3f; p = %.3f)\n", $avgAcc, $avgPop, $t, $lowerP if ($debug);
						}
				}
			printf "\n\tTime taken to identify $badAlignCount poorly aligning sequences: %s seconds\n", time - $badZero;
		}

# Read in clustal output and back-translate using framed sequences
	print "\nReading in ClustalW aligned AA data and back-translating to DNA ...\n";
		if ($seqOrder eq "clustal")
			{
			undef (@accNum);
			undef %clustalSpeciesCount;
			}
	setLineBreak("clustalAA_fasta.aln");
	open (ALIGN, "<clustalAA_fasta.aln") or die "Cannot open file containing aligned AA data, clustalAA_fasta.aln\n";
		while (<ALIGN>)
			{
			chomp;
			next unless ($_);
			next unless ($_ =~ /\S/);	# ClustalW will occasionally put in a line of all spaces
			next if ($_ =~ /PileUp/ or $_ =~ /Check:/ or $_ =~ /\\\\/);	# Fixes for additional lines in GCG format
			next if ($_ =~ /\s+\d+\s+\d+/);	# Fix for header line in PHYLIP format
			
			my (@alignedLine) = split(/\s+/);
				my $species = shift(@alignedLine);
				my $sequence = join('', @alignedLine);
					$sequence =~ s/\./\-/g;	# Fix for GCG format, where gaps are periods
			if (defined $accPresent{$species} and $deletedSeq{$species} == 0)	# Only read in recognized sequences that have not been deleted for any reason
				{
				$alignedAA{$species} .= $sequence;
				if ($seqOrder eq "clustal")
					{
					$clustalSpeciesCount{$species}++;
					push @accNum, $species if ($clustalSpeciesCount{$species} == 1);
					}
				}
			}
	close ALIGN;
	
	# Quickly parse aligned sequences to remove any constant gaps introduced by removing poorly aligning sequences
		if ($badAlignCount)
			{
			my $gapZero = time;
			print "\nChecking for constant gaps introduced by removal of poorly aligning species ...\n";
			my $delPosition = -1;
			my $totalLength = length($AAseqs[0]);
			for (my $bp = $totalLength - 1; $bp >= 0; $bp--)	# Checked backwards so that constant sites can be deleted "simultaneously"
				{
				my $gapCount = 0;
				foreach my $accession (@AAseqs)
					{
					next if ($deletedSeq{$accession});
					if (substr($alignedAA{$accession}, $bp, 1) eq "-")	# Increment gap count as appropriate
						{
						$gapCount++;
						}
					else
						{
						last unless ($delPosition >= 0);
						}
					if ($delPosition >= 0)	# Delete following position if all gaps
						{
						substr($alignedAA{$accession}, $delPosition, 1) = "";
						}
					}
				if ($gapCount == (scalar(@AAseqs) - $badAlignCount))	# All sequences have a gap at position; mark for deletion
					{
					$delPosition = $bp;
					}
				else
					{
					$delPosition = -1;
					}
				}
			if ($delPosition == 0)	# First sites is all gaps, but cannot be cleared in above loop
				{
				foreach my $accession (@AAseqs)
					{
					substr($alignedAA{$accession}, $delPosition, 1) = "";
					}
				}
			printf "\n\tTime taken to remove any constant gaps: %s seconds)\n", time - $gapZero;
			}

	# Relate back to framed sequences
		my $transZero = time;
		foreach my $accession (@AAseqs)
			{
			next if ($deletedSeq{$accession});
			if ($debug)
				{
				print "\tSequence: $nameLabel{$accession}";
					print " ($accession)" unless ($accession =~ /^tAlign/);
					print "\n";
				}
			my $alignStart = 0;
			for (my $framedPos = 0; $framedPos < length($framedSeq{$accession}); $framedPos += 3)	# Sequentially read in codons from DNA sequence
				{
				my $framedCodon = substr($framedSeq{$accession}, $framedPos, 3);
				my $framedAA;
					if ($framedCodon =~ /-/)
						{
						$framedAA = "-";
						}
					else
						{
						$framedAA = $DNAtoAA{$geneticCode{$accession}}{$framedCodon};
						$framedAA = "-" if ($framedAA eq "?" or $framedAA eq "*");
						}
				print "\t\tPosition $framedPos (sequence codon: $framedCodon --> $framedAA)\n" if ($debug);
				
				if ($framedAA eq "-")
					{
					my $gapCodon .= $framedCodon;
					
					# Find end of ambiguous stretch in framed DNA
						unless ($framedPos + 3 > length($framedSeq{$accession}) - 1)	# Already at end of sequence
							{
							until ($framedAA ne "-" or $framedPos > length($framedSeq{$accession}) - 1)
								{
								$framedPos += 3;
								my $nextCodon = substr($framedSeq{$accession}, $framedPos, 3);
								unless ($framedPos > length($framedSeq{$accession}) - 1)
									{
									if ($nextCodon =~ /-/)
										{
										$framedAA = "-";
										}
									else
										{
										$framedAA = $DNAtoAA{$geneticCode{$accession}}{$nextCodon};
										$framedAA = "-" if ($framedAA eq "?" or $framedAA eq "*");
										}
									print "\t\t\tPosition $framedPos (sequence codon: $nextCodon --> $framedAA)\n" if ($debug);
									$gapCodon .= $nextCodon if ($framedAA eq "-");
									}
								}
							$framedPos -= 3;
							}
						print "\n\t\t\tAmbiguous DNA sequence: $gapCodon\n\n" if ($debug);
					
					# Find end of gap in clustal alignment
						my $clustalAA = substr($alignedAA{$accession}, $alignStart, 1);
							my $gapClustal .= $clustalAA;
							
						print "\t\t\tClustal AA: $clustalAA\n" if ($debug);
						until ($clustalAA ne "-" or $alignStart == length($alignedAA{$accession}))
							{
							$alignStart++;
							$clustalAA = substr($alignedAA{$accession}, $alignStart, 1);
							$gapClustal .= $clustalAA if ($clustalAA eq "-");
							print "\t\t\tClustal AA: $clustalAA\n" if ($debug);
							}
							
						$alignStart--;
						print "\n\t\t\tGap sequence: $gapClustal\n" if ($debug);
#							print "\t\t\t\tEnd of gap ...\n" if ($debug);
					
					# Determine whether ambiguous stretch fits better to start or end of gap
						my $gapSequence;
						if ($gapCodon =~ /^-/ or not defined $finalSeq{$accession})	# Stretch begins with gap or at start of sequence; fits better to end
							{
							$gapSequence .= "-" x ((length($gapClustal) * 3) - length($gapCodon));
							$gapSequence .= $gapCodon;
							print "\t\t\t\tSequence $gapSequence fitted to end of gap\n" if ($debug);
							}
						else	# Fits better to start
							{
							$gapSequence .= $gapCodon;
							$gapSequence .= "-" x ((length($gapClustal) * 3) - length($gapCodon));
							print "\t\t\t\tSequence $gapSequence fitted to start of gap\n" if ($debug);
							}
						$finalSeq{$accession} .= $gapSequence;
					}
				else
					{
					my $clustalAA = substr($alignedAA{$accession}, $alignStart, 1);
						print "\t\t\tClustal AA: $clustalAA\n" if ($debug);
						if ($clustalAA eq "-")	# Alignment gap inserted by clustal
							{
							until ($clustalAA ne "-" or $alignStart == length($alignedAA{$accession}))	# Find end of gap
								{
								$alignStart++;
								$finalSeq{$accession} .= "---";	# Only ClustalW can add three gaps to sequence; transAlign cannot
								$clustalAA = substr($alignedAA{$accession}, $alignStart, 1);
								print "\t\t\tClustal AA: $clustalAA\n" if ($debug);
								}
							print "\t\t\t\tEnd of gap ...\n" if ($debug);
							$finalSeq{$accession} .= $framedCodon;
							}
						else
							{
							$finalSeq{$accession} .= $framedCodon;
							}
					}
				$alignStart++;
				}
			if (length($alignedAA{$accession}) * 3 > length($finalSeq{$accession}))	# Sequence has ended cleanly; need to add trailing gaps
				{
				$finalSeq{$accession} .= "-" x ((length($alignedAA{$accession}) * 3) - length($finalSeq{$accession}));
				}
			$maxLength = length($finalSeq{$accession}) if (length($finalSeq{$accession}) > $maxLength);
			}
		
		# Check and correct for any constant leading gaps
			# Make a quick dirty consensus
				my $leadingConsensus = "--";
				foreach my $accession (@AAseqs)
					{
					last if ($leadingConsensus eq "NN");
					substr($leadingConsensus, 0, 1, "N") if (substr($finalSeq{$accession}, 0, 1) ne "-");
					substr($leadingConsensus, 1, 1, "N") if (substr($finalSeq{$accession}, 1, 1) ne "-");
					}

			# Correct for any constant leading gaps				
				my $leadingGaps = 0;
					$leadingGaps = ($leadingConsensus =~ tr/\-//);

				if ($leadingGaps)
					{
					print "\n\tStripping $leadingGaps constant leading gap(s) from all sequences ...\n" if ($verbose);
					substr($finalSeq{$_}, 0, $leadingGaps, "") foreach (@AAseqs);
					}

		printf "\n\tTime taken to back-translate to DNA: %s seconds\n", time - $transZero;

# Profile align any frame shifted sequences (if desired) and infer frame shift locations
	if (@shiftedSeqs and $frameShift eq "DNA")
		{
		# Print both DNA data sets to separate files
			open (P1, ">clustalP1_fasta.txt") or die "Cannot write aligned DNA sequences to file clustalP1_fasta.txt\n";
				foreach my $accession (@AAseqs)
					{
					my $fastaSeq = $finalSeq{$accession};
						my $breakPoint = 79;
						until ($breakPoint > length($fastaSeq))
							{
							my $replaceString = "\n" . substr($fastaSeq, $breakPoint, 1);
							substr($fastaSeq, $breakPoint, 1) = $replaceString;
							$breakPoint += 80;
							}
					print P1 ">$accession\n$fastaSeq\n";
					}
			close P1;

			open (P2, ">clustalP2_fasta.txt") or die "Cannot write frame shifted DNA sequences to file clustalP2_fasta.txt\n";
				foreach my $accession (@shiftedSeqs)
					{
					my $fastaSeq = $sequence{$accession};
						my $breakPoint = 79;
						until ($breakPoint > length($fastaSeq))
							{
							my $replaceString = "\n" . substr($fastaSeq, $breakPoint, 1);
							substr($fastaSeq, $breakPoint, 1) = $replaceString;
							$breakPoint += 80;
							}
					print P2 ">$accession\n$fastaSeq\n";
					}
			close P2;

		# Call up ClustalW to align as profiles; number of comparisons is (aligned * shifted) + (shifted choose 2)
			print "\nPassing DNA sequence data to ClustalW for alignment as profiles ...\n";
			print "\tRunning status from ClustalW saved to file clustal_profile_log.txt\n";
			my $clustalZero = time;
				my $clustalString = "-profile1=clustalP1_fasta.txt -profile2=clustalP2_fasta.txt -sequences -outfile=clustalP2_fasta.aln > clustal_profile_log.txt";
				system("$path $clustalString");

			printf "\n\tTime taken to align DNA sequences: %s seconds\n", time - $clustalZero;

		# Infer frame shifts locations
			print "\nReading in ClustalW aligned DNA data and inferring indel positions ...\n";
				undef %finalSeq;
				if ($seqOrder eq "clustal")
					{
					undef (@accNum);
					undef %clustalSpeciesCount;
					}
			setLineBreak("clustalP2_fasta.aln");
			open (ALIGN, "<clustalP2_fasta.aln") or die "Cannot open file contained aligned DNA profile data, clustalP2_fasta.aln\n";
				while (<ALIGN>)
					{
					chomp;
					next unless ($_);
					next unless ($_ =~ /\S/);	# ClustalW will occasionally put in a line of all spaces
					next if ($_ =~ /PileUp/ or $_ =~ /Check:/ or $_ =~ /\\\\/);	# Fixes for additional lines in GCG format
					next if ($_ =~ /\s+\d+\s+\d+/);	# Fix for header line in PHYLIP format

					my (@alignedLine) = split(/\s+/);
						my $species = shift(@alignedLine);
						my $sequence = join('', @alignedLine);
							$sequence =~ s/\./\-/g;	# Fix for GCG format, where gaps are periods
					unless (not defined $accPresent{$species})
						{
						$finalSeq{$species} .= $sequence;
						}
					if ($seqOrder eq "clustal" and defined $accPresent{$species})
						{
						$clustalSpeciesCount{$species}++;
						push @accNum, $species if ($clustalSpeciesCount{$species} == 1);
						}
					}
			close ALIGN;
			
			my $inferZero = time;

			# Infer locations of insertions; will be manifested in consensus of "good" sequences as constant "N-N" or "N--N" (not possible to spot in regions not shared with shifted sequences)
				print "\tChecking for insertions ...\n" if ($debug);
				# Quickly infer dirty consensus sequence
					my $consensusSeq;
					for (my $position = 0; $position < length($finalSeq{$AAseqs[0]}); $position++)
						{
						my $gapCount = 0;
						foreach my $accession (@AAseqs)
							{
							last unless (substr($finalSeq{$accession}, $position, 1) eq "-");
							$gapCount++;
							}
	
						if ($gapCount == scalar(@AAseqs))
							{
							$consensusSeq .= "-";
							}
						else
							{
							$consensusSeq .= "N";
							}
						}
					print "\t\tConsensus sequence: $consensusSeq\n" if ($debug);
				
				# Only process if dirty consensus has constant gaps, indicating putative insertions
					if ((my $numConst = ($consensusSeq =~ tr/-//)) > 0)
						{
						print "\n\t\t\tChecking if constant sites match N-N or N--N motifs\n" if ($debug);
						for (my $position = 0; $position < (length($consensusSeq) - 3); $position++)
							{
							next unless (substr($consensusSeq, $position, 3) eq "N-N" or substr($consensusSeq, $position, 4) eq "N--N");
							print "\t\t\t\tMotif found at position $position; putative insertion in\n" if ($debug);
							
							foreach my $shift (@shiftedSeqs)
								{
								if (substr($finalSeq{$shift}, $position + 1, 1) ne "-")	# Check first gap position
									{
									my $location = $position + 2;
									push @{ $insertions{$shift}}, $location;
									$insertionPoint{$location} = 1;
									print "\t\t\t\t\t$nameLabel{$shift} (first position)\n" if ($debug);
									}
								if (substr($finalSeq{$shift}, $position + 2, 1) ne "-" and substr($consensusSeq, $position, 4) eq "N--N")	# Check second gap position if it exists
									{
									my $location = $position + 3;
									push @{ $insertions{$shift}}, $location;
									$insertionPoint{$location} = 1;
									print "\t\t\t\t\t$nameLabel{$shift} (second position)\n" if ($debug);
									}
								}
							}
						}

			# Infer locations of deletions; will be manifested in each "shifted" sequence as a single or double gap compared to raw sequence; use similar procedure as for back-translation
				print "\n\tChecking for deletions ...\n" if ($debug);
				foreach my $accession (@shiftedSeqs)
					{
					next if (length($finalSeq{$accession}) == length($sequence{$accession}));	# No gaps inserted by Clustal; unlikely, but also means no deletions
					if ($debug)
						{
						print "\t\tSequence: $nameLabel{$accession}";
							print " ($accession)" unless ($accession =~ /^tAlign/);
							print "\n";
						}
					my $alignStart = 0;
					my $codingStart;	# Determine start of coding portion for deletion checking
					for (my $alignedPos = 0; $alignedPos < length($sequence{$accession}); $alignedPos++)	# Sequentially read in bps from raw DNA sequence
						{
						my $rawBP = substr($sequence{$accession}, $alignedPos, 1);
						print "\t\t\tPosition $alignedPos (BP: $rawBP)\n" if ($debug);
						
						if ($rawBP eq "-")
							{
							my $gapDNA .= $rawBP;
							
							# Find end of ambiguous stretch in raw DNA
								until ($rawBP ne "-" or $alignedPos > length($sequence{$accession}))
									{
									$alignedPos++;
									my $nextBP = substr($sequence{$accession}, $alignedPos, 1);
									print "\t\t\t\tPosition $alignedPos (BP: $nextBP)\n" if ($debug);
									$gapDNA .= $nextBP if ($nextBP eq "-");
									}
									
								$alignedPos -= 1;
								print "\n\tt\t\tAmbiguous DNA sequence: $gapDNA\n\n" if ($debug);
							
							# Find end of gap in aligned DNA
								my $clustalDNA = substr($finalSeq{$accession}, $alignStart, 1);
									$codingStart = $alignStart if ($clustalDNA ne "-" and not defined $codingStart);
									my $gapClustal .= $clustalDNA;
									
								print "\t\t\t\tAligned DNA: $clustalDNA (position: $alignStart)\n" if ($debug);
								until ($clustalDNA ne "-" or $alignStart == length($finalSeq{$accession}))
									{
									$alignStart++;
									$clustalDNA = substr($finalSeq{$accession}, $alignStart, 1);
										$codingStart = $alignStart if ($clustalDNA ne "-" and not defined $codingStart);
									$gapClustal .= $clustalDNA if ($clustalDNA eq "-");
									print "\t\t\t\tAligned DNA: $clustalDNA (position: $alignStart)\n" if ($debug);
									}
									
								$alignStart--;
								print "\n\t\t\t\tGap sequence: $gapClustal\n" if ($debug);
		#							print "\t\t\t\tEnd of gap ...\n" if ($debug);
							}
						else
							{
							my $clustalDNA = substr($finalSeq{$accession}, $alignStart, 1);
								$codingStart = $alignStart if ($clustalDNA ne "-" and not defined $codingStart);
							my $gapLength;
								print "\t\t\t\tAligned DNA: $clustalDNA (position: $alignStart)\n" if ($debug);
								if ($clustalDNA eq "-")	# Alignment gap inserted by clustal
									{
									$gapLength++;
									until ($clustalDNA ne "-" or $alignStart == length($finalSeq{$accession}))	# Find end of gap
										{
										$alignStart++;
										$clustalDNA = substr($finalSeq{$accession}, $alignStart, 1);
											$codingStart = $alignStart if ($clustalDNA ne "-" and not defined $codingStart);
										$gapLength++ if ($clustalDNA eq "-");
										print "\t\t\t\tAligned DNA: $clustalDNA (position: $alignStart)\n" if ($debug);
										}
									print "\t\t\t\t\tEnd of gap (length: $gapLength) ...\n" if ($debug);
									if ($alignStart - 1 > $codingStart)	# Check for deletions in non-leading gaps
										{
										print "\t\t\t\t\t\tChecking for deletions ...\n" if ($debug);	# Check all gaps; not only those of lengths in multiples of three because other insertions might mask deletions
											my @delPositions;
											for (my $gapPosition = $alignStart - $gapLength; $gapPosition < $alignStart; $gapPosition++)	# Store potential deletion positions
												{
												next if ($gapPosition < $codingStart);
												next if (defined $insertionPoint{$gapPosition});	# Deletion actually due to an insertion in anopther shifted sequence
												next if (substr($consensusSeq, $gapPosition, 1) eq "-");
												
												print "\t\t\t\t\t\t\tStoring position $gapPosition\n" if ($debug);
												push @delPositions, $gapPosition + 1;
												}
											if (scalar(@delPositions) > 1)	# Store those deletion positions that do not occur in groups of multiples of three
												{
												my $count = 1;
												while ($count < scalar(@delPositions))
													{
													if ($delPositions[$count] - $delPositions[$count-1] == 1)
														{
														my $delLength = 1;
														my $lower = $count - 1;
														while ($count < scalar(@delPositions) and $delPositions[$count] - $delPositions[$count-1] == 1)
															{
															$delLength++;
															$count++;
															}
														push @{ $deletions{$accession}}, $delPositions[$lower]."-".$delPositions[$count - 1] unless (int($delLength / 3) * 3 == $delLength);	# Store those deletion positions that do not occur in groups of multiples of three
														}
													else
														{
														push @{ $deletions{$accession}}, $delPositions[$count-1];
														$count++;
														}
													}
												}
											else
												{
												push @{ $deletions{$accession}}, $delPositions[0] if (defined $delPositions[0]);
												}
										}
#									if ($gapLength == 1 and not defined $insertionPoint{$alignStart})	# Insertions in other shifted sequences will appear as putative single deletion
#										{
#										print "\t\t\t\t\tPutative deletion\n" if ($debug);
#										push @{ $deletions{$accession}}, $alignStart;
#										}
#									if ($gapLength == 2)
#										{
#										my $startStart = $alignStart - 1;
#										print "\t\t\t\t\tPutative deletion\n" if ($debug);
#										if (defined $insertionPoint{$alignStart})	# Insertions in other shifted sequences will make single deletion appear as putative multiple deletion
#											{
#											push @{ $deletions{$accession}}, $startStart;
#											}
#										else
#											{
#											push @{ $deletions{$accession}}, "$startStart-$alignStart";
#											}
#										}
									}
							}
						$alignStart++;
						}
					}
			
			# Check for frame shifting indels in any unique leading stretches in "shifted" sequences
				print "\n\tChecking for leading indels ...\n" if ($debug);

				# Determine number of gaps added to 5' end of "good" sequences
					my $rawSequence = $framedSeq{$AAseqs[0]};
						my $rawGapCount = 0;
						while ($rawSequence =~ s/^\-//)
							{
							$rawGapCount++;
							}
					my $alignedSequence = $finalSeq{$AAseqs[0]};
						my $alignedGapCount = 0;
						while ($alignedSequence =~ s/^\-//)
							{
							$alignedGapCount++;
							}

					my $addedGaps = $alignedGapCount - $rawGapCount;
						if (int($addedGaps/3) * 3 != $addedGaps)	# Number of added gaps not a multiple of three; potential frame shift
							{
							print "\t\tPotential frame shifting indel detected; checking shifted sequences\n" if ($debug);
							foreach my $shift (@shiftedSeqs)
								{
								if ($debug)
									{
									print "\t\t\tSequence: $nameLabel{$shift}";
										print " ($shift)" unless ($shift =~ /^tAlign/);
										print "\n";
									}
								unless ($finalSeq{$shift} =~ /^\-/)
									{
									print "\t\t\t\tPutative leading indel\n" if ($debug);
									push @{ $indel{$shift}}, "1-$addedGaps";
									}
								}
							}

			# Print out positions of putative indels
				print "\n\tLocation of putative indels causing frame shifts:\n";
				foreach my $shift (@shiftedSeqs)
					{
					print "\t\tSequence: $nameLabel{$shift}";
						print " ($shift)" unless ($shift =~ /^tAlign/);
						print "\n";
					
					print "\t\t\tInsertions: ";
						if (defined @{ $insertions{$shift}})
							{
							print join(" ", @{ $insertions{$shift}})
							}
						else
							{
							print "none";
							}
						print "\n";
					print "\t\t\tDeletions: ";
						if (defined @{ $deletions{$shift}})
							{
							print join(" ", @{ $deletions{$shift}})
							}
						else
							{
							print "none";
							}
						print "\n";
					print "\t\t\tUnspecified indels: ";
						if (defined @{ $indel{$shift}})
							{
							print join(" ", @{ $indel{$shift}})
							}
						else
							{
							print "none";
							}
						print "\n";
					}
			
			printf "\n\tTime taken to infer locations of frame shifts: %s seconds\n", time - $inferZero;
		}

# Print results!
	my $ntax = $seqCount - $frameDelCount - $nDelCount - $badAlignCount;
	@accNum = sort { $nameLabel{$a} cmp $nameLabel{$b} } keys %nameLabel if ($seqOrder eq "alphabetical");

	print "\nPrinting results ...\n";
		seqPrint();

exit(0);

### Subroutines used in the program

sub setLineBreak	# Check line breaks of input files and set input record separator accordingly
	{
	my $inFile = shift;
	$/ ="\n";
	open (IN, "<$inFile") or die "Cannot open $inFile to check form of line breaks.\n";
		while (<IN>)
			{
			if ($_ =~ /\r\n/)
				{
				print "\tDOS line breaks detected ...\n" if ($verbose);
				$/ ="\r\n";
				last;
				}
			elsif ($_ =~ /\r/)
				{
				print "\tMac line breaks detected ...\n" if ($verbose);
				$/ ="\r";
				last;
				}
			else
				{
				print "\tUnix line breaks detected ...\n" if ($verbose);
				$/ ="\n";
				last;
				}
			}
	close IN;
	}

sub geneticCoder	# Create translation tables for all genetic codes
	{
	my %geneticCode = ('1' => 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
					   '2' => 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG',
					   '3' => 'FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
					   '4' => 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
					   '5' => 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG',
					   '6' => 'FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
					   '9' => 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
					   '10' => 'FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
					   '11' => 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
					   '12' => 'FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
					   '13' => 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG',
					   '14' => 'FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
					   '15' => 'FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
					   '16' => 'FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
					   '21' => 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
					   '22' => 'FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
					   '23' => 'FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG');
	
	foreach my $code (qw(1 2 3 4 5 6 9 10 11 12 13 14 15 16 21 22 23))
		{
		# Establish basic translation table for each genetic code
#			print "\nEstablishing \"$transTable{$code}\" genetic code ...\n" if ($debug);
			my $position = 0;
			foreach my $base1 (qw (T C A G))
				{
				foreach my $base2 (qw (T C A G))
					{
					foreach my $base3 (qw (T C A G))
						{
						my $codon = $base1.$base2.$base3;
						$DNAtoAA{$code}{$codon} = substr($geneticCode{$code}, $position, 1);
#							print "\t$codon = $DNAtoAA{$code}{$codon}\n" if ($debug);
						$position++;
						}
					}
				}
	
		# Extend translation table to account for ambiguity codes (note: does not account for gaps)
#			print "\nExtending translation table to account for ambiguity codes ...\n" if ($debug);
			foreach my $firstPos (@ambigList)
				{
				foreach my $secondPos (@ambigList)
					{
					foreach my $thirdPos (@ambigList)
						{
						my $codon = $firstPos.$secondPos.$thirdPos;
						next if (defined $DNAtoAA{$code}{$codon});
						my $refAA = "";
						foreach my $firstNT (@ {$constitNTlist{$firstPos} })
							{
							last if (defined $DNAtoAA{$code}{$codon});
							foreach my $secondNT (@ {$constitNTlist{$secondPos} })
								{
								last if (defined $DNAtoAA{$code}{$codon});
								foreach my $thirdNT (@ {$constitNTlist{$thirdPos} })
									{
									my $testCodon = $firstNT.$secondNT.$thirdNT;
									if (not $refAA)
										{
										$refAA = $DNAtoAA{$code}{$testCodon};
										}
									else
										{
										if ($DNAtoAA{$code}{$testCodon} ne $refAA)
											{
											$DNAtoAA{$code}{$codon} = "?";
											last;
											}
										}
									}
								}
							}
						$DNAtoAA{$code}{$codon} = $refAA if (not defined $DNAtoAA{$code}{$codon});
#						print "\t$codon = $DNAtoAA{$code}{$codon}\n" if ($debug);
						}
					}
				}
		}
	return;
	}
	
sub translate	# Translate a DNA sequence to an AA sequence (note: does not account for gaps)
	{
	my $DNAseq = shift;
	my $userCode = shift;
	
	my $protSeq;
	for (my $codonStart = 0; $codonStart < length($DNAseq); $codonStart += 3)
		{
		if (length($DNAseq) - $codonStart >= 3)	# Codon is complete; translate
			{
			my $codon = substr($DNAseq, $codonStart, 3);
			if ($codon =~ /-/ or $codon =~ /\./)
				{
				$protSeq .= "?";
				}
			else
				{
				$protSeq .= $DNAtoAA{$userCode}{$codon};
				}
			}
		else	# Incomplete codon; automatically translates as ?
			{
			$protSeq .= "?";
			}
		}

	return $protSeq;
	}

sub seqRead
	{
	my $seqFile = shift;

	print "\nReading in sequence data from file $seqFile (type is $inputType) ...\n" if ($inputType);
	setLineBreak($seqFile);
	open (SEQ, "<$seqFile") or die "Cannot open file containing sequences, $seqFile\n";
		my ($header, $tempAcc, $tempName, $tempSeq);
		my $fastaAcc;
		my (%nexusSpecies, %nexusAcc, $nexusRead);
		my ($phylipLineCount, $phylipTaxa, $phylipChars, %phylipSeq);
		my $sealCode;
		my ($sealDelFlag, $owner) = (0, 0);

		while (<SEQ>)
			{
			chomp;
			my $lineRead = $_;
			next unless ($lineRead);
			
			# Autodetect sequence format
				if (not $inputType)
					{
					$inputType = "fasta" if ($lineRead =~ /^>/);
					$inputType = "nexus" if ($lineRead =~ /\#nexus/i);
					$inputType = "phylip" if ($lineRead =~ /^\s*\d+\s+\d+/);
					$inputType = "Se-Al" if ($lineRead =~ /^\s*Database=\{/i);
					print "\nReading in sequence data from file $seqFile (type determined to be $inputType) ...\n" if ($inputType);
					}
			
			if ($inputType eq "nexus")
				{
				# Only read in data lines
					if ($lineRead =~ /^\s*matrix/i)
						{
						$nexusRead = 1;
						next;
						}
					$nexusRead = 0 if ($lineRead =~ /;\s*$/);
					next unless ($nexusRead);
					next unless ($lineRead =~ /a/i or $lineRead =~ /c/i or $lineRead =~ /g/i or $lineRead =~ /t/i);
				# Clean up input line
					$lineRead =~ s/^\s+//;
					$lineRead =~ s/\'//g;
				my ($species, $seq) = split(/\s+/, $lineRead);
					$species =~ s/\s+/_/g;
				if (not defined $nexusSpecies{$species})
					{
					$nexusSpecies{$species} = 1;
					$seqCount++;
					$nexusAcc{$species} = "tAlign_".$seqCount;
					push @accNum, $nexusAcc{$species};
						$nameLabel{$nexusAcc{$species}} = $species;
						$sequence{$nexusAcc{$species}} = uc($seq);
						$geneticCode{$nexusAcc{$species}} = $globalGenCode;
					}
				else	# Sequences are in interleaved format; append sequence
					{
					$sequence{$nexusAcc{$species}} .= uc($seq);
					}
				}

			if ($inputType eq "fasta")
				{
				if ($lineRead =~/^\s*>/)
					{
					my $species;
					$seqCount++;
					(my $tempSpecies = $lineRead) =~ s/^\s*>//;
					
						if ($tempSpecies =~ /^Mit\.\s+/)	# Entry comes from European RNA project
							{
							$tempSpecies =~ s/^Mit\.\s+//i;	# To fix entries from European RNA project
							my @speciesInfo = split(/\s+/, $tempSpecies);
								$species = join('_', $speciesInfo[0], $speciesInfo[1]);
							if (defined $speciesInfo[2])
								{
								$fastaAcc = $speciesInfo[2];
								}
							else
								{
								$fastaAcc = "tAlign_".$seqCount;
								}
							}
						else
							{
							my @speciesLine = split(/\s+/, $tempSpecies);
							if ($speciesLine[$#speciesLine] =~ /^\(?[A-Z]+\d+\)?$/ and scalar(@speciesLine) > 1)	# Check whether last entry is an accession number
								{
								$fastaAcc = pop (@speciesLine);
								$fastaAcc =~ s/^\(//g;
								$fastaAcc =~ s/\)$//g;
								}
							else
								{
								$fastaAcc = "tAlign_".$seqCount;
								}
							$species = join('_', @speciesLine);
								$species = "Sequence_".$seqCount if ($species eq "");
							}
					push @accNum, $fastaAcc;
						$geneticCode{$fastaAcc} = $globalGenCode;
					$nameLabel{$fastaAcc} = $species;
					}
				else
					{
					$sequence{$fastaAcc} .= uc($lineRead);
					}
				}

			if ($inputType eq "Se-Al")
				{
				my $header;
				$sealDelFlag = 1 if ($lineRead =~/MCoL/);	# Se-Al sometimes places deleted species at end of file; do not read in remainder of file
					next if ($sealDelFlag == 1);
				next unless ($lineRead =~/NumSites/i or $lineRead =~/Owner/i or $lineRead =~/Name/i or $lineRead =~/Accession/i or $lineRead =~/Sequence/i or $lineRead =~/GeneticCode/i);
				if ($lineRead =~/Owner\s*\=\s*(\d+)/i)
					{
					$owner = $1;
					}
				if ($lineRead =~/Accession/i and $owner == 2)
					{
					$seqCount++;
					if ($lineRead =~ /null/ or $lineRead =~ /\"\"/)
						{
						$tempAcc = "tAlign_".$seqCount;
						}
					else
						{
						($header, $tempAcc) = split (/=/, $lineRead);
							$tempAcc =~ s/\"//g;
							$tempAcc =~ s/;//g;
						}
					push @accNum, $tempAcc;
					}
				if ($lineRead =~/Name/i and $owner == 2)
					{
					($header, $tempName) = split (/=/, $lineRead);
						$tempName =~ s/\"//g;
						$tempName =~ s/\s*;//g;
					}
				if ($lineRead =~/GeneticCode/i and $owner == 2)
					{
					($header, $sealCode) = split (/=/, $lineRead);
						$sealCode =~ s/\"//g;
						$sealCode =~ s/\s*;//g;
						$geneticCode{$tempAcc} = $sealCode + 1;
					}
				if ($lineRead =~/Sequence/i and $owner == 2)
					{
					($header, $tempSeq) = split (/=/, $lineRead);
						$tempSeq =~ s/\"//g;
						$tempSeq =~ s/;//g;
					$nameLabel{$tempAcc} = $tempName;
					$sequence{$tempAcc} = uc($tempSeq);
					}
				}

			if ($inputType eq "phylip")
				{
				if ($lineRead =~ /^\s*(\d+)\s+(\d+)/)
					{
					$phylipTaxa = $1;
					$phylipChars = $2;
					$phylipLineCount = 0;
					}
				else
					{
					$phylipLineCount++;
					
					$lineRead =~ s/\s//g;
					
					$phylipSeq{$phylipLineCount} .= $lineRead;
					
					$phylipLineCount = 0 if ($phylipLineCount == $phylipTaxa);
					}
				}
			}
	close SEQ;
	
	if ($inputType eq "phylip")	# Postprocess input to derive taxon names and sequence; accounts for both sequential and extended formatting
		{
		for (my $i = 1; $i <= $phylipTaxa; $i++)
			{
			my $phylipAcc = "tAlign_" . $i;
			
			push @accNum, $phylipAcc;
			$geneticCode{$phylipAcc} = $globalGenCode;
			
			# Derive taxon name and sequence
				$sequence{$phylipAcc} = uc(substr($phylipSeq{$i}, 0 - $phylipChars));
				$nameLabel{$phylipAcc} = substr($phylipSeq{$i}, 0, length($phylipSeq{$i}) - $phylipChars);
					$nameLabel{$phylipAcc} =~ s/\s+//g;
			}
		}
	}
	
sub seqPrint
	{
	# Print fasta-formatted file (always)
		if ($fastaPrint)
			{
			print "\tWriting to fasta-formatted file $fastaOut ...\n";
			open (FASTA, ">$fastaOut") or die "Cannot open fasta file for aligned DNA sequences, $fastaOut";
				foreach my $entry (@accNum)
					{
					next if ($deletedSeq{$entry});
					my $fastaSeq = $finalSeq{$entry};
						my $breakPoint = 79;
						until ($breakPoint > length($fastaSeq))
							{
							my $replaceString = "\n" . substr($fastaSeq, $breakPoint, 1);
							substr($fastaSeq, $breakPoint, 1) = $replaceString;
							$breakPoint += 80;
							}
					print FASTA ">$nameLabel{$entry}";
						print FASTA "\t($entry)" unless ($entry =~ /^tAlign/);
					print FASTA "\n$fastaSeq\n";
					}
			close FASTA;
			}

	# Print nexus-formatted file (on demand)
		if ($nexusPrint)
			{
			print "\tWriting to nexus file $nexusOut ...\n";
			open (NEX, ">$nexusOut") or die "Cannot open nexus file for aligned DNA sequences, $nexusOut";
				print NEX "#nexus\n\n";
				print NEX "[File created from $seqFile using transAlign.pl v1.1 on ".localtime()."]\n\n";
				print NEX "begin data;\n";
				print NEX "\tdimensions ntax = $ntax nchar = $maxLength;\n";
				print NEX "\tformat datatype = DNA gap = - missing = ?;\n\n";
				print NEX "\tmatrix\n\n";
				foreach my $entry (@accNum)
					{
					next if ($deletedSeq{$entry});
					if ($nameLabel{$entry} =~ /\W/)
						{
						print NEX "'$nameLabel{$entry}'";
						}
					else
						{
						print NEX "$nameLabel{$entry}";
						}
					print NEX "\t$finalSeq{$entry}\n";
					}
				print NEX "\t;\nend;\n";
			close NEX;
			}

	# Print phylip-formatted file (on demand)
		if ($phylipTradPrint or $phylipExtPrint)
			{
			my $maxTaxLength = 50;
				$maxTaxLength = 10 if ($phylipTradPrint);
			my %shortNameCount;	
				
			print "\tWriting to phylip file $phylipOut ...\n";
			open (PHYLIP, ">$phylipOut") or die "Cannot open phylip file for aligned DNA sequences, $phylipOut";
				print PHYLIP "\t$ntax\t$maxLength\n";
				foreach my $entry (@accNum)
					{
					next if ($deletedSeq{$entry});
					
					my $phylipName = $nameLabel{$entry};

					# Check name label and adjust to proper length if needed
						if (length($phylipName) < $maxTaxLength)
							{
							$shortNameCount{$phylipName}++;
							$phylipName .= " " x ($maxTaxLength - length($phylipName));	# Pad end of name with spaces as needed
							}
						else
							{
							my $trimmedName = substr($phylipName, 0 , $maxTaxLength);
							$shortNameCount{$trimmedName}++;
							if ($shortNameCount{$trimmedName} > 1)	# Check for duplicates among shortened names and make unique by adding numbers
								{
								$phylipName = substr($phylipName, 0, $maxTaxLength - length($shortNameCount{$trimmedName}));
									$phylipName .= $shortNameCount{$trimmedName};
									$phylipName .= " " x ($maxTaxLength - length($phylipName));	# Pad end of name with spaces as needed
								}
							else
								{
								$phylipName = $trimmedName;
								}
							}
						
					print PHYLIP "$phylipName";
						print PHYLIP " " if ($phylipExtPrint);
					print PHYLIP "$finalSeq{$entry}\n";
					}
			close PHYLIP;
			}

	# Print Se-Al-formatted file (on demand)
		if ($sealPrint)
			{
			print "\tWriting to Se_Al file $sealOut ...\n";
			open (SEAL, ">$sealOut") or die "Cannot open Se-Al file for aligned DNA sequences, $sealOut\n";
				print SEAL "Database={\n";
				print SEAL "\tID='MLst';\n";
				print SEAL "\tOwner=null;\n";
				print SEAL "\tName=null;\n";
				print SEAL "\tDescription=null;\n";
				print SEAL "\tFlags=0;\n";
				print SEAL "\tCount=2;\n";
				print SEAL "\t{\n\t\t{\n";
				
				print SEAL "\t\t\tID='PAli';\n";
				print SEAL "\t\t\tOwner=1;\n";
				print SEAL "\t\t\tName=\"$seqFile\";\n";
				print SEAL "\t\t\tDescription=null;\n";
				print SEAL "\t\t\tFlags=0;\n";
				print SEAL "\t\t\tNumSites=$maxLength;\n";
				print SEAL "\t\t\tType=\"Nucleotide\";\n";
				print SEAL "\t\t\tFeatures=null;\n";
				print SEAL "\t\t\tColourMode=1;\n";
				print SEAL "\t\t\tLabelMode=0;\n";
				print SEAL "\t\t\ttriplets=false;\n";
				print SEAL "\t\t\tinverse=true;\n";
				print SEAL "\t\t\tCount=$ntax;\n";
				print SEAL "\t\t\t{\n";
				
				my $i = 0;
				foreach my $sequence (@accNum)
					{
					next if ($deletedSeq{$sequence});
					$i++;
					print SEAL "\t\t\t\t{\n";
					print SEAL "\t\t\t\t\tID='PSeq';\n";
					print SEAL "\t\t\t\t\tOwner=2;\n";
					print SEAL "\t\t\t\t\tName=\"$nameLabel{$sequence}\";\n";
					print SEAL "\t\t\t\t\tDescription=null;\n";
					print SEAL "\t\t\t\t\tFlags=0;\n";
					print SEAL "\t\t\t\t\tAccession=";
						if ($sequence =~/^tAlign_/)
							{
							print SEAL "null;\n";
							}
						else
							{
							print SEAL "$sequence;\n";
							}
							
					print SEAL "\t\t\t\t\tType=\"DNA\";\n";
					print SEAL "\t\t\t\t\tLength=".length($finalSeq{$sequence}).";\n";
					print SEAL "\t\t\t\t\tSequence=\"$finalSeq{$sequence}\";\n";
					my $sealCode = $geneticCode{$sequence} - 1;
					print SEAL "\t\t\t\t\tGeneticCode=$sealCode;\n";
					print SEAL "\t\t\t\t\tCodeTable=null;\n";
					print SEAL "\t\t\t\t\tFrame=1;\n";
					print SEAL "\t\t\t\t\tFeatures=null;\n";
					print SEAL "\t\t\t\t\tParent=null;\n";
					print SEAL "\t\t\t\t\tComplemented=false;\n";
					print SEAL "\t\t\t\t\tReversed=false;\n";
					print SEAL "\t\t\t\t}";
					print SEAL "," unless ($i == $ntax);
					print SEAL "\n";
					}
				
				print SEAL "\t\t\t};\n";
				print SEAL "\t\t},\n";
				print SEAL "\t\t{\n";
				print SEAL "\t\t\tID='MCoL';\n";
				print SEAL "\t\t\tOwner=1;\n";
				print SEAL "\t\t\tName=\"Genetic Codes\";\n";
				print SEAL "\t\t\tDescription=\"Custom Genetic Codes\";\n";
				print SEAL "\t\t\tFlags=0;\n";
				print SEAL "\t\t\tCount=0;\n";
				print SEAL "\t\t}\n";
				print SEAL "\t};\n";
				print SEAL "};\n";
			close SEAL;
			}
	}

sub _subtprob	# Calculates upper probability from Student's t-distribution; taken from the CPAN Perl module Statistics::Distributions by Michael Kospach
	{
	my ($n, $x) = @_;

	my ($a,$b);
	my $w = atan2($x / sqrt($n), 1);
	my $z = cos($w) ** 2;
	my $y = 1;

	for (my $i = $n-2; $i >= 2; $i -= 2)
		{
		$y = 1 + ($i-1) / $i * $z * $y;
		} 

	if ($n % 2 == 0)
		{
		$a = sin($w)/2;
		$b = .5;
		}
	else
		{
		$a = ($n == 1) ? 0 : sin($w)*cos($w)/PI;
		$b= .5 + $w/PI;
		}
	return max(0, 1 - $b - $a * $y);
	}
	
sub max	# Finds maximum value in an array; taken from the CPAN Perl module Statistics::Distributions by Michael Kospach
	{
	my $max = shift;
	my $next;
	while (@_)
		{
		$next = shift;
		$max = $next if ($next > $max);
		}	
	return $max;
	}

# Version history
#
#	v1.2 (June 22, 2005)
#		- major bug fix to calculation of p-values for Student's t-distribution (no longer uses gamma function and Lanczos approximation of it)
#		- minor bug fix to seqRead subroutine
#		- uses File::Basename module to parse filenames
#
#	v1.1 (April 12, 2005)
#		- added features
#			- added ability to automatically delete pooorly aligning sequences as determined by
#			  a one-tailed two-sample t-test of ClustalW alignment scores
#			- added ability for stop codon threshold to be a percentage
#			- added user option regarding stripping of gaps in raw sequences
#			- added user option as to maximum percentage of Ns allowed in a sequence
#			- added user option regarding output order of sequences
#			- added phylip format to input and output files
#			- added autodetection of recognized input formats
#		- forced ClustalW to output alignments to a standard name because some versions
#		  will not output clustal format (thanks to Bjorn Ostman for pointing this out)
#		- ClustalW output can now be in GCG and PHYLIP formats as well
#		- more robust nexus parser
#		- sequence reading and printing shunted to subroutines
#		- more efficient checking whether or not files exist in user-input mode
#		- (other) minor bug fixes
#
#	v1.0 (September 24, 2004)
#		- initial release
