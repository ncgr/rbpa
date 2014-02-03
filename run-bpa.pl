#!/usr/bin/env perl

use strict;

use vars qw/%ops/;

&use_get_opts;

&help_main if !$ops{m};

my ($output, $mode) = ("BPA", $ops{m});

if ($mode =~ /unitigs/i) {
	&help_unitigs if $ops{h};
	&run_unitigs;
} elsif ($mode =~ /olc/i) {
	&help_OLC if $ops{h};
	&run_OLC;
} elsif ($mode =~ /scaffold/i) {
	&help_scaffold if $ops{h};
	&run_scaffold;
} elsif ($mode =~ /Gapclose/i) {
	&help_gapclose if $ops{h};
	&run_gapclose;
} elsif ($mode =~ /Realign/i) {
	&help_realign if $ops{h};
	&run_realign;
} elsif ($mode =~ /BLAST/i) {
	&help_BLAST if $ops{h};
	&run_BLAST;
} elsif ($mode =~ /ESTscan/i) {
	&help_ESTscan if $ops{h};
	&run_ESTscan;
} elsif ($mode =~ /HMMscan/i) {
	&help_HMMscan if $ops{h};
	&run_HMMscan;
} elsif ($mode =~ /assemble/i) {
	&help_assemble if $ops{h};
	&run_OLC;
	&run_scaffold;
	&run_gapclose;
	&run_realign;
} elsif ($mode =~ /annotate/i) {
	&help_annotate if $ops{h};
	&run_BLAST;
	&help_ESTscan;
	&help_HMMscan;
} else {
	warn "\nERROR:$mode is not an option\n";
	&help_main;
}

sub run_unitigs {
	die "qsub command did not execute properly: exit value=$?\n" if system("qsub -help >/dev/null") != 0;
	die "mpirun did not execute properly: exit value=$?\n" if system("mpirun -h 2>/dev/null") != 0;
	die "ABYSS-P did not execture properly, please ensure that you exported the binary directory into your PATH variable\n" if system("ABYSS-P --help >/dev/null 2>/dev/null") != 0;
	die "Please specify at least one fastq read file with -i\n" if !$ops{i};
	die "Please specify the number of cores to use on your cluster with -c\n" if !$ops{c};
	die "Please provide a file with kmer values, one per line with -k\nexamples are shown in the Readme and the usage\n" if !$ops{k};
	my ($reads, $cores, $kmers, $coverage) = ($ops{i}, $ops{c}, $ops{k}, ($ops{v}) ? $ops{v} : "");
	system ("run-abyss.bash $kmers \"$reads\" $cores $coverage");
	die "run-abyss.bash failed, check logs or STDERR captures for information\n" if $? != 0;
}

sub run_OLC {
	die "Cap3 did not execture properly, please ensure that you exported the binary directory into your PATH variable\n" if system("cap3 2>/dev/null") != 256;
	die "Please specify at least one fasta file with -i\n" if !$ops{i};
	my $threads = ($ops{t}) ? $ops{t} : 0;
	my ($fasta, $cap3, $cdhit) = ($ops{i}, ($ops{u}) ? "$ops{u}" : "-o100 -h50", ($ops{e}) ? "-T $threads -c $ops{e}" : "");
	system ("run-OLC.bash -f \"$fasta\" -c \"$cap3\" -d \"$cdhit\"") if $ops{e};
	system ("run-OLC.bash -f \"$fasta\" -c \"$cap3\"") unless $ops{e};
	die "run-OLC.bash failed, check logs or STDERR captures for more information\n" if $? != 0;
}

sub run_scaffold {
	die "abyss did not execute properly, please ensure that you exported the binary directory into your PATH variable\n" if system("abyss-scaffold --help >/dev/null") != 0;
	die "bwa did not execute properly, please ensure that you exported the binary directory into your PATH variable\n" if system("bwa aln 2>/dev/null") != 0;
	die "Please specify at least one set of paired fastq files -r, or a bam file with b\n" if !$ops{r} || !$ops{b};
	die "Please specify at least one fasta file with -R, if -b provided reference must be alignment reference\n" if !$ops{R};
	die "Please specify a kmer value to use with abyss -k\n" if !$ops{k};
	my $params = "-r \"$ops{r}\" -R $ops{R} -k $ops{k}";
	$params=$params." -m $ops{i}" if $ops{i};
	$params=$params." -p $ops{p}" if $ops{p};
	$params=$params." -o $ops{o}" if $ops{o};
	$params=$params." -s $ops{s}" if $ops{s};
	$params=$params." -t $ops{t}" if $ops{t};
	system("run-scaffold.bash $params");
	die "run-scaffold.bash failed, check logs or STDERR captures for more information\n" if $? != 0;
}

sub run_gapclose {
	die "GapCloser did not execute properly, please ensure correct export in PATH variable\n" if system("GapCloser -h >/dev/null") != 256;
	die "Please specify at least one file for mate one of paired fastq files, forward -f\n" if !$ops{f};
	die "Please specify at least one file for mate two of paired fastq files, reverse -r\n" if !$ops{r};
	die "Please provide reference fasta file -R\n" if !$ops{R};
	die "Please specify a maximum read length -l\n" if !$ops{l};
	my ($r1, $r2, $ref, $ins, $mlen, $threads) = ($ops{f}, $ops{r}, $ops{R}, ($ops{i}) ? $ops{i} : "", $ops{l}, ($ops{t}) ? $ops{t} : 1);
	system("run-gapclose.bash \"$r1\" \"$r2\" $ref $mlen $threads $ins");
	die "run-gapclose.bash failed, check logs or STDERR captures for more information\n" if $? != 0;
}

sub run_realign {
	die "samtools did not execute properly, please ensure that you exported the binary directory into your PATH variable\n" if system("samtools view 2>/dev/null") != 256;
	die "bwa did not execute properly, please ensure that you exported the binary directory into your PATH variable\n" if system("bwa aln 2>/dev/null") != 0;
	die "Please specify at least one fastq file to realign -f\n" if !$ops{f};
#       die "Please specify at least one file for mate two of paired fastq files, reverse -r\n" if !$ops{r};
	die "Please specify fasta input file with -R\n" if !$ops{R};

	my ($r1, $ref, $threads) = ($ops{f}, $ops{R}, ($ops{t}) ? $ops{t} : 1);
	my $mism = ($ops{m}) ? $ops{m} : "";
	$output = $ops{p} if $ops{p};
	system("run-bwa-index.bash $ref");
	die "index did not complete properly, please see logs or stderr\n" if $? != 0;
	if ($ops{r}){
		my $r2 = $ops{r};
		system("run-bwa-paired.bash \"$r1\" \"$r2\" $ref $output $threads $mism");
	} else {
		system("run-bwa-single.bash \"$r1\" $ref $output $threads $mism");
	}
	die "realign did not finish properly, please see logs or stderr\n" if $? != 0;
}

sub run_BLAST {
	die "blast did not execute properly, please ensure that you exported the binary directory into your PATH variable\n" if system("blastn -h >/dev/null") != 256;
#	die "Please specify a fasta input file with -i or a directory for cluster searching with -d\n" if !$ops{i};
	die "Please specify a fasta input file with -i\n" if !$ops{i};
	die "Please specify BLAST formatted reference basename -R\n" if !$ops{R};
	my ($clust, $query, $ref, $out, $threads, $shell, $evalue, $outfmt, $gencode, $que) = (($ops{c}) ? $ops{c} : 0, $ops{i}, $ops{R}, ($ops{o}) ? $ops{o} : $output, ($ops{t}) ? $ops{t} : 1, ($ops{s}) ? $ops{s} : "blastx.bash", ($ops{e}) ? $ops{e} : "1e-5", ($ops{f}) ? $ops{f} : 7, ($ops{g}) ? $ops{g} : 1, ($ops{q}) ? $ops{q} : "");
	if ($clust == 0) {
		system("$shell $query $ref $out $threads $evalue $outfmt $gencode");
	#} else {
	#	die "qsub command did not execute properly: exit value=$?\n" if system("qsub -help >/dev/null") != 0;
	#	die "Please specify a directory for cluster searching with -d, full path required\n" if !$ops{d};
	#	my $query_dir = $ops{d};
	#	system("launch-cluster-BLAST.bash run-BLAST.bash $shell $query_dir $ref $out $threads $evalue $outfmt $gencode $que");
	}
	die "blast did not execute properly, please see logs or stderr\n" if $? != 0;

}

sub run_ESTscan {
	die "ESTscan did not execute properly, please ensure that you exported the binary directory into your PATH variable\n" if system("estscan -h 2>/dev/null") != 256;
	die "Please specify a fasta input file with -i\n" if !$ops{i};
	die "Please specify a scoring matrix with -m\n" if !$ops{m};
	my ($in, $smat, $out, $minlen) = ($ops{i}, $ops{m}, ($ops{o}) ? $ops{o} : $output, ($ops{l}) ? $ops{l} : "");
	system("run-ESTscan.bash $in $smat $out $minlen");
	die "run-ESTscan.bash did not execute properly, please see logs of stderr\n" if $? != 0;
}

sub run_HMMscan {
	die "hmmscan did not execute properly, please ensure that you exported the binary directory into your PATH variable\n" if system("hmmscan -h >/dev/null") != 0;
	die "Please specify a protein fasta input file with -i\n" if !$ops{i};
	die "Please specify a formatted database to search against using hmmer3 -d\n" if !$ops{d};
	my ($in, $db, $out, $threads) = ($ops{i}, $ops{d}, ($ops{o}) ? $ops{o} : $output, ($ops{t}) ? $ops{t} : 1);
	system("run-hmmer.bash $in $db $out $threads");
	die "run-hmmer.bash did not execute properly, please see logs of stderr\n" if $? != 0;
}

sub help_main {

print STDERR << "_HELP";

run-bpa.pl:

This wrapper controls a BPA pipeline run.  Please see README for more 
information.
 
        usage:$0 <-m mode>

	mode    : Unitigs, OLC, Scaffold, Gapclose, Realign, BLAST, ESTscan, HMMscan

	Program mode details:

_HELP
	&help_unitigs;
	&help_OLC;
	&help_scaffold;
	&help_gapclose;
	&help_realign;
	&help_BLAST;
	&help_ESTscan;
	&help_HMMscan;
        exit 1;
}

sub help_unitigs {

print STDERR << "_HELP";

Unitigs:

	Usage: $0 -m $ops{m} <-i> <-c> <-k> [-v]

	-i) fastq input list. If multiple fastq files are provided please utilize quotations (can use wildcard).
	-c) number of cores to use on SGE cluster. This implies PE ORTE
	-k) kmer list, text file with one line per value
	-v) minimum kmer coverage [5]

_HELP

}

sub help_OLC {

print STDERR << "_HELP";

OLC:

	Usage: $0 -m $ops{m} <-i "fasta(s)"> [-t threads] [-u "cap3 options"] [-e <= 1.0]

	-i) fasta input file. Multiple files may be specified with quotations
	-t) threads to be used in cd-hit-est
	-u) CAP3 command [-o 100]. Use quotations for multiple options.
	-e) percent identity for cd-hit-est. If not provided reduction reduction will not be performed

_HELP

}
sub help_scaffold {

print STDERR << "_HELP";

Scaffold:

	Usage: $0 -m $ops{m} <-r "read file(s)"> <-R reference.fa> <-k kmer> [-m maximum mismatch bwa] [-p min pairs] [-o output prefix] [-s seed length] [-t threads]

	-r) fastq read files. Expected orientation is forward reverse.
	-R) fasta reference to scaffold. Sequences will be renamed sequentially
	-k) kmer value to use in abyss.
	-m) maximum mismatch for Burrows . Wheeler aligner
	-p) minimum pairing evidence to use for scaffolds
	-o) output file prefix
	-s) minimum scaffold seed length
	-t) threads

_HELP

}

sub help_gapclose {

print STDERR << "_HELP";

Gapclose:

	Usage: $0 -m $ops{m} <-f forward reads> <-r reverse reads> <-R fasta file to gap close> <-l maximum read length> [-i insert size] [-t threads]

	-f) mate one of fastq paired reads
	-r) mate two of fastq paired reads
	-R) fasta reference to gap close
	-l) maximum fastq read length
	-i) insert size. Must be provided if scaffolding directory is not contained in parent.
	-t) threads

_HELP

}

sub help_realign {

print STDERR << "_HELP";

Realign:

	Usage: $0 -m $ops{m} <-f forward reads> [-r reverse reads] <-R reference> [-o output prefix] [-t threads] [-m mismatch]

	-f) mate one of fastq read files. Multiple files may be included with quotations order must match mate two
	-r) mate two of fastq read files. Multiple files may be included with quotations order must match mate one
	-R) fasta reference file
	-o) output file prefix
	-t) threads
	-m) maximum mismatch in Burrows . Wheeler aligner

_HELP

}

sub help_BLAST {

print STDERR << "_HELP";

BLAST:

	Usage: $0 -m $ops{m} <-i input file> <-R blast formatted reference> [-o output prefix] [-t threads] [-s shell] [-e e-value] [-f out format] [-g genetic code]

	-i) fasta input file to be used as query
	-R) blast formatted reference prefix
	-o) output prefix
	-t) threads
	-s) shell or blast command use quotations for command
	-e) e- value
	-f) blast output format 
	-g) genetic code to be used for translations

_HELP

}

sub help_ESTscan {

print STDERR << "_HELP";

ESTscan:

	Usage: $0 -m $ops{m} <-i input file> <-m scoring matrix> [-o output prefix] [-l minimum result length]

	-i) fasta input file
	-m) scoring matrix. See documentation for ESTscan for information on how to create these
	-o) output prefix
	-l) minimum result length

_HELP

}

sub help_HMMscan {

print STDERR << "_HELP";

HMMscan:

	Usage: $0 -m $ops{m} <-i input file> <-d HMMER3 database> [-o output prefix] [-t threads]

	-i) peptide input file in fasta format
	-d) HMMER3 database to use with hmmscan
	-o) output file prefix
	-t) threads

_HELP

}

sub use_get_opts {
        use Getopt::Std;
        my $opts = 'hm:p:i:c:k:v:u:t:e:r:R:o:s:f:l:g:q:d:';
        getopts ("$opts", \%ops) or &help_main;
}

