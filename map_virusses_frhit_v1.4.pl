#!/usr/bin/perl -w

#mapping non-human fragments
# v1.2 added read low complexity filter
# v1.2.2 added option fasta-ref
# v1.4 added change_fastq_toN.pl for some samples

use strict;
use Bio::SearchIO;
use Cwd;
use Data::Dumper;

my $fastq_dir  = $ARGV[0];
my $work_dir   = $ARGV[1];
my $fastaDB=$ARGV[2];
my $temp_dir   = $work_dir . "/temp";
my $result_dir = $work_dir . "/fr-hit_results";

my @samples;

getSamples();

mkdir $result_dir;

foreach my $sample (@samples)
{
	mkdir $temp_dir;

	print "start processing sample ", $sample, "\n";

## gunzip extraction

	my $startTime = time;
	print "extracting";
	system( "gunzip -c $fastq_dir/$sample" . "_*.gz > $temp_dir/$sample.fq" );
	print "\telapsed time: ", time - $startTime, " s\n";

## fastq to fasta

	$startTime = time;
	print "converting to fasta";

system(
"perl ~/Perl_scripts/virusScan/FR-HIT_pipeline/change_fastq_toN.pl $temp_dir/$sample"
	);
	system(
"perl ~/progs/prinseq-lite-0.20.3/prinseq-lite.pl -fastq $temp_dir/$sample.fq -out_format 1 -lc_method dust -lc_threshold 4 -out_good $temp_dir/$sample"
	);
	print "\telapsed time: ", time - $startTime, " s\n";

## mapping with fr-hit

	$startTime = time;
	print "mapping on virus";
	frhitSamples($sample);
	print "\telapsed time: ", time - $startTime, " s\n\n";

	system("rm -rf $temp_dir");
}

sub getSamples
{
	open( LS, "ls $fastq_dir|" ) || die "could not make ls command";

	while ( my $file = <LS> )
	{
		chomp $file;
		if ( $file =~ /(.+)_1\.fq\.gz/ )
		{
			unless ( -e "$result_dir/$1.psl" )
			{
				push @samples, $1;
			}
		}

	}
	close LS;
}

sub frhitSamples
{
	my $sample = $_[0];
	system(
"/data/klaas/progs/fr-hit-v0.7-2011-12-16/fr-hit -a $temp_dir/$sample.fasta -d $fastaDB -o $result_dir/$sample.psl -g 1 -T 0 -r 1 -u 0 -e 0.00001"
	);
}

