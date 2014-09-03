#!/usr/bin/perl -w
use strict;
use Data::Dumper;

my $dir_table_file = $ARGV[0];
my $finaldir       = $ARGV[1];
my $workdir        = $ARGV[2];
my $tempdir;

my $fst_map             = "/1st_map";
my $fst_multiple_mapped = "/1st_multiple_mapped";
my $fst_unmapped        = "/1st_unmapped";

my $index = "/storage2/bowtie/indexes/h_sapiens_37_asm";
my $samples;

mkdir $finaldir;
mkdir $finaldir . $fst_map;
mkdir $finaldir . $fst_multiple_mapped;
mkdir $finaldir . $fst_unmapped;

checkNewSamples();
print scalar( keys %$samples ), " samples to be done\n";

foreach my $sample ( sort keys %$samples )
{

	$tempdir = $workdir . "/temp$sample";
	mkdir $tempdir;

	print "start processing sample ", $sample, "\n";

	decompress($sample);

	general_bowtie($sample);
	print "\n\n";

	system("rm -r $tempdir");

}

removeShared();

sub removeShared
{
	open( IPCS, "ipcs |" ) || die "could not make ipcs command";

	while ( my $line = <IPCS> )
	{

		if ( $line =~ /^(\S+)\s.*klaas / )
		{
			system "ipcrm -M $1";
		}
	}

}

sub checkNewSamples
{
	open( NEWSAMPLES, "<$dir_table_file" );
	while ( my $line = <NEWSAMPLES> )
	{
		chomp $line;
		my @sample = split( /\t/, $line );

		if (   -e $finaldir . $fst_unmapped . "/" . $sample[0] . "_1.fq.gz"
			&& -e $finaldir . $fst_unmapped . "/" . $sample[0] . "_2.fq.gz" )
		{
			print( "sample " . $sample[0] . " already done\n" );
		}
		else
		{
			$sample[1] =~ tr/,/ /;
			$samples->{ $sample[0] }->{"input1"} = $sample[1];

			$sample[2] =~ tr/,/ /;
			$samples->{ $sample[0] }->{"input2"} = $sample[2];

		}
	}
	close NEWSAMPLES;
}

sub decompress
{
	my $sample = $_[0];

	$samples->{$sample}->{"fastq1"} = "$tempdir/$sample" . "_1.fastq";
	$samples->{$sample}->{"fastq2"} = "$tempdir/$sample" . "_2.fastq";

my @samplesdir = split(/ /, $samples->{$sample}->{"input1"});

	if ( -e $samplesdir[0] )
	{
		system( "pigz -d -c "
			  . $samples->{$sample}->{"input1"} . " > "
			  . $samples->{$sample}->{"fastq1"} );

		system( "pigz -d -c "
			  . $samples->{$sample}->{"input2"} . " > "
			  . $samples->{$sample}->{"fastq2"} );
	}
	else
	{

		$samples->{$sample}->{"input1"} =~ s/fastq.gz/fastq/;
		system( "cat " . $samples->{$sample}->{"input1"} . " > " . $samples->{$sample}->{"fastq1"} );

		$samples->{$sample}->{"input2"} =~ s/fastq.gz/fastq/;
		system( "cat " . $samples->{$sample}->{"input2"} . " > " . $samples->{$sample}->{"fastq2"} );

	}

}

sub general_bowtie
{
	my $sample = $_[0];
	my $input1 = $samples->{$sample}->{"fastq1"};
	my $input2 = $samples->{$sample}->{"fastq2"};

	my $mapped_output       = "$tempdir/$sample.bowtie";
	my $unmapped_output     = "$tempdir$fst_unmapped/$sample.fq";
	my $multi_mapped_output = "$tempdir$fst_multiple_mapped/$sample.fq";

	mkdir "$tempdir$fst_unmapped";
	mkdir "$tempdir$fst_multiple_mapped";

	my $command =
"bowtie -a -X400 -n3 -m1 -p32 --shmem --un $unmapped_output --max $multi_mapped_output --chunkmbs 1024 $index -1 $input1 -2 $input2 $mapped_output";

	system($command);

	system( "pigz -c $tempdir$fst_unmapped/$sample"
		  . "_1.fq > $finaldir$fst_unmapped/$sample"
		  . "_1.fq.gz" );
	system( "pigz -c $tempdir$fst_unmapped/$sample"
		  . "_2.fq > $finaldir$fst_unmapped/$sample"
		  . "_2.fq.gz" );

#	system( "pigz -c $tempdir$fst_multiple_mapped/$sample"
#		  . "_1.fq > $finaldir$fst_multiple_mapped/$sample"
#		  . "_1.fq.gz" );
#	system( "pigz -c $tempdir$fst_multiple_mapped/$sample"
#		  . "_2.fq > $finaldir$fst_multiple_mapped/$sample"
#		  . "_2.fq.gz" );

}
