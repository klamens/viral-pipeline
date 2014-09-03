#!/usr/bin/perl -w
#add removal duplicate reads
#revision duplicate reads removal

use strict;
use Bio::SearchIO;
use Data::Dumper;

my $samples;
my $VirusID;
my $virus_name;

my $work_dir      = $ARGV[0];
my $result_dir    = $ARGV[1];
my $virus_csv_dir = $ARGV[2];

print "load virus_csv\n";
open( VIRUS, "<$virus_csv_dir" ) || die "could not load virus.csv";
while ( my $line = <VIRUS> )
{
	chomp $line;
	my @virus_acc = split( /\t/, $line );
	$virus_name->{ $virus_acc[0] } = $virus_acc[1];
}

print "check old report\n";
my $newReportInt = checkOldReport();

print "start checkSamples\n";
getSamples();

print "start openOldReport\n";
openOldReport();

print "start report\n";
summarize();

#print Dumper($samples);
report();

sub checkOldReport
{
	my $reportNR = 1;
	until ( !-e "$work_dir/fr-hit_report_1.4_$reportNR.csv" )
	{
		$reportNR++;
	}
	return $reportNR;
}

sub getSamples
{
	open( LS, "ls $result_dir|" ) || die "could not make ls command";

	while ( my $file = <LS> )
	{
		chomp $file;
		if ( $file =~ /(.+)\.psl/ )
		{
			$samples->{$1} = "empty";
		}

	}
	close LS;
}

sub openOldReport
{

	open( OLD, "<$work_dir/fr-hit_report_1.4_" . ( $newReportInt - 1 ) . ".csv" )
	  || die "$work_dir/fr-hit_report_1.4_" . ( $newReportInt - 1 ) . ".csv",
	  return;

	my $firstline = <OLD>;
	my @oldSamples = split /\t/, $firstline;
	shift(@oldSamples);
	pop(@oldSamples);

	while ( defined( my $line = <OLD> ) )
	{
		chomp $line;
		my @lines = split /\t/, $line;
		my $virus = shift(@lines);

		for ( my $i = 0 ; $i < @oldSamples ; $i++ )
		{
			$VirusID->{$virus}->{ $oldSamples[$i] } = $lines[$i];
		}
	}
	close OLD;
}

sub summarize
{
	foreach my $sample ( keys %$samples )
	{
		my $reads;
		unless ( $samples->{$sample} eq "empty" ) { next; }

		open( RESULTS, "<$result_dir/$sample.psl" );
		my $duplicates;
		while ( my $line = <RESULTS> )
		{
			
			chomp $line;

			my ( $ReadName, $ReadLength, $Eval, $AlignmentLength, $qBegin, $qEnd, $Strand, $Identity,
				$ReferenceSequenceName, $Begin, $End )
			  = split( /\t/, $line );

			$ReferenceSequenceName =~ /\|([^\|]*?)\|/;
			my $virus_id = $1;

			$ReadName =~ /(.+:\d+:\d+:\d+:\d+)(\/\d)?/;    #Asclepios verschil reverse and forward name
			$ReadName = $1;                                #Asclepios verschil reverse and forward name

			if ( $reads->{$ReadName} && $reads->{$ReadName}->[0] eq $virus_id )
			{
				my @ends =
				  sort { $a <=> $b }
				  ( $Begin, $End, $reads->{$ReadName}->[1], $reads->{$ReadName}->[2] );    #get most extreme positions

				if ( !$duplicates->{$virus_id}->{ $ends[0] } || $duplicates->{$virus_id}->{ $ends[0] } != $ends[3] )
				{
					$VirusID->{$virus_id}->{$sample}++;
					$duplicates->{$virus_id}->{ $ends[0] } = $ends[3];

				}
			}
			else
			{
				$reads->{$ReadName} = [ $virus_id, $Begin, $End ];
			}
		}

	}
}

sub report
{

	open( REPORT, ">$work_dir/fr-hit_report_1.4_" . ($newReportInt) . ".csv" );

	print REPORT "ACCESIONNR";
	foreach my $sample ( sort { $a <=> $b } keys %$samples )
	{
		print REPORT "\t$sample";
	}
	print REPORT "\tDESCRIPTION\n";

	foreach my $found_virus ( keys %$VirusID )
	{

		print REPORT "$found_virus";
		foreach my $sample ( sort { $a <=> $b } keys %$samples )
		{
			if ( exists $VirusID->{$found_virus}->{$sample} )
			{
				print REPORT "\t" . $VirusID->{$found_virus}->{$sample};
			}
			else { print REPORT "\t0"; }

		}
		print REPORT "\t" . $virus_name->{$found_virus} . "\n";
	}
	close REPORT;
}
