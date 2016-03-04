#! /usr/bin/perl -w
use strict;
use 5.0.10;
die "Usage: perl $0 <intron.gff>\n" if (@ARGV!=1);
my (%poses_gene,$gene,%type_gene,%score_gene,@genes,%gene_seq);
open IN,$ARGV[0];
$_=<IN>;
chomp $_;
my ($scf,$type,$s,$e,$score,$gene_attr)=(split,$_)[0,2,3,4,5,-1];
$gene=(split /=/,$gene_attr)[1];
push @genes,$gene;
$poses_gene{$gene}="$s\t$e";
$type_gene{$gene}=$type;
$score_gene{$gene}=$score;
my $test_line=1;
while (<IN>){
	chomp;
	$test_line++;
	($type,$s,$e,$gene_attr)=(split)[2,3,4,-1];
	if ($type=~/gene|tRNA/){
		$gene=(split /=/,$gene_attr)[1];
		if (exists $poses_gene{$gene}){
			#print "$test_line\t$gene\t$poses_gene{$gene}\n";
			if (not exists $gene_seq{$gene}){
				$gene_seq{$gene}=1; ###assignment first is important!!
				$gene.="_1";
			}else{
				$gene_seq{$gene}++;  ##plus first is important!!
				$gene.="_".$gene_seq{$gene};
			}
		}
		push @genes,$gene;
		$poses_gene{$gene}="$s\t$e";
		$type_gene{$gene}=$type;
		$score_gene{$gene}=(split)[5];
	}elsif($type=~/intron/){
		my $tp_se="\t$s\t$e";
		$poses_gene{$gene}="$poses_gene{$gene}$tp_se"
	}
}
close IN;

foreach my $tp_gene (sort keys %poses_gene){
	my @poses=sort {$a<=>$b} (split /\t/,$poses_gene{$tp_gene});
	my $source=($type_gene{$tp_gene} eq "tRNA")?"RNAweasel":"blastn";
	print "$scf\t$source\t$type_gene{$tp_gene}\t$poses[0]\t$poses[-1]\t$score_gene{$tp_gene}\t+\t.\tgene=$tp_gene\n";
	if (@poses>2){
		my ($exon_s,$exon_e);
		for (my $i=0;$i<=($#poses);$i+=2){
			$exon_s=($i==0)?$poses[0]:($poses[$i]+1);
			$exon_e=($i==($#poses-1))?$poses[-1]:($poses[$i+1]-1);
			print "$scf\tRNAweasel\texon\t$exon_s\t$exon_e\t.\t+\t.\tParent=$tp_gene\n";
		}
	}
}
