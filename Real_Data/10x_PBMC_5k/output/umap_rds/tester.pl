my @methods = qw(BROCKMAN
chromVAR_kmers2
chromVAR_kmers
chromVAR_motifs2
chromVAR_motifs
Cicero2
Cicero
cisTopic
Control
Cusanovich2018
GeneScoring2
GeneScoring
scABC
Scasat
SCRAT2
SCRAT
SnapATAC
);

my @genes = ('S100A12', 'MS4A1', 'GAPDH');

foreach my $method (@methods)
{
	my $methodName = $method;
	$methodName =~ s/2$/ after PCA/;
	foreach my $gene (@genes)
	{
		print "plot_$method"."_$gene = plotStuff('$method','$methodName','$gene') + theme(legend.position = 'none')\n";
	}
}
print "============\n";

foreach my $gene (@genes)
{
	foreach my $method (@methods)
	{
		print "plot_$method"."_$gene,\n";
	}
}

