#!/usr/bin/perl -w

#chmod +x compare_GeneInGroup.pl 
#./compare_GeneInGroup.pl
#output: shows fusions which matches in all samples in one group

use strict;
use warnings;


my $writer = ">/home/yhoang/workspace_Mariani/DeFuse_Perl/compare_GeneInGroup.txt";
my @counts = 0; # nxn table with relative matches
my @save;	# save sample number
my $count = 0;
my $j = 0;
my $it = 0;
my $x = 1;	#group count

my @grouparray;

open WRITER, $writer or die $!;

for (my $g = 1; $g < 5; $g++){

	my $GROUP = "</home/yhoang/workspace_Mariani/DeFuse_Perl/group".$g;
#	print "GROUP: $GROUP\n";
	open GROUP, $GROUP or die "no open";
	print WRITER "Group $x: ";++$x;
	while (<GROUP>) {
		chomp();
		$grouparray[$g][$it] = $_;				# $g group count
		print WRITER $grouparray[$g][$it]."\t";	# $it sample count of one group
		++$it;
	}
	close(GROUP);

	print WRITER "\n\n";	

	my @list;
	for (my $c = 0; $c <$it;$c++){
		my $where = "</home/yhoang/workspace_Mariani/DeFuse_Perl/compare_myself/".$grouparray[$g][$c];	# sample $it
#		print "\tWHERE: $where\n";
		open LISTi, $where or die "no open1";
		while (<LISTi>) {
			chomp();
			my @l =split( /\t/, $_ );
			$list[$c][$j] = join ("\t",$l[0],$l[1],$l[2],$l[3]);	
		#	print WRITER "\t$c $j $list[$c][$j]\n";
			++$j;
		}
		$j = 0;
	}	
	for (my $i = 0; $i < $#{$list[0]}; $i ++){		#$sample1
		print WRITER "$list[0][$i]\t";
		for (my $n =1; $n < $it ; $n++){				#each sample>1
			$counts[$n]= 0;
			for (my $m = 0 ; $m < $#{$list[$n]};$m++){	#each fusion
				if ($list[0][$i] eq $list[$n][$m]){
				#	print WRITER "\t\tsample".$grouparray[$g][$n]." $list[$n][$m]\n";
					$counts[$n]=1;
				}		
			}
		}
		for (my $t = 1;$t<$it;$t++){
			if ($counts[$t] == 1) {
	#			$save[$count] = "sample $grouparray[$g][$t]\t";
				++$count;
			}
		}
		for (my $t = 1;$t<$it;$t++){
			$counts[$t] = 0;
		}
		if ($count == 8){
			++$count;
	#	print WRITER "@save\n";
		print WRITER "<--- $count samples have the same event.";
		}
		print WRITER "\n";	
		@save=0;
		$count=0;
	}
	$it = 0;
	
	print WRITER "\n";		
}
close(WRITER);
	
