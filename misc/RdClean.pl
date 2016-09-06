#!/usr/bin/perl
use strict;
use warnings;

chdir("./pkg/man");

for my $file (<*.Rd>){

	#print $file, "\n";
	open(IN,"<$file");
	my @lines = <IN>;
	close(IN);

	open(OUT,">$file");
	for my $l (@lines){
		if( ($l !~ /^ \\tab \\code\{\.\S+\} \\tab  -\\cr\n$/) and
		    ($l !~ /\\keyword\{.+\}/)	
		   ){

			$l =~ s/\[R_GlobalEnv\]//g;
			$l =~ s/R_GlobalEnv//g;
			print OUT $l;
		}
	}
	close(OUT);

}
