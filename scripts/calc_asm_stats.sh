#!/bin/bash
for x in fermi raymeta sga masurca abyss;
	do assemblathon_stats.pl *$x'_contigs.fa';
done >> asm_stats.txt
