load taxGIFile='TAX2GI';
update;
import blastFile='diamond/PREFIX_contigs_diamond_nr.sam.gz', 'diamond/PREFIX_pe_diamond_nr.sam.gz', 'diamond/PREFIX_single_diamond_nr.sam.gz', 'diamond/PREFIX_merged_diamond_nr.sam.gz' meganFile='OUT_FILE' maxMatches=10 minScore=50.0 maxExpected=0.01 topPercent=10.0 minSupport=1 minComplexity=-1.0 useMinimalCoverageHeuristic=false useSeed=false useCOG=false useKegg=false paired=true suffix1='/1' suffix2='/2' useIdentityFilter=false textStoragePolicy=InRMAZ blastFormat=SAM mapping='Taxonomy:GI_MAP=true';
update;
quit;
