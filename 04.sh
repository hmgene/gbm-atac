

#biocyc.txt                chromosome.txt            gwas.txt                  lipidmaps.txt             pfam.txt                  smart.txt                 
#biological_process.txt    cosmic.txt                interactions.txt          molecular_function.txt    prints.txt                smpdb.txt                 
#C16_NA_LUT_Up.tsv         gene3d.txt                interpro.txt              msigdb.txt                prosite.txt               wikipathways.txt          
#cellular_component.txt    geneOntology.html         kegg.txt                  pathwayInteractionDB.txt  reactome.txt              

rm results/fnc_long.txt
#rm results/fnc_gene_long.txt
for f in results/*/{biological_process,kegg,reactome}.txt ;do
    c=`echo "$f" | cut -d"/" -f 2`
    p=${f##*/};p=${p%.txt};
#TermID	Term	Enrichment	logP	Genes in Term	Target Genes in Term	Fraction of Targets in Term	Total Target Genes	Total Genes	Entrez Gene IDs	Gene Symbols
    hm cutn $f '^Term\$',logP,Gene_Symbols | awk '$2 < log(0.0001)' |\
    perl -ne 'chomp;my ($t,$p,$g)=split/\t/,$_; 
        my $c="'$c'";
        my @gg=split/,/,$g;
        next if $#gg < 10 || $#gg > 20;
        my @tt=split/_/,$t;
        my $t1=join("_",@tt[0..3]);
        print "$t1\t$c\t1\n";
        #map { print "$t1\t$_\t1\n"; } @gg;
    ' >> results/fnc_long.txt 
    #' >> results/fnc_gene_long.txt 
done 

ca rc2mat results/fnc_long.txt | hm heatmap - sample_and_pathway.pdf 7 7 
