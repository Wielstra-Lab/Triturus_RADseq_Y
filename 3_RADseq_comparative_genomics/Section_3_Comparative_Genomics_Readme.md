## Section 3: Comparative Genomics with P. waltl (and L. Vulgaris)


## BLAST against _P. waltl_ 

> a fasta file was made containing the sequences placed on the linkage map (and the candidate Y-linked markers)

```sh
cat Trit_Y_map_mapped.txt| cut -f 1 | while read line; do grep -A 1 ">$line " ~/data1/RADseq/Linkage/catalog.fa; done > Trit_RADmap_markers_1.fa

cat Trit_RADmap_markers_1.fa Trit_Y_candidates_1.fa > Trit_RADmap_markers_All_1.fa
```

> A BLAST database was created for the _P. waltl_ genome

```sh
makeblastdb -dbtype nucl -in aPleWal1.pri.20220803.fasta
```

> The sequences placed on the linkage map were BLASTed against the _P. waltl_ database

```sh
blastn -query Trit_RADmap_markers_All_1.fa -db aPleWal1.pri.20220803.fasta -outfmt 6 -evalue 1e-20 -word_size 11 -num_threads 8 > Trit_RADmap_blast_1_All_raw.txt

sort -k1,1 -k12,12nr Trit_RADmap_blast_1_All_raw.txt > Trit_RADmap_blast_1_All_sorted.txt
```

> The sorted results from the BLAST were then filtered with Filter_blast.R
```sh
Rscript Filter_blast.R Trit_RADmap_blast_1_All_sorted.txt Trit_RADmap_blast_1_All_filtered.txt
```

> The linkage map and filtered BLAST results were used to create an oxford plot, this was combined with the data from a previous study on Lissotriton vulgaris to create duel oxford plots for the Triturus and Lissotriton against _P. waltl_ with the following script:

Plot_map_v_genome_twin_charts.R