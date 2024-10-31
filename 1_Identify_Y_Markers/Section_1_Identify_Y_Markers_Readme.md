## Section 1: Identification of candidate Triturus Y-linked markers from RADseq data from individuals of known-sex

### Process RADseq data with Stacks to produce initial candidate markers

> Use the denovo_map.pl script of the Stacks package to assign reads to RADtags:

```sh
denovo_map.pl -T 16 --samples ~/data1/RADseq/samples --popmap ~/data1/RADseq/popmap_sexed_1.txt -o ~/data1/RADseq/sexed -M 10 -n 10 --paired

> With the output stacks files in directory …/RADseq/sexed run the script RADcov.pl 
> This loops over all .bam files in directory and uses the depth function of SAMtools return the depth of sequencing at each basepair
> The script RADcoverage.sh then produces a file for each sample listing the coverage of each RADtag

```sh
perl /RADcov.pl -d ~/data1/RADseq/sexed
```

> The script RADsort_sexed_1.R then ingests the coverage files for all samples to produce an initial table of candidate male-linked markers

```sh
Rscript RADsort_sexed_1.R ~/data1/RADseq/popmap_sexed_1.txt ~/data1/RADseq/sexed/Coverage Trit_Y_candidates_1.txt 
```

### Filter and rank the candidates

> Produce a fasta file of the sequences of the candidate markers from the catalog.fa file produced by Stacks (this must be unzipped first) 

```sh
gunzip ~/data1/RADseq/sexed/catalog.fa.gz

tail -n +2 Trit_Y_candidates_1.txt | cut -f 1 | while read line; do grep -A 1 ">$line " ~/data1/RADseq/sexed/catalog.fa; done > Trit_Y_candidates_1.fa
```

> BLAST the candidate Y_linked sequences against the catalog to get number of paralogs

```sh
blastn -query Trit_Y_candidates_1.fa -db ~/data1/RADseq/sexed/catalog.fa -perc_identity 90 -qcov_hsp_perc 25 -outfmt 6 > Trit_Y_blast_tables.txt
```

> The script RADsort_sexed_2.R ranks the candidate markers by absence in females, paralogy, presence in males and coverage in males

```sh
Rscript RADsort_sexed_2.R Trit_Y_candidates_1.txt Trit_Y_blast_tables.txt Trit_Y_candidates_ordered
```

> Primers can then be designed for selected candidate markers 

