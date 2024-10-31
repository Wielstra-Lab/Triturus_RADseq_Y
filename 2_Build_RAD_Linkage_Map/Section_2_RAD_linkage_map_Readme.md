## Section 2: Construct linkage map based on RADseq data from Triturus family and locate previously identified candidate Y-linked markers

> Use the denovo_map.pl script of the Stacks package to assign reads to RADtags, followed by the populations function to produce a multisample .vcf file for the RADtags:

```sh
denovo_map.pl -T 16 --samples ~/data1/RADseq/samples --popmap ~/data1/RADseq/popmap_linkage_1.txt -o ~/data1/RADseq/Linkage --paired

populations -P ~/data1/RADseq/Linkage -O ~/data1/RADseq/Linkage/VCF --popmap ~/data1/RADseq/popmap_linkage_1.txt -t 16 --vcf
````

> The multisample .vcf file is then strictly filtered for minor allele frequency, mean depth and missing data. Indels were removed and 1 SNP per marker was selected:

```sh
vcftools --vcf ~/data1/RADseq/Linkage/VCF/populations.snps.vcf --recode --recode-INFO-all --out Trit_linkage_joint.filtered --maf 0.2 --min-meanDP 10 --max-missing 0.95 --remove-indels --thin 500
````

### Create genotype calls for candidate Y-linked presence/absence markers 

> A BLAST database was created from the catalog.fa file created by denovo_map.pl

````sh
makeblastdb -dbtype nucl -in ~/data1/RADseq/Linkage/catalog.fa
````

> The candidate Y-linked markers identified in the know-sex adults were blasted against the linkage map catalog

```sh
blastn -qcov_hsp_perc 30 -outfmt 6 -perc_identity 95 -query Trit_Y_candidates_1.fa -db ~/data1/RADseq/Linkage/catalog.fa | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge > Y_candidates_in_linkage_map.blast
````

> The coverage of every marker in each sample in the linkage family was calculated using RADcov.pl (which calls RADcoverage.sh)

```sh
perl RADcov.pl -d ~/data1/RADseq/Linkage
````

> Coverage_sort_1.R was then used to produce a table of coverage in each sample for the candidate Y-linked markers

```sh
Rscript Coverage_sort_1.R ~/data1/RADseq/popmap_linkage_1.txt ~/data1/RADseq/Linkage/Coverage Y_candidates_in_linkage_map.blast sexed_Y_candidates_raw.cov
```

> This was then filtered with Coverage_filter_1.R, to eliminate markers with very low coverage, or coverage in the mother

```sh
Rscript Coverage_filter_1.R sexed_Y_candidates_raw.cov sexed_Y_candidates_filtered.cov 1 BW_0009all
```

> Add_to_call_table_sexed.R was used to transform the coverage table into Psuedo-SNP calls in a format compatable with LepMAP 3

```sh
Rscript Add_to_call_table_sexed.R sexed_Y_candidates_raw.cov Y_marker_calls_1.txt Y
```

### Construct Linkage maps with The LepMAP 3 package

> Starting with identifying recombination's with the ParentCall 2 program and the filtered multisample .vcf file

```sh
java -cp ~/LepMap/bin/ParentCall2 data = Trit_RAD_samples.ped vcfFile = Trit_linkage_joint.filtered.recode.vcf > Triturus_RAD_parent.call
```

> The resulting .call file was then concatenated with the output of the Add_to_call_table_sexed.R script to include the candiate Y-linked presence/absence markers

```sh
cat Triturus_RAD_parent.call Y_marker_calls_1.txt > Triturus_RAD_parent_with_Y.call
```

> The next 2 steps of the the LepMAP 3 pipeline (SeperateChromosomes 2 and JoinSingles 2) were run

```sh
java -cp ~/LepMap/bin/ SeparateChromosomes2  data = Triturus_RAD_parent_with_Y.call lodLimit = 20 distortionLod = 1 > Triturus_RAD_Y_map_1.txt
java -cp ~/LepMap/bin/ JoinSingles2All  data = Triturus_RAD_parent_with_Y.call  map = Triturus_RAD_Y_map_1.txt lodLimit = 15 > Triturus_RAD_Y_map_with_singles_1.txt
```

> Markers were ordered with the following loop running the OrderMarkers 2 function of LepMAP 3

```sh

input=Triturus_RAD_Y_map_with_singles_1.txt
output_stem=Trit_Y_map
outdir=Ordered_maps_1

mkdir $outdir

for i in {1..12}

do
  output=$outdir"/"$output_stem"_"$i".txt"
  output_SA=$outdir"/"$output_stem"_"$i"_SA.txt"

  echo "working on "$output

  java -cp ~/LepMap/bin/ OrderMarkers2 map = $input data = Triturus_RAD_parent_with_Y.call  chromosome = $i  numMergeIterations = 12 numPolishIterations = 6 minError=0.02 scale = M/N 2 numThreads = 8 sexAveraged=1 > $output
  tail -n +4 $output | while read marker pat_position mat_position line; do echo -e $marker"\t"$pat_position"\t"$mat_position"\t"$i; done >> $outdir"/"$output_stem"_ordered.txt"


  echo "done "$output
done

cat Lissotriton_RAD_1_Y_marker_parent.call|cut -f 1,2|awk '(NR>=7)' > snps.txt

awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' snps.txt $outdir"/"$output_stem"_ordered.txt" | cut -f 1,3,4,5 > $outdir"/"$output_stem"_mapped.txt"
```

> Resulting in the map file Trit_Y_map_mapped.txt



