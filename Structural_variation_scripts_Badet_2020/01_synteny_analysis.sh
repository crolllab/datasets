##############################
### software requirements ###
##############################

## faSplit >(https://anaconda.org/bioconda/ucsc-fasplit)
## MUMmer (here version 3) > http://mummer.sourceforge.net/
## ngmlr (here version 0.2.7) > https://github.com/philres/ngmlr
## samtools (here version 1.9) > https://github.com/samtools/samtools
## bedtools (here version v2.26.0) > https://bedtools.readthedocs.io/en/latest/
## syri (here version 1.0) > https://github.com/schneebergerlab/syri
## R (here version 3.5.0) > https://www.r-project.org/

##############################
##### required datasets #####
##############################

## 1) fasta whole-genome assembly files (here $isolate.assembly.fa) > available at https://github.com/crolllab/datasets/tree/master/19-isolate_pangenome_BMC_Biology_2020
## 2) syri requires equal number of sequences for the pairwise annotation of rearrangements > process the 13 core chromosomes together and each accessory chromosome individually

## bash pipeline for synteny analysis
for i in *.genome.fa
do
## split genome assemblies by chromosome
faSplit byname $i outSplit/
done

# cd outSplit/
for i in *.chr_1.fa
do
isolate=$( echo $i | sed 's/.chr_1.fa//' )
# echo $isolate
chromosomes=$( ls *.fa )
core=$( echo $chromosomes | tr ' ' '\n' | egrep -v 'chr_14|chr15|chr_15|chr_16|chr_17|chr_18|chr_19|chr_20|chr_21' | grep "$isolate" )
## concatenante core chromosomes by isolate
cat $core > ${isolate}
done

## pipeline for core chromosomes analysis
# mkdir core
# mv *_core.fa core/
# cd core/

# split reference genome into 10kb bins
samtools faidx IPO323_core.fa
awk '{print $1"\t"$2}' IPO323_core.fa.fai > IPO323_core_genome.bed
bedtools makewindows -g IPO323_core_genome.bed -w 10000 > IPO323_10kb_windows_core_genome.bed

# proceed with query genome assemblies
for i in *_core.fa
do
query=$( echo $i | sed 's/_core.fa//')
reference="IPO323_core.fa"
# perform whole-genome alignement
nucmer --maxmatch -p IPO323.$query $reference $query
# filter based on minimal alignement identity and length
delta-filter -m IPO323.$query.delta -i 90 -l 100 > IPO323.$query.delta.filter
# extract coordinates
show-coords -THrd IPO323.$query.delta.filter > IPO323.$query.coords
# annotate syntenic/rearranged regions
syri -c IPO323.$query.coords --nc 4 -k -d IPO323.$query.delta -r $reference -q $query --all --prefix IPO323.$query
## reorder start and stop for bedtools and remove regions missing in the reference
awk '$2 > $3 { temp = $3; $3 = $2; $2 = temp } 1' OFS='\t' IPO323.${query}syri.out | grep -ve "^-" > IPO323.$query.syri.bed
# focus on local alignments (remove alignements with non-focal chromosome)
cat IPO323.$query.syri.sorted | egrep "^(chr_[0-9]+).*\1\t" > IPO323.$query.syri.local.bed
# sort the bed file
bedtools sort -i IPO323.$query.syri.local.bed | awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$11}' - > IPO323.$query.syri.local.sorted
# get the syri intersect in those 10kb windows
bedtools intersect -a IPO323_10kb_windows_core_genome.bed -b IPO323.$query.syri.local.sorted -loj > IPO323_${query}_10kb_core_syri.local.intersect
done

# concatenate all resulting files
cat *syri.local.intersect > IPO323_all_isolates_10kb_core_syri.local.intersect

# compute synteny ratios in R
#library(dplyr)
intersect <- read.table("IPO323_YEQ92_4_10kb_core_syri.intersect")
names(intersect) <- c("chr", "start", "end", "rchr", "rstart", "rend", "query", "sv")
## count the occurrences of each type of rearrangement by 10kb window
counts <- as.data.frame(intersect %>%
group_by(chr, start, end, sv, add=TRUE) %>%
summarise(n = dplyr::n()) )
## sum all rearrangement by 10kb window
stats <- plyr::ddply(counts, c("chr","start","end"), dplyr::summarise, syri_sum = sum(n))
## subset the syntenic regions
syntenic <- subset(counts, sv %in% c("SYN", "SYNAL"))
## count the syntenic regions by 10kb window
syntenic <- plyr::ddply(syntenic, c("chr","start","end"), dplyr::summarise, synteny_sum = sum(n))
## merge the stats for syntenic regions together with all regions
table <- merge(stats, syntenic, by=c("chr","start","end"))
## compute the synteny ratio in each 10kb window (given by the ratio of the number syntenic regions over the total number of regions) 
table$synteny_ratio <- table$synteny_sum / table$syri_sum
#write.table(table, file="IPO323_all_isolates_10kb_core_syri_synteny_ratios.txt", quote = F, sep = "\t", row.names = F, col.names=T)

