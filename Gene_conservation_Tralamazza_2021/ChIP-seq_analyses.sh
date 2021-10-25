### ChIP-Seq analysis

#SRA accession in Supplementary Table 1.
prefetch --option-file SRR_acession_list | fastq-dump

#sequence trimming
trimmomatic SE $PathToFastq/$fastq.file $PathToTrimmo/$fastq.trim ILLUMINACLIP:Trueseq3_SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
#alignment
bowtie2 -p 1 --local -x Genome_reference -U $PathToTrimmo/$fastq.trim -S $PathToSam/$sample.sam
#convert and sort file
samtools view -bS $PathToBam/$sample.sam | samtools sort - > $PathToBam/$sample.bam
#convert to bed
bamToBed -i $PathToBam/$sample.bam > $PathToBed/$sample.raw.bed
#sort
sort -k1,1 -k2,2n -k3,3n -k6,6 -o $PathToBed/$sample.bed $PathToBed/$sample_sort.bed

### HOMER

#create tag directory
makeTagDirectory /chip_seq/$folder_name -format bed -forceBED -mapq 1 $PathToBed/$sample_sort.bed
#peak calling
findPeaks /chip_seq/$ folder_name / -style histone -region -size 1000 -minDist 500 -C 0 > /findPeak_output/$Peak
#convert to bed
pos2bed.pl /chip_seq/findPeak_output/$Peak > /findPeak_output/$bed
#annotate peak based on gene coordinates
bedtools annotate -i reference.genes.bed -files $bed > $annot.bed

### TSS analysis +/- 1kb

#convert c hipseq bam to bw format
bamCoverage--bam $sample.bam -o $sample.bw –normalizeTo1x$genomeSize --binSize 10
#Deeptools for TSS analysis

computeMatrix reference-point --referencePoint TSS -b 1000 -a 1000 -R FgramR_fgsc_variable_H3k27.bed FgramR_highly_conserved_H3k27.bed -S chipseq/coverage_output_bw/wt_low_H3K27 -o matrix_output_TSS.gz --outFileSortedRegions ouput_k27me3_genes.bed
#plot heatmap
plotHeatmap -m matrix_output_TSS.gz --colorMap 'summer' -out heatmap_matrix_output_TSS.pdf
