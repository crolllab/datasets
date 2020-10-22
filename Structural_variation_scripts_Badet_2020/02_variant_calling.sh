##############################
### software requirements ###
##############################

## bam2fastq (here version 1.3.0) > https://github.com/PacificBiosciences/bam2fastx
## ngmlr (here version 0.2.7) > https://github.com/philres/ngmlr
## samtools (here version 1.9) > https://github.com/samtools/samtools
## sniffles (here version 1.0.11) > https://github.com/fritzsedlazeck/Sniffles
## SURVIVOR (here version 1.0.7) > https://github.com/fritzsedlazeck/SURVIVOR

##############################
##### required datasets #####
##############################

## 1) fasta whole-genome assembly files (here $isolate.assembly.fa) > available at https://github.com/crolllab/datasets/tree/master/19-isolate_pangenome_BMC_Biology_2020
## 2) long-reads sequencing files (here $isolate.subreads.bam)

## bash pipeline
for i in *.bam
do
isolate=$( echo $i | sed '/.subreads.bam//')
# convert pacbio subreads into fastq files
bam2fastq -o $isolate $i
# align fastq files to the reference genome assembly
ngmlr -t 4 -r IPO323.assembly.fa -q $isolate.fastq.gz -o $isolate.ngmlr.sam
# convert sam alignments to bam
samtools view -Sb $isolate.ngmlr.sam | samtools sort -o $isolate.ngmlr.sorted.bam -@ 4
# call structural variants
sniffles -t 2 -m $isolate.ngmlr.sorted.bam -d 5000 -f 0.8 -l 10 -q 30 -v $isolate.sniffles.vcf
grep -v "IMPRECISE" $isolate.sniffles.vcf > $isolate.sniffles.precise.vcf
done
# merge variants individually mapped to in isolate (similar 1kb distant variants are merged into a single event)
ls *.precise.vcf > all_vcf_calls.txt
SURVIVOR merge all_vcf_calls.txt 1000 -1 1 -1 -1 -1 all_sniffles_precise_SURVIVOR_1kb.vcf

# NB1 > the all_sniffles_precise_SURVIVOR_1kb.vcf file was further filtered to remove variants larger than 100'000 bp in size.
# NB2 > the same pipeline was applied to the progeny isolates assemblies using 1A5 as the reference genome