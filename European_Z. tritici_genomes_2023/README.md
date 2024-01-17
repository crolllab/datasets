### ### ### ### ### ### ### ### ### 
# HiFi assemblies annotation 
### ### ### ### ### ### ### ### ### 

## `Briefly :`
``As a firest step, single-end RNA-seq reads generated upon minimal media condition for 14 pangenome isolates were mapped to each of the four HiFi genome assemblies using hisat2 with default parameters (RNA-seq data from Arg00, Aus01, CH95, CNR93, CRI10, I93, IPO323, IR01_48b, ISY92, KE94, OregS90, TN09, UR95 and YEQ92 isolates).
The resulting alignment files were then used together with IPO323 proteins (Zymoseptoria_tritici.MG2.Grandaubert2015.proteins.fa) to train intron / start / stop hints during gene annotation using braker.``

## `Full pipeline here :`

### `-- read trimming using trimmomatic version 0.39 --`
##### ``trimmomatic SE -threads 4 raw_reads/$isolate.fastq.gz trim/${isolate}_trim.fastq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:5:10 MINLEN:50``

### `-- mapping RNA reads iteratively to each of the HiFi assemblies using HISAT2 version 2.2.1 --`
##### ``hisat2 -p 4 -x $assembly -U trim/${isolate}_trim.fastq.gz --rna-strandness R -S bam/$isolate.$rna.sam``

### `-- generate bam format alignment files using samtools version 1.18 --`
##### ``samtools view -bS bam/$isolate.$rna.sam -@ 4 | samtools sort -@ 4 - > bam/$isolate.$rna.bam``

### `-- gene annotation using RNA (given the 14 RNA-seq datasets) and protein (given IPO323 proteins) evidence for each assembly using braker version 2.1.6 --`
##### ``braker.pl --genome=$assembly --bam=$bam_list --alternatives-from-evidence=false --prot_seq=$pep --etpmode --gff3 --cores 8 --fungus --species=$isolate``
