### RNA-Seq analysis

#trimming
trimmomatic PE -threads $n -basein $PathToFastq/$fastq.file -baseout $PathToTrimmo/$fastq.trim ILLUMINACLIP:Trueseq3_PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
#alignment
hisat2 -p $n -q -x /reference_IND -1 $PathToTrimmo/$1P.fastq.trim -2 $PathToTrimmo/$2P.fastq.trim  -S $PathToBam/$isolate.sam
#convert sam to bam
samtools view -bS $PathToBam/$isolate.sam | samtools sort - > $PathToBam/$isolate.bam
#read count
htseq-count -i transcript_id -f bam -a 10 -s reverse -m union $PathToBam/$isolate.bam $REF_GFF > $PathToCount/$isolate.count.txt
