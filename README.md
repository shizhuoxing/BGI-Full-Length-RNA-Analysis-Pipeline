# Full-Length-RNA-Analysis-Best-Practice
*This is a Full-Length RNA Analysis pipeline developted by BGI RD group.*

As we all knows, with the progress of single molecule sequencing technology, full-length transcript sequencing will become more popular. Compared to the second generation sequencing technology, the three generations sequencing technology can detect full-length transcript from 5-end to polyA tail, this enables us to take the more accurate way to quantifying gene and isoform expression, and can take more accurate way to research isoform structure, such as alternative splicing(AS), alternative polyadenylation(APA), allele specific expression(ASE), fusion gene, UTR length and UTR secondary structure, etc.   
Here, we provide a bioinformatics pipeline for PacBio IsoSeq data analysis from raw subreads.bam. This pipeline contain quality control, basic statistics, full-length transcripts identification, clustering, error correction and isoform quantification.   

# Dependencies   
SMRTlink 6.0 or later  
blast   
R

# Usage
set `smrtlink/smrtcmds/bin` `blast` to you path.

## Step1 raw data chunking
Chunk and parallel processing of the raw data can significantly reduce computing time.
```
dataset create --type SubreadSet raw.subreadset.xml *.subreads.bam
dataset split --zmws --chunks 3 raw.subreadset.xml
```
## Step2 CCS for each chunk
```
perl creat_chunk_rtc.pl raw.chunk1.subreadset.xml ./ > resolved-tool-contract-1.json   
ccs --resolved-tool-contract resolved-tool-contract-1.json   
perl creat_chunk_rtc.pl raw.chunk2.subreadset.xml ./ > resolved-tool-contract-1.json   
ccs --resolved-tool-contract resolved-tool-contract-1.json  
perl creat_chunk_rtc.pl raw.chunk3.subreadset.xml ./ > resolved-tool-contract-3.json   
ccs --resolved-tool-contract resolved-tool-contract-3.json  
```
## Step3 classify ccs by primer blast
```
bamtools convert -format fastq -in ccs.bam -out ccs.fq 
awk 'NR%4 == 1 {print ">" substr($0,2)} NR%4 == 2 {print}' ccs.fq > ccs.fa 
samtools view ccs.bam | awk '{print $1"\t"length($11)"\t"$13"\t"$14}' | sed 's/np:i://' | sed 's/rq:f://' > ccs.stat 
blastn -query ccs.fa -db ./blastdb/gi.primer.fa -outfmt 7 -word_size 5 > mapped.m7 
perl classify_by_primer.pl mapped.m7 ccs.fa ./ 
samtools view ccs.bam > ccs.sam 
perl fl_to_sam.pl ccs.sam isoseq_flnc.fasta > isoseq_flnc.sam 
```
## Step4 isoform cluster

## Step5 isoform expression quntify


# Contact
If you have any questions, encounter problems or potential bugs, don’t hesitate to contact us. Either report issues on github or write an email to:

Zhuoxing Shi - shizhuoxing@bgi.com
