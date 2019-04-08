# Full-Length-RNA-Analysis-Best-Practice
*This is a Full-Length RNA Analysis pipeline developted by BGI RD group.*

As we all knows, with the progress of single molecule sequencing technology, full-length transcript sequencing will become more popular. Compared to the second generation sequencing technology, the three generations sequencing technology can detect full-length transcript from 5-end to polyA tail, this enables us to take the more accurate way to quantifying gene and isoform expression, and can take more accurate way to research isoform structure, such as alternative splicing(AS), alternative polyadenylation(APA), Allele specific expression(ASE), fusion gene, UTR length and UTR secondary structure, etc.   
Here, we provide a bioinformatics pipeline for PacBio IsoSeq data analysis from raw subreads.bam. This pipeline contain quality control, basic statistics, full-length transcripts identification, clustering, error correction and isoform quantification.   

`SMRTlink   
Chunk CCS   
Classify by primer   
IsoSeq3   
Merge and quantify`

Dependency   
SMRTlink 6.0   
blast   
R

## Step1 chunking raw data
If we
```
*/smrtlink/smrtcmds/bin/dataset create --type SubreadSet */raw.subreadset.xml */m54269_190219_090012.subreads.bam
*/smrtlink/smrtcmds/bin/dataset split --zmws --chunks 3 */raw.subreadset.xml
```
## Step2 CCS
```
perl creat_chunk_rtc.pl */raw.chunk1.subreadset.xml */CHUNK1 > *CHUNK1/resolved-tool-contract-1.json   
*/smrtlink/smrtcmds/bin/ccs --resolved-tool-contract */CHUNK1/resolved-tool-contract-1.json   
perl creat_chunk_rtc.pl */raw.chunk2.subreadset.xml */CHUNK2 > *CHUNK2/resolved-tool-contract-1.json   
*/smrtlink/smrtcmds/bin/ccs --resolved-tool-contract */CHUNK2/resolved-tool-contract-1.json  
perl creat_chunk_rtc.pl */raw.chunk3.subreadset.xml */CHUNK3 > *CHUNK96/resolved-tool-contract-3.json   
*/smrtlink/smrtcmds/bin/ccs --resolved-tool-contract */CHUNK3/resolved-tool-contract-3.json  
```

*/smrtlink/smrtcmds/bin/bamtools convert -format fastq -in */CHUNK96/ccs.bam -out */CHUNK96/ccs.fq && 

awk 'NR%4 == 1 {print ">" substr($0,2)} NR%4 == 2 {print}' */CHUNK96/ccs.fq > */CHUNK96/ccs.fa && 

*/smrtlink/smrtcmds/bin/samtools view */CHUNK96/ccs.bam | awk '{print $1"\t"length($11)"\t"$13"\t"$14}' | sed 's/np:i://' | sed 's/rq:f://' > */CHUNK96/ccs.stat && 

*/blastn -query */CHUNK96/ccs.fa -db */blastdb/gi.primer.fa -outfmt 7 -word_size 5 > */CHUNK96/mapped.m7 && 

perl classify_by_primer.pl */CHUNK96/mapped.m7 */CHUNK96/ccs.fa */CHUNK96/ && 

*/smrtlink/smrtcmds/bin/samtools view */CHUNK96/ccs.bam > */CHUNK96/ccs.sam && 

perl fl_to_sam.pl */CHUNK96/ccs.sam */CHUNK96/isoseq_flnc.fasta > */CHUNK96/isoseq_flnc.sam && 
```
