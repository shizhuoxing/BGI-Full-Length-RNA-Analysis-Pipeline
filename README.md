# Full-Length-RNA-Analysis-Best-Practice
*This is a Full-Length RNA Analysis pipeline developted by BGI RD group.*

As we all knows, with the progress of single molecule sequencing technology, full-length transcript sequencing will become more popular. Compared to the second generation sequencing technology, the third generations sequencing technology can detect full-length transcript from 5-end to polyA tail, this enables us to take the more accurate way to quantifying gene and isoform expression, and can take more accurate way to research isoform structure, such as alternative splicing(AS), alternative polyadenylation(APA), allele specific expression(ASE), transcription start site(TSS), fusion gene, UTR length and UTR secondary structure, etc.   
Here, we provide a command line's version bioinformatics pipeline for PacBio IsoSeq data analysis from raw `subreads.bam`, this pipeline was work well in both PacBio official IsoSeq library construction protocol and **BGI patented** `multi-transcripts in one ZMW library (MTZL) construction protocol` and `full-length polyA tail detection library construction protocol`. This pipeline contain quality control, basic statistics, full-length transcripts identification, isoform clustering, error correction and isoform quantification, and it is very easy to install and use.   

More about the library construction protocol detail and performance can find in this wiki：https://github.com/shizhuoxing/Full-Length-RNA-Analysis-Best-Practice/wiki

# Dependencies   
* SMRTlink 6.0 or later `you can install it in light way: smrtlink-*.run --rootdir smrtlink --smrttools-only`   
* ncbi-blast-2.2.26+ or later   
* R-3.4.1 or later with ggplot2| gridExtra | grid

# Usage
export`smrtlink` `blast` `R` to you path first.
```
export PATH=$PATH:/smrtlink/smrtcmds/bin
export PATH=$PATH:/ncbi-blast-2.2.28+/bin
export PATH=$PATH:/R-3.1.1/bin
```

## Step1 raw data statistics and chunking
### 1.1) get raw data statistics
```
samtools view *.subreads.bam | awk '{print $1"\t"length($10)}' > tmp.len
sed '1i Subreads\tLength' tmp.len > subreads.len
perl PolymerraseReads.stat.pl subreads.len ./
perl SubReads.stat.pl subreads.len ./
```
### 1.2) raw data chunking
Chunk and parallel processing of the raw data can significantly reduce computing time.   
If the compute nodes of your computing cluster allow it, the `--chunks` set up to >100 will have more significant speedup, `--chunks` set up to 100 can complete CCS analysis in few hours.
```
dataset create --type SubreadSet raw.subreadset.xml *.subreads.bam
dataset split --zmws --chunks 3 raw.subreadset.xml
```
## Step2 run CCS for each chunk
```
mkdir CHUNK1 && cd CHUNK1 && perl creat_chunk_rtc.pl raw.chunk1.subreadset.xml ./ > resolved-tool-contract-1.json && ccs --resolved-tool-contract resolved-tool-contract-1.json   
mkdir CHUNK2 && cd CHUNK2 && perl creat_chunk_rtc.pl raw.chunk2.subreadset.xml ./ > resolved-tool-contract-2.json && ccs --resolved-tool-contract resolved-tool-contract-2.json  
mkdir CHUNK3 && cd CHUNK3 && perl creat_chunk_rtc.pl raw.chunk3.subreadset.xml ./ > resolved-tool-contract-3.json && ccs --resolved-tool-contract resolved-tool-contract-3.json  
```
Here, the configure file `resolved-tool-contract.json` you can download in this repository must be in the same directory as the `creat_chunk_rtc.pl`.   
The configure file `resolved-tool-contract.json` contain CCS parameter as follow, you can easily modify it:
```
"options": {
              "pbccs.task_options.by_strand": false,
              "pbccs.task_options.max_drop_fraction": 0.8,
              "pbccs.task_options.max_length": 20000,
              "pbccs.task_options.min_length": 100,
              "pbccs.task_options.min_passes": 0,
              "pbccs.task_options.min_predicted_accuracy": 0.75,
              "pbccs.task_options.min_read_score": 0.65,
              "pbccs.task_options.min_snr": 3.75,
              "pbccs.task_options.min_zscore": -9999.0,
              "pbccs.task_options.model_path": "",
              "pbccs.task_options.model_spec": "",
              "pbccs.task_options.polish": true,
              "pbccs.task_options.report_file": "ccs_report.txt",
              "pbccs.task_options.num_threads": 8
},
```
## Step3 classify CCS by primer blast
### 3.1) cat ccs result in bam format from each chunk
```
ls CHUNK*/ccs.bam > ccs.bam.list && bamtools merge -list ccs.bam.list -out ccs.bam  
samtools view ccs.bam > ccs.sam
bamtools convert -format fasta -in ccs.bam -out ccs.fa 
samtools view ccs.bam | awk '{print $1"\t"length($11)"\t"$13"\t"$14}' | sed 's/np:i://' | sed 's/rq:f://' > ccs_stat.xls 
```
### 3.2) make primer blast to CCS
```
makeblastdb -in primer.fa -dbtype nucl
blastn -query ccs.fa -db primer.fa -outfmt 7 -word_size 5 > mapped.m7 
```
Here is a PacBio official IsoSeq library construction protocol and `BGI patented multi-isoforms in one ZMW library construction protocol` commonly used primer sequence.
```
$ cat primer.fa
>primer_F
AAGCAGTGGTATCAACGCAGAGTACGGGGGGGG
>primer_S
GTACTCTGCGTTGATACCACTGCTTACTAGT
```
### 3.3) classify CCS by primer
Here is an example for classify CCS generate from PacBio official IsoSeq library construction protocol and `BGI patented multi-isoforms in one ZMW library construction protocol`.
```
perl classify_by_primer.pl -blastm7 mapped.m7 -ccsfa ccs.fa -umilen 6 -min_primerlen 13 -outdir ./ 
```
`classify_by_primer.pl` wraps a tool to detection full-length transcript from CCS base on PacBio official IsoSeq library construction protocol and `BGI patented multi-isoforms in one ZMW library construction protocol`.
```
$ perl classify_by_primer.pl

Despriprion: BGI version's full-length transcript detection algorithm for PacBio official IsoSeq library construction protocol and BGI patented multi-isoforms in one ZMW library construction protocol.
Usage: perl classify_by_primer.pl -blastm7 mapped.m7 -ccsfa ccs.fa -umilen 8 -min_primerlen 15 -min_isolen 200 -outdir ./

Options:
        -blastm7*:              result of primer blast to ccs.fa in blast -outfmt 7 format
        -ccsfa*:                the ccs.fa you want to classify to get full-length transcript
        -umilen*:               the UMI length in your library, if set to 0 means nonUMI for library construction
        -min_primerlen*:        the minimum primer alignment length in ccs.fa
        -min_isolen*:           the minimum output's transcript length whithout polyA tail
        -outdir*:               output directory
```
## Step4 isoform clustering and polishing the consensus
### 4.1) make isoseq_flnc.sam base on ccs.sam and isoseq_flnc.fasta
```
perl fl_to_sam.pl ccs.sam isoseq_flnc.fasta > isoseq_flnc.sam   
samtools view -bS isoseq_flnc.sam > isoseq_flnc.bam
```
### 4.2) run isoseq3 cluster and split cluster result for multi-chunks
```
isoseq3 cluster isoseq_flnc.bam unpolished.bam --split-bam 3
```
### 4.3) run isoseq3 polish for each chunk of cluster result
Chunk and parallel processing of the data can significantly reduce polishing computing time.   
If the compute nodes of your computing cluster allow it, the `--split-bam` set up to ~50 will have more significant speedup, `--split-bam` set up to 50 can complete polish analysis under 1 hours.
```
isoseq3 polish unpolished.0.bam *.subreads.bam polished.0.bam --verbose
isoseq3 polish unpolished.1.bam *.subreads.bam polished.1.bam --verbose
isoseq3 polish unpolished.2.bam *.subreads.bam polished.2.bam --verbose
```
## Step5 isoform expression quantify and plotting
```
ls polished.*.bam > polished.bam.list && bamtools merge -list polished.bam.list -out polish.bam
samtools view polish.bam > polish.sam
perl classify_stat.pl ccs_stat.xls isoseq_flnc.polyAstat.xls > classify_stat.xls
perl get_polish_fl.pl ccs_stat.xls polish.sam isoseq_flnc.fasta ./ SampleName LibraryName
perl get_ccs_stat_polt.pl ./ ./
perl get_classify_stat_polt.pl ./ ./
perl get_cluster_stat_polt.pl ./ ./
```

# Contact
If you have any questions, encounter problems or potential bugs, don’t hesitate to contact us. Either report issues on github or write an email to:

Zhuoxing Shi - shizhuoxing@bgi.com
