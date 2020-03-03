# Full-Length-RNA-Analysis-Best-Practice
*This is a Full-Length RNA Analysis pipeline developted by BGI RD group.*

As we all know, with the progress of single molecule sequencing technology, full-length transcript sequencing will become more popular. Compared to the second generation sequencing technology, the third generation sequencing technology can detect full-length transcript from 5-end to polyA tail, this enables us to take the more accurate way to quantifying gene and isoform expression, and can take more accurate way to research isoform structure, such as alternative splicing(AS), alternative polyadenylation(APA), allele specific expression(ASE), transcription start site(TSS), fusion gene, UTR length and UTR secondary structure, etc.
Here, we provide a command line's version bioinformatics pipeline for PacBio IsoSeq data analysis from raw `subreads.bam`, this pipeline works well in both PacBio official IsoSeq library construction protocol and **BGI patented** `multi-transcripts in one ZMW library (MTZL) construction protocol` and `full-length polyA tail detection library construction protocol`. This pipeline contains quality control, basic statistics, full-length transcripts identification, UMI detection, isoform clustering, error correction and isoform quantification, which is free of compilation and very easy to use.

More about the library construction protocol detail and performance can find in this wiki：https://github.com/shizhuoxing/BGI-Full-Length-RNA-Analysis-Pipeline/wiki

# Dependencies
* SMRTlink 8.0 or later `you can install it in light way: smrtlink_*.run --rootdir smrtlink --smrttools-only`
* ncbi-blast-2.2.26+ or later
* R-3.4.1 or later with ggplot2| gridExtra | grid

# Usage
export`smrtlink` `blast` `R` to you path first.
```
export PATH=$PATH:/smrtlink/smrtcmds/bin
export PATH=$PATH:/ncbi-blast-2.2.28+/bin
export PATH=$PATH:/R-3.1.1/bin
```

## Step1 raw data statistics (optional)
```
samtools view *.subreads.bam | awk '{print $1"\t"length($10)}' > tmp.len
sed '1i Subreads\tLength' tmp.len > subreads.len
perl PolymeraseReads.stat.pl subreads.len ./
perl SubReads.stat.pl subreads.len ./
```
## Step2 run CCS
```
ccs *.subreads.bam ccs.bam --min-passes 0 --min-length 50 --max-length 21000 --min-rq 0.75 -j 4
```
Start from SMRTlink8.0, CCS4.0 significantly speeds up the analysis and can be easily parallelized by using `--chunk`.

## Step3 classify CCS by primer blast
### 3.1) cat ccs result in bam format from each chunk
```
samtools view ccs.bam > ccs.sam
samtools view ccs.bam | awk '{print ">"$1"\n"$10}' > ccs.fa
```
### 3.2) make primer blast to CCS
```
makeblastdb -in primer.fa -dbtype nucl
blastn -query ccs.fa -db primer.fa -outfmt 7 -word_size 5 > mapped.m7
```
The following primer sequence is commonly used by PacBio official IsoSeq library construction protocol and `BGI patented multi-transcripts in one ZMW library construction protocol`.
```
$ cat primer.fa
>primer_F
AAGCAGTGGTATCAACGCAGAGTACATGGGGGGGG
>primer_S
GTACTCTGCGTTGATACCACTGCTTACTAGT
```
The following primer sequence is used by `BGI patented full-length polyA tail detection library construction protocol`.
```
$ cat primer.fa
>primer_F
AAGCAGTGGTATCAACGCAGAGTAC
>primer_S
AAGCAGTGGTATCAACGCAGAGTACATCGATCCCCCCCCCCCCTTT
```
### 3.3) classify CCS by primer
Here is an example for classifying CCS generate from PacBio official IsoSeq library construction protocol and `BGI patented multi-transcripts in one ZMW library construction protocol`.
```
classify_by_primer -blastm7 mapped.m7 -ccsfa ccs.fa -umilen 8 -min_primerlen 16 -min_isolen 200 -outdir ./
```
`classify_by_primer` wraps a tool to detect full-length transcript from CCS base on PacBio official IsoSeq library construction protocol and `BGI patented multi-transcripts in one ZMW library construction protocol`.
```
$ classify_by_primer -h

Despriprion: BGI version's full-length transcript detection algorithm for PacBio official IsoSeq library construction protocol and BGI patented multi-transcripts in one ZMW library construction protocol.
Usage: classify_by_primer -blastm7 mapped.m7 -ccsfa ccs.fa -umilen 8 -min_primerlen 15 -min_isolen 200 -outdir ./

Options:
        -blastm7*:              result of primer blast to ccs.fa in blast -outfmt 7 format
        -ccsfa*:                the ccs.fa you want to classify to get full-length transcript
        -umilen*:               the UMI length in your library, if set to 0 means nonUMI for library construction
        -min_primerlen*:        the minimum primer alignment length in ccs.fa
        -min_isolen*:           the minimum output's transcript length whithout polyA tail
        -outdir*:               output directory
```
Here is an example for classifying CCS generate from `BGI patented full-length polyA tail detection library construction protocol`, the parameters and usage are the same as in `classify_by_primer`.
```
classify_by_primer.fullpa -blastm7 mapped.m7 -ccsfa ccs.fa -umilen 8 -min_primerlen 16 -min_isolen 200 -outdir ./
```
## Step4 isoform clustering and polish the consensus (optional)
### 4.1) make isoseq_flnc.sam based on ccs.sam and isoseq_flnc.fasta
```
flnc2sam ccs.sam isoseq_flnc.fasta > isoseq_flnc.sam
samtools view -bS isoseq_flnc.sam > isoseq_flnc.bam
```
### 4.2) run `isoseq3 cluster` and split cluster result for multi-chunks
```
isoseq3 cluster isoseq_flnc.bam unpolished.bam --split-bam 3
```
### 4.3) run `isoseq3 polish` for each chunk of cluster result
Chunk and parallel processing of the data can significantly reduce polishing computing time.
If the compute nodes of your computing cluster allow it, the `--split-bam` set up to ~50 will have more significant speedup, `--split-bam` set up to 50 can complete polish analysis under 1 hours.
```
isoseq3 polish unpolished.0.bam *.subreads.bam polished.0.bam --verbose
isoseq3 polish unpolished.1.bam *.subreads.bam polished.1.bam --verbose
isoseq3 polish unpolished.2.bam *.subreads.bam polished.2.bam --verbose
```

***We recommend making isoforms Cluster and Polish as the optional step. Cluster and Polish can removed allele specific sequence variation, because this step takes the FLNC and clusters the isoforms by similarity and makes a multiple alignment of each cluster and performs error correction using this alignment.***

***As the Sequel and Sequel II system have improve the polymerase read length, the proportion of HiFi reads (CCS QV>0.99) have significant improve than before. In fact, in the case of that have a reference genome, we can directly map isforms to reference genome and then filtered according to mapping the quality to instead of polishing low quality isoforms.***

# Contact
If you have any questions, encounter problems or potential bugs, don’t hesitate to contact us. Either report issues on github or write an email to:

Zhuoxing Shi - shizhuoxing@bgi.com
