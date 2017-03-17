# About

Realign Tumor/Normal Bam files  simultaneously.

dualRealign is a post alignment re-aligner useing SRMA (https://sourceforge.net/projects/srma/).
Especially, could be used for Tumor/Normal set of BAM files.   

# Licence

 GPL2.0


# Install

- Download dualRealign.jar from web sites.
- Install Java Runtime 1.5 or later.

# Required files

- Bam files for normal reads Alignment.
- Bam files for tumor reads Alignment.
- 2.bit reference file
- CaptureTarget Information (bed format,optional)


# Run


```
command

usage: dualRealign.jar realign -n <arg> -t <arg> -r <arg> -o <arg> [-ct <arg>]  [-nt <arg>]
 
 -n,--normalBam <arg>        normal bam file
 -t,--tumorBam <arg>         tumor bam file
 -r,--reference <arg>        2 bit genome reference file
 -o,--outdir <arg>           output directory
 -ct,--captureTarget <arg>   Capture target regions(bed format)
 -nt,--num threads <arg>     number of threads
