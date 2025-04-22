# bisulfite-align-NF: Modification

Modified to improve use with UB HPC and to allow use of NuGEN deduplication step. Includes updated config file with singularity profile option. We are working on including NuGEN deduplication step. Although deduplication is generally not recommended for RRBS data, the NuGEN Ovation assay specifically applies unique molecular identifiers for this purpose, contained in read 2. 

      -profile                          Choice of: 'singularity', 'awsbatch', 'docker', 'conda' **NOTE THE LACK OF DOUBLE DASH!**
---- 

# bisulfite-align-NF: Nextflow Bisulfite/MethylC-Seq Pipeline

Modified expansion of [nf-core/methylseq](https://github.com/nf-core/methylseq). Alignment is performed using [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and [Bismark](https://github.com/FelixKrueger/Bismark) with added support for [NuGEN Ovation RRBS](https://github.com/nugentechnologies/NuMetRRBS) libraries.

# Overview

This total [Nextflow](https://www.nextflow.io/) pipeline encompasses all major data processing steps for WGBS, RRBS, scWGBS, and scRRBS. Based upon [nf-core/methylseq](https://github.com/nf-core/methylseq), this pipeline has also been expanded to correct for diversity bases used in the [NuGEN Ovation RRBS](https://github.com/nugentechnologies/NuMetRRBS) kit. Additionally, the pipeline logic has been re-hashed to allow for running of specific processing steps or inputing of intermediate files for re-processing.

Steps include:

1. QC of raw reads with FastQC
2. Read trimming with Trim Galore!
3. Nugen read trimming (*if applicable*)
4. Alignment using Bowtie 2 & Bismark
5. Deduplication of alignments (*if applicable*)
6. Extraction of methylation calls using Bismark
7. Sample report for Bismark (*if applicable*) *
8. BAM QC using Qualimap (*if applicable*) *
9. Read saturation analysis with Preseq (*if applicable*) *
10. Overall sumamry report using MultiQC (*if applicable*) *

The default Docker image used for the pipeline is [rfchan/bisulfite_align](https://hub.docker.com/repository/docker/rfchan/bisulfite_align).

**Each step of the pipeline can be skipped** and the appropriate intermediate files can be used to start a given step as needed. 

***NOTE**: when using skip parameters steps 7-10 will not be run as the necessary output/report files are no longer contained with the pipeline flow.

## Usage
The typical command for running the pipeline:

    nextflow run bisulfite_align/main.nf -profile 'docker' --reads '*_R{1,2}.fastq.gz' --bismark_index </path/to/bismark/index/dir> --outdir <path/to/output/dir>

The typical command for running the pipeline within AWS Batch:

    Command: ['s3://path/to/staging/dir',
             '--reads', 's3://path/to/fastq/files/*_R{1,2}.fastq.gz',
             '--outdir', 's3://path/to/output/dir',
             '--bismark_index', 's3://path/to/bismark/genome/index/dir',
             '-profile', 'awsbatch']

## Arguments

Basic Arguments:

      -profile                          Choice of: 'awsbatch', 'docker', 'conda' **NOTE THE LACK OF DOUBLE DASH!**
     --reads [file]                     Path to input data; not needed if skipping FastQC & Trim Galore!
     --bismark_index [dir]              Path to Bismark genome reference dir; not needed if skipping alignment and extraction
     --outdir [dir]                     Directory to save results; must point to an S3 Bucket if on AWS Batch; default ''./results'

 Optional Arguments:

     --rrbs                             Set if MspI digested template DNA; will auto skip deduplication; default false
     --nugen                            Set if processing data generated with NuGEN Ovation RRBS Methyl-Seq libraries
     --custom_container [uri]           Input URI path for custom container; default image is rfchan/bisulfite_align

Intermediate Arguments:

    Use these when skipping certain steps is necessary or desired. 
    Note, skipping any step will also break the downstream generation of a MultiQC report.

    --skip_fastqc                       Skips FastQC of raw reads
    --skip_trim                         Skips read trimming; will automatically skip raw read FastQC step
    --skip_align                        Skips Bismark Alignment; automatically invoked when "--bams" provided 
    --skip_dedup                        Skips BAM deduplication; automatically invoked with "--rrbs"
    --skip_extract                      Skips Bismark methylation call extraction

    --trimmed_reads [file]              Use in place of --reads to align trimmed read fastq files; automatically invokes "--skip_trim"
    --bams [file]                       Use when skipping alignment and/or deduplication. If running deduplication input raw 
                                        aligned BAMs. If performing methylation call extraction from BAMs "--skip_align" is 
                                        automatically invoked. Invoke "--skip_dedup" if deduplication is not desired for input 
                                        BAMs. For methylation call extraction from RRBS BAM invoking only "--rrbs" is needed

Trim Galore! Options:

     --adapter1 [str]                   Use custom adapter 1 sequence for trim; quoted string; default = false
     --adapter2 [str]                   Use custom adapter 2 sequence for trim; quoted string; default = false
     --clip_r1 [int]                    Trim the specified number of bases from the 5' end of read 1 (or single-end reads); default = 0
     --clip_r2 [int]                    Trim the specified number of bases from the 5' end of read 2 (paired-end only); default = 0
     --three_prime_clip_r1 [int]        Trim the specified number of bases from the 3' end of read 1 AFTER adapter/quality trimming; default  = 0
     --three_prime_clip_r2 [int]        Trim the specified number of bases from the 3' end of read 2 AFTER adapter/quality trimming; default = 0
     --min_insert [int]                 Set the minimum length of post-trim read pairs; default 20

    
    WARNING: the following overwrite command line settings for: --clip_r1 --clip_r2 --three_prime_clip_r1 --three_prime_clip_r2 !!

     --truseq_epi                       Presets for WGBS applications using Illumina TruSeq Epigenome
     --single_cell                      Presets for increased stringency in scWGBS and scRRBS; ignored when "--nugen" supplied

Bismark Alignment Options:

     --comprehensive                    Output information for all cytosine contexts; default = false
     --cytosine_report                  Output stranded cytosine report during Bismark's bismark_methylation_extractor step; default = false
     --non_directional                  Run alignment against all four possible strands; default = false
     --unmapped                         Save unmapped reads to fastq files; default = false
     --num_mismatches [float]           E.g. 0.6 will allow a penalty of bp * -0.6 - for 100bp reads (default is 0.2)
     --local_alignment                  Allow soft-clipping of reads (potentially useful for single-cell experiments); default = false
     --bismark_align_cpu_per_multicore [int]    Specify how many CPUs are required per --multicore for bismark align; default = 3
     --bismark_align_mem_per_multicore [str]    Specify how much memory is required per --multicore for bismark align; default = 15.GB

# Output

Results for a full pipeline run will be saved to **--outdir** as:

    - MultiQC/
        - multiqc_data/
            - multiqc.log
            - multiqc_bismark_alignment.txt
            - multiqc_bismark_methextract.txt
            - multiqc_cutadapt.txt
            - multiqc_data.json
            - multiqc_fastqc.txt
            - multiqc_general_stats.txt
            - multiqc_qualimap_bamqc_genome_results.txt
            - multiqc_sources.txt
        - multiqc_report.html
    - bismark_alignments/
        - *.bam
    - bismark_deduplicated/
        - *.deduplicated.bam
    - bismark_methylation_calls/
        - bedGraph/
            - *.deduplicated.bedGraph.gz
        - m-bias/
            - *.deduplicated.M-bias.txt
        - methylation_calls/
            - CpG_OT_*.deduplicated.txt.gz
            - CpG_OB_*.deduplicated.txt.gz
            - CHH_OT_*.deduplicated.txt.gz
            - CHH_OB_*.deduplicated.txt.gz
            - CHG_OT_*.deduplicated.txt.gz
            - CHG_OB_*.deduplicated.txt.gz
        - methylation_coverage/
            - *.deduplicated.bismark.cov.gz
        - reports/
            - *.deduplicated_splitting_report.txt
        - stranded_CpG_report/
            - *.deduplicated.CpG_report.txt.gz
    - bismark_summary/
        - bismark_summary_report.html
        - bismark_summary_report.txt
    - fastqc/
        - *_fastqc.html
    - nugen_trim/
        - *.fq_trimmed.fq.gz
    - preseq/
        - *.ccurve.txt
    - qualimap/
        - <sample>_qualimap/
            - css/
            - images_qualimapReport/
            - raw_data_qualimapReport/
            - genome_results.txt
            - qualimapReport.html
    - trim_galore/
        - FastQC/
            - *_fastqc.html
            - *_fastqc.zip
        - reports/
            - *_fastq.gz_trimming_report.txt
        - trimmed_reads/
            - *.fq.gz
