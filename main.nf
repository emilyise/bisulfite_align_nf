#!/usr/bin/env nextflow

/*
 *=============================================================================
 *               WGBS & RRBS Alignment Pipeline for Nextflow
 *=============================================================================
 *              https://github.com/robinfchan/bisulfite_align
 *-----------------------------------------------------------------------------
 */

def helpMessage() {
    log.info"""

    Usage:

    The typical command for running the pipeline:

       nextflow run bisulfite_align/main.nf --profile 'docker' --reads '*_R{1,2}.fastq.gz' --bismark_index </path/to/index/dir/>

    
    The typical command for running the pipeline within AWS Batch:

       Command: ['s3://path/to/staging/dir',
                 '--reads', 's3://path/to/fastq/files/*_R{1,2}.fastq.gz',
                 '--outdir', 's3://path/to/output/dir',
                 '--bismark_index', 's3://path/to/bismark/genome/index/dir',
                 '--profile', 'awsbatch']

    Basic Arguments:
     --profile                          Choice of: 'awsbatch', 'docker', 'conda'
     --reads [file]                     Path to input data; not needed if skipping FastQC & Trim Galore!
     --bismark_index [dir]              Path to Bismark genome reference dir; not needed if skipping alignment and extraction
     --outdir [dir]                     Directory to save results; must point to an S3 Bucket if on AWS Batch; default ''./results'
    
    Optional Arguments:
     --rrbs                             Set if MspI digested template DNA; will auto skip deduplication; default false
     --nugen                            Set if processing data generated with NuGEN Ovation RRBS Methyl-Seq libraries
     --custom_container [uri]           Input URI path for custom container; default image is rfchan/bisulfite_align

    Intermediate Arguments:
        Use these when skipping certain steps is necessary or desired. 
        Note, skipping any step will also break the downstream generation of a MultiQC report:
    --skip_fastqc                       Skips FastQC of raw reads
    --skip_trim                         Skips read trimming; will automatically skip raw read FastQC step
    --skip_align                        Skips Bismark Alignment; automatically invoked when "--bams" provided 
    --skip_dedup                        Skips BAM deduplication; automatically invoked with "--rrbs"
    --skip_extract                      Skips Bismark methylation call extraction

    --trimmed_reads [file]              Use in place of --reads to align trimmed read fastq files; automatically invokes "--skip_trim" 
    --bams [file]                       Use when skipping alignment and/or deduplication. If running deduplication input raw 
                                        aligned BAMs. If performing methylation call extraction from BAMs "--skip_align" is 
                                        automatically invoked. Invoke "--skip_dedup" if deduplication is not desired for input 
                                        BAMs. For methylation call extraction from RRBS BAMs invoking only "--rrbs" is needed

    Trim Galore! Options:
     --clip_r1 [int]                    Trim the specified number of bases from the 5' end of read 1 (or single-end reads); default 0
     --clip_r2 [int]                    Trim the specified number of bases from the 5' end of read 2 (paired-end only); default 0
     --three_prime_clip_r1 [int]        Trim the specified number of bases from the 3' end of read 1 AFTER adapter/quality trimming; default 0
     --three_prime_clip_r2 [int]        Trim the specified number of bases from the 3' end of read 2 AFTER adapter/quality trimming; default 0
       
        
    !! WARNING: the following overwrite command line settings for --clip_r1 --clip_r2 --three_prime_clip_r1 --three_prime_clip_r2 !!
     
     --truseq_epi                       Presets for WGBS applications using Illumina TruSeq Epigenome
     --single_cell                      Presets for increased stringency in scWGBS and scRRBS; ignored when "--nugen" supplied

    Bismark Alignment Options:
     --comprehensive                    Output information for all cytosine contexts; default = false
     --cytosine_report                  Output stranded cytosine report during Bismark's bismark_methylation_extractor step; default = false
     --non_directional                  Run alignment against all four possible strands; default = false
     --unmapped                         Save unmapped reads to fastq files; default = false
     --relax_mismatches                 Turn on to relax stringency for alignment (set allowed penalty with --num_mismatches); default = false
     --num_mismatches [float]           default value = 0.6 will allow a penalty of bp * -0.6 - for 100bp reads (bismark default is 0.2)
     --local_alignment                  Allow soft-clipping of reads (potentially useful for single-cell experiments); default = false
     --bismark_align_cpu_per_multicore [int]    Specify how many CPUs are required per --multicore for bismark align; default = 3
     --bismark_align_mem_per_multicore [str]    Specify how much memory is required per --multicore for bismark align; default = 15.GB

    """.stripIndent()
}

// Show help message
if ( params.help ) {
    helpMessage()
    exit 0
}

/*
 * Initialize Variables & Channels
 */

// Bismark genome ref
if( params.bismark_index ){
    Channel
        .fromPath(params.bismark_index, checkIfExists: true)
        .ifEmpty { exit 1, "Bismark index files not found: ${params.bismark_index}" }
        .into { ch_bismark_index_for_bismark_align; ch_bismark_index_for_bismark_methXtract }
    } else if (!params.skip_align || !params.skip_extract && !params.bismark_index) { 
        exit 1, "No reference genome index path specified!"
    } else {
    Channel
        .empty()
        .into { ch_bismark_index_for_bismark_align; ch_bismark_index_for_bismark_methXtract }
        }

// Trimming presets
clip_r1 = params.clip_r1
clip_r2 = params.clip_r2
three_prime_clip_r1 = params.three_prime_clip_r1
three_prime_clip_r2 = params.three_prime_clip_r2
if( params.truseq_epi ){
    clip_r1 = 8
    clip_r2 = 8
    three_prime_clip_r1 = 8
    three_prime_clip_r2 = 8
}
else if( params.single_cell && !params.nugen ){
    clip_r1 = 6
    clip_r2 = 6
    three_prime_clip_r1 = 6
    three_prime_clip_r2 = 6
}

if( workflow.profile == 'awsbatch') {
  // Check outdir paths to be S3 buckets if running on AWSBatch
  if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
  // Prevent trace files to be stored on S3 since S3 does not support rolling files
  if (params.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// Read files channels
if ( params.reads && !params.skip_fastqc && !params.skip_trim ) {
    Channel
        .fromFilePairs( params.reads, size: 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}" }
        .into { ch_read_files_for_fastqc; ch_read_files_for_trim_galore }
    } else if ( params.reads && params.skip_fastqc && !params.skip_trim ) {
        Channel
            .fromFilePairs( params.reads, size: 2 )
            .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}" }
            .set { ch_read_files_for_trim_galore }
        ch_read_files_for_fastqc = Channel.empty()
    } else if  ( params.reads && !params.skip_fastqc && params.skip_trim ) {
        Channel
            .fromFilePairs( params.reads, size: 2 )
            .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}" }
            .set { ch_read_files_for_fastqc }
        ch_read_files_for_trim_galore = Channel.empty()
    } else if ( !params.skip_fastqc && !params.skip_trim && !params.reads && !params.trimmed_reads && !params.bams ) {
        exit 1, "Must provide reads or BAMs!"
    } else {
        ch_read_files_for_fastqc = Channel.empty()
        ch_read_files_for_trim_galore = Channel.empty()
    }

// Assign pre-trimmed reads to channel if skipping Trim Galore!
if ( params.skip_trim || params.trimmed_reads && !params.skip_align ) {
    Channel
        .fromFilePairs( params.trimmed_reads, size: 2 )
        .ifEmpty { exit 1, "Cannot find any trimmed reads in dir: ${params.trimmed_reads}" }
        .set { ch_trimmed_reads_for_alignment_preproc }
    } else if ( params.skip_trim && !params.reads && !params.trimmed_reads && !params.skip_align ) {
        exit 1, "Must provide reads to align!"
    } else { ch_trimmed_reads_for_alignment_preproc = Channel.empty() }

// Assign pre-aligned bams for use in deduplication and/or methylation call extraction
if( !params.skip_dedup || !params.skip_extract ) {
    need_bams = true
} else { need_bams = false }

if ( params.bams ) {
    Channel
        .fromFilePairs( params.bams, size: 1 )
        .ifEmpty { exit 1, "Cannot find any bams in: ${params.bams}" }
        .set { ch_bam_for_bismark_deduplicate_preproc }
    } else if ( need_bams && !params.reads && !params.trimmed_reads && !params.bams ) {
        exit 1, "Must provide reads or BAMs to deduplicate and/or extract methylation calls!"
    } else { ch_bam_for_bismark_deduplicate_preproc = Channel.empty() }


// Parameters Summary
def summary = [:]
summary['Run Name']         = workflow.runName
if(params.reads)            summary['Reads'] = params.reads
if(params.trimmed_reads)    summary['Trimmed Reads'] = params.trimmed_reads
if(params.bams)             summary['BAMs'] = params.bams
if(params.bismark_index)    summary['Bismark Index'] = params.bismark_index
if(params.skip_fastqc || params.skip_trim || params.trimmed_reads)    summary['Skip FastQC'] = "Yes"
if(params.skip_trim || params.trimmed_reads) summary['Skip Trim Galore!'] = "Yes"
if(params.skip_align || params.bams)  summary['Skip Align'] = "Yes"
if(params.rrbs || params.skip_dedup) summary['Skip Deduplication'] = "Yes"
if(params.skip_extract) summary['Skip mC Extract'] = "Yes"
if(params.rrbs)             summary['RRBS Mode'] = 'On'
if(params.nugen)            summary['Nugen Trim Mode'] = 'On'
if(params.truseq_epi)       summary['Trimming Profile'] = 'TruSeq Epigenome'
if(params.single_cell && !params.nugen)     summary['Trimming Profile'] = 'Single Cell'
if(params.single_cell && params.nugen)      summary['Trimming Profile'] = 'Nugen Trim + Single Cell'
summary['Trimming']         = "5'R1: $clip_r1 / 5'R2: $clip_r2 / 3'R1: $three_prime_clip_r1 / 3'R2: $three_prime_clip_r2"
summary['All C Contexts']   = params.comprehensive ? 'Yes' : 'No'
summary['Cytosine report']  = params.cytosine_report ? 'Yes' : 'No'
if(params.unmapped)         summary['Save Unmapped Reads'] = 'Yes'
if(params.bismark_align_cpu_per_multicore) summary['Bismark align CPUs per --multicore'] = params.bismark_align_cpu_per_multicore
if(params.bismark_align_mem_per_multicore) summary['Bismark align memory per --multicore'] = params.bismark_align_mem_per_multicore
if(params.custom_container) summary['Custom Container'] = params.custom_container
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Pipeline dir']     = workflow.projectDir
summary['User']             = workflow.userName
summary['Config Profile']   = workflow.profile
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"

log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"


/*
 * -------------------------------
 *        PIPELINE START
 * -------------------------------
 */

/*
 * Pre-trim FastQC
 */
 if ( params.skip_fastqc || params.skip_trim || params.trimmed_reads){
     ch_fastqc_results_for_multiqc = Channel.from(false)
 } else {
    process fastqc {
        if (params.custom_container) container "${params.custom_container}"

        tag "$name"
        publishDir "${params.outdir}/fastqc", mode: 'copy', overwrite: true,
            saveAs: { filename ->
                        filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"
                    }

        input:
        set val(name), file(reads) from ch_read_files_for_fastqc

        output:
        file '*_fastqc.{zip,html}' into ch_fastqc_results_for_multiqc

        script:
        """
        fastqc --quiet --threads $task.cpus $reads
        """
    }
 }

/*
 * Trim Galore!
 */
 if ( params.skip_trim || params.trimmed_reads ){
     ch_trimmed_reads_for_alignment = Channel.empty()
     ch_trim_galore_results_for_multiqc = Channel.from(false)
 } else {
    process trim_galore {
        if (params.custom_container) container "${params.custom_container}"
        
        tag "$name"
        publishDir "${params.outdir}/trim_galore", mode: 'copy', overwrite: true,
            saveAs: {filename ->
                if( filename.indexOf("_fastqc") > 0 ) "FastQC/$filename"
                else if( filename.indexOf("trimming_report.txt" ) > 0) "reports/$filename"
                else if( filename.indexOf("fq.gz") > 0 ) "trimmed_reads/$filename"
                else null
            }

        input:
        set val(name), file(reads) from ch_read_files_for_trim_galore

        output:
        set val(name), file('*fq.gz') into ch_trimmed_reads_for_alignment
        file "*trimming_report.txt" into ch_trim_galore_results_for_multiqc
        file "*_fastqc.{zip,html}"

        script:
        def c_r1 = clip_r1 > 0 ? "--clip_r1 $clip_r1" : ''
        def c_r2 = clip_r2 > 0 ? "--clip_r2 $clip_r2" : ''
        def tpc_r1 = three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 $three_prime_clip_r1" : ''
        def tpc_r2 = three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 $three_prime_clip_r2" : ''
        def rrbs = params.rrbs && !params.nugen ? "--rrbs" : ''
        def nugen = params.nugen ? "-a AGATCGGAAGAGC -a2 AAATCAAAAAAAC" : ''
        def cores = 1
        if(task.cpus){
            cores = (task.cpus as int) - 4
            if (cores < 1) cores = 1
            if (cores > 1) cores = 4
        }
        """
        trim_galore \
        --fastqc \
        --gzip \
        --paired $reads \
        $rrbs \
        $nugen \
        $c_r1 \
        $c_r2 \
        $tpc_r1 \
        $tpc_r2 \
        --cores $cores
        """

    }
 }

/*
 * Nugen MspI site filter & Diversity Trimming
 */
if ( !params.nugen || params.trimmed_reads ){
    ch_nugen_trimmed_reads_for_alignment = Channel.empty()
} else if ( params.nugen ){
    process nugen_trim{
        if (params.custom_container) container "${params.custom_container}"

        tag "$name"
        publishDir "${params.outdir}/nugen_trim", mode: 'copy', overwrite: true

        input:
        set val(name), file(reads) from ch_trimmed_reads_for_alignment

        output:
        set val(name), file('*_trimmed.fq.gz') into ch_nugen_trimmed_reads_for_alignment

        script:
        reads_chunk = "-1 ${reads[0]} -2 ${reads[1]}"

        """
        python2 /trimRRBSdiversityAdaptCustomers.py $reads_chunk
        """
    }
}

/*
 * Bismark alignment
 */
if ( params.skip_align || params.bams ){
    ch_bam_for_bismark_deduplicate = Channel.empty()
    ch_bam_for_bismark_summary = Channel.empty()
    ch_bam_for_preseq = Channel.empty()
    ch_bismark_align_log_for_bismark_report = Channel.from(false)
    ch_bismark_align_log_for_bismark_summary = Channel.from(false)
    ch_bismark_align_log_for_multiqc = Channel.from(false)
} else {
    process bismark_align { 
    if (params.custom_container) container "${params.custom_container}"

    if ( params.nugen ) {
        ch_final_trimmed_reads_for_alignment = Channel.empty().mix(ch_nugen_trimmed_reads_for_alignment)
        } else { 
            ch_final_trimmed_reads_for_alignment = Channel.empty().mix(ch_trimmed_reads_for_alignment) 
        }

    tag "$name"
    publishDir "${params.outdir}/bismark_alignments", mode: 'copy', overwrite: true,
        saveAs: {filename ->
            if( filename.indexOf(".fq.gz") > 0 ) "unmapped/$filename"
            else if( filename.indexOf("report.txt") > 0 ) "reports/$filename"
            else if( filename.indexOf(".bam") > 0 ) "$filename"
            else null
        }

    input:
    set val(name), file(reads) from ch_final_trimmed_reads_for_alignment.mix(ch_trimmed_reads_for_alignment_preproc)
    file(index) from ch_bismark_index_for_bismark_align.collect()

    output:
    set val(name), file("*.bam") into ch_bam_for_bismark_deduplicate, ch_bam_for_bismark_summary, ch_bam_for_preseq
    set val(name), file("*report.txt") into ch_bismark_align_log_for_bismark_report, ch_bismark_align_log_for_bismark_summary, ch_bismark_align_log_for_multiqc
    file "*.fq.gz" optional true

    script:
    reads_chunk = "-1 ${reads[0]} -2 ${reads[1]}"

    // Optional Bismark parameters
    non_directional = params.non_directional ? "--non_directional" : ''
    unmapped = params.unmapped ? "--unmapped" : ''
    mismatches = params.relax_mismatches ? "--score_min L,0,-${params.num_mismatches}" : ''
    soft_clipping = params.local_alignment ? "--local" : ''

    // Set Bismark cores & memory according to the given task
    multicore = ''
    if( task.cpus ){
        // Numbers suggested for mouse genome adjusted for human genome (+14%)
        if( params.single_cell || params.non_directional ){
            cpu_per_multicore = 5
            mem_per_multicore = (21.GB).toBytes()
        } else {
            cpu_per_multicore = 3
            mem_per_multicore = (15.GB).toBytes()
        }
        // Check if the user has specified this and overwrite if so
        if(params.bismark_align_cpu_per_multicore) {
            cpu_per_multicore = (params.bismark_align_cpu_per_multicore as int)
        }
        if(params.bismark_align_mem_per_multicore) {
            mem_per_multicore = (params.bismark_align_mem_per_multicore as nextflow.util.MemoryUnit).toBytes()
        }
        // How many multicore splits can we afford with the cpus we have?
        ccore = ((task.cpus as int) / cpu_per_multicore ) as int
        // Check that we have enough memory, assuming 15GB memory per instance
        try {
            tmem = (task.memory as nextflow.util.MemoryUnit).toBytes()
            mcore = (tmem / mem_per_multicore) as int
            ccore = Math.min(ccore, mcore)
        } catch (all) {
            log.debug "Warning: Not able to set bismark align multicore based on available resources"
        }
        if( ccore > 1 ){
            multicore = "--multicore $ccore"
        }
    }

    // Main command
    """
    bismark $reads_chunk \
        --genome $index \
        $multicore \
        $non_directional \
        $unmapped \
        $mismatches \
        $soft_clipping
    """
    }
}

/*
 * Bismark deduplicate
 */
if ( params.rrbs || params.skip_dedup ){
    run_dedup = false
    } else {run_dedup = true}

if( !run_dedup && !params.skip_align ) {
    ch_bam_for_bismark_deduplicate.into { ch_bam_dedup_for_bismark_methXtract; ch_bam_dedup_for_qualimap }
    ch_bismark_dedup_log_for_bismark_report = Channel.from(false)
    ch_bismark_dedup_log_for_bismark_summary = Channel.from(false)
    ch_bismark_dedup_log_for_multiqc  = Channel.from(false)
} else if ( !run_dedup && params.skip_align ){
    ch_bam_for_bismark_deduplicate_preproc.into { ch_bam_dedup_for_bismark_methXtract; ch_bam_dedup_for_qualimap }
    ch_bismark_dedup_log_for_bismark_report = Channel.from(false)
    ch_bismark_dedup_log_for_bismark_summary = Channel.from(false)
    ch_bismark_dedup_log_for_multiqc  = Channel.from(false)
} else {
    process bismark_deduplicate {
        if (params.custom_container) container "${params.custom_container}"
        
        tag "$name"
        publishDir "${params.outdir}/bismark_deduplicated", mode: 'copy', overwrite: true,
            saveAs: {filename -> filename.indexOf("_report.txt") > 0 ? "reports/$filename" : "$filename"}

        input:
        set val(name), file(bam) from ch_bam_for_bismark_deduplicate.mix(ch_bam_for_bismark_deduplicate_preproc)

        output:
        set val(name), file("*.deduplicated.bam") into ch_bam_dedup_for_bismark_methXtract, ch_bam_dedup_for_qualimap
        set val(name), file("*.deduplication_report.txt") into ch_bismark_dedup_log_for_bismark_report, ch_bismark_dedup_log_for_bismark_summary, ch_bismark_dedup_log_for_multiqc

        script:
        fq_type = '-p' // change to '-s' for single end; will add feature later
        """
        deduplicate_bismark $fq_type --bam $bam
        """
    }
}

/*
 * Bismark methylation extraction
 */
 if ( params.skip_extract ){
     // Do nothing
 } else {
    process bismark_methXtract {
        if (params.custom_container) container "${params.custom_container}"

        tag "$name"
        publishDir "${params.outdir}/bismark_methylation_calls", mode: 'copy', overwrite: true,
            saveAs: {filename ->
                if( filename.indexOf("splitting_report.txt" ) > 0 ) "reports/$filename"
                else if( filename.indexOf("M-bias" ) > 0) "m-bias/$filename"
                else if( filename.indexOf(".cov" ) > 0 ) "methylation_coverage/$filename"
                else if( filename.indexOf("bedGraph" ) > 0 ) "bedGraph/$filename"
                else if( filename.indexOf("CpG_report" ) > 0 ) "stranded_CpG_report/$filename"
                else "methylation_calls/$filename"
            }

        input:
        set val(name), file(bam) from ch_bam_dedup_for_bismark_methXtract
        file index from ch_bismark_index_for_bismark_methXtract.collect()

        output:
        set val(name), file("*splitting_report.txt") into ch_bismark_splitting_report_for_bismark_report, ch_bismark_splitting_report_for_bismark_summary, ch_bismark_splitting_report_for_multiqc
        set val(name), file("*.M-bias.txt") into ch_bismark_mbias_for_bismark_report, ch_bismark_mbias_for_bismark_summary, ch_bismark_mbias_for_multiqc
        file '*.{png,gz}'

        script:
        comprehensive = params.comprehensive ? '--comprehensive --merge_non_CpG' : ''
        cytosine_report = params.cytosine_report ? "--cytosine_report --genome_folder ${index} " : ''
        
        // Set multicore params
        multicore = ''
        if( task.cpus ){
            // Numbers based on Bismark docs
            ccore = ((task.cpus as int) / 3) as int
            if( ccore > 1 ){
                multicore = "--multicore $ccore"
            }
        }
        buffer = ''
        if( task.memory ){
            mbuffer = (task.memory as nextflow.util.MemoryUnit) - 2.GB
            // Only set if more than 6GB available
            if( mbuffer.compareTo(4.GB) == 1 ){
                buffer = "--buffer_size ${mbuffer.toGiga()}G"
            }
        }

        """
        bismark_methylation_extractor \
            $comprehensive \
            $cytosine_report \
            $multicore \
            $buffer \
            --ignore_r2 2 \
            --ignore_3prime_r2 2 \
            --bedGraph \
            --counts \
            --gzip \
            -p \
            --no_overlap \
            --report \
            $bam
        """
    }
 }

if (params.skip_fastqc || params.skip_trim || params.skip_align || params.skip_dedup || params.skip_extract || params.trimmed_reads || params.bams ){
    // Do nothing; skip Bismark Sample Report & Summary Report as now certain reports are outside the current session
} else {

    // Combine reports from Bismark steps above
    ch_bismark_align_log_for_bismark_report
        .join(ch_bismark_dedup_log_for_bismark_report)
        .join(ch_bismark_splitting_report_for_bismark_report)
        .join(ch_bismark_mbias_for_bismark_report)
        .set{ ch_bismark_logs_for_bismark_report }


    /*
    * Bismark Sample Report
    */
    process bismark_report {
        if (params.custom_container) container "${params.custom_container}"

        tag "$name"
        publishDir "${params.outdir}/bismark_reports", mode: 'copy', overwrite: true

        input:
        set val(name), file(align_log), file(dedup_log), file(splitting_report), file(mbias) from ch_bismark_logs_for_bismark_report

        output:
        file '*{html,txt}' into ch_bismark_reports_results_for_multiqc

        script:
        """
        bismark2report \
            --alignment_report $align_log \
            --dedup_report $dedup_log \
            --splitting_report $splitting_report \
            --mbias_report $mbias
        """
    }

    /*
    * Bismark Summary Report
    */
    process bismark_summary {
        if (params.custom_container) container "${params.custom_container}"

        publishDir "${params.outdir}/bismark_summary", mode: 'copy', overwrite: true

        input:
        file ('*') from ch_bam_for_bismark_summary.collect()
        file ('*') from ch_bismark_align_log_for_bismark_summary.collect()
        file ('*') from ch_bismark_dedup_log_for_bismark_summary.collect()
        file ('*') from ch_bismark_splitting_report_for_bismark_summary.collect()
        file ('*') from ch_bismark_mbias_for_bismark_summary.collect()

        output:
        file '*{html,txt}' into ch_bismark_summary_results_for_multiqc

        script:
        """
        bismark2summary
        """
    }
}

if (params.skip_fastqc || params.skip_trim || params.skip_align || params.skip_dedup || params.skip_extract || params.trimmed_reads || params.bams ){
    // Do nothing
} else {
    /*
    * BAM Sorting & Qualimap
    */
    process qualimap {
        if (params.custom_container) container "${params.custom_container}"

        tag "$name"
        publishDir "${params.outdir}/qualimap", mode: 'copy', overwrite: true

        input:
        set val(name), file(bam) from ch_bam_dedup_for_qualimap

        output:
        file "${bam.baseName}_qualimap" into ch_qualimap_results_for_multiqc

        script:
        def avail_mem = task.memory ? ((task.memory.toGiga() - 6) / task.cpus).trunc() : false
        def sort_mem = avail_mem && avail_mem > 2 ? "-m ${avail_mem}G" : ''
        """
        samtools sort $bam \
            -@ ${task.cpus} $sort_mem \
            -o ${bam.baseName}.sorted.bam
        qualimap bamqc -bam ${bam.baseName}.sorted.bam \
            -outdir ${bam.baseName}_qualimap \
            --collect-overlap-pairs \
            --java-mem-size=${task.memory.toGiga()}G \
            -nt ${task.cpus}
        """
    }
}

if (params.skip_fastqc || params.skip_trim || params.skip_align || params.skip_dedup || params.skip_extract || params.trimmed_reads || params.bams ){
    // Do nothing
} else {
    /*
    * preseq
    */
    process preseq {
        if (params.custom_container) container "${params.custom_container}"

        tag "$name"
        publishDir "${params.outdir}/preseq", mode: 'copy', overwrite: true

        input:
        set val(name), file(bam) from ch_bam_for_preseq

        output:
        file "${bam.baseName}.ccurve.txt" into preseq_results

        script:
        def avail_mem = task.memory ? ((task.memory.toGiga() - 6) / task.cpus).trunc() : false
        def sort_mem = avail_mem && avail_mem > 2 ? "-m ${avail_mem}G" : ''
        """
        samtools sort $bam \\
            -@ ${task.cpus} $sort_mem \\
            -o ${bam.baseName}.sorted.bam
        preseq lc_extrap -v -B ${bam.baseName}.sorted.bam -o ${bam.baseName}.ccurve.txt
        """
    }
}

if (params.skip_fastqc || params.skip_trim || params.skip_align || params.skip_dedup || params.skip_extract || params.trimmed_reads || params.bams ){
    // Do nothing; skipping steps breaks MultiQC input channel structure -- need to run this manually in local env
} else {
    /*
    * MultiQC
    */
    process multiqc {
        if (params.custom_container) container "${params.custom_container}"

        publishDir "${params.outdir}/MultiQC", mode: 'copy', overwrite: true

        input:
        file ('fastqc/*') from ch_fastqc_results_for_multiqc.collect().ifEmpty([])
        file ('trim_galore/*') from ch_trim_galore_results_for_multiqc.collect().ifEmpty([])
        file ('bismark/*') from ch_bismark_align_log_for_multiqc.collect().ifEmpty([])
        file ('bismark/*') from ch_bismark_dedup_log_for_multiqc.collect().ifEmpty([])
        file ('bismark/*') from ch_bismark_splitting_report_for_multiqc.collect().ifEmpty([])
        file ('bismark/*') from ch_bismark_mbias_for_multiqc.collect().ifEmpty([])
        file ('bismark/*') from ch_bismark_reports_results_for_multiqc.collect().ifEmpty([])
        file ('bismark/*') from ch_bismark_summary_results_for_multiqc.collect().ifEmpty([])
        file ('qualimap/*') from ch_qualimap_results_for_multiqc.collect().ifEmpty([])
        file ('preseq/*') from preseq_results.collect().ifEmpty([])

        output:
        file "*multiqc_report.html" into ch_multiqc_report
        file "*_data"

        script:
        """
        multiqc . -m picard -m qualimap -m bismark -m samtools -m preseq -m cutadapt -m fastqc
        """
    }
}
