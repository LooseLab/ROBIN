#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// set up and create an output directory
outdir = file(params.outdir)
outdir.mkdir()

nextflow_version="v.0.8.2"

// default behaviour is to NOT:
// run sturgeon classifier 
// run NanoPlot

// default the various optional parameters to false
params.sturgeon = false
params.nanoplot = false
params.seq = false

log.info """\

        ================================================================
        Human-variation -> Rapid_CNS2 - Nextflow P I P E L I N E ${nextflow_version}
        ================================================================

        INPUTS
        ================================================================
        sample                  : ${params.sample}
        patient                 : ${params.patient}
        input_bam               : ${params.bam}
        outdir                  : ${params.outdir}
        reference               : ${params.reference}
        annotations             : ${params.annotations}
        threads                 : ${params.threads}
        bam_min_coverage        : ${params.bam_min_coverage}
        min_mgmt_coverage       : ${params.minimum_mgmt_cov}
        sturgeon                : ${params.sturgeon}
        nanoplot                : ${params.nanoplot}

        ================================================================
        To run with SLURM, add -process.executor='slurm' to your nextflow command.

        To ALSO run the sturgeon classifier, add the --sturgeon flag. (Default behaviour is to not run sturgeon)
        To generate a NanoPlot report, add the --nanoplot flag. (Default behaviour is to not generate QC report)
        ================================================================

        ================================================================
        Latest Changes in v.0.8.1:
        - wf-human-variation "--mod" module uses the BAM file with the merged mods
        - wf-human-variation "--snp" uses the raw BAM subsetted to just the target regions
        - wf-human-variation "--sv" use the subsetted BAM, which has all associated supplementary alignments retained
        - wf-human-variation "--cnv" uses the raw input BAM
        - two versions of the output report are created; a lite v with a simple mutations table and a full v with interactive igvreport
        - hv_rapidCNS2 workflow version number now included in the final report
        - sequencer information derived from input BAM and included in report (can be specified with --seq)
        Changes in v0.8.2:
        - MGMT coverage, methylation status, and methylation plot now correctly displayed in lite and full reports       
         ================================================================

        """
        .stripIndent()

process check_bam_has_meth_data {
    input:
        path(input_bam)
        val(threads)
    
    output:
        stdout emit: meth_check

    script:
        """
        samtools \
        view \
        -@${threads} \
        ${input_bam} \
        | grep -m 1 MM:Z
        """
}

process modkit_adjust_mods {
    input:
        path(input_bam)
        val(sample)
        val(threads)

    publishDir("${params.outdir}")

    output:
        path "*_modkit_merge.bam", emit: modkit_merged_bam

    script:
        """
        /modkit \
        adjust-mods \
        --convert h m \
        ${input_bam} \
        ${sample}_modkit_merge.bam \
        --threads ${threads}
        """
}

process index_input_bam {
    input:
        path(input_bam)
        val(threads)

    output:
        path "*.bai", emit: indexed_bam

    script:
        """
        samtools \
        index \
        -@${threads} \
        ${input_bam}
        """
}

process index_merged_bam {
    input:
        path(input_bam)
        val(threads)

    output:
        path "*.bai", emit: indexed_bam

    script:
        """
        samtools \
        index \
        -@${threads} \
        ${input_bam}
        """
}

process nanoplot {
    input:
        path(input_bam)
        val(threads)
        val(sample)
        path(indexed_bam)

    publishDir("${params.outdir}")

    output:
        path "*NanoPlot-report.html", emit: nanoplot

    script:
        """
        NanoPlot \
        -t ${threads} \
        --bam ${input_bam} \
        -p ${sample}
        """
}

process check_mgmt_coverage {
    input:
        path(input_bam)
	path(mgmt_bed)
	val(minimum_mgmt_coverage)
        path(indexed_bam)
        val(threads)

    publishDir("${params.outdir}")

    output:
	val true
	path "*_cov.txt", emit: mgmt_avg_cov_file
        path "mgmt_cov.mosdepth.summary.txt"	
        stdout emit: mgmt_avg_cov

    script:
        """
        /mosdepth \
        -t ${threads} \
        -n \
        --by ${mgmt_bed} \
        mgmt_cov ${input_bam}
        
        cov="\$(grep "^chr10_region" mgmt_cov.mosdepth.summary.txt | awk '{ print \$4 }')"
        
        echo \${cov}
        if awk 'BEGIN{exit ARGV[1]>ARGV[2]}' "\$cov" ${minimum_mgmt_coverage}
        then
            echo \${cov} > mgmt_below_thresh_cov.txt
        else
            echo \${cov} > mgmt_avg_cov.txt
        fi
	"""
}

process draw_mgmt_methylartist {
    maxRetries 5
    errorStrategy { (task.attempt <= maxRetries) ? 'retry' : 'ignore' }

    input:
        path(indexed_bam)
        path(bam)
        path(reference)
        val(params.outdir)
	val ready // mgmt_coverage has run

    publishDir("${params.outdir}")
    
    output:
        val true
        path "*.png", emit: mgmt_plot optional true
    
    script:
        cov_file = file("${params.outdir}/mgmt_avg_cov.txt")
        if ( cov_file.exists() == true )    
            """
            python3 /methylartist/methylartist \
            locus \
            -i chr10:129466536-129467536 \
            -b ${bam} \
            --ref ${reference} \
            --motif CG \
            --mods m
            """
        else

            """
            exit 1
            """
}

process mosdepth {
    input:
        val(threads)
        path(targets)
        path(input_bam)
        val(sample)
        path(indexed_bam)
    
    publishDir("${params.outdir}")
	
    output:
        path "*.mosdepth.summary.txt", emit: mosdepth_out

    script:
        """
        /mosdepth \
        -t ${threads} \
        -n \
        --by ${targets} \
        --fast-mode \
        ${sample} \
        ${input_bam}
        """
}

process index_subsetted_bam {
    input:
        path(input_bam)
        val(threads)
    
    output:
        val true

    script:
        """
        samtools \
        index \
        -@${threads} \
        ${input_bam}
        """
}

process human_variation_mods {
    input:
        path(input_bam)
        path(targets_bed)
        path(reference)
        val(sample)
        val(outdir)
        val(bam_min_coverage)
        path(indexed_bam)
        val(threads)
    
    publishDir("${params.outdir}")

    output:
        val true
    
    script:
        """
        nextflow run epi2me-labs/wf-human-variation \
        -with-report ${PWD}/${params.outdir}/human_variation_mods_nextflow_report.html \
        -profile standard \
        -w ${PWD}/${params.outdir}/workspace \
        --ref ${reference} \
        --mod \
        --bam ${input_bam} \
        --bed ${targets_bed} \
        --sample_name ${sample} \
        --bam_min_coverage ${bam_min_coverage} \
        --out_dir ${PWD}/${params.outdir} \
        --threads ${threads}
        """ 
}

process human_variation_cnv {
    input:
        path(input_bam)
        path(targets_bed)
        path(reference)
        val(sample)
        val(outdir)
        val(bam_min_coverage)
        path(indexed_bam)
        val(threads)

    publishDir("${params.outdir}")

    output:
        val true

    script:
        """
        nextflow run epi2me-labs/wf-human-variation \
        -with-report ${PWD}/${params.outdir}/human_variation_cv_nextflow_report.html \
        -profile standard \
        -w ${PWD}/${params.outdir}/workspace \
        --ref ${reference} \
        --cnv \
        --bam ${input_bam} \
        --bed ${targets_bed} \
        --sample_name ${sample} \
        --bam_min_coverage ${bam_min_coverage} \
        --out_dir ${PWD}/${params.outdir} \
        --threads ${threads}
        """
}

process human_variation_sv {
    input:
        path(input_bam)
        path(targets_bed)
        path(reference)
        val(sample)
        val(outdir)
        val(bam_min_coverage)
        path(indexed_bam)
        val(threads)

    publishDir("${params.outdir}")

    output:
        val true   

    script:
        """
        nextflow run epi2me-labs/wf-human-variation \
        -with-report ${PWD}/${params.outdir}/human_variation_sv_cnv_nextflow_report.html \
        -profile standard \
        -w ${PWD}/${params.outdir}/workspace \
        --ref ${reference} \
        --sv \
        --bam ${input_bam} \
        --bed ${targets_bed} \
        --sample_name ${sample} \
        --bam_min_coverage ${bam_min_coverage} \
        --out_dir ${PWD}/${params.outdir} \
        --threads ${threads} \
        --sniffles_args="--non-germline"
        """
}

process subset_bam_by_bed {
    input:
        path(input_bam)
        path(input_bed)
        val(sample)
        val(threads)

    publishDir("${params.outdir}")

    output:
        path "*subset.bam", emit: subsetted_bam

    script:
        """
        samtools view -@{threads} \
        -b \
        -h \
        -L ${input_bed} \
        ${input_bam} \
        > ${sample}_subset.bam
        """
}

process get_wanted_read_ids {
    input:
        path(subset_bam)
        val(threads)

    output:
        path "temp_ids.txt", emit: wanted_ids

    script:
        """
        samtools \
        view \
        -@${threads} \
        ${subset_bam} \
        | cut -f 1 \
        > temp_ids.txt
        """

}

process add_suppl_alignments_to_subset {
    input:
        path(subset_bam)
        path(input_bam)
        val(threads)
        val(sample)
        path(list_of_ids)

    publishDir("${params.outdir}")

    output:
        path "*_subset.suppl.bam", emit: suppl_subset_bam

    script:
        """
        samtools \
        view \
        -@${threads} \
        -h \
        -f 2048 \
        -N ${list_of_ids} \
        ${input_bam} \
        | samtools \
        merge \
        -@${threads} \
        -o ${sample}_subset.suppl.bam \
        ${subset_bam} \
        -
        """
}

process index_suppl_subset_bam {
    input:
        path(input_bam)
        val(threads)

    output:
        path "*.bai", emit: suppl_subset_indexed_bam

    script:
        """
        samtools \
        index \
        -@${threads} \
        ${input_bam}   
        """ 
}

process human_variation_snp {
    input:
        path(input_bam)
        path(targets_bed)
        path(reference)
        val(sample)
        val(outdir)
        val(bam_min_coverage)
        val(ready)
        val(threads)

    output:
        val true

    script:
        """
        nextflow run epi2me-labs/wf-human-variation \
        -with-report ${params.outdir}/human_variation_snp_nextflow_report.html \
        -profile standard \
        -w ${PWD}/${params.outdir}/workspace \
        --ref ${reference} \
        --snp \
        --bam ${input_bam} \
        --bed ${targets_bed} \
        --sample_name ${sample} \
        --bam_min_coverage ${bam_min_coverage} \
        --out_dir ${PWD}/${params.outdir} \
        --threads ${threads}
        """
 }

 process igv_reports {
    input:
        val(ready)  // filter-report ready
        val(sample)
        path(reference)
        path(input_bam)
        path(indexed_bam)
        path(annotations)

    publishDir("${params.outdir}")

    output:
        val true
 
    script:
        """
        sed -e 's/,/\t/g' -e 's/\"//g' \
        ${PWD}/${params.outdir}/${sample}_clair3_report.csv > ${PWD}/${params.outdir}/${sample}_clair3_report.fmt.csv 
        create_report ${PWD}/${params.outdir}/${sample}_clair3_report.fmt.csv \
        --fasta ${reference} \
        --sequence 1 \
        --begin 2 \
        --end 3 \
        --flanking 1000 \
        --info-columns Chr Start End Func Gene ExonicFunc AAChange cytoBand 1000g_EUR COSMIC \
        --output ${PWD}/${params.outdir}/${sample}_igv-report.html \
        --standalone \
        --tracks ${input_bam} ${annotations}
        """
}

process cnvpytor {
    input:
        val(sample)
        path(input_bam)
        val(threads)

    output:
        val true
    
    publishDir("${params.outdir}")

    script:
        """
        cnvpytor -root ${sample}_CNV.pytor -rd ${input_bam.toRealPath()} -j ${threads}
        cnvpytor -root ${sample}_CNV.pytor -his 1000 10000 100000 -j ${threads} # 4 mins
        cnvpytor -root ${sample}_CNV.pytor -partition 1000 10000 100000 -j ${threads} # SLOW
        cnvpytor -root ${sample}_CNV.pytor -call 1000 -j ${threads} > ${sample}.cnvpytor.calls.1000.tsv
        cnvpytor -root ${sample}_CNV.pytor -call 10000 -j ${threads} > ${sample}.cnvpytor.calls.10000.tsv
        cnvpytor -root ${sample}_CNV.pytor -call 100000 -j ${threads} > ${sample}.cnvpytor.calls.100000.tsv
        cnvpytor -root ${sample}_CNV.pytor -plot manhattan 100000 -chrom chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY -o ${PWD}/${params.outdir}/${sample}_cnvpytor_100k.png
        """
}

process vcftools {
    maxRetries 5
    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }

    input:
        val(ready) // human-variation SNP module
        val(sample)

    output:
        path "*.vcf", emit: variants_out

    script:
        """
        vcftools \
        --gzvcf ${PWD}/${params.outdir}/${sample}.wf_snp.vcf.gz \
        --remove-filtered-all \
        --recode \
        --out ${sample}
        """
}

process gzip {
    input:
        path(input_file)

    output:
        path "*.gz", emit: compressed_out

    script:
        """
        pigz \
        -1 \
        -c \
        ${input_file} \
        > ${input_file}.gz
        """
}

process bedtools_intersect {
    input:
        path(input1)
        path(input2)
        val(output_file)
        val(ext)

    output:   
        path "*.vcf", emit: intersect_vcf
        
    script:
        """
        bedtools \
        intersect \
        -a ${input1} \
        -b ${input2} > ${output_file}${ext}
        """
}

process convert2annovar{
    input:
        path(input)
        val(output_file)
        val(ext)

    output:
        path "*.avinput", emit: annovar_input

    script:
        """
        /annovar/convert2annovar.pl \
        -format vcf4 ${input} \
        -withfreq \
        -includeinfo \
        > ${output_file}${ext}
        """
}

process table_annovar {
    input:
        path(input)
        val(annovar_ver)
        val(output_file)
        val(ext)
        val(threads)
    
    output:
        path "*_multianno.csv", emit: clair3_output
      
    script:
        """
        /annovar/table_annovar.pl ${input} \
        /annovar/humandb/ \
        -buildver ${annovar_ver} \
        -out ${output_file}${ext} \
        -protocol refGene,cytoBand,avsnp147,dbnsfp30a,1000g2015aug_eur,cosmic70 \
        -operation gx,r,f,f,f,f \
        -nastring . \
        -csvout \
        -polish \
        -otherinfo \
        -thread ${threads}
        """
}

process bedtools_intersect2 {
    input:
        val(ready) // filtered bedmethyl file from filter to just 5mC
        path(input2)
        val(output_file)
        val(ext)
        val(sample)

    publishDir("${params.outdir}")

    output:   
        path "*.bed", emit: intersect_bed
        
    script:
        """
        bedtools \
        intersect \
        -a ${PWD}/${params.outdir}/${sample}.wf_mods.bedmethyl.gz \
        -b ${input2} > ${output_file}.${ext}
        """
}

process mgmt_pred {
    input:
        path(mgmt_pred)
        file(intersect_bed)
        path(probes)
        path(model)
        val(sample)
        val(params.outdir)
        val(threads)

    output:
        val true

    script:
        cov_file = file("${params.outdir}/mgmt_avg_cov.txt")
        if ( cov_file.exists() == true )
            """
            Rscript ${mgmt_pred} \
            --input ${intersect_bed} \
            --probes ${probes} \
            --model ${model} \
            --out_dir ${PWD}/${params.outdir} \
            --sample ${sample} \
            --threads ${threads}
            """
        else
            """
            """
}

process meth_classification {
    input:
        path(meth_class)
        val(sample)
        val(params.outdir)
        path(topprobes)
        path(trainingdata)
        path(arrayfile)
        val(threads)
        val(ready) // filter to just 5mC

    publishDir("${params.outdir}")

    output:
        val true

    script:
        """
        Rscript ${meth_class} \
        --sample ${sample} \
        --out_dir ${PWD}/${params.outdir} \
        --in_file ${PWD}/${params.outdir}/${sample}.wf_mods.bedmethyl.gz \
        --probes ${topprobes} \
        --training_data ${trainingdata} \
        --array_file ${arrayfile} \
        --threads ${threads}
        """
}

process filter_report {
    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries 5

    input:
        path(filterreport)
        path(clair3_multianno)
        val(sample)
        val(params.outdir)
    
    publishDir("${params.outdir}")

    output:
        val true
    
    script:
        """
        Rscript ${filterreport} \
        --input ${clair3_multianno} \
        --output ${PWD}/${params.outdir}/${sample}_clair3_report.csv \
        --sample ${sample}
        """        
}

process make_report {

    // cache false forces it to regenerate the report each time and not use the cache
    cache false

    input:
        path(makereport)
        val(ready) // cnvpytor
        val(ready) // mgmt_pred
        val(ready) // meth_classification
        val(ready) // filter_report
        val(sample)
        path(report_UKHD)
        path(mosdepth_plot_data) // mosdepth
        val(mgmt_cov)
        val(sturgeon) // sturgeon flag
        val(sturgeon) // sturgeon channel indicator
        val(igv_reports) //igv_reports has run
        val(nextflow_version)
        path(input_bam)
        val(seq)
    
    output:
	val true

    script:
        """
        # check for a specified sequencer, this overrides the checks below
        if [ "${seq}" != "false" ]
        then
            seq=${seq}
        else
            # try to get the sequencer model from the @RG group (if it exists)
            RG_seq=\$(samtools view -@4 -H ${input_bam} | grep ^@RG | grep -Po "PM:.*?\t" | awk '{print substr(\$NF,4,3)}')
            # if found it, save as seq
            if [ "\$RG_seq" ]
            then 
                seq=\$RG_seq
            else
            # if didn't find @RG, try for the fn:Z tag
                FN_seq=\$(samtools view ${input_bam} | grep -Po "fn:Z:[F,P]" | head -n 1 | awk '{print substr(\$NF,6,6)}')        
                # if found it, save as seq
                if [ "\$FN_seq" ]
                then
                    seq=\$FN_seq
                else
                    # try the f5:Z tag
                    F5_seq=\$(samtools view ${input_bam} | grep -Po "f5:Z:[F,P]" | head -n 1 | awk '{print substr(\$NF,6,6)}')
                    if [ "\$F5_seq" ]
                    then
                        seq=\$F5_seq
                    fi
                fi
            fi
        fi
        ## if all else fails, set it as unknown
        if [ "${seq}" != "false" ]
        then
            seq="Unknown"
        fi
       
        Rscript ${makereport} \
        --prefix ${sample} \
        --mutations ${PWD}/${params.outdir}/${sample}_clair3_report.csv \
        --cnv_plot ${PWD}/${params.outdir}/${sample}_cnvpytor_100k.global.0000.png \
        --rf_details ${PWD}/${params.outdir}/${sample}_rf_details.tsv \
        --votes ${PWD}/${params.outdir}/${sample}_votes.tsv \
        --output_dir ${PWD}/${params.outdir} \
        --patient ${sample} \
        --coverage ${PWD}/${params.outdir}/${sample}.mosdepth.summary.txt \
        --sample ${sample} \
        --report_UKHD ${report_UKHD} \
        --methylartist ${PWD}/${params.outdir}/*.chr10_129466536_129467536.m.ms1.smw20.locus.meth.png \
        --mgmt ${PWD}/${params.outdir}/${sample}_mgmt_status.csv \
        --sturgeon ${sturgeon} \
        --sturgeon_csv ${PWD}/${params.outdir}/${sample}_sturgeon_scores.csv \
        --sturgeon_pdf ${PWD}/${params.outdir}/${sample}_sturgeon_classification.pdf \
        --igv_report ${PWD}/${params.outdir}/${sample}_igv-report.html \
        --nextflow_ver ${nextflow_version} \
        --seq \${seq} \
        --promoter_mgmt_coverage ${mgmt_cov}
        """
}

process STURGEON_modkit_extract {
    input:
        path(mod_merged_bam)
        val(sample)
        val(threads)

    output:
        path "*_modkit_output.txt", emit: modkit_extract_output

    script:
        """
        /modkit \
        extract \
        ${mod_merged_bam} \
        ${sample}_modkit_output.txt \
        --threads ${threads}
        """
}

process STURGEON_inputtobed {
    input:
        path(modkit_out)
        val(params.outdir)

    output:
        path "${params.outdir}/merged_probes_methyl_calls.bed", emit: sturgeon_bed_convert

    script:
        """
        source /sturgeon-0.4.2/venv/bin/activate
        sturgeon \
        inputtobed \
        --reference-genome hg38 \
        -i ${modkit_out} \
        -o ${params.outdir} \
        -s modkit
        deactivate
        """
}

process STURGEON_predict {
    input:
        path(input_bed) // sturgeon inputtobed output

    output:
        val(true)

    script:
        """
        source /sturgeon-0.4.2/venv/bin/activate
        sturgeon \
        predict \
        -i  ${input_bed} \
        -o ${PWD}/${params.outdir} \
        --model-files /sturgeon-0.4.2/sturgeon/include/models/CAPPER_MODEL.zip \
        --plot-results        
        """
}

process STURGEON_rename {
    input:
        val(sturgeon) //  sturgeon predict has run
        val(sample)

    publishDir("${params.outdir}")

    output:

    script:
        """
        mv ${PWD}/${params.outdir}/merged_probes_methyl_calls_CAPPER_MODEL.pdf \
        ${PWD}/${params.outdir}/${sample}_sturgeon_classification.pdf
        mv ${PWD}/${params.outdir}/merged_probes_methyl_calls_CAPPER_MODEL.csv \
        ${PWD}/${params.outdir}/${sample}_sturgeon_scores.csv
        """
}

///////////////////////////
// MAIN WORKFLOW SECTION //
///////////////////////////

workflow {

    // Read all input parameters

    // User defined parameters

    Channel.fromPath(params.bam, checkIfExists: true)
    .set {input_bam}

    Channel.fromPath(params.reference, checkIfExists: true)
    .set {reference}

    Channel.from(params.sample)
    .set {sample}

    Channel.from(params.outdir)
    .set {outdir}

    Channel.from(params.patient)
    .set {patient}

    Channel.from(params.threads)
    .set {threads}

    Channel.from(params.minimum_mgmt_cov)
    .set {minimum_mgmt_cov}

    Channel.fromPath(params.annotations, checkIfExists: true)
    .set {annotations}

    // Collect variables and scripts from bin

    Channel.fromPath("${projectDir}/bin/NPHD_panel_hg38_clean.bed", checkIfExists: true)
    .set {targets}

    Channel.fromPath("${projectDir}/bin/mgmt_hg38.bed", checkIfExists: true)
    .set {mgmt_bed}

    Channel.fromPath("${projectDir}/bin/mgmt_probes.Rdata", checkIfExists: true)
    .set {probes}

    Channel.fromPath("${projectDir}/bin/mgmt_137sites_mean_model.Rdata", checkIfExists: true)
    .set {model}

    Channel.fromPath("${projectDir}/bin/mgmt_pred_v0.3.R", checkIfExists: true)
    .set {mgmt_pred}

    Channel.fromPath("${projectDir}/bin/methylation_classification_nanodx_v0.1.R", checkIfExists: true)
    .set {meth_class}

    Channel.fromPath("${projectDir}/bin/top_probes_hm450.Rdata", checkIfExists: true)
    .set {topprobes}
    
    Channel.fromPath("${projectDir}/bin/capper_top_100k_betas_binarised.Rdata", checkIfExists: true)
    .set {trainingdata}

    Channel.fromPath("${projectDir}/bin/HM450.hg38.manifest.gencode.v22.Rdata", checkIfExists: true)
    .set {arrayfile}

    Channel.fromPath("${projectDir}/bin/filter_report_v0.1.R", checkIfExists: true)
    .set {filterreport}

    Channel.fromPath("${projectDir}/bin/make_report_v0.4.R", checkIfExists: true)
    .set {makereport}

    Channel.fromPath("${projectDir}/bin/Rapid_CNS2_report_UKHD_v0.5.Rmd", checkIfExists: true)
    .set {report_UKHD}

    // WORKFLOW

    // check input bam file for methylation tags
    check_ch = check_bam_has_meth_data(input_bam, threads)

    // use modkit to merge 'm' and 'h' mods into a single score
    modkit_adjust_ch = modkit_adjust_mods(input_bam, sample, threads)

    // index the input bam
    index_ch = index_input_bam(input_bam, threads)

    // index the merged bam file 
    merged_index_ch = index_merged_bam(modkit_adjust_ch.modkit_merged_bam, threads)

    if ( params.nanoplot != false) {
        // generate NanoPlot QC Report
        nanoplot(input_bam, threads, sample, index_ch.indexed_bam)
    }

    // check the coverage of the mgmt region in the sequence data
    mgmt_coverage_ch = check_mgmt_coverage(input_bam, mgmt_bed, minimum_mgmt_cov, index_ch.indexed_bam, threads)

    // methylartist mgmt plot
    // methylartist only runs if the coverage is detected to be high enough (> minimum_mgmt_cov)
    methyl_artist_ch = draw_mgmt_methylartist(merged_index_ch.indexed_bam, modkit_adjust_ch.modkit_merged_bam, reference, outdir, check_mgmt_coverage.out[0])

    // mosdepth coverage plots
    mosdepth_ch = mosdepth(threads, targets, input_bam, sample, index_ch.indexed_bam)

    // subset the input bam file by the bed file - to speed up the SNP portion of human-variation
    subsetted_bam_ch = subset_bam_by_bed(modkit_adjust_ch.modkit_merged_bam, targets, sample, threads)

    // index the subsetted bam generated above
    index_subsetted_bam(subsetted_bam_ch.subsetted_bam, threads)

    // get the read IDs out of the subsetted bam, find all alignmnets associated with IDs and merge them back in
    wanted_ids_ch = get_wanted_read_ids(subsetted_bam_ch.subsetted_bam, threads)

    // using the wanted IDs, pull out the supplementary alignments from the input bam and merge them into the subsetted bam
    subset_suppl_bam_ch = add_suppl_alignments_to_subset(subsetted_bam_ch.subsetted_bam, input_bam, threads, sample, wanted_ids_ch.wanted_ids)

    // index the subsetted bam with the suppl reads added back in
    index_suppl_subset_bam_ch=index_suppl_subset_bam(subset_suppl_bam_ch.suppl_subset_bam, threads)

    // call and run the epi2me-labs/wf-human-variation : mods
    human_variation_mods(modkit_adjust_ch.modkit_merged_bam, targets, reference, sample, outdir, 1, merged_index_ch.indexed_bam, threads)

    // call and run the epi2me-labs/wf-human-variation : sv
    human_variation_sv(subset_suppl_bam_ch.suppl_subset_bam, targets, reference, sample, outdir, 1, index_suppl_subset_bam_ch.suppl_subset_indexed_bam, threads)

    // call and run the epi2me-labs/wf-human-variation : cnv
    human_variation_cnv(input_bam, targets, reference, sample, outdir, 1, index_ch.indexed_bam, threads)

    // run the SNP human variation workflow on the subsetted bam
    human_variation_snp(subsetted_bam_ch.subsetted_bam, targets, reference, sample, outdir, 1, index_subsetted_bam.out, threads)

    // run cnvPYTOR
    cnvpytor(sample, input_bam, threads)

    // run vcftools
    variants_ch = vcftools(human_variation_snp.out, sample)

    // compress output from vcftools
    compressed_variants_ch = gzip(variants_ch.variants_out)

    // run bedtools intersect on output of vcftools
    vcf_intersect_ch = bedtools_intersect(compressed_variants_ch.compressed_out, targets, sample, '_clair3_panel.vcf')

    // convert vcf to annovar input
    converted_annovar_ch = convert2annovar(vcf_intersect_ch.intersect_vcf, sample, '_clair3_panel.avinput')

    // run table_annovar on the converted annovar input
    clair3_annovar_ch = table_annovar(converted_annovar_ch.annovar_input, 'hg38', sample, '_clair3_panel', threads)

    // run bedtools intersect
    intersect_bed_ch = bedtools_intersect2(human_variation_mods.out, mgmt_bed, 'mgmt_5mC.hg38', 'bed', sample)

    // run the mgmt_pred script
    mgmt_pred_ch = mgmt_pred(mgmt_pred, intersect_bed_ch.intersect_bed, probes, model, sample, params.outdir, threads)
   
    // run the meth classification script
    meth_classification(meth_class, sample, params.outdir, topprobes, trainingdata, arrayfile, threads, human_variation_mods.out)
    
    // collect report data and generate report
    filter_report(filterreport, clair3_annovar_ch.clair3_output, sample, params.outdir)  

    //  produce igv_report for each SNP in the clair3 report
    igv_reports(filter_report.out, sample, reference, input_bam, index_ch.indexed_bam, annotations)

   // STURGEON section

    if ( params.sturgeon != false) {

        // modkit extract values
        modkit_extract_ch = STURGEON_modkit_extract(modkit_adjust_ch.modkit_merged_bam, sample, threads)

        // sturgeon convert input to bed file
        sturgeon_inputtobed_ch = STURGEON_inputtobed(modkit_extract_ch.modkit_extract_output, params.outdir)

        // sturgeon predict
        STURGEON_predict(sturgeon_inputtobed_ch.sturgeon_bed_convert)

        // sturgeon rename
        STURGEON_rename(STURGEON_predict.out, sample)

        // generate report inc. sturgeon
        make_report(makereport, cnvpytor.out, mgmt_pred.out, meth_classification.out, filter_report.out, sample, report_UKHD, mosdepth_ch.mosdepth_out, mgmt_coverage_ch.mgmt_avg_cov, params.sturgeon, STURGEON_predict.out, igv_reports.out, nextflow_version, input_bam, params.seq)

    } else {
        // collect data and generate final report - no sturgeon
        make_report(makereport, cnvpytor.out, mgmt_pred.out, meth_classification.out, filter_report.out, sample, report_UKHD, mosdepth_ch.mosdepth_out, mgmt_coverage_ch.mgmt_avg_cov, params.sturgeon, "NULL", igv_reports.out, nextflow_version, input_bam, params.seq)
    }
}

