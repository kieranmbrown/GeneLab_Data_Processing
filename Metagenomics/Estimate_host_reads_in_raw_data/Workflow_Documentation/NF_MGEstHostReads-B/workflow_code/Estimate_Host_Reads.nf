

log.info """\
    ESTIMATE HOST READS
    ===================================
    Download DB:        ${params.DL_kraken}
    Single end reads:   ${params.single_end}
    projectDir:         ${projectDir}
    GLDS number:        ${params.gldsAccession}
    Runsheet:           ${params.use_runsheet}
    Specify reads:      ${params.specify_reads}"""
    .stripIndent()



process MAP_TO_GENOME {
    tag "Mapping species to reference genome for $species_name using ${params.ref_genomes_csv_path}"
    
    input:
    val species_name
    
    output:
    env genome_link 

    script:
    """
    genome_link=\$(awk -F',' -v species="$species_name" 'species == \$1 {print \$2}' ${params.ref_genomes_csv_path})
    
    """
}


process BUILD_KRAKEN2_DB {
    container params.kraken2container
    tag "Downloading host reads database to $db_path\n"
    publishDir path: "$projectDir" 

    input:
    val link
    path db_path
    

    output:
    path "${params.host_db_name}/"

    script:
    """
    kraken2-build --download-taxonomy --db $db_path
    
    curl -LO $link
    gunzip *_genomic.fna.gz

    kraken2-build --add-to-library *_genomic.fna.gz --no-masking --db $db_path

    kraken2-build --build --db $db_path --threads 30 --no-masking

    kraken2-build --clean --db $db_path
    """

}


process RUNSHEET_FROM_GLDS {
  // Downloads isa Archive and creates runsheet using GeneLab API
    container = params.dp_tools_container
    tag "$glds"
    publishDir "${projectDir}/$glds/Metadata",
    pattern: "*.{zip,csv}"

    input:
    val glds

    output:
    path("${glds}_*_v?_runsheet.csv"), emit: runsheet
    path("*.zip"), emit: isaArchive

    script:
    """

    dpt-get-isa-archive --accession $glds
    ls ${params.dp_tools_plugin}

    dpt-isa-to-runsheet --accession $glds \
      --plugin-dir ${params.dp_tools_plugin} \
      --config-type estimatehost \
      --isa-archive *.zip 

    """
}


process DOWNLOAD_GLDS_READS {
    container = params.genelab_utils_container
    tag "$glds"

    input:
    val glds

    output:
    path("*.gz"), emit: raw_reads

    script:
    """
    if [ ${params.RawFilePattern} == null ];then
        
        GL-download-GLDS-data -f -g $glds -p ".fastq.gz" -o ${params.reads_dir}
        
    else

        GL-download-GLDS-data -f -g $glds -p  ${params.RawFilePattern} -o ${params.reads_dir}

    fi

    """

}


process EXTRACT_PARAMS {
    input:
        tuple val(host_organism), val(sample_id), val(paired), val(r1_suffix), val(r2_suffix), val(readPaths)
        //will be angry about lacking r2 info with SE data, make optional?

    output:
        tuple val(host_organism), val(sample_id), val(paired), val(r1_suffix), val(r2_suffix), val(readPaths)

    script:
        """
        echo "Host organism: $host_organism"
        echo "Sample ID: $sample_id"
        echo "Paired-end: $paired"
        echo "R1 Suffix: $r1_suffix"
        echo "R2 Suffix: $r2_suffix"
        echo "Read Paths: $readPaths"
        """
}


process PE_kraken2 {
    container params.kraken2container 
    tag "$sample_id"
    publishDir "$params.kraken_output_dir", pattern: "*.{txt,tsv}"

    input:
    path database 
    tuple val(sample_id), path(reads_ch)
    

    output:
    path "${sample_id}-kraken2-output.txt"
    path "${sample_id}-kraken2-report.tsv"
    path "${sample_id}-host-read-count-summary.tsv"

    script:
    """
    kraken2 --db $database --gzip-compressed \
            --threads 2 --use-names --paired  \
            --output ${sample_id}-kraken2-output.txt \
            --report ${sample_id}-kraken2-report.tsv \
            ${reads_ch[0]} ${reads_ch[1]}

    total_fragments=\$(wc -l ${sample_id}-kraken2-output.txt | sed 's/^ *//' | cut -f 1 -d " ")
    fragments_classified=\$(grep -w -c "^C" ${sample_id}-kraken2-output.txt)

    # Calculate percentage of host fragments
    perc_host=\$(awk 'BEGIN {printf \"%.2f\", ('"\${fragments_classified}"' / '"\${total_fragments}"') * 100}')

    # Generate summary file
    printf "Sample_ID\\tTotal_fragments\\tTotal_host_fragments\\tPercent_host\\n" > ${sample_id}-host-read-count-summary.tsv
    printf "${sample_id}\\t\${total_fragments}\\t\${fragments_classified}\\t\${perc_host}\\n" >> ${sample_id}-host-read-count-summary.tsv
    

    """
}


process SE_kraken2 {
    
    container params.kraken2container 
    tag "$sample_id"
    publishDir "$params.kraken_output_dir", pattern: "*.{txt,tsv}"
    publishDir "$params.kraken_output_dir/reads", pattern: "*.fastq.gz"

    input:
    path database 
    tuple val(sample_id), path(reads_ch)

    output:
    path "${sample_id}-kraken2-output.txt"
    path "${sample_id}-kraken2-report.tsv"
    path "${sample_id}${params.SE_reads_out_suffix}.gz"

    script:
    """
    kraken2 --db $database --gzip-compressed --threads 2 --use-names  \
            --output ${sample_id}-kraken2-output.txt \
            --report ${sample_id}-kraken2-report.tsv \
            ${reads_ch[0]}

    """
}


process COMBINE_SUMMARIES {
    tag "Joining summary files"
    publishDir "$params.kraken_output_dir"

    input:
    file removal_info_files

    output:
    file "Host-read-count-summary.tsv"

    script:
    """
    {
        printf "Sample_ID\\tTotal_fragments\\tTotal_host_fragments\\tPercent_host\\n"
        for f in ${removal_info_files.join(' ')}; do
            tail -n +2 "\$f"
        done
    } > Host-read-count-summary.tsv

    rm ${removal_info_files.join(' ')}
    
    """
}



workflow {


if (params.DL_kraken == true) {
    log.info "\nPreparing to download new host reads database"
    database_ch =  BUILD_KRAKEN2_DB()
    database_ch.view{"database path: ${it}"}
}

else {
    log.info "\nAccessing previous human reads database"
    database_ch = Channel.value(params.host_db_path)
    database_ch.view{"database path: ${it}"}
}


if (params.gldsAccession != false) {
    log.info "Preparing to download available reads for ${params.gldsAccession}"
    DOWNLOAD_GLDS_READS(params.gldsAccession)
    runsheet = RUNSHEET_FROM_GLDS(params.gldsAccession).out.runsheet
}
else {
    log.info "Processing local data (no GLDS accession number provided)"
    runsheet = params.runsheet_path
}


if (params.use_runsheet == true) {

    parse_runsheet = Channel
                .fromPath(params.runsheet_path)
                .splitCsv(header: true, sep: ',', strip: true)
                .map { row ->
                    def host_organism = row['host organism']
                    def sample_id = row['Sample Name']
                    def paired = row.paired_end == "True" ? true : false
                    def r1_suffix = row.raw_R1_suffix
                    def r2_suffix = row.raw_R2_suffix
                    def readPaths = paired ? [row.read1_path, row.read2_path] : [row.read1]
                    return [host_organism, sample_id, paired, r1_suffix, r2_suffix, readPaths]
                    }
                
    extracted_runsheet = EXTRACT_PARAMS(parse_runsheet)
    extracted_params = extracted_runsheet.first() //holds metadata in a value channel


        if (extracted_params.view{it[2]} == false) {  //if runsheet states single-end
            log.info "\nReading Single-end data from ${params.reads_dir}\n"

            reads_ch = extracted_runsheet.map{[it[1], [params.reads_dir + it[1] + it[3]]]}
        
            reads_ch.view{"reads: ${it}"}
            output_ch = SE_kraken2(database_ch, reads_ch)
        
        }
        else {      //if runsheet states paired-end
            log.info "\nReading Paired-end data from ${params.reads_dir}\n"
            reads_ch = extracted_runsheet.map{[it[1], [params.reads_dir + it[5][0], params.reads_dir + it[5][1]]]}

            reads_ch.view{"reads: ${it}"}
            output_ch = PE_kraken2(database_ch, reads_ch)
        }

}



else {

    if (params.single_end == true) {
        log.info "\nReading Single-end data from ${params.reads_dir}\n"

        if (params.specify_reads) {
            reads_ch = Channel
            .fromPath("${params.sample_id_list}")
            .splitText()
            .map { it.trim() }
            .map { sample_id ->
                def files = file("${params.reads_dir}${sample_id}${params.SE_reads_suffix}")
                return [sample_id, files]
            }
        }
        else {
        reads_ch = Channel
            .fromPath("${params.reads_dir}/*${params.SE_reads_suffix}", checkIfExists: true)
            .map { readfile ->
                def sampleId = readfile.name.replaceAll("${params.SE_reads_suffix}\$", "")
                return tuple(sampleId, readfile)
            }
        }
        reads_ch.view{"reads: ${it}"}
        output_ch = SE_kraken2(database_ch, reads_ch)
        
    }
    else {
        log.info "\nReading Paired-end data from ${params.reads_dir}\n"

        if (params.specify_reads) {
            reads_ch = Channel
            .fromPath("${params.sample_id_list}")
            .splitText()
            .map { it.trim() }
            .map { sample_id ->
                def files = file("${params.reads_dir}${sample_id}${params.PE_reads_suffix}").toList().sort()
                return [sample_id, files]
            }
        }
        else {
            reads_ch = Channel.fromFilePairs(params.reads_dir + "*" + params.PE_reads_suffix, checkIfExists: true)
        }
        reads_ch.view{"reads: ${it}"}
        output_ch = PE_kraken2(database_ch, reads_ch)
    }
}


removal_info_files_ch = output_ch[2] 
COMBINE_SUMMARIES(removal_info_files_ch.collect())


final_percent = output_ch[1]
    .collect{1-(it.text[0..5])
    .toFloat()}
    .average()
    .trunc(2)
    .view{"\nRESULT: ${it}% of input reads were associated with the host genome. Output files are available in ${params.kraken_output_dir}"}

}
