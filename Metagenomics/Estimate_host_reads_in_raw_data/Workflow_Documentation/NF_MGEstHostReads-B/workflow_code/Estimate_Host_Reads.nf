

log.info """\
    ESTIMATE HOST READS
    ===================================
    Download DB:        ${params.DL_kraken}
    Single end reads:   ${params.single_end}
    projectDir:         ${projectDir}"""
    .stripIndent()






process RUNSHEET_FROM_GLDS {
  // Downloads isa Archive and creates runsheet using GeneLab API
  tag "${ gldsAccession }"
  publishDir "${ params.outputDir }/${ gldsAccession }/Metadata",
    pattern: "*.{zip,csv}",
    mode: params.publish_dir_mode

  input:
    // TEMP: RESTORE ONCE OSD SUPPORT ADDED val(osdAccession)
    val(gldsAccession)
    path(dp_tools_plugin)

  output:
    path("${ gldsAccession }_*_v?_runsheet.csv"), emit: runsheet
    path("*.zip"), emit: isaArchive

  script:
    def injects = params.biomart_attribute ? "--inject biomart_attribute='${ params.biomart_attribute }'" : ''
    """

    dpt-get-isa-archive --accession ${ gldsAccession }
    ls ${dp_tools_plugin}

    dpt-isa-to-runsheet --accession ${ gldsAccession } \
      --plugin-dir ${dp_tools_plugin} \
      --isa-archive *.zip ${ injects }

    cat > versions.yaml <<EOF
    - name: "dp_tools"
      version: "\$(python -c 'import dp_tools; print(dp_tools.__version__)')"
      homepage: "https://github.com/J-81/dp_tools"
      workflow task: "${task.process}"
    EOF
    """
}
/*
process reads_from_runsheet{
    output:
    reads_ch
    
    reads_ch = Channel
                .fromPath(params.runsheet_path)
                .splitCsv(header: true, sep: ',', strip: true)
                .map { row ->
                    def sampleID = row.sample_name
                    def paired = row.paired_end == "True" ? true : false
                    def readPaths = paired ? [row.read1_path, row.read2_path] : [row.read1]
                    return [sampleID, readPaths]}
                .view{it}
}
*/


process map_species_to_genome_link {
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

process set_up_host_db {
    tag "Downloading human reads database to ${params.host_db_path}\n"
    publishDir path: "$projectDir" 

    output:
    path "${params.host_db_name}/"

    script:
    """
    kraken2-build --download-taxonomy --db ${params.kraken2-host-db}
    
    curl -LO ${params.species_link}
    gunzip *_genomic.fna.gz

    kraken2-build --add-to-library *_genomic.fna.gz --no-masking --db ${params.kraken2-host-db}

    kraken2-build --build --db ${params.kraken2-host-db} --threads 30 --no-masking

    kraken2-build --clean --db ${params.kraken2-host-db}
    """

}
/*
process generate_summary_report {
    tag "Generating summary report in ${params.kraken_output_dir}"

    script:
    """
    total_fragments=$(wc -l sample-1-kraken2-output.txt | sed 's/^ *      //' | cut -f 1 -d " ")

    fragments_classified=$(grep -w -c "^C" sample-1-kraken2-output.txt)

    perc_host=$(printf "%.2f\n" $(echo "scale=4; ${fragments_classified} / ${total_fragments} * 100" | bc -l))

    cat <( printf "Sample_ID\tTotal_fragments\tTotal_host_fragments\tPercent_host\n" ) \
        <( printf "Sample-1\t${total_fragments}\t${fragments_classified}\t${perc_host}\n" ) > Host-read-count-summary.tsv
    
    """
}
*/

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
    path "${sample_id}-kraken2-summary.txt"

    script:
    """
    kraken2 --db $database --gzip-compressed \
            --threads 2 --use-names --paired  \
            --output ${sample_id}-kraken2-output.txt \
            --report ${sample_id}-kraken2-report.tsv \
            ${reads_ch[0]} ${reads_ch[1]}

    # Compute total and classified fragments
    total_fragments=\$(wc -l ${sample_id}-kraken2-output.txt | sed 's/^ *//' | cut -f 1 -d " ")
    fragments_classified=\$(grep -w -c "^C" ${sample_id}-kraken2-output.txt)
    
    # Calculate percentage of classified fragments
    perc_host=\$(printf "%.2f\n" \$(echo "scale=4; \$fragments_classified / \$total_fragments * 100" | bc -l))
    
    # Write summary info
    printf "${sample_id}\t\$total_fragments\t\$fragments_classified\t\$perc_host\n" > ${sample_id}-kraken2-summary.txt

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
    path "${sample_id}-kraken2-summary.txt"
    

    script:
    """
    kraken2 --db $database --gzip-compressed --threads 2 --use-names  \
            --output ${sample_id}-kraken2-output.txt \
            --report ${sample_id}-kraken2-report.tsv \
            ${reads_ch[0]}

    """
}
/*
workflow {
    map_species_to_genome_link("blah").view{it}
    

}
*/

workflow {
if(params.DL_kraken == true){
    log.info "\nPreparing to download new host reads database"
    database_ch =  set_up_host_db()
    database_ch.view{"database path: ${it}"}
}

else {
    log.info "\nAccessing previous human reads database"
    database_ch = Channel.value(params.host_db_path)
    database_ch.view{"database path: ${it}"}
}

if(params.single_end == true) {
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

final_percent = output_ch[1]
    .collect{(it.text[0..5]).toFloat()}
    .average().trunc(2)
    .view{"\nRESULT: ${it}% of input reads were unclassified, available in ${params.kraken_output_dir}/reads "}
}

