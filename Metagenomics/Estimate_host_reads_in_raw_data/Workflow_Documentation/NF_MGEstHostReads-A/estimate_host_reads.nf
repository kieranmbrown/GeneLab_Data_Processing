log.info """\
    ESTIMATE HOST READS
    ===================================
    Download DB:        ${params.DL_kraken}
    Single end reads:   ${params.single_end}
    projectDir:         ${projectDir}"""
    .stripIndent()


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

    script:
    """
    kraken2 --db $database --gzip-compressed \
            --threads 2 --use-names --paired  \
            --output ${sample_id}-kraken2-output.txt \
            --report ${sample_id}-kraken2-report.tsv \
            ${reads_ch[0]} ${reads_ch[1]}

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
    .view{"\nRESULT: ${it}% of input reads were unclassified. Output analysis files are available in ${params.kraken_output_dir}"}
}
