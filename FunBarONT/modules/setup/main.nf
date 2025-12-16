process setup_processing_folder {
    input:
        tuple val(fastq_path_absolute), val(barcode_name), path(fastq_file), path(BLASTDB_PATH)

    output:
        tuple val(fastq_path_absolute), val(barcode_name), path(fastq_file), path(BLASTDB_PATH), path("processing_${barcode_name}"), emit: data_tuple

    script:
    """
    mkdir -p processing_$barcode_name
    touch processing_${barcode_name}/processing.log
    echo "\$(date '+%Y-%m-%d %H:%M:%S') ✅ Created processing directory: processing" | tee -a processing_$barcode_name/processing.log
    """
}

process unzip_merge_fastq {
    input:
        tuple val(fastq_path_absolute), val(barcode_name), path(fastq_file), path(BLASTDB_PATH), path(processing_dir)

    output:
        tuple val(fastq_path_absolute), val(barcode_name), path(fastq_file), path(BLASTDB_PATH), path(processing_dir), path("$processing_dir/${barcode_name}.fastq.gz"), emit: data_tuple

    script:
    """
    echo "\$(date '+%Y-%m-%d %H:%M:%S') 📂 Processing: ${processing_dir}" | tee -a $processing_dir/processing.log
    # Check if file is already gzipped
    if [[ ${fastq_file} == *.gz ]]; then
        cp ${fastq_file} $processing_dir/${barcode_name}.fastq.gz
    else
        gzip -c ${fastq_file} > $processing_dir/${barcode_name}.fastq.gz
    fi
    echo "\$(date '+%Y-%m-%d %H:%M:%S') ✅ Prepared: $processing_dir/${barcode_name}.fastq.gz" | tee -a $processing_dir/processing.log
    """
}