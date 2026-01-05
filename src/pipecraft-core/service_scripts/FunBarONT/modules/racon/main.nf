// Polish sequences with Racon
process polish_with_racon {

    input:
    val(run_id)
    tuple val(barcode_dir_absolute), val(barcode_name), path(barcode_dir), path(BLASTDB_PATH), path(processing_dir), path(fastq_file), path(filtlong_file), path(centroids_file), path(minimap_file)
    val cpu_threads

    output:
    tuple val(barcode_dir_absolute), val(barcode_name), path(barcode_dir), path(BLASTDB_PATH), path(processing_dir), path(fastq_file), path(filtlong_file), path(centroids_file), path(minimap_file), path("$processing_dir/combined.${barcode_name}.racon.fasta"), emit: data_tuple
    path("${barcode_name}.racon.fasta"), emit: racon_fasta

    publishDir "${run_id}_results/04_polished_sequences/", pattern: "${barcode_name}.racon.fasta", mode: 'copy'

    script:
    """
    echo "\$(date '+%Y-%m-%d %H:%M:%S') 🛠️ Running Racon polishing" | tee -a $processing_dir/processing.log
    racon $processing_dir/$fastq_file -q 20 -w 100 -t $cpu_threads $processing_dir/$minimap_file $centroids_file > $processing_dir/combined.${barcode_name}.racon.fasta 2>> $processing_dir/processing.log
    # Copy for publishing
    cp $processing_dir/combined.${barcode_name}.racon.fasta ${barcode_name}.racon.fasta
    echo "\$(date '+%Y-%m-%d %H:%M:%S') ✅ Polishing complete with Racon" | tee -a $processing_dir/processing.log
    """
}