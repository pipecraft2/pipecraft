
// Second round of polishing - Medaka
process polish_with_medaka {

    input:
    val(run_id)
    tuple val(barcode_dir_absolute), val(barcode_name), path(barcode_dir), path(BLASTDB_PATH), path(processing_dir), path(fastq_file), path(filtlong_file), path(centroids_file), path(minimap_file), path(racon_file)
    val medaka_model
    val cpu_threads

    output:
    tuple val(barcode_dir_absolute), val(barcode_name), path(barcode_dir), path(BLASTDB_PATH), path(processing_dir), path(fastq_file), path(filtlong_file), path(centroids_file), path(minimap_file), path("$processing_dir/${barcode_name}_medaka_output"), emit: data_tuple
    path("${barcode_name}.medaka.consensus.fasta"), emit: medaka_fasta

    publishDir "${run_id}_results/04_polished_sequences/", pattern: "${barcode_name}.medaka.consensus.fasta", mode: 'copy'

    script:
    """
    echo "\$(date '+%Y-%m-%d %H:%M:%S') 🛠️ Running Medaka polishing" | tee -a $processing_dir/processing.log
    medaka_consensus \
        -i $processing_dir/$fastq_file \
        -d $processing_dir/$racon_file \
        -o $processing_dir/${barcode_name}_medaka_output \
        -t $cpu_threads \
        -m "r1041_e82_260bps_sup_g632" 2>> $processing_dir/processing.log
    # Copy for publishing
    cp $processing_dir/${barcode_name}_medaka_output/consensus.fasta ${barcode_name}.medaka.consensus.fasta
    echo "\$(date '+%Y-%m-%d %H:%M:%S') ✅ Polishing complete with Medaka" | tee -a $processing_dir/processing.log
    """
}