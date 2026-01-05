
process its_extraction {
    input:
    val(run_id)
    tuple val(barcode_dir_absolute), val(barcode_name), path(barcode_dir), path(BLASTDB_PATH), path(processing_dir), path(fastq_file), path(filtlong_file), path(centroids_file), path(minimap_file), path(medaka_file)
    val(use_itsx)
    val(cpu_threads)

    output:
    tuple val(barcode_dir_absolute), val(barcode_name), path(barcode_dir), path(BLASTDB_PATH), path(processing_dir), path(fastq_file), path(filtlong_file), path(centroids_file), path(minimap_file), path(medaka_file), path("${barcode_name}.after_itsx.fasta"), emit: data_tuple
    path("${barcode_name}.its.fasta"), emit: its_fasta

    publishDir "${run_id}_results/05_its_extracted/", pattern: "${barcode_name}.its.fasta", mode: 'copy'


    script:
    """
    echo "\$(date '+%Y-%m-%d %H:%M:%S') 🧬 Running ITSx extraction" | tee -a $processing_dir/processing.log
    if [ $use_itsx = 1 ]; then
        mkdir -p ${barcode_name}_itsx_output
        ITSx -i $medaka_file/consensus.fasta -o ${barcode_name}_itsx_output/itsx_output --cpu $cpu_threads
        # clean fasta headers
        sed '/^>/ s/|.*//' ${barcode_name}_itsx_output/itsx_output.full.fasta > ${barcode_name}.after_itsx.fasta
    else
        cp $medaka_file/consensus.fasta ${barcode_name}.after_itsx.fasta
    fi
    # Copy for publishing
    cp ${barcode_name}.after_itsx.fasta ${barcode_name}.its.fasta
    echo "\$(date '+%Y-%m-%d %H:%M:%S') ✅ ITSx extraction complete" | tee -a $processing_dir/processing.log
    """
}