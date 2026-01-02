nextflow.enable.dsl = 2

params.help = 0
if ( params.help != 0 ) {
    help = """FUNGAL BARCODING WITH ONT: This pipeline streamlines the conversion of Oxford Nanopore Technologies (ONT) basecaller output into high-quality Internal Transcribed Spacer (ITS) sequences.
             |
             |Required arguments:
             |
             |  --FASTQ_DIRECTORY  Location of the directory containing demultiplexed .fastq or .fq files.
             |
             |  --BLASTDB_PATH  Location of the BLAST database directory.
             |
             |  --RUN_ID  Run identifier for this analysis.
             |
             |Optional arguments:
             |
             |  --MEDAKA_MODEL  Medaka inference model. [default: r1041_e82_400bps_hac_variant_v4.3.0]
             |
             |  --USE_ITSX  Set to 0 if you want to ommit extraction of full ITS region using ITSx. [default: 1]
             |
             |  --CHOPPER_MIN_READ_LENGTH Reads shorter than this value wont be used for clusters generation. [default: 150]
             |
             |  --CHOPPER_MAX_READ_LENGTH  Reads longer than this value wont be used for clusters generation. [default: 1000]
             |
             |  --OUTPUT_ALL_POLISHED_SEQS Output all polished sequences even those without UNITE hist (if your are working witn non-ITS sequences). [default: 0]
             |
             |  --REL_ABU_THRESHOLD  Output only clusters with barcode-wise relative abundance above this value. [default: 10]
             |
             |  --CPU_THREADS  Number of CPU threads. [default: 8]
             |
""".stripMargin()

    println(help)
    exit(0)
}
params.help

// Define input barcode directory as a parameter
params.FASTQ_DIRECTORY = null
params.BLASTDB_PATH = null
params.RUN_ID = null
params.MEDAKA_MODEL = 'r1041_e82_400bps_hac_variant_v4.3.0'
params.USE_ITSX = 1
params.OUTPUT_ALL_POLISHED_SEQS = 0

// Chopper min and max length of reads
params.CHOPPER_MIN_READ_LENGTH = 150
params.CHOPPER_MAX_READ_LENGTH = 1000

// Rel abu threshold on final table [0-100]
params.REL_ABU_THRESHOLD = 10

// CPU threads
params.CPU_THREADS = 8

include { 
    ont_barcode_workflow
} from './subworkflows/ont_barcode_workflow.nf'

include {
    create_final_table
} from './modules/create_final_table/main.nf'

workflow {
    
    // Validate required parameters
    if (!params.FASTQ_DIRECTORY) {
        error "ERROR: --FASTQ_DIRECTORY parameter is required!"
    }
    if (!params.BLASTDB_PATH) {
        error "ERROR: --BLASTDB_PATH parameter is required!"
    }
    if (!params.RUN_ID) {
        error "ERROR: --RUN_ID parameter is required!"
    }
    
    def fastq_dir = params.FASTQ_DIRECTORY
    
    // Find all .fastq and .fq files in the directory
    def cmd = ["bash", "-c", "find ${fastq_dir} -maxdepth 1 -type f \\( -name '*.fastq' -o -name '*.fq' \\)"]
    def file_list = cmd.execute().text.readLines()
    file_list.each { println "📂 Found: $it" }
    
    def ch_barcodes = Channel.from(file_list).map { fastq_file ->
        def file_obj = file(fastq_file)
        def barcode_name = file_obj.getSimpleName()  // Get filename without extension
        def fastq_fullpath = file_obj.toAbsolutePath()

        tuple(
            fastq_fullpath, barcode_name, file(fastq_file), file(params.BLASTDB_PATH)
        )
    }

    barcode_results = ont_barcode_workflow(
        params.RUN_ID,
        params.MEDAKA_MODEL,
        params.USE_ITSX,
        params.OUTPUT_ALL_POLISHED_SEQS,
        params.CHOPPER_MIN_READ_LENGTH,
        params.CHOPPER_MAX_READ_LENGTH,
        params.CPU_THREADS,
        ch_barcodes
    ).collect()

    // Use absolute path to the script since we're running from /sequences
    script_path = file("$projectDir/results_aggregation_script.py")
    create_final_table(script_path, barcode_results, params.RUN_ID, params.USE_ITSX, params.REL_ABU_THRESHOLD)

}

