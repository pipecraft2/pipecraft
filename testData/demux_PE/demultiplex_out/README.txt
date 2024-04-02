# Demultiplexing was performed using cutadapt (see 'Core command' below for the used settings).

Files in 'demultiplex_out' directory represent per sample sequence files, that were generated based on the specified indexes file (indexFile_paired.fasta).
index_*fasta file(s) = indexFile_paired.fasta but with added search window size for cutadapt.

Paired-end data, has been demultiplexed taken into account that some sequences
may be also in reverse complementary orientation (two rounds of cutadapt runs, see below).
Output R1 and R2 reads are synchronized for merging paired-end data. 

Files in 'unnamed_index_combinations' directory [if present; only when using dual indexes] represent
index combinations that do not correspond to combinations used in indexes file (indexFile_paired.fasta).

IF SEQUENCE YIELD PER SAMPLE IS LOW (OR ZERO), DOUBLE-CHECK THE INDEXES FORMATTING.

Core commands -> 
Round1: cutadapt -g file:/input/demultiplex_out/index_fwd.fasta -G file:/input/demultiplex_out/index_rev.fasta -e 0 --no-indels --overlap 8 --minimum-length 32 outR1 outR2 inputR1 inputR2
Round2 (RC; R1 and R2 position switched!): cutadapt -g file:index_fwd.fasta -G file:index_rev.fasta -e 0 --no-indels --overlap 8 --minimum-length 32 outR2_round2 outR1_round2 input_for_round2_R2 input_for_round2_R1

Summary of sequence counts in 'seq_count_summary.txt'

Total run time was 3 sec.

##########################################################
###Third-party applications used for this process [PLEASE CITE]:
#cutadapt v4.4 for demultiplexing
    #citation: Martin, Marcel (2011) Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, 17(1), 10-12.
    #https://cutadapt.readthedocs.io/en/stable/index.html
#seqkit v2.3.0 for validating indexes file and adjusting sample names
    #citation: Shen W, Le S, Li Y, Hu F (2016) SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE 11(10): e0163962. https://doi.org/10.1371/journal.pone.0163962
    #https://bioinf.shenwei.me/seqkit/
##################################################################