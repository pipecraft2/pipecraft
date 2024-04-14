import Vue from "vue";
import Vuex from "vuex";
import router from "../router/index.js";
var _ = require("lodash");

Vue.use(Vuex);

export default new Vuex.Store({
  state: {
    dockerStatus: "",
    Qcheck: {
      fileExtension: "",
      folderPath: "",
      reportReady: false,
      reportLoading: false,
      dockerActive: false,
    },
    runInfo: {
      active: false,
      type: null,
      step: null,
      nrOfSteps: null,
      containerID: null,
    },
    pullLoader: {
      active: false,
    },
    pullLoader2: {
      active: false,
    },
    workingDir: "/input",
    inputDir: "",
    data: {
      readType: "",
      fileFormat: "",
      dada2mode: "FORWARD",
      debugger: false,
    },
    env_variables: ["FOO=bar", "BAZ=quux"],
    selectedSteps: [],
    steps: [
      {
        stepName: "demultiplex",
        disabled: "demultiplexed",
        services: [
          {
            tooltip:
              "demultiplex data to per-sample files based on specified index file. Note that for read1 and read2 will get .R1 and .R2 identifiers when demultiplexing paired-end data",
            scriptName: "demux_paired_end_data.sh",
            imageName: "pipecraft/cutadapt:4.4",
            serviceName: "demultiplex",
            selected: false,
            showExtra: false,
            extraInputs: [
              {
                name: "cores",
                value: 1,
                disabled: "never",
                tooltip: "number of cores to use",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
              {
                name: "no_indels",
                value: true,
                disabled: "never",
                tooltip:
                  "do not allow insertions or deletions in the index sequence",
                type: "bool",
              },
            ],
            Inputs: [
              {
                name: "index_file",
                value: "undefined",
                btnName: "select fasta",
                disabled: "never",
                tooltip:
                  "select your fasta formatted indexes file for demultiplexing, where fasta headers are sample names, and sequences are sample specific index or index combination",
                type: "file",
              },
              {
                name: "index_file_example",
                value:
                  "https://pipecraft2-manual.readthedocs.io/en/latest/quicktools.html#indexes-file-example-fasta-formatted",
                disabled: "never",
                type: "link",
                tooltip: "link to PipeCraft2 manual page, index file examples",
              },
              {
                name: "index_mismatch",
                value: 0,
                disabled: "never",
                tooltip: "allowed mismatches during the index search",
                type: "numeric",
                rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
              },
              {
                name: "overlap",
                value: 8,
                disabled: "never",
                tooltip:
                  "number of overlap bases with the index. Recommended overlap is the max length of the index for confident sequence assignments to samples in the indexes file",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
              {
                name: "search_window",
                value: 35,
                disabled: "never",
                tooltip:
                  "the index search window size. The default 35 means that the forward index is searched among the first 35 bp and the reverse index among the last 35 bp. This search restriction prevents random index matches in the middle of the sequence",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
            ],
          },
        ],
      },
      {
        stepName: "reorient",
        disabled: "never",
        services: [
          {
            tooltip: "reorient reads based on specified primer sequences",
            scriptName: "reorient_paired_end_reads.sh",
            imageName: "pipecraft/reorient:1",
            serviceName: "reorient",
            selected: false,
            showExtra: false,
            extraInputs: [],
            Inputs: [
              {
                name: "mismatches",
                value: 1,
                disabled: "never",
                tooltip: "allowed mismatches in the primer search",
                type: "numeric",
                rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
              },
              {
                name: "forward_primers",
                value: [],
                disabled: "never",
                tooltip: "specify forward primer (5'-3'); add up to 13 primers",
                type: "chip",
                iupac: true,
                rules: [(v) => v.length <= 13 || "TOO MANY PRIMERS"],
              },
              {
                name: "reverse_primers",
                value: [],
                disabled: "never",
                tooltip: "specify reverse primer (3'-5'); add up to 13 primers",
                type: "chip",
                iupac: true,
                rules: [(v) => v.length <= 13 || "TOO MANY PRIMERS"],
              },
            ],
          },
        ],
      },
      {
        stepName: "cut primers",
        disabled: "never",
        services: [
          {
            tooltip: "remove primers sequences from the reads",
            scriptName: "cut_primers_paired_end_reads.sh",
            imageName: "pipecraft/cutadapt:4.4",
            serviceName: "cutadapt",
            selected: false,
            showExtra: false,
            extraInputs: [
              {
                name: "cores",
                value: 1,
                disabled: "never",
                tooltip:
                  "number of cores to use. For paired-end data in fasta format, set to 1 [default]. For fastq formats you may set the value to 0 to use all cores",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
              {
                name: "no_indels",
                value: true,
                disabled: "never",
                tooltip:
                  "do not allow insertions or deletions is primer search. Mismatches are the only type of errors accounted in the error rate parameter. ",
                type: "bool",
              },
            ],
            Inputs: [
              {
                name: "forward_primers",
                value: [],
                disabled: "never",
                tooltip: "specify forward primer (5'-3'); add up to 13 primers",
                type: "chip",
                iupac: true,
                rules: [(v) => v.length <= 13 || "TOO MANY PRIMERS"],
              },
              {
                name: "reverse_primers",
                value: [],
                disabled: "never",
                tooltip: "specify reverse primer (3'-5'); add up to 13 primers",
                type: "chip",
                iupac: true,
                rules: [(v) => v.length <= 13 || "TOO MANY PRIMERS"],
              },
              {
                name: "mismatches",
                value: 1,
                disabled: "never",
                tooltip: "allowed mismatches in the primer search",
                type: "numeric",
                rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
              },
              {
                name: "min_overlap",
                value: 21,
                disabled: "never",
                tooltip:
                  "number of overlap bases with the primer sequence. Partial matches are allowed, but short matches may occur by chance, leading to erroneously clipped bases. Specifying higher overlap than the length of primer sequnce will still clip the primer (e.g. primer length is 22 bp, but overlap is specified as 25 - this does not affect the identification and clipping of the primer as long as the match is in the specified mismatch error range)",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
              {
                name: "seqs_to_keep",
                items: ["keep_all", "keep_only_linked"],
                value: "keep_all",
                disabled: "never",
                tooltip:
                  "'keep_all' keeps the sequences where at least one primer was found (fwd or rev); recommended when cutting primers from paired-end data (unassembled), where individual R1 or R2 read lenghts are shorther than the expected amplicon length. 'keep_only_linked' = keep sequences if primers are found in both ends (fwd…rev); discards the read if both primers were not found in this read; maybe useful for single-end data",
                type: "select",
              },
              {
                name: "pair_filter",
                items: ["both", "any"],
                value: "both",
                disabled: "single_end",
                tooltip:
                  "applies only for paired-end data. Option 'both' discards a read-pair when both reads do not contain a primer sequence. Option 'any' discards the read-pair when one of the reads does not contain a primer sequence",
                type: "select",
              },
            ],
          },
        ],
      },
      {
        stepName: "quality filtering",
        disabled: "never",
        services: [
          {
            tooltip: "quality filtering with vsearch",
            scriptName: "quality_filtering_paired_end_vsearch.sh",
            imageName: "pipecraft/vsearch_dada2:2",
            serviceName: "vsearch",
            selected: false,
            showExtra: false,
            extraInputs: [
              {
                name: "max_length",
                value: null,
                disabled: "never",
                tooltip:
                  "discard sequences with more than the specified number of bases. Note that if 'trunc length' setting is specified, then 'max length' SHOULD NOT be lower than 'trunc lenght' (otherwise all reads are discared) [empty field = no action taken]",
                type: "numeric",
                rules: [
                  (v) =>
                    (v >= 1) | (v == "") ||
                    "ERROR: specify values >= 1 or leave it empty (=no action taken)",
                ],
              },
              {
                name: "qmin",
                value: 0,
                disabled: "never",
                tooltip:
                  "the minimum quality score accepted for FASTQ files. The default is 0, which is usual for recent Sanger/Illumina 1.8+ files. Older formats may use scores between -5 and 2",
                type: "numeric",
                rules: [(v) => v >= -5 || "ERROR: specify values >= -5"],
              },
              {
                name: "maxee_rate",
                value: null,
                disabled: "never",
                tooltip:
                  "discard sequences with more than the specified number of expected errors per base (empty field = no action taken)",
                type: "numeric",
                rules: [
                  (v) =>
                    (v >= 0.001) | (v == "") ||
                    "ERROR: specify values >=0.001 or leave it empty (= no action taken)",
                ],
              },
              {
                name: "truncqual",
                value: null,
                disabled: "never",
                tooltip:
                  "tuncate sequences starting from the first base with the specified base quality score value or lower (empty field = no action taken)",
                type: "numeric",
                rules: [
                  (v) =>
                    (v >= 1) | (v == "") ||
                    "ERROR: specify values >=1 or leave it empty (= no action taken)",
                ],
              },
              {
                name: "truncee",
                value: null,
                disabled: "never",
                tooltip:
                  "runcate sequences so that their total expected error is not higher than the specified value (empty field = no action taken)",
                type: "numeric",
                rules: [
                  (v) =>
                    (v >= 0.001) | (v == "") ||
                    "ERROR: specify values >=0.001 or leave it empty (= no action taken)",
                ],
              },
              {
                name: "cores",
                value: 1,
                disabled: "never",
                tooltip: "number of cores to use",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
            ],
            Inputs: [
              {
                name: "maxee",
                value: 1,
                disabled: "never",
                tooltip:
                  "maximum number of expected errors per sequence. Sequences with higher error rates will be discarded",
                type: "numeric",
                rules: [(v) => v >= 0.001 || "ERROR: specify values >=0.001"],
              },
              {
                name: "maxNs",
                value: 0,
                disabled: "never",
                tooltip:
                  "discard sequences with more than the specified number of Ns",
                type: "numeric",
                rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
              },
              {
                name: "min_length",
                value: 32,
                disabled: "never",
                tooltip:
                  "minimum length of the filtered output sequence. Note that if 'trunc length' setting is specified, then 'min length' SHOULD BE lower than 'trunc lenght' (otherwise all reads are discared)",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
              {
                name: "trunc_length",
                value: null,
                disabled: "never",
                tooltip:
                  "truncate sequences to the specified length. Shorter sequences are discarded; thus if specified, check that 'min length' setting is lower than 'trunc length' ('min lenght' therefore has basically no effect) [empty field = no action taken]",
                type: "numeric",
                rules: [
                  (v) =>
                    (v >= 5) | (v == "") ||
                    "ERROR: specify values >= 5 or leave it empty (=no action taken)",
                ],
              },
              {
                name: "qmax",
                value: 41,
                disabled: "never",
                tooltip:
                  "specify the maximum quality score accepted when reading FASTQ files. The default is 41, which is usual for recent Sanger/Illumina 1.8+ files. For PacBio data use 93",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
            ],
          },
          {
            tooltip: "quality filtering with trimmomatic",
            scriptName: "quality_filtering_paired_end_trimmomatic.sh",
            imageName: "pipecraft/trimmomatic:0.39",
            serviceName: "trimmomatic",
            selected: false,
            showExtra: false,
            extraInputs: [
              {
                name: "leading_qual_threshold",
                value: null,
                disabled: "never",
                tooltip:
                  "quality score threshold to remove low quality bases from the beginning of the read. As long as a base has a value below this threshold the base is removed and the next base will be investigated (empty field = no action taken)",
                type: "numeric",
                rules: [
                  (v) =>
                    (v >= 1) | (v == "") ||
                    "ERROR: specify values >= 1 or leave it empty (=no action taken)",
                ],
              },
              {
                name: "trailing_qual_threshold",
                value: null,
                disabled: "never",
                tooltip:
                  "quality score threshold to remove low quality bases from the end of the read. As long as a base has a value below this threshold the base is removed and the next base will be investigated (empty field = no action taken)",
                type: "numeric",
                rules: [
                  (v) =>
                    (v >= 1) | (v == "") ||
                    "ERROR: specify values >= 1 or leave it empty (=no action taken)",
                ],
              },
              {
                name: "cores",
                value: 1,
                disabled: "never",
                tooltip: "number of cores to use",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
              {
                name: "phred",
                items: [33, 64],
                value: 33,
                disabled: "never",
                tooltip:
                  "phred quality scored encoding. Use phred64 if working with data from older Illumina (Solexa) machines",
                type: "select",
              },
            ],
            Inputs: [
              {
                name: "window_size",
                value: 5,
                disabled: "never",
                tooltip:
                  "the number of bases to average base qualities. Starts scanning at the 5'-end of a sequence and trimms the read once the average required quality (required_qual) within the window size falls below the threshold",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
              {
                name: "required_quality",
                value: 27,
                disabled: "never",
                tooltip:
                  "the average quality required for selected window size",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
              {
                name: "min_length",
                value: 32,
                disabled: "never",
                tooltip: "minimum length of the filtered output sequence",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
            ],
          },
          {
            tooltip: "quality filtering with fastp",
            scriptName: "quality_filtering_paired_end_fastp.sh",
            imageName: "pipecraft/fastp:0.23.2",
            serviceName: "fastp",
            selected: false,
            showExtra: false,
            extraInputs: [
              {
                name: "trim_polyG",
                value: 12,
                disabled: "never",
                tooltip:
                  "specify the minimum length for polyG tail trimming. Useful for Illumina NextSeq/NovaSeq data, where read tails contain multiple Gs, which means no signal in the Illumina two-color systems (empty field = disable_trim_poly_g)",
                type: "numeric",
                rules: [
                  (v) =>
                    (v >= 1) | (v == "") ||
                    "ERROR: specify values >= 1 or leave it empty (=> --disable_trim_poly_g)",
                ],
              },
              {
                name: "trim_polyX",
                value: null,
                disabled: "never",
                tooltip:
                  "specify the minimum length to trim polyX (i.e. multiple A, T or Cs) in the read tail (empty field = no action taken)",
                type: "numeric",
                rules: [
                  (v) =>
                    (v >= 1) | (v == "") ||
                    "ERROR: specify values >= 1 or leave it empty (= no action)",
                ],
              },
              {
                name: "max_length",
                value: 0,
                disabled: "never",
                tooltip:
                  "reads longer than 'max length' will be discarded, default 0 means no limitation",
                type: "numeric",
                rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
              },
              {
                name: "trunc_length",
                value: 0,
                disabled: "never",
                tooltip:
                  "truncate sequences to specified length. Shorter sequences are discarded; thus check that 'min length' setting is lower than 'trunc length'",
                type: "numeric",
                rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
              },
              {
                name: "aver_qual",
                value: 0,
                disabled: "never",
                tooltip:
                  "if one read's average quality score <'aver_qual', then this read/pair is discarded. Default 0 means no requirement",
                type: "numeric",
                rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
              },
              {
                name: "low_complexity_filter",
                value: null,
                disabled: "never",
                tooltip:
                  "enable low complexity filter and specify the threshold for low complexity filter. The complexity is defined as the percentage of base that is different from its next base (base[i] != base[i+1]). E.g. vaule 30 means then 30% complexity is required. Not specified = filter not applied",
                type: "numeric",
                rules: [
                  (v) =>
                    (v >= 1) | (v == "") ||
                    "ERROR: specify values >= 1 or leave it empty (=no action taken)",
                ],
              },
              {
                name: "cores",
                value: 4,
                disabled: "never",
                tooltip: "number of cores to use",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
            ],
            Inputs: [
              {
                name: "window_size",
                value: 4,
                disabled: "never",
                tooltip: "the window size for calculating mean quality",
                type: "numeric",
                rules: [
                  (v) =>
                    (v >= 1) | (v <= 1000) ||
                    "ERROR: specify values in range 1-1000",
                ],
              },
              {
                name: "required_qual",
                value: 27,
                disabled: "never",
                tooltip:
                  "the mean quality requirement per sliding window (window_size)",
                type: "numeric",
                rules: [
                  (v) =>
                    (v >= 1) | (v <= 36) ||
                    "ERROR: specify values in range 1-36",
                ],
              },
              {
                name: "min_qual",
                value: 15,
                disabled: "never",
                tooltip:
                  "the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
              {
                name: "min_qual_thresh",
                value: 40,
                disabled: "never",
                tooltip:
                  "how many percents of bases are allowed to be unqualified (0-100)",
                type: "numeric",
                rules: [
                  (v) =>
                    (v >= 0) | (v <= 100) ||
                    "ERROR: specify values in range 0-100",
                ],
              },
              {
                name: "maxNs",
                value: 0,
                disabled: "never",
                tooltip:
                  "discard sequences with more than the specified number of Ns",
                type: "numeric",
                rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
              },
              {
                name: "min_length",
                value: 32,
                disabled: "never",
                tooltip:
                  "minimum length of the filtered output sequence. Shorter sequences are discarded",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
            ],
          },
          {
            tooltip: "quality filtering with DADA2 'filterAndTrim' function",
            scriptName: "quality_filtering_paired_end_dada2.sh",
            imageName: "pipecraft/vsearch_dada2:2",
            serviceName: "DADA2",
            selected: false,
            showExtra: false,
            extraInputs: [
              {
                name: "truncQ",
                value: 2,
                disabled: "never",
                tooltip:
                  "truncate reads at the first instance of a quality score less than or equal to truncQ",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
              {
                name: "truncLen",
                value: 0,
                disabled: "never",
                tooltip:
                  "truncate reads after truncLen bases (applies to R1 reads when working with paired-end data). Reads shorter than this are discarded. Explore quality profiles (with QualityCheck module) and see whether poor quality ends needs to be truncated",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
              {
                name: "truncLen_R2",
                value: 0,
                disabled: "single_end",
                tooltip:
                  "applies only for paired-end data. Truncate R2 reads after truncLen bases. Reads shorter than this are discarded. Explore quality profiles (with QualityCheck module) and see whether poor quality ends needs to be truncated",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
              {
                name: "maxLen",
                value: 9999,
                disabled: "never",
                tooltip:
                  "remove reads with length greater than maxLen. maxLen is enforced on the raw reads. In dada2, the default = Inf, but here set as 9999",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
              {
                name: "minQ",
                value: 0,
                disabled: "never",
                tooltip:
                  "after truncation, reads contain a quality score below minQ will be discarded",
                type: "numeric",
                rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
              },
              {
                name: "matchIDs",
                value: true,
                disabled: "single_end",
                tooltip:
                  "applies only for paired-end data. If TRUE, then double-checking (with seqkit pair) that only paired reads that share ids are outputted",
                type: "bool",
              },
            ],
            Inputs: [
              {
                name: "read_R1",
                value: ["_R1"],
                disabled: "single_end",
                tooltip:
                  "applies only for paired-end data. Identifyer string that is common for all R1 reads. R1 reads must contain R1 strings in the file names; sample names must not contain R1! (e.g. when all R1 files have '.R1' string, then enter '\\.R1'. Note that backslash is only needed to escape dot regex; e.g. when all R1 files have '_R1' string, then enter '_R1'.)'",
                type: "chip",
                rules: [(v) => v.length <= 1 || "ADD ONLY ONE IDENTIFIER"],
              },
              {
                name: "read_R2",
                value: ["_R2"],
                disabled: "single_end",
                tooltip:
                  "applies only for paired-end data. Identifyer string that is common for all R2 reads. R2 reads must contain R2 strings in the file names; sample names must not contain R2! (e.g. when all R2 files have '.R2' string, then enter '\\.R2'. Note that backslash is only needed to escape dot regex; e.g. when all R2 files have '_R1' string, then enter '_R2'.)",
                type: "chip",
                rules: [(v) => v.length <= 1 || "ADD ONLY ONE IDENTIFIER"],
              },
              // {
              //   name: "samp_ID",
              //   value: ["_"],
              //   disabled: "single_end",
              //   tooltip:
              //     "applies only for paired-end data. Identifyer string that separates the sample name from redundant charachters (e.g. file name = sample1.R1.fastq, then '\\.' would be the 'identifier string' (sample name = sample1)); note that backslash is only needed to escape dot regex (e.g. when file name = sample1_R1.fastq then specify as '_')",
              //   type: "chip",
              //   rules: [(v) => v.length <= 1 || "ADD ONLY ONE IDENTIFIER"],
              // },
              {
                name: "maxEE",
                value: 2,
                disabled: "never",
                tooltip:
                  "discard sequences with more than the specified number of expected errors",
                type: "numeric",
                rules: [(v) => v >= 0.001 || "ERROR: specify values >=0.001"],
              },
              {
                name: "maxN",
                value: 0,
                disabled: "never",
                tooltip:
                  "discard sequences with more than the specified number of N’s (ambiguous bases)",
                type: "numeric",
                rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
              },
              {
                name: "minLen",
                value: 20,
                disabled: "never",
                tooltip:
                  "remove reads with length less than minLen. minLen is enforced after all other trimming and truncation",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
            ],
          },
        ],
      },
      {
        stepName: "assemble paired-end",
        disabled: "single_end",
        services: [
          {
            tooltip: "assemble paired-end reads with vsearch",
            scriptName: "assemble_paired_end_data_vsearch.sh",
            imageName: "pipecraft/vsearch_dada2:2",
            serviceName: "vsearch",
            disabled: "single_end",
            selected: false,
            showExtra: false,
            extraInputs: [
              {
                name: "max_diffs",
                value: 20,
                disabled: "never",
                tooltip:
                  "the maximum number of non-matching nucleotides allowed in the overlap region",
                type: "numeric",
                rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
              },
              {
                name: "max_Ns",
                value: 0,
                disabled: "never",
                tooltip:
                  "discard sequences with more than the specified number of N’s",
                type: "numeric",
                rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
              },
              {
                name: "max_len",
                value: 600,
                disabled: "never",
                tooltip: "maximum length of the merged sequence",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
              {
                name: "keep_disjointed",
                value: false,
                disabled: "never",
                tooltip:
                  "output reads that were not merged into separate FASTQ files",
                type: "bool",
              },
              {
                name: "fastq_qmax",
                value: 41,
                disabled: "never",
                tooltip:
                  "maximum quality score accepted when reading FASTQ files. The default is 41, which is usual for recent Sanger/Illumina 1.8+ files",
                type: "numeric",
                rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
              },
            ],
            Inputs: [
              {
                name: "read_R1",
                value: ["_R1"],
                disabled: "single_end",
                tooltip:
                  "identifyer string that is common for all R1 reads. R1/R2 reads must contain R1/R2 strings in the file names; sample names must not contain R1/R2! (e.g. when all R1 files have '.R1' string, then enter '\\.R1'. Note that backslash is only needed to escape dot regex; e.g. when all R1 files have '_R1' string, then enter '_R1'.)",
                type: "chip",
                rules: [(v) => v.length <= 1 || "ADD ONLY ONE IDENTIFIER"],
              },
              {
                name: "min_overlap",
                value: 12,
                disabled: "never",
                tooltip: "minimum overlap between the merged reads",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
              {
                name: "min_lenght",
                value: 32,
                disabled: "never",
                tooltip: "minimum length of the merged sequence",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
              {
                name: "allow_merge_stagger",
                value: true,
                disabled: "never",
                tooltip:
                  "allow to merge staggered read pairs. Staggered pairs are pairs where the 3' end of the reverse read has an overhang to the left of the 5' end of the forward read. This situation can occur when a very short fragment is sequenced",
                type: "bool",
              },
              {
                name: "include_only_R1",
                value: false,
                disabled: "never",
                tooltip:
                  "Include unassembled R1 reads to the set of assembled reads per sample. This may be relevant when working with e.g. ITS2 sequences, because the ITS2 region in some taxa is too long for assembly, therefore discarded completely after assembly process. Thus, including also unassembled R1 reads, partial ITS2 sequences for these taxa will be represented in the final output",
                type: "bool",
              },
            ],
          },
          {
            tooltip:
              "denoise and assemble paired-end reads with DADA2 'mergePairs' and 'dada' functions. Note that only FASTA is outputted!",
            scriptName: "assemble_paired_end_data_dada2.sh",
            imageName: "pipecraft/vsearch_dada2:2",
            serviceName: "DADA2_denoise_and_merge",
            selected: false,
            disabled: "single_end",
            showExtra: false,
            extraInputs: [
              {
                name: "pool",
                items: ["FALSE", "TRUE", "psuedo"],
                value: "FALSE",
                disabled: "never",
                tooltip:
                  "Denoising option. If pool = TRUE, the algorithm will pool together all samples prior to sample inference. Pooling improves the detection of rare variants, but is computationally more expensive. If pool = 'pseudo', the algorithm will perform pseudo-pooling between individually processed samples. This argument has no effect if only 1 sample is provided, and pool does not affect error rates, which are always estimated from pooled observations across samples",
                type: "select",
              },
              // {
              //   name: "selfConsist",
              //   disabled: "never",
              //   value: false,
              //   tooltip:
              //     "Denoising option. If selfConsist = TRUE, the algorithm will alternate between sample inference and error rate estimation until convergence",
              //   type: "bool",
              // },
              {
                name: "qualityType",
                items: ["Auto", "FastqQuality"],
                value: "Auto",
                disabled: "never",
                tooltip:
                  "Auto means to attempt to auto-detect the fastq quality encoding. This may fail for PacBio files with uniformly high quality scores, in which case use 'FastqQuality'",
                type: "select",
              },
              {
                name: "BAND_SIZE",
                value: 16,
                disabled: "never",
                tooltip:
                  "default = 16. Banding for Needleman-Wunsch alignments. Default sets here = 16 for Illumina and 32 for PacBio",
                type: "numeric",
                rules: [(v) => v >= -1 || "ERROR: specify values >= -1"],
              },
              {
                name: "OMEGA_A",
                value:
                  (0.0000000000000000000000000000000000000001).toExponential(),
                disabled: "never",
                tooltip:
                  "default = 1e-40. Denoising setting; see DADA2 'setDadaOpt()' for detalis. Default value  is a conservative setting to avoid making false positive inferences, but comes at the cost of reducing the ability to identify some rare variants",
                type: "numeric",
                rules: [(v) => v >= 0 || "ERROR: specify values > 0"],
              },
              {
                name: "OMEGA_P",
                value: (0.0001).toExponential(),
                disabled: "never",
                tooltip:
                  "default = 1e-4 (0.0001). Denoising setting; see DADA2 'setDadaOpt()' for detalis",
                type: "numeric",
                rules: [(v) => v >= 0 || "ERROR: specify values > 0"],
              },
              {
                name: "OMEGA_C",
                value:
                  (0.0000000000000000000000000000000000000001).toExponential(),
                disabled: "never",
                tooltip:
                  "default = 1e-40. Denoising setting; see DADA2 'setDadaOpt()' for detalis",
                type: "numeric",
                rules: [(v) => v >= 0 || "ERROR: specify values > 0"],
              },
              {
                name: "DETECT_SINGLETONS",
                disabled: "never",
                value: false,
                tooltip:
                  "Default = FALSE. Denoising setting; see DADA2 'setDadaOpt()' for detalis. If set to TRUE, this removes the requirement for at least two reads with the same sequences to exist in order for a new ASV to be detected",
                type: "bool",
              },
            ],
            Inputs: [
              {
                name: "read_R1",
                value: ["_R1"],
                disabled: "single_end",
                tooltip:
                  "identifyer string that is common for all R1 reads. R1/R2 reads must contain R1/R2 strings in the file names; sample names must not contain R1/R2! (e.g. when all R1 files have '.R1' string, then enter '\\.R1'. Note that backslash is only needed to escape dot regex; e.g. when all R1 files have '_R1' string, then enter '_R1'.)",
                type: "chip",
                rules: [(v) => v.length <= 1 || "ADD ONLY ONE IDENTIFIER"],
              },
              {
                name: "read_R2",
                value: ["_R2"],
                disabled: "single_end",
                tooltip:
                  "identifyer string that is common for all R2 reads (e.g. when all R2 files have '.R2' string, then enter '\\.R2'. Note that backslash is only needed to escape dot regex; e.g. when all R2 files have '_R1' string, then enter '_R2'.)",
                type: "chip",
                rules: [(v) => v.length <= 1 || "ADD ONLY ONE IDENTIFIER"],
              },
              // {
              //   name: "samp_ID",
              //   value: ["_"],
              //   disabled: "never",
              //   tooltip:
              //     "identifyer string that separates the sample name from redundant charachters (e.g. file name = sample1.R1.fastq, then '\\.' would be the 'identifier string' (sample name = sample1)); note that backslash is only needed to escape dot regex (e.g. when file name = sample1_R1.fastq then specify as '_')",
              //   type: "chip",
              //   rules: [(v) => v.length <= 1 || "ADD ONLY ONE IDENTIFIER"],
              // },
              {
                name: "minOverlap",
                value: 12,
                disabled: "never",
                tooltip:
                  "the minimum length of the overlap required for merging the forward and reverse reads",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
              {
                name: "maxMismatch",
                value: 0,
                disabled: "never",
                tooltip: "the maximum mismatches allowed in the overlap region",
                type: "numeric",
                rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
              },
              {
                name: "trimOverhang",
                value: false,
                disabled: "never",
                tooltip:
                  "if TRUE, overhangs in the alignment between the forwards and reverse read are trimmed off. Overhangs are when the reverse read extends past the start of the forward read, and vice-versa, as can happen when reads are longer than the amplicon and read into the other-direction primer region",
                type: "bool",
              },
              {
                name: "justConcatenate",
                value: false,
                disabled: "never",
                tooltip:
                  "if TRUE, the forward and reverse-complemented reverse read are concatenated rather than merged, with a NNNNNNNNNN (10 Ns) spacer inserted between them",
                type: "bool",
              },
            ],
          },
        ],
      },
      {
        stepName: "chimera filtering",
        disabled: "never",
        services: [
          {
            tooltip:
              "tick the checkbox to filter chimeras with vsearch --uchime_denovo",
            scriptName: "chimera_filtering_vsearch.sh",
            imageName: "pipecraft/vsearch_dada2:2",
            serviceName: "uchime_denovo",
            selected: false,
            showExtra: false,
            extraInputs: [
              {
                name: "cores",
                value: 4,
                disabled: "never",
                tooltip:
                  "Number of cores to use (only for reference based chimera filtering)",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
              {
                name: "abundance_skew",
                value: 2,
                disabled: "never",
                tooltip:
                  "The abundance skew is used to distinguish in a threeway alignment which sequence is the chimera and which are the parents. The assumption is that chimeras appear later in the PCR amplification process and are therefore less abundant than their parents. The default value is 2.0, which means that the parents should be at least 2 times more abundant than their chimera. Any positive value equal or greater than 1.0 can be used",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
              {
                name: "min_h",
                value: 0.28,
                disabled: "never",
                tooltip:
                  "Minimum score (h). Increasing this value tends to reduce the number of false positives and to decrease sensitivity. Values ranging from 0.0 to 1.0 included are accepted",
                max: 1,
                min: 0,
                step: 0.01,
                type: "slide",
              },
            ],
            Inputs: [
              {
                name: "pre_cluster",
                value: 0.98,
                disabled: "never",
                tooltip:
                  "Identity percentage when performing 'pre-clustering' with --cluster_size for denovo chimera filtering with --uchime_denovo",
                max: 1,
                min: 0,
                step: 0.01,
                type: "slide",
              },
              {
                name: "min_unique_size",
                value: 1,
                disabled: "never",
                tooltip:
                  "Minimum amount of a unique sequences in a fasta file. If value = 1, then no sequences are discarded after dereplication; if value = 2, then sequences, which are represented only once in a given file are discarded; and so on",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
              {
                name: "denovo",
                value: true,
                disabled: "never",
                tooltip:
                  "Perform denovo chimera filtering with --uchime_denovo",
                type: "bool",
              },
              {
                name: "reference_based",
                active: false,
                btnName: "select file",
                value: "undefined",
                disabled: "never",
                tooltip:
                  "Perform reference database based chimera filtering with --uchime_ref. If denovo = TRUE, then reference based chimera filtering will be performed after denovo",
                type: "boolfile",
              },
            ],
          },
          {
            tooltip:
              "tick the checkbox to filter chimeras with vsearch --uchime3_denovo [for denoised sequences]",
            scriptName: "chimera_filtering_vsearch_uchime3.sh",
            imageName: "pipecraft/vsearch_dada2:2",
            serviceName: "uchime3_denovo",
            selected: false,
            showExtra: false,
            extraInputs: [
              {
                name: "cores",
                value: 4,
                disabled: "never",
                tooltip:
                  "Number of cores to use (only for reference based chimera filtering)",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
              {
                name: "abundance_skew",
                value: 16,
                disabled: "never",
                tooltip:
                  "The abundance skew is used to distinguish in a threeway alignment which sequence is the chimera and which are the parents. The assumption is that chimeras appear later in the PCR amplification process and are therefore less abundant than their parents",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
            ],
            Inputs: [
              {
                name: "pre_cluster",
                value: 0.98,
                disabled: "never",
                tooltip:
                  "Identity percentage when performing 'pre-clustering' with --cluster_size for denovo chimera filtering with --uchime_denovo",
                max: 1,
                min: 0,
                step: 0.01,
                type: "slide",
              },
              {
                name: "min_unique_size",
                value: 1,
                disabled: "never",
                tooltip:
                  "Minimum amount of a unique sequences in a fasta file. If value = 1, then no sequences are discarded after dereplication; if value = 2, then sequences, which are represented only once in a given file are discarded; and so on",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
              {
                name: "denovo",
                value: true,
                disabled: "never",
                tooltip:
                  "Perform denovo chimera filtering with --uchime_denovo",
                type: "bool",
              },
              {
                name: "reference_based",
                active: false,
                btnName: "select file",
                value: "undefined",
                disabled: "never",
                tooltip:
                  "Perform reference database based chimera filtering with --uchime_ref. If denovo = TRUE, then reference based chimera filtering will be performed after denovo",
                type: "boolfile",
              },
            ],
          },
        ],
      },
      {
        stepName: "ITS Extractor",
        disabled: "never",
        services: [
          {
            tooltip:
              "if data set consists of ITS sequences; identify and extract the ITS regions using ITSx",
            scriptName: "ITS_extractor.sh",
            imageName: "pipecraft/itsx:1.1.3",
            serviceName: "itsx",
            selected: false,
            showExtra: false,
            extraInputs: [
              {
                name: "cores",
                value: 6,
                disabled: "never",
                tooltip: "number of cores to use",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
              {
                name: "e_value",
                value: (0.01).toExponential(),
                disabled: "never",
                tooltip:
                  "domain E-value cutoff a sequence must obtain in the HMMER-based step to be included in the output. Here, the defaul 1e-2 = 0.01 (more relaxed compared with the default ITSx [1e-5])",
                type: "numeric",
                rules: [(v) => v >= 0 || "ERROR: specify only values > 0"],
              },
              {
                name: "scores",
                value: 0,
                disabled: "never",
                tooltip:
                  "domain score cutoff that a sequence must obtain in the HMMER-based step to be included in the output. Leave as default if unsure how to set",
                type: "numeric",
                rules: [(v) => v >= 0 || "ERROR: specify only values > 0"],
              },
              {
                name: "domains",
                value: 2,
                disabled: "never",
                tooltip:
                  "the minimum number of domains (different HMM gene profiles) that must match a sequence for it to be included in the output (detected as an ITS sequence). Setting the value lower than two will increase the number of false positives, while increasing it above two will decrease ITSx detection abilities on fragmentary data",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
              {
                name: "complement",
                value: true,
                disabled: "never",
                tooltip:
                  "if TRUE, then ITSx checks also reverse complementary strands for matches to HMM-profiles",
                type: "bool",
              },
              {
                name: "only_full",
                value: false,
                disabled: "never",
                tooltip:
                  "if TRUE, the output is limited to full-length ITS1 and ITS2 regions only",
                type: "bool",
              },
              {
                name: "truncate",
                value: true,
                disabled: "never",
                tooltip:
                  "if TRUE, ITSx removes ends of ITS sequences if they are outside of the ITS region. If off, the whole input sequence is saved when ITS region is detected",
                type: "bool",
              },
            ],
            Inputs: [
              {
                name: "organisms",
                items: [
                  "all",
                  "alveolata",
                  "bryophyta",
                  "bacillariophyta",
                  "amoebozoa",
                  "euglenozoa",
                  "fungi",
                  "chlorophyta",
                  "rhodophyta",
                  "phaeophyceae",
                  "marchantiophyta",
                  "metazoa",
                  "oomycota",
                  "haptophyceae",
                  "raphidophyceae",
                  "rhizaria",
                  "synurophyceae",
                  "tracheophyta",
                  "eustigmatophyceae",
                  "apusozoa",
                  "parabasalia",
                ],
                value: ["fungi"],
                disabled: "never",
                tooltip:
                  "set of profiles to use for the search. Can be used to restrict the search to only a few organism groups types to save time, if one or more of the origins are not relevant to the dataset under study",
                type: "combobox",
              },
              {
                name: "regions",
                items: ["all", "SSU", "ITS1", "5.8S", "ITS2", "LSU"],
                value: ["all"],
                disabled: "never",
                tooltip:
                  "ITS regions to output (note that 'all' will output also full ITS region [ITS1-5.8S-ITS2])",
                type: "combobox",
              },
              {
                name: "partial",
                value: 50,
                disabled: "never",
                tooltip:
                  "if larger than 0, ITSx will save additional FASTA-files for full and partial ITS sequences longer than the specified cutoff value. If his setting is left to 0 (zero), it means OFF",
                type: "numeric",
                rules: [(v) => v >= 0 || "ERROR: specify only values >= 0"],
              },
            ],
          },
        ],
      },
      {
        stepName: "clustering",
        disabled: "never",
        services: [
          {
            scriptName: "clustering_vsearch.sh",
            tooltip: "tick the checkbox to cluster reads with vsearch",
            imageName: "pipecraft/vsearch_dada2:2",
            serviceName: "vsearch",
            selected: false,
            showExtra: false,
            extraInputs: [
              {
                name: "similarity_type",
                items: ["0", "1", "2", "3", "4"],
                value: "2",
                disabled: "never",
                tooltip: "pairwise sequence identity definition (--iddef)",
                type: "select",
              },
              {
                name: "sequence_sorting",
                items: ["size", "length", "no"],
                value: "size",
                disabled: "never",
                tooltip:
                  'size = sort the sequences by decreasing abundance; "length" = sort the sequences by decreasing length (--cluster_fast); "no" = do not sort sequences (--cluster_smallmem --usersort)',
                type: "select",
              },
              {
                name: "centroid_type",
                items: ["similarity", "abundance"],
                value: "similarity",
                disabled: "never",
                tooltip:
                  '"similarity" = assign representative sequence to the closest (most similar) centroid (distance-based greedy clustering); "abundance" = assign representative sequence to the most abundant centroid (abundance-based greedy clustering; --sizeorder), --maxaccepts should be > 1',
                type: "select",
              },
              {
                name: "maxaccepts",
                value: 1,
                disabled: "never",
                tooltip:
                  "maximum number of hits to accept before stopping the search (should be > 1 for abundance-based selection of centroids [centroid type])",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
              {
                name: "mask",
                items: ["dust", "none"],
                value: "dust",
                disabled: "never",
                tooltip:
                  'mask regions in sequences using the "dust" method, or do not mask ("none").',
                type: "select",
              },
              {
                name: "cores",
                value: 4,
                disabled: "never",
                tooltip: "number of cores to use",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
            ],
            Inputs: [
              {
                name: "OTU_type",
                items: ["centroid", "consout"],
                disabled: "never",
                tooltip:
                  '"centroid" = output centroid sequences; "consout" = output consensus sequences',
                value: "centroid",
                type: "select",
              },
              {
                name: "similarity_threshold",
                value: 0.97,
                disabled: "never",
                tooltip:
                  "define OTUs based on the sequence similarity threshold; 0.97 = 97% similarity threshold",
                max: 1,
                min: 0,
                step: 0.01,
                type: "slide",
              },
              {
                name: "strands",
                items: ["both", "plus"],
                disabled: "never",
                tooltip:
                  "when comparing sequences with the cluster seed, check both strands (forward and reverse complementary) or the plus strand only",
                value: "both",
                type: "select",
              },
              {
                name: "remove_singletons",
                value: false,
                disabled: "never",
                tooltip:
                  "remove singleton OTUs (e.g., if TRUE, then OTUs with only one sequence will be discarded)",
                type: "bool",
              },
            ],
          },
          {
            scriptName: "clustering_unoise.sh",
            tooltip:
              "tick the checkbox to cluster reads with vsearch --cluster_unoise (and optionally remove chimeras with --uchime3_denovo)",
            imageName: "pipecraft/vsearch_dada2:2",
            serviceName: "unoise3",
            selected: false,
            showExtra: false,
            extraInputs: [
              {
                name: "unoise_alpha",
                value: 2.0,
                disabled: "never",
                tooltip:
                  "default = 2.0. alpha parameter to the vsearch --cluster_unoise command",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
              {
                name: "denoise_level",
                items: ["global", "individual"],
                value: "global",
                disabled: "never",
                tooltip:
                  "at which level to perform denoising; global = by pooling samples, individual = independently for each sample (if samples are denoised individually, reducing minsize to 4 may be more reasonable for higher sensitivity)",
                type: "select",
              },
              {
                name: "remove_chimeras",
                value: true,
                disabled: "never",
                tooltip:
                  "perform chimera removal with UCHIME3 de novo algoritm",
                type: "bool",
              },
              {
                name: "abskew",
                value: 16,
                disabled: "never",
                tooltip:
                  "the abundance skew of chimeric sequences in comparsion with parental sequences (by default, parents should be at least 16 times more abundant than their chimera)",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
              {
                name: "similarity_type",
                items: ["0", "1", "2", "3", "4"],
                value: "2",
                disabled: "never",
                tooltip: "pairwise sequence identity definition (--iddef)",
                type: "select",
              },
              {
                name: "maxaccepts",
                value: 1,
                disabled: "never",
                tooltip:
                  "maximum number of hits to accept before stopping the search",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
              {
                name: "maxrejects",
                value: 32,
                disabled: "never",
                tooltip:
                  "maximum number of non-matching target sequences to consider before stopping the search",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
              {
                name: "mask",
                items: ["dust", "none"],
                value: "dust",
                disabled: "never",
                tooltip:
                  'mask regions in sequences using the "dust" method, or do not mask ("none").',
                type: "select",
              },
              {
                name: "cores",
                value: 4,
                disabled: "never",
                tooltip: "number of cores to use",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
            ],
            Inputs: [
              {
                name: "zOTUs_thresh",
                value: 1,
                disabled: "never",
                tooltip:
                  "sequence similarity threshold for zOTU table creation; 1 = 100% similarity threshold for zOTUs",
                max: 1,
                min: 0,
                step: 0.01,
                type: "slide",
              },
              {
                name: "similarity_threshold",
                value: 1,
                disabled: "never",
                tooltip:
                  "cluster zOTUs to OTUs based on the sequence similarity threshold; if id = 1, no OTU clustering will be performed",
                max: 1,
                min: 0,
                step: 0.01,
                type: "slide",
              },
              {
                name: "strands",
                items: ["both", "plus"],
                disabled: "never",
                tooltip:
                  "when comparing sequences with the cluster seed, check both strands (forward and reverse complementary) or the plus strand only",
                value: "both",
                type: "select",
              },
              {
                name: "minsize",
                value: 8,
                disabled: "never",
                tooltip: "minimum abundance of sequences for denoising",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
            ],
          },
        ],
      },
      {
        stepName: "postclustering",
        disabled: "never",
        services: [
          {
            tooltip: "postclustering with LULU algorithm",
            scriptName: "lulu.sh",
            imageName: "pipecraft/vsearch_dada2:2",
            serviceName: "LULU",
            selected: false,
            showExtra: false,
            extraInputs: [
              {
                name: "match_list_soft",
                items: ["vsearch", "BLAST"],
                value: "vsearch",
                disabled: "never",
                tooltip:
                  "use either 'blastn' or 'vsearch' to generate match list for LULU. Default is 'vsearch' (much faster)",
                type: "select",
              },
              {
                name: "vsearch_similarity_type",
                items: ["0", "1", "2", "3", "4"],
                value: "2",
                disabled: "never",
                tooltip:
                  "applies only when 'vsearch' is used as 'match_list_soft'. Pairwise sequence identity definition (--iddef)",
                type: "select",
              },
              {
                name: "perc_identity",
                value: 84,
                disabled: "never",
                tooltip:
                  "percent identity cutoff for match list. Excluding pairwise comparisons with lower sequence identity percentage than specified threshold",
                type: "numeric",
                rules: [
                  (v) => v >= 1 || "ERROR: specify values >= 1",
                  (v) => v <= 100 || "ERROR: specify values <= 100",
                ],
              },
              {
                name: "coverage_perc",
                value: 80,
                disabled: "never",
                tooltip:
                  "percent query coverage per hit. Excluding pairwise comparisons with lower sequence coverage than specified threshold",
                type: "numeric",
                rules: [
                  (v) => v >= 1 || "ERROR: specify values >= 1",
                  (v) => v <= 100 || "ERROR: specify values <= 100",
                ],
              },
              {
                name: "strands",
                items: ["plus", "both"],
                value: "both",
                disabled: "never",
                tooltip:
                  "query strand to search against database. Both = search also reverse complement",
                type: "select",
              },
              {
                name: "cores",
                value: 4,
                disabled: "never",
                tooltip:
                  "number of cores to use for generating match list for LULU",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
            ],
            Inputs: [
              {
                name: "table",
                active: false,
                btnName: "select file",
                value: "undefined",
                disabled: "never",
                tooltip:
                  "select OTU/ASV table. If no file is selected, then PipeCraft will look OTU_table.txt or ASV_table.txt in the working directory",
                type: "boolfile",
              },
              {
                name: "rep_seqs",
                active: false,
                btnName: "select file",
                value: "undefined",
                disabled: "never",
                tooltip:
                  "select fasta formatted sequence file containing your OTU/ASV reads",
                type: "boolfile",
              },
              {
                name: "min_ratio_type",
                items: ["min", "avg"],
                value: "min",
                disabled: "never",
                tooltip:
                  "sets whether a potential error must have lower abundance than the parent in all samples 'min' (default), or if an error just needs to have lower abundance on average 'avg'",
                type: "select",
              },
              {
                name: "min_ratio",
                value: 1,
                disabled: "never",
                tooltip:
                  "default = 1. Sets the minimim abundance ratio between a potential error and a potential parent to be identified as an error",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
              {
                name: "min_match",
                value: 90,
                disabled: "never",
                tooltip:
                  "default = 90%. specify minimum threshold of sequence similarity for considering any OTU as an error of another",
                type: "numeric",
                rules: [
                  (v) => v >= 1 || "ERROR: specify values >= 1",
                  (v) => v <= 100 || "ERROR: specify values <= 100",
                ],
              },
              {
                name: "min_rel_cooccurence",
                value: 0.95,
                disabled: "never",
                tooltip:
                  "minimum co-occurrence rate. Default = 0.95 (meaning that 1 in 20 samples are allowed to have no parent presence)",
                max: 1,
                min: 0,
                step: 0.01,
                type: "slide",
              },
            ],
          },

          {
            tooltip:
              "applies to DADA2 output ASV table (rds). Collaplse identical ASVs or/and filter ASVs by length [SELECT WORKDIR (data format, extension and read types are irrelevant here)]",
            scriptName: "table_filtering_dada2.sh",
            imageName: "pipecraft/vsearch_dada2:2",
            serviceName: "DADA2 collapse ASVs",
            disabled: "never",
            selected: false,
            showExtra: false,
            extraInputs: [
              {
                name: "minOverlap",
                value: 20,
                disabled: "never",
                tooltip:
                  "collapseNoMismatch setting. Default = 20. The minimum overlap of base pairs between ASV sequences required to collapse them together",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
              {
                name: "vec",
                value: true,
                disabled: "never",
                tooltip:
                  "collapseNoMismatch setting. Default = TRUE. Use the vectorized aligner. Should be turned off if sequences exceed 2kb in length",
                type: "bool",
              },
            ],
            Inputs: [
              {
                name: "DADA2_table",
                value: "undefined",
                btnName: "select rds",
                disabled: "never",
                tooltip:
                  "select the RDS file (ASV table), output from DADA2 workflow; usually in ASVs_out.dada2/ASVs_table.denoised-merged.rds",
                type: "file",
              },
              {
                name: "collapseNoMismatch",
                value: true,
                disabled: "never",
                tooltip:
                  "collapses ASVs in an ASV table that are identical up to shifts or length variation, i.e. that have no mismatches or internal indels (dada2 'collapseNoMismatch')",
                type: "bool",
              },
              {
                name: "by_length",
                value: 249,
                disabled: "never",
                tooltip:
                  "discard ASVs from the ASV table that are shorter than specified value (in base pairs). Value 0 means OFF; no filtering by length.",
                type: "numeric",
                rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
              },
            ],
          },
        ],
      },

      {
        stepName: "assign taxonomy",
        disabled: "never",
        services: [
          {
            tooltip:
              "assign taxonomy with BLAST against a selected database [SELECT WORKDIR that contains only ONE fasta file for the process]",
            scriptName: "taxonomy_BLAST.sh",
            imageName: "pipecraft/blast:2.14",
            serviceName: "BLAST",
            selected: false,
            showExtra: false,
            extraInputs: [
              {
                name: "e_value",
                value: 10,
                disabled: "never",
                tooltip:
                  "a parameter that describes the number of hits one can expect to see by chance when searching a database of a particular size. The lower the e-value the more 'significant' the match is",
                type: "numeric",
                default: 10,
                rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
              },
              {
                name: "word_size",
                value: 11,
                disabled: "never",
                tooltip:
                  "the size of the initial word that must be matched between the database and the query sequence",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
              {
                name: "reward",
                value: 2,
                disabled: "never",
                tooltip: "reward for a match",
                type: "numeric",
                rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
              },
              {
                name: "penalty",
                value: -3,
                disabled: "never",
                tooltip: "penalty for a mismatch",
                type: "numeric",
                rules: [(v) => v <= 0 || "ERROR: specify values <= 0"],
              },
              {
                name: "gap_open",
                value: 5,
                disabled: "never",
                tooltip: "cost to open a gap",
                type: "numeric",
                rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
              },
              {
                name: "gap_extend",
                value: 2,
                disabled: "never",
                tooltip: "cost to extend a gap",
                type: "numeric",
                // rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
              },
              {
                name: "cores",
                value: 4,
                disabled: "never",
                tooltip: "number of cores to use",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
            ],
            Inputs: [
              {
                name: "database_file",
                btnName: "select file",
                value: "undefined",
                disabled: "never",
                tooltip:
                  "database file (may be fasta formated - automatically will convert to BLAST database format)",
                type: "file",
              },
              {
                name: "task",
                items: ["blastn", "megablast"],
                value: "blastn",
                disabled: "never",
                tooltip: "task (blastn or megablast)",
                type: "select",
              },
              {
                name: "strands",
                items: ["plus", "both"],
                value: "both",
                disabled: "never",
                tooltip:
                  "query strand to search against database. Both = search also reverse complement",
                type: "select",
              },
            ],
          },
          {
            tooltip:
              "assign taxonomy with RDP Classifier against a selected database [SELECT WORKDIR that contains only ONE fasta file for the process]",
            scriptName: "taxonomy_RDP.sh",
            imageName: "pipecraft/metaworks:1.12.0",
            serviceName: "RDP_classifier",
            selected: false,
            showExtra: false,
            extraInputs: [
              {
                name: "mem",
                value: 10,
                disabled: "never",
                tooltip:
                  "default is 10GB. The amount of memory to allocate to the RDP classifier",
                type: "numeric",
                rules: [
                  (v) => v >= 0 || "ERROR: specify values >0",
                  (v) => v <= 300 || "ERROR: specify values <= 300",
                ],
              },
            ],
            Inputs: [
              {
                name: "database",
                btnName: "select RDP db",
                value: "undefined",
                disabled: "never",
                tooltip:
                  "RDP-trained reference database for the RDP classifier. Click on the header to download trained reference databases the RDP classifier, link MetaWorks user guide: https://terrimporter.github.io/MetaWorksSite/#classifier_table",
                type: "file",
              },
              {
                name: "confidence",
                value: 0.8,
                disabled: "never",
                tooltip:
                  "default is 0.8. Assignment confidence cutoff used to determine the assignment count for each taxon. Range [0-1]",
                type: "numeric",
                rules: [
                  (v) => v >= 0 || "ERROR: specify values >0",
                  (v) => v <= 1 || "ERROR: specify values <= 1",
                ],
              },
            ],
          },
          {
            tooltip:
              "assign taxonomy with DADA2 'assignTaxonomy' function [SELECT WORKDIR that contains only ONE fasta file for the process]",
            scriptName: "taxonomy_dada2.sh",
            imageName: "pipecraft/vsearch_dada2:2",
            serviceName: "DADA2 classifier",
            disabled: "never",
            selected: false,
            showExtra: false,
            extraInputs: [],
            Inputs: [
              {
                name: "dada2_database",
                btnName: "select fasta",
                value: "undefined",
                disabled: "never",
                tooltip:
                  "Select a reference database fasta file for taxonomy annotation. Click on the header to download DADA2-formatted reference databases https://benjjneb.github.io/dada2/training.html",
                type: "file",
              },
              {
                name: "minBoot",
                value: 50,
                disabled: "never",
                tooltip:
                  "the minimum bootstrap confidence for assigning a taxonomic level",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
              {
                name: "tryRC",
                value: false,
                disabled: "never",
                tooltip:
                  "the reverse-complement of each sequences will be used for classification if it is a better match to the reference sequences than the forward sequence",
                type: "bool",
              },
            ],
          },
        ],
      },
      {
        stepName: "postprocessing",
        disabled: "never",
        services: [
          {
            scriptName: "clustering_vsearch_ASVs2OTUs.sh",
            tooltip:
              "clustering ASVs to OTUs with vsearch; and making an OTU table",
            imageName: "pipecraft/vsearch_dada2:2",
            serviceName: "ASV_to_OTU",
            selected: false,
            showExtra: false,
            extraInputs: [
              {
                name: "OTU_type",
                items: ["centroid", "consout"],
                disabled: "never",
                tooltip:
                  '"centroid" = output centroid sequences; "consout" = output consensus sequences',
                value: "centroid",
                type: "select",
              },
              {
                name: "strands",
                items: ["both", "plus"],
                disabled: "never",
                tooltip:
                  "when comparing sequences with the cluster seed, check both strands (forward and reverse complementary) or the plus strand only",
                value: "both",
                type: "select",
              },
              {
                name: "remove_singletons",
                value: false,
                disabled: "never",
                tooltip:
                  "remove singleton OTUs (e.g., if TRUE, then OTUs with only one sequence will be discarded)",
                type: "bool",
              },
              {
                name: "similarity_type",
                items: ["0", "1", "2", "3", "4"],
                value: "2",
                disabled: "never",
                tooltip: "pairwise sequence identity definition (--iddef)",
                type: "select",
              },
              {
                name: "sequence_sorting",
                items: ["size", "length", "no"],
                value: "size",
                disabled: "never",
                tooltip:
                  'size = sort the sequences by decreasing abundance; "length" = sort the sequences by decreasing length (--cluster_fast); "no" = do not sort sequences (--cluster_smallmem --usersort)',
                type: "select",
              },
              {
                name: "centroid_type",
                items: ["similarity", "abundance"],
                value: "similarity",
                disabled: "never",
                tooltip:
                  '"similarity" = assign representative sequence to the closest (most similar) centroid (distance-based greedy clustering); "abundance" = assign representative sequence to the most abundant centroid (abundance-based greedy clustering; --sizeorder), --maxaccepts should be > 1',
                type: "select",
              },
              {
                name: "maxaccepts",
                value: 1,
                disabled: "never",
                tooltip:
                  "maximum number of hits to accept before stopping the search (should be > 1 for abundance-based selection of centroids [centroid type])",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
              {
                name: "mask",
                items: ["dust", "none"],
                value: "dust",
                disabled: "never",
                tooltip:
                  'mask regions in sequences using the "dust" method, or do not mask ("none").',
                type: "select",
              },
              {
                name: "cores",
                value: 4,
                disabled: "never",
                tooltip: "number of cores to use",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
            ],
            Inputs: [
              {
                name: "ASV_fasta",
                active: false,
                btnName: "select file",
                value: "undefined",
                disabled: "never",
                tooltip:
                  "select fasta formatted ASVs sequence file (ASV IDs must match with the ones in the ASVs table) [output will be in the directory as specified under 'SELECT WORKDIR']",
                type: "file",
              },
              {
                name: "ASV_table",
                active: true,
                btnName: "select file",
                value: "undefined",
                disabled: "never",
                tooltip:
                  "select ASVs_table file [1st col is ASVs ID, 2nd col must be 'Sequences' (default PipeCraft's output)]",
                type: "file",
              },
              {
                name: "similarity_threshold",
                value: 0.97,
                disabled: "never",
                tooltip:
                  "define OTUs based on the sequence similarity threshold; 0.97 = 97% similarity threshold",
                max: 1,
                min: 0,
                step: 0.01,
                type: "slide",
              },
            ],
          },
          {
            scriptName: "tag_jump_removal.sh",
            tooltip: "using UNCROSS2 to filter out putative tag-jumps",
            imageName: "vmikk/nextits:0.5.0",
            serviceName: "filter_tag-jumps",
            selected: false,
            showExtra: false,
            extraInputs: [],
            Inputs: [
              {
                name: "OTU_table",
                active: true,
                btnName: "select file",
                value: "undefined",
                disabled: "never",
                tooltip:
                  "select TAB-DELIMITED OTU/ASV table, where the 1st column is the OTU/ASV IDs and the following columns represent samples; 2nd column may be Sequence column, with the colName 'Sequence' [output will be in the directory as specified under 'SELECT WORKDIR']",
                type: "file",
              },
              {
                name: "f_value",
                value: 0.01,
                max: 1,
                min: 0,
                step: 0.01,
                disabled: "never",
                tooltip:
                  "f-parameter of UNCROSS2, which defines the expected cross-talk rate. Default is 0.01 (equivalent to 1%). A higher value enforces stricter filtering",
                type: "slide",
                rules: [(v) => v > 0 || "ERROR: specify values > 0"],
              },
              {
                name: "p_value",
                value: 1,
                disabled: "never",
                tooltip:
                  "p-parameter, which controls the severity of tag-jump removal. It adjusts the exponent in the UNCROSS formula. Default is 1. Opt for 0.5 or 0.3 to steepen the curve",
                type: "numeric",
                rules: [(v) => v > 0 || "ERROR: specify values > 0"],
              },
            ],
          },
          {
            tooltip:
              "pseudogene fintering with ORFfinder (search open reading frames) and/or HMMs",
            scriptName: "ORFfinder_hmm.sh",
            imageName: "pipecraft/metaworks:1.12.0",
            serviceName: "filter_numts",
            selected: false,
            showExtra: false,
            extraInputs: [
              {
                name: "ignore_nested",
                value: true,
                disabled: "never",
                tooltip:
                  "ignore nested open reading frames (completely placed within another)",
                type: "bool",
              },
              {
                name: "strand",
                items: ["plus", "minus", "both"],
                value: "both",
                disabled: "never",
                tooltip:
                  "output open reading frames (ORFs) on specified strand only",
                type: "select",
              },
            ],
            Inputs: [
              {
                name: "fasta_file",
                active: false,
                btnName: "select fasta file",
                value: "undefined",
                disabled: "never",
                tooltip:
                  "select fasta formatted sequence file containing your OTU/ASV reads. Sequence IDs cannot contaon underlines '_' [output will be in the directory as specified under 'SELECT WORKDIR']",
                type: "file",
              },
              {
                name: "min_len",
                value: 309,
                disabled: "never",
                tooltip: "minimum length of an output sequence",
                type: "numeric",
                rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              },
              // {
              //   name: "max_len",
              //   value: 317,
              //   disabled: "never",
              //   tooltip: "maximum length of an output sequence",
              //   type: "numeric",
              //   rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
              // },
              {
                name: "genetic_code",
                value: 5,
                disabled: "never",
                tooltip:
                  "5 = invertebrate mitochondrial code. Use 1 for rbcL. Specify values from 1 to 33; see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi",
                type: "numeric",
                rules: [
                  (v) => v >= 1 || "ERROR: specify values >= 1",
                  (v) => v <= 34 || "ERROR: specify values <= 34",
                ],
              },
              {
                name: "start_codon",
                value: 2,
                disabled: "never",
                tooltip:
                  "0 = ATG only; 1 = ATG and alternative initation codons; 2 = any sense codon",
                type: "numeric",
                rules: [
                  (v) => v >= 1 || "ERROR: specify values >= 1",
                  (v) => v <= 34 || "ERROR: specify values <= 34",
                ],
              },
              {
                name: "arthropod_hmm",
                value: true,
                disabled: "never",
                tooltip:
                  "on top of ORFfinder, remove sequences that don't match well to a profile HMM based on arthropod COI barcode sequences",
                type: "bool",
              },
            ],
          },
          {
            tooltip:
              "DEICODE (Robust Aitchison PCA on sparse compositional metabarcoding data)",
            scriptName: "DEICODE.sh",
            imageName: "pipecraft/deicode:0.2.4",
            serviceName: "DEICODE",
            selected: false,
            showExtra: false,
            extraInputs: [],
            Inputs: [
              {
                name: "table",
                active: false,
                btnName: "select file",
                value: "undefined",
                disabled: "never",
                tooltip:
                  "select OTU/ASV table. If no file is selected, then PipeCraft will look OTU_table.txt or ASV_table.txt in the WORKDIR [output will be in the directory as specified under 'SELECT WORKDIR']",
                type: "boolfile",
              },
              {
                name: "subset_IDs",
                active: false,
                btnName: "select file",
                value: "undefined",
                disabled: "never",
                tooltip:
                  "select list of OTU/ASV IDs for analysing a subset from the full table",
                type: "boolfile",
              },
              {
                name: "min_otu_reads",
                value: 10,
                disabled: "never",
                tooltip:
                  "cutoff for reads per OTU/ASV. OTUs/ASVs with lower reads then specified cutoff will be excluded from the analysis",
                type: "numeric",
                rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
              },
              {
                name: "min_sample_reads",
                value: 500,
                disabled: "never",
                tooltip:
                  "cutoff for reads per sample. Samples with lower reads then specified cutoff will be excluded from the analysis",
                type: "numeric",
                rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
              },
            ],
          },
        ],
      },
    ],
    vsearch_OTUs: [
      {
        tooltip: "remove primers sequences from the reads",
        scriptName: "cut_primers_paired_end_reads.sh",
        imageName: "pipecraft/cutadapt:4.4",
        serviceName: "cut primers",
        manualLink:
          "https://pipecraft2-manual.readthedocs.io/en/latest/quicktools.html#cut-primers",
        disabled: "never",
        selected: false,
        showExtra: false,
        extraInputs: [
          {
            name: "cores",
            value: 1,
            disabled: "never",
            tooltip:
              "number of cores to use. For paired-end data in fasta format, set to 1 [default]. For fastq formats you may set the value to 0 to use all cores",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "no_indels",
            value: true,
            disabled: "never",
            tooltip:
              "do not allow insertions or deletions is primer search. Mismatches are the only type of errors accounted in the error rate parameter",
            type: "bool",
          },
        ],
        Inputs: [
          {
            name: "forward_primers",
            value: [],
            disabled: "never",
            tooltip: "specify forward primer (5'-3'); add up to 13 primers",
            type: "chip",
            iupac: true,
            rules: [(v) => v.length <= 13 || "TOO MANY PRIMERS"],
          },
          {
            name: "reverse_primers",
            value: [],
            disabled: "never",
            tooltip: "specify reverse primer (3'-5'); add up to 13 primers",
            type: "chip",
            iupac: true,
            rules: [(v) => v.length <= 13 || "TOO MANY PRIMERS"],
          },
          {
            name: "mismatches",
            value: 1,
            disabled: "never",
            tooltip: "allowed mismatches in the primer search",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
          },
          {
            name: "min_overlap",
            value: 21,
            disabled: "never",
            tooltip:
              "number of overlap bases with the primer sequence. Partial matches are allowed, but short matches may occur by chance, leading to erroneously clipped bases. Specifying higher overlap than the length of primer sequnce will still clip the primer (e.g. primer length is 22 bp, but overlap is specified as 25 - this does not affect the identification and clipping of the primer as long as the match is in the specified mismatch error range)",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "seqs_to_keep",
            items: ["keep_all", "keep_only_linked"],
            value: "keep_all",
            disabled: "never",
            tooltip:
              "'keep_all' keeps the sequences where at least one primer was found (fwd or rev); recommended when cutting primers from paired-end data (unassembled), where individual R1 or R2 read lenghts are shorther than the expected amplicon length. 'keep_only_linked' = keep sequences if primers are found in both ends (fwd…rev); discards the read if both primers were not found in this read; maybe useful for single-end data",
            type: "select",
          },
          {
            name: "pair_filter",
            items: ["both", "any"],
            value: "both",
            disabled: "single_end",
            tooltip:
              "applies only for paired-end data. Option 'both' discards a read-pair when both reads do not contain a primer sequence. Option 'any' discards the read-pair when one of the reads does not contain a primer sequence",
            type: "select",
          },
        ],
      },
      {
        tooltip: "assemble paired-end reads with vsearch",
        scriptName: "assemble_paired_end_data_vsearch.sh",
        imageName: "pipecraft/vsearch_dada2:2",
        serviceName: "merge reads",
        manualLink:
          "https://pipecraft2-manual.readthedocs.io/en/latest/quicktools.html#merge-vsearch",
        selected: "always",
        disabled: "single_end",
        showExtra: false,
        extraInputs: [
          {
            name: "max_diffs",
            value: 20,
            disabled: "never",
            tooltip:
              "the maximum number of non-matching nucleotides allowed in the overlap region",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
          },
          {
            name: "max_Ns",
            value: 0,
            disabled: "never",
            tooltip:
              "discard sequences with more than the specified number of Ns",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
          },
          {
            name: "max_len",
            value: 600,
            disabled: "never",
            tooltip: "maximum length of the merged sequence",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "keep_disjointed",
            value: false,
            disabled: "never",
            tooltip:
              "output reads that were not merged into separate FASTQ files",
            type: "bool",
          },
          {
            name: "fastq_qmax",
            value: 41,
            disabled: "never",
            tooltip:
              "maximum quality score accepted when reading FASTQ files. The default is 41, which is usual for recent Sanger/Illumina 1.8+ files",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
        ],
        Inputs: [
          {
            name: "read_R1",
            value: ["_R1"],
            disabled: "single_end",
            tooltip:
              "identifyer string that is common for all R1 reads. R1/R2 reads must contain R1/R2 strings in the file names; sample names must not contain R1/R2! (e.g. when all R1 files have '.R1' string, then enter '\\.R1'. Note that backslash is only needed to escape dot regex; e.g. when all R1 files have '_R1' string, then enter '_R1'.)",
            type: "chip",
            rules: [(v) => v.length <= 1 || "ADD ONLY ONE IDENTIFIER"],
          },
          {
            name: "min_overlap",
            value: 12,
            disabled: "never",
            tooltip: "minimum overlap between the merged reads",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "min_lenght",
            value: 32,
            disabled: "never",
            tooltip: "minimum length of the merged sequence",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "allow_merge_stagger",
            value: true,
            disabled: "never",
            tooltip:
              "allow to merge staggered read pairs. Staggered pairs are pairs where the 3' end of the reverse read has an overhang to the left of the 5' end of the forward read. This situation can occur when a very short fragment is sequenced",
            type: "bool",
          },
          {
            name: "include_only_R1",
            value: false,
            disabled: "never",
            tooltip:
              "include unassembled R1 reads to the set of assembled reads per sample. This may be relevant when working with e.g. ITS2 sequences, because the ITS2 region in some taxa is too long for assembly, therefore discarded completely after assembly process. Thus, including also unassembled R1 reads, partial ITS2 sequences for these taxa will be represented in the final output",
            type: "bool",
          },
        ],
      },
      {
        tooltip: "quality filtering with vsearch",
        scriptName: "quality_filtering_paired_end_vsearch.sh",
        imageName: "pipecraft/vsearch_dada2:2",
        serviceName: "quality filtering",
        manualLink:
          "https://pipecraft2-manual.readthedocs.io/en/latest/quicktools.html#qfilt-vsearch",
        disabled: "never",
        selected: "always",
        showExtra: false,
        extraInputs: [
          {
            name: "cores",
            value: 1,
            disabled: "never",
            tooltip: "number of cores to use",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "max_length",
            value: null,
            disabled: "never",
            tooltip:
              "discard sequences with more than the specified number of bases. Note that if 'trunc length' setting is specified, then 'max length' SHOULD NOT be lower than 'trunc lenght' (otherwise all reads are discared)",
            type: "numeric",
            rules: [
              (v) =>
                (v >= 1) | (v == "") ||
                "ERROR: specify values >= 1 or leave it empty (=no action taken)",
            ],
          },
          {
            name: "qmax",
            value: 41,
            disabled: "never",
            tooltip:
              "specify the maximum quality score accepted when reading FASTQ files. The default is 41, which is usual for recent Sanger/Illumina 1.8+ files. For PacBio data use 93",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "qmin",
            value: 0,
            disabled: "never",
            tooltip:
              "the minimum quality score accepted for FASTQ files. The default is 0, which is usual for recent Sanger/Illumina 1.8+ files. Older formats may use scores between -5 and 2",
            type: "numeric",
            rules: [(v) => v >= -5 || "ERROR: specify values >= -5"],
          },
          {
            name: "maxee_rate",
            value: null,
            disabled: "never",
            tooltip:
              "discard sequences with more than the specified number of expected errors per base",
            type: "numeric",
            rules: [
              (v) =>
                (v >= 0.001) | (v == "") ||
                "ERROR: specify values >=0.001 or leave it empty (= no action taken)",
            ],
          },
        ],
        Inputs: [
          {
            name: "maxee",
            value: 1,
            disabled: "never",
            tooltip:
              "maximum number of expected errors per sequence. Sequences with higher error rates will be discarded",
            type: "numeric",
            rules: [(v) => v >= 0.001 || "ERROR: specify values >= 0.001"],
          },
          {
            name: "maxNs",
            value: 0,
            disabled: "never",
            tooltip:
              "discard sequences with more than the specified number of Ns",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
          },
          {
            name: "min_length",
            value: 32,
            disabled: "never",
            tooltip:
              "minimum length of the filtered output sequence. Note that if 'trunc length' setting is specified, then 'min length' SHOULD BE lower than 'trunc lenght' (otherwise all reads are discared)",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "trunc_length",
            value: null,
            disabled: "never",
            tooltip:
              "truncate sequences to the specified length. Shorter sequences are discarded; thus if specified, check that 'min length' setting is lower than 'trunc length' ('min lenght' therefore has basically no effect)",
            type: "numeric",
            rules: [
              (v) =>
                (v >= 5) | (v == "") ||
                "ERROR: specify values >= 5 or leave it empty (=no action taken)",
            ],
          },
        ],
      },
      {
        tooltip:
          "chimera filtering with vsearch. Untick the checkbox to skip this step",
        scriptName: "chimera_filtering_vsearch.sh",
        imageName: "pipecraft/vsearch_dada2:2",
        serviceName: "chimera filtering",
        manualLink:
          "https://pipecraft2-manual.readthedocs.io/en/latest/quicktools.html#chimera-filtering",
        disabled: "never",
        selected: "always",
        showExtra: false,
        extraInputs: [
          {
            name: "cores",
            value: 4,
            disabled: "never",
            tooltip:
              "Number of cores to use (only for reference based chimera filtering)",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "abundance_skew",
            value: 2,
            disabled: "never",
            tooltip:
              "the abundance skew is used to distinguish in a threeway alignment which sequence is the chimera and which are the parents. The assumption is that chimeras appear later in the PCR amplification process and are therefore less abundant than their parents. The default value is 2.0, which means that the parents should be at least 2 times more abundant than their chimera. Any positive value equal or greater than 1.0 can be used",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "min_h",
            value: 0.28,
            disabled: "never",
            tooltip:
              "minimum score (h). Increasing this value tends to reduce the number of false positives and to decrease sensitivity. Values ranging from 0.0 to 1.0 included are accepted",
            max: 1,
            min: 0,
            step: 0.01,
            type: "slide",
          },
        ],
        Inputs: [
          {
            name: "pre_cluster",
            value: 0.97,
            disabled: "never",
            tooltip:
              "identity percentage when performing 'pre-clustering' with --cluster_size for denovo chimera filtering with --uchime_denovo",
            max: 1,
            min: 0,
            step: 0.01,
            type: "slide",
          },
          {
            name: "min_unique_size",
            value: 1,
            disabled: "never",
            tooltip:
              "minimum amount of a unique sequences in a fasta file. If value = 1, then no sequences are discarded after dereplication; if value = 2, then sequences, which are represented only once in a given file are discarded; and so on",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "denovo",
            value: true,
            disabled: "never",
            tooltip:
              "if TRUE, then perform denovo chimera filtering with --uchime_denovo",
            type: "bool",
          },
          {
            name: "reference_based",
            active: false,
            btnName: "select file",
            value: "undefined",
            disabled: "never",
            tooltip:
              "perform reference database based chimera filtering with --uchime_ref. If denovo = TRUE, then reference based chimera filtering will be performed after denovo",
            type: "boolfile",
          },
        ],
      },
      {
        tooltip:
          "if data set consists of ITS sequences; identify and extract the ITS regions using ITSx. Select appropriate ITSx output region for CLUSTERING after the process is finished",
        scriptName: "ITS_extractor.sh",
        imageName: "pipecraft/itsx:1.1.3",
        serviceName: "itsx",
        manualLink:
          "https://pipecraft2-manual.readthedocs.io/en/latest/quicktools.html#id17",
        disabled: "never",
        selected: false,
        showExtra: false,
        extraInputs: [
          {
            name: "cores",
            value: 6,
            disabled: "never",
            tooltip: "number of cores to use",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "e_value",
            value: (0.01).toExponential(),
            disabled: "never",
            tooltip:
              "domain E-value cutoff a sequence must obtain in the HMMER-based step to be included in the output. Here, the defaul 1e-2 = 0.01 (more relaxed compared with the default ITSx [1e-5])",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
          },
          {
            name: "scores",
            value: 0,
            disabled: "never",
            tooltip:
              "domain score cutoff that a sequence must obtain in the HMMER-based step to be included in the output",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
          },
          {
            name: "domains",
            value: 2,
            disabled: "never",
            tooltip:
              "the minimum number of domains (different HMM gene profiles) that must match a sequence for it to be included in the output (detected as an ITS sequence). Setting the value lower than two will increase the number of false positives, while increasing it above two will decrease ITSx detection abilities on fragmentary data",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
          },
          {
            name: "complement",
            value: true,
            disabled: "never",
            tooltip:
              "if TRUE, ITSx checks both DNA strands for matches to HMM-profiles",
            type: "bool",
          },
          {
            name: "only_full",
            value: false,
            disabled: "never",
            tooltip:
              "If TRUE, the output is limited to full-length ITS1 and ITS2 regions only",
            type: "bool",
          },
          {
            name: "truncate",
            value: true,
            disabled: "never",
            tooltip:
              "if TRUE, ITSx removes ends of ITS sequences if they are outside of the ITS region. If off, the whole input sequence is saved when ITS region is detected",
            type: "bool",
          },
        ],
        Inputs: [
          {
            name: "organisms",
            items: [
              "all",
              "alveolata",
              "bryophyta",
              "bacillariophyta",
              "amoebozoa",
              "euglenozoa",
              "fungi",
              "chlorophyta",
              "rhodophyta",
              "phaeophyceae",
              "marchantiophyta",
              "metazoa",
              "oomycota",
              "haptophyceae",
              "raphidophyceae",
              "rhizaria",
              "synurophyceae",
              "tracheophyta",
              "eustigmatophyceae",
              "apusozoa",
              "parabasalia",
            ],
            value: ["fungi"],
            disabled: "never",
            tooltip:
              "set of profiles to use for the search. Can be used to restrict the search to only a few organism groups types to save time, if one or more of the origins are not relevant to the dataset under study",
            type: "combobox",
          },
          {
            name: "regions",
            items: ["all", "SSU", "ITS1", "5.8S", "ITS2", "LSU"],
            value: ["all"],
            disabled: "never",
            tooltip:
              "ITS regions to output (note that 'all' will output also full ITS region [ITS1-5.8S-ITS2])",
            type: "combobox",
          },
          {
            name: "partial",
            value: 50,
            disabled: "never",
            tooltip:
              "if larger than 0, ITSx will save additional FASTA-files for full and partial ITS sequences longer than the specified cutoff value. If his setting is left to 0 (zero), it means OFF",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
          },
          {
            name: "region_for_clustering",
            items: ["full_ITS", "SSU", "ITS1", "5.8S", "ITS2", "LSU"],
            value: "ITS2",
            disabled: "never",
            tooltip:
              "specify region for clustering (because multiple output folders are generated during this process)",
            type: "select",
          },
          {
            name: "cluster_full_and_partial",
            value: true,
            disabled: "never",
            tooltip:
              "if setting 'partial' is not 0, then at the NEXT STEP cluster 'full and partial' (e.g.) ITS2 reads (dir /ITS2/full_and_partial). If OFF, then cluster only full ITS2 reads",
            type: "bool",
          },
        ],
      },
      {
        tooltip: "cluster reads to OTUs with vsearch",
        scriptName: "clustering_vsearch.sh",
        imageName: "pipecraft/vsearch_dada2:2",
        serviceName: "clustering",
        manualLink:
          "https://pipecraft2-manual.readthedocs.io/en/latest/quicktools.html#clustering-vsearch",
        disabled: "never",
        selected: "always",
        showExtra: false,
        extraInputs: [
          {
            name: "similarity_type",
            items: ["0", "1", "2", "3", "4"],
            value: "2",
            disabled: "never",
            tooltip: "pairwise sequence identity definition (--iddef)",
            type: "select",
          },
          {
            name: "sequence_sorting",
            items: ["size", "length", "no"],
            value: "size",
            disabled: "never",
            tooltip:
              'size = sort the sequences by decreasing abundance; "length" = sort the sequences by decreasing length (--cluster_fast); "no" = do not sort sequences (--cluster_smallmem --usersort)',
            type: "select",
          },
          {
            name: "centroid_type",
            items: ["similarity", "abundance"],
            value: "similarity",
            disabled: "never",
            tooltip:
              '"similarity" = assign representative sequence to the closest (most similar) centroid (distance-based greedy clustering); "abundance" = assign representative sequence to the most abundant centroid (abundance-based greedy clustering; --sizeorder), --maxaccepts should be > 1',
            type: "select",
          },
          {
            name: "maxaccepts",
            value: 1,
            disabled: "never",
            tooltip:
              "maximum number of hits to accept before stopping the search (should be > 1 for abundance-based selection of centroids [centroid type])",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "mask",
            items: ["dust", "none"],
            value: "dust",
            disabled: "never",
            tooltip:
              'mask regions in sequences using the "dust" method, or do not mask ("none")',
            type: "select",
          },
          {
            name: "cores",
            value: 4,
            disabled: "never",
            tooltip: "number of cores to use",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
        ],
        Inputs: [
          {
            name: "OTU_type",
            items: ["centroid", "consout"],
            disabled: "never",
            tooltip:
              '"centroid" = output centroid sequences; "consout" = output consensus sequences',
            value: "centroid",
            type: "select",
          },
          {
            name: "similarity_threshold",
            value: 0.97,
            disabled: "never",
            tooltip:
              "define OTUs based on the sequence similarity threshold; 0.97 = 97% similarity threshold",
            max: 1,
            min: 0,
            step: 0.01,
            type: "slide",
          },
          {
            name: "strands",
            items: ["both", "plus"],
            disabled: "never",
            tooltip:
              "when comparing sequences with the cluster seed, check both strands (forward and reverse complementary) or the plus strand only",
            value: "both",
            type: "select",
          },
          {
            name: "remove_singletons",
            value: false,
            disabled: "never",
            tooltip:
              "remove singleton OTUs (e.g., if TRUE, then OTUs with only one sequence will be discarded)",
            type: "bool",
          },
        ],
      },
      {
        tooltip: "assign taxonomy with BLAST against a selected database",
        scriptName: "taxonomy_BLAST.sh",
        imageName: "pipecraft/blast:2.14",
        serviceName: "assign taxonomy",
        manualLink:
          "https://pipecraft2-manual.readthedocs.io/en/latest/quicktools.html#blast-camacho-et-al-2009",
        selected: false,
        showExtra: false,
        extraInputs: [
          {
            name: "e_value",
            value: 10,
            disabled: "never",
            tooltip:
              "a parameter that describes the number of hits one can expect to see by chance when searching a database of a particular size. The lower the e-value the more 'significant' the match is",
            type: "numeric",
            default: 10,
            rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
          },
          {
            name: "word_size",
            value: 11,
            disabled: "never",
            tooltip:
              "the size of the initial word that must be matched between the database and the query sequence",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "reward",
            value: 2,
            disabled: "never",
            tooltip: "reward for a match",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
          },
          {
            name: "penalty",
            value: -3,
            disabled: "never",
            tooltip: "penalty for a mismatch",
            type: "numeric",
            rules: [(v) => v <= 0 || "ERROR: specify values <= 0"],
          },
          {
            name: "gap_open",
            value: 5,
            disabled: "never",
            tooltip: "cost to open a gap",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
          },
          {
            name: "gap_extend",
            value: 2,
            disabled: "never",
            tooltip: "cost to extend a gap",
            type: "numeric",
            // rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
          },
          {
            name: "cores",
            value: 4,
            disabled: "never",
            tooltip: "number of cores to use",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
        ],
        Inputs: [
          {
            name: "database_file",
            btnName: "select file",
            value: "undefined",
            disabled: "never",
            tooltip:
              "database file (may be fasta formated - automatically will convert to BLAST database format)",
            type: "file",
          },
          {
            name: "task",
            items: ["blastn", "megablast"],
            value: "blastn",
            disabled: "never",
            tooltip: "task (blastn or megablast)",
            type: "select",
          },
          {
            name: "strands",
            items: ["plus", "both"],
            value: "both",
            disabled: "never",
            tooltip:
              "query strand to search against database. Both = search also reverse complement",
            type: "select",
          },
        ],
      },
    ],
    UNOISE_ASVs: [
      {
        tooltip: "remove primers sequences from the reads",
        scriptName: "cut_primers_paired_end_reads.sh",
        imageName: "pipecraft/cutadapt:4.4",
        serviceName: "cut primers",
        manualLink:
          "https://pipecraft2-manual.readthedocs.io/en/latest/quicktools.html#cut-primers",
        disabled: "never",
        selected: false,
        showExtra: false,
        extraInputs: [
          {
            name: "cores",
            value: 1,
            disabled: "never",
            tooltip:
              "number of cores to use. For paired-end data in fasta format, set to 1 [default]. For fastq formats you may set the value to 0 to use all cores",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "no_indels",
            value: true,
            disabled: "never",
            tooltip:
              "do not allow insertions or deletions is primer search. Mismatches are the only type of errors accounted in the error rate parameter",
            type: "bool",
          },
        ],
        Inputs: [
          {
            name: "forward_primers",
            value: [],
            disabled: "never",
            tooltip: "specify forward primer (5'-3'); add up to 13 primers",
            type: "chip",
            iupac: true,
            rules: [(v) => v.length <= 13 || "TOO MANY PRIMERS"],
          },
          {
            name: "reverse_primers",
            value: [],
            disabled: "never",
            tooltip: "specify reverse primer (3'-5'); add up to 13 primers",
            type: "chip",
            iupac: true,
            rules: [(v) => v.length <= 13 || "TOO MANY PRIMERS"],
          },
          {
            name: "mismatches",
            value: 1,
            disabled: "never",
            tooltip: "allowed mismatches in the primer search",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
          },
          {
            name: "min_overlap",
            value: 21,
            disabled: "never",
            tooltip:
              "number of overlap bases with the primer sequence. Partial matches are allowed, but short matches may occur by chance, leading to erroneously clipped bases. Specifying higher overlap than the length of primer sequnce will still clip the primer (e.g. primer length is 22 bp, but overlap is specified as 25 - this does not affect the identification and clipping of the primer as long as the match is in the specified mismatch error range)",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "seqs_to_keep",
            items: ["keep_all", "keep_only_linked"],
            value: "keep_all",
            disabled: "never",
            tooltip:
              "'keep_all' keeps the sequences where at least one primer was found (fwd or rev); recommended when cutting primers from paired-end data (unassembled), where individual R1 or R2 read lenghts are shorther than the expected amplicon length. 'keep_only_linked' = keep sequences if primers are found in both ends (fwd…rev); discards the read if both primers were not found in this read; maybe useful for single-end data",
            type: "select",
          },
          {
            name: "pair_filter",
            items: ["both", "any"],
            value: "both",
            disabled: "single_end",
            tooltip:
              "applies only for paired-end data. Option 'both' discards a read-pair when both reads do not contain a primer sequence. Option 'any' discards the read-pair when one of the reads does not contain a primer sequence",
            type: "select",
          },
        ],
      },
      {
        tooltip: "assemble paired-end reads with vsearch",
        scriptName: "assemble_paired_end_data_vsearch.sh",
        imageName: "pipecraft/vsearch_dada2:2",
        serviceName: "merge reads",
        manualLink:
          "https://pipecraft2-manual.readthedocs.io/en/latest/quicktools.html#merge-vsearch",
        selected: "always",
        disabled: "single_end",
        showExtra: false,
        extraInputs: [
          {
            name: "max_diffs",
            value: 20,
            disabled: "never",
            tooltip:
              "the maximum number of non-matching nucleotides allowed in the overlap region",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
          },
          {
            name: "max_Ns",
            value: 0,
            disabled: "never",
            tooltip:
              "discard sequences with more than the specified number of Ns",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
          },
          {
            name: "max_len",
            value: 600,
            disabled: "never",
            tooltip: "maximum length of the merged sequence",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "keep_disjointed",
            value: false,
            disabled: "never",
            tooltip:
              "output reads that were not merged into separate FASTQ files",
            type: "bool",
          },
          {
            name: "fastq_qmax",
            value: 41,
            disabled: "never",
            tooltip:
              "maximum quality score accepted when reading FASTQ files. The default is 41, which is usual for recent Sanger/Illumina 1.8+ files",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
        ],
        Inputs: [
          {
            name: "read_R1",
            value: ["_R1"],
            disabled: "single_end",
            tooltip:
              "identifyer string that is common for all R1 reads. R1/R2 reads must contain R1/R2 strings in the file names; sample names must not contain R1/R2! (e.g. when all R1 files have '.R1' string, then enter '\\.R1'. Note that backslash is only needed to escape dot regex; e.g. when all R1 files have '_R1' string, then enter '_R1'.)",
            type: "chip",
            rules: [(v) => v.length <= 1 || "ADD ONLY ONE IDENTIFIER"],
          },
          {
            name: "min_overlap",
            value: 12,
            disabled: "never",
            tooltip: "minimum overlap between the merged reads",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "min_lenght",
            value: 32,
            disabled: "never",
            tooltip: "minimum length of the merged sequence",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "allow_merge_stagger",
            value: true,
            disabled: "never",
            tooltip:
              "allow to merge staggered read pairs. Staggered pairs are pairs where the 3' end of the reverse read has an overhang to the left of the 5' end of the forward read. This situation can occur when a very short fragment is sequenced",
            type: "bool",
          },
          {
            name: "include_only_R1",
            value: false,
            disabled: "never",
            tooltip:
              "include unassembled R1 reads to the set of assembled reads per sample. This may be relevant when working with e.g. ITS2 sequences, because the ITS2 region in some taxa is too long for assembly, therefore discarded completely after assembly process. Thus, including also unassembled R1 reads, partial ITS2 sequences for these taxa will be represented in the final output",
            type: "bool",
          },
        ],
      },
      {
        tooltip: "quality filtering with vsearch",
        scriptName: "quality_filtering_paired_end_vsearch.sh",
        imageName: "pipecraft/vsearch_dada2:2",
        serviceName: "quality filtering",
        manualLink:
          "https://pipecraft2-manual.readthedocs.io/en/latest/quicktools.html#qfilt-vsearch",
        disabled: "never",
        selected: "always",
        showExtra: false,
        extraInputs: [
          {
            name: "cores",
            value: 1,
            disabled: "never",
            tooltip: "number of cores to use",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "max_length",
            value: null,
            disabled: "never",
            tooltip:
              "discard sequences with more than the specified number of bases. Note that if 'trunc length' setting is specified, then 'max length' SHOULD NOT be lower than 'trunc lenght' (otherwise all reads are discared)",
            type: "numeric",
            rules: [
              (v) =>
                (v >= 1) | (v == "") ||
                "ERROR: specify values >= 1 or leave it empty (=no action taken)",
            ],
          },
          {
            name: "qmax",
            value: 41,
            disabled: "never",
            tooltip:
              "specify the maximum quality score accepted when reading FASTQ files. The default is 41, which is usual for recent Sanger/Illumina 1.8+ files. For PacBio data use 93",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "qmin",
            value: 0,
            disabled: "never",
            tooltip:
              "the minimum quality score accepted for FASTQ files. The default is 0, which is usual for recent Sanger/Illumina 1.8+ files. Older formats may use scores between -5 and 2",
            type: "numeric",
            rules: [(v) => v >= -5 || "ERROR: specify values >= -5"],
          },
          {
            name: "maxee_rate",
            value: null,
            disabled: "never",
            tooltip:
              "discard sequences with more than the specified number of expected errors per base",
            type: "numeric",
            rules: [
              (v) =>
                (v >= 0.001) | (v == "") ||
                "ERROR: specify values >=0.001 or leave it empty (= no action taken)",
            ],
          },
        ],
        Inputs: [
          {
            name: "maxee",
            value: 1,
            disabled: "never",
            tooltip:
              "maximum number of expected errors per sequence. Sequences with higher error rates will be discarded",
            type: "numeric",
            rules: [(v) => v >= 0.001 || "ERROR: specify values >= 0.001"],
          },
          {
            name: "maxNs",
            value: 0,
            disabled: "never",
            tooltip:
              "discard sequences with more than the specified number of Ns",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
          },
          {
            name: "min_length",
            value: 32,
            disabled: "never",
            tooltip:
              "minimum length of the filtered output sequence. Note that if 'trunc length' setting is specified, then 'min length' SHOULD BE lower than 'trunc lenght' (otherwise all reads are discared)",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "trunc_length",
            value: null,
            disabled: "never",
            tooltip:
              "truncate sequences to the specified length. Shorter sequences are discarded; thus if specified, check that 'min length' setting is lower than 'trunc length' ('min lenght' therefore has basically no effect)",
            type: "numeric",
            rules: [
              (v) =>
                (v >= 5) | (v == "") ||
                "ERROR: specify values >= 5 or leave it empty (=no action taken)",
            ],
          },
        ],
      },
      {
        tooltip:
          "if data set consists of ITS sequences; identify and extract the ITS regions using ITSx. Select appropriate ITSx output region for UNOISE after the process is finished",
        scriptName: "ITS_extractor.sh",
        imageName: "pipecraft/itsx:1.1.3",
        serviceName: "itsx",
        manualLink:
          "https://pipecraft2-manual.readthedocs.io/en/latest/quicktools.html#id17",
        disabled: "never",
        selected: false,
        showExtra: false,
        extraInputs: [
          {
            name: "cores",
            value: 6,
            disabled: "never",
            tooltip: "number of cores to use",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "e_value",
            value: (0.01).toExponential(),
            disabled: "never",
            tooltip:
              "domain E-value cutoff a sequence must obtain in the HMMER-based step to be included in the output. Here, the defaul 1e-2 = 0.01 (more relaxed compared with the default ITSx [1e-5])",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
          },
          {
            name: "scores",
            value: 0,
            disabled: "never",
            tooltip:
              "domain score cutoff that a sequence must obtain in the HMMER-based step to be included in the output",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
          },
          {
            name: "domains",
            value: 2,
            disabled: "never",
            tooltip:
              "the minimum number of domains (different HMM gene profiles) that must match a sequence for it to be included in the output (detected as an ITS sequence). Setting the value lower than two will increase the number of false positives, while increasing it above two will decrease ITSx detection abilities on fragmentary data",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
          },
          {
            name: "complement",
            value: true,
            disabled: "never",
            tooltip:
              "if TRUE, ITSx checks both DNA strands for matches to HMM-profiles",
            type: "bool",
          },
          {
            name: "only_full",
            value: false,
            disabled: "never",
            tooltip:
              "If TRUE, the output is limited to full-length ITS1 and ITS2 regions only",
            type: "bool",
          },
          {
            name: "truncate",
            value: true,
            disabled: "never",
            tooltip:
              "if TRUE, ITSx removes ends of ITS sequences if they are outside of the ITS region. If off, the whole input sequence is saved when ITS region is detected",
            type: "bool",
          },
        ],
        Inputs: [
          {
            name: "organisms",
            items: [
              "all",
              "alveolata",
              "bryophyta",
              "bacillariophyta",
              "amoebozoa",
              "euglenozoa",
              "fungi",
              "chlorophyta",
              "rhodophyta",
              "phaeophyceae",
              "marchantiophyta",
              "metazoa",
              "oomycota",
              "haptophyceae",
              "raphidophyceae",
              "rhizaria",
              "synurophyceae",
              "tracheophyta",
              "eustigmatophyceae",
              "apusozoa",
              "parabasalia",
            ],
            value: ["fungi"],
            disabled: "never",
            tooltip:
              "set of profiles to use for the search. Can be used to restrict the search to only a few organism groups types to save time, if one or more of the origins are not relevant to the dataset under study",
            type: "combobox",
          },
          {
            name: "regions",
            items: ["all", "SSU", "ITS1", "5.8S", "ITS2", "LSU"],
            value: ["all"],
            disabled: "never",
            tooltip:
              "ITS regions to output (note that 'all' will output also full ITS region [ITS1-5.8S-ITS2])",
            type: "combobox",
          },
          {
            name: "partial",
            value: 50,
            disabled: "never",
            tooltip:
              "if larger than 0, ITSx will save additional FASTA-files for full and partial ITS sequences longer than the specified cutoff value. If his setting is left to 0 (zero), it means OFF",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
          },
          {
            name: "region_for_clustering",
            items: ["full_ITS", "SSU", "ITS1", "5.8S", "ITS2", "LSU"],
            value: "ITS2",
            disabled: "never",
            tooltip:
              "specify region for next step, UNOISE (because multiple output folders are generated during this process)",
            type: "select",
          },
          {
            name: "cluster_full_and_partial",
            value: true,
            disabled: "never",
            tooltip:
              "if setting 'partial' is not 0, then at the NEXT STEP cluster 'full and partial' (e.g.) ITS2 reads (dir /ITS2/full_and_partial). If OFF, then cluster only full ITS2 reads",
            type: "bool",
          },
        ],
      },
      {
        scriptName: "clustering_unoise.sh",
        tooltip:
          "tick the checkbox to cluster reads with vsearch --cluster_unoise (and optionally remove chimeras with --uchime3_denovo)",
        imageName: "pipecraft/vsearch_dada2:2",
        serviceName: "unoise3",
        selected: false,
        showExtra: false,
        extraInputs: [
          {
            name: "unoise_alpha",
            value: 2.0,
            disabled: "never",
            tooltip:
              "default = 2.0. alpha parameter to the vsearch --cluster_unoise command",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "denoise_level",
            items: ["global", "individual"],
            value: "global",
            disabled: "never",
            tooltip:
              "at which level to perform denoising; global = by pooling samples, individual = independently for each sample (if samples are denoised individually, reducing minsize to 4 may be more reasonable for higher sensitivity)",
            type: "select",
          },
          {
            name: "remove_chimeras",
            value: true,
            disabled: "never",
            tooltip: "perform chimera removal with UCHIME3 de novo algoritm",
            type: "bool",
          },
          {
            name: "abskew",
            value: 16,
            disabled: "never",
            tooltip:
              "the abundance skew of chimeric sequences in comparsion with parental sequences (by default, parents should be at least 16 times more abundant than their chimera)",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "similarity_type",
            items: ["0", "1", "2", "3", "4"],
            value: "2",
            disabled: "never",
            tooltip: "pairwise sequence identity definition (--iddef)",
            type: "select",
          },
          {
            name: "maxaccepts",
            value: 1,
            disabled: "never",
            tooltip:
              "maximum number of hits to accept before stopping the search",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "maxrejects",
            value: 32,
            disabled: "never",
            tooltip:
              "maximum number of non-matching target sequences to consider before stopping the search",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "mask",
            items: ["dust", "none"],
            value: "dust",
            disabled: "never",
            tooltip:
              'mask regions in sequences using the "dust" method, or do not mask ("none").',
            type: "select",
          },
          {
            name: "cores",
            value: 4,
            disabled: "never",
            tooltip: "number of cores to use",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
        ],
        Inputs: [
          {
            name: "zOTUs_thresh",
            value: 1,
            disabled: "never",
            tooltip:
              "sequence similarity threshold for zOTU table creation; 1 = 100% similarity threshold for zOTUs",
            max: 1,
            min: 0,
            step: 0.01,
            type: "slide",
          },
          {
            name: "similarity_threshold",
            value: 1,
            disabled: "never",
            tooltip:
              "cluster zOTUs to OTUs based on the sequence similarity threshold; if id = 1, no OTU clustering will be performed",
            max: 1,
            min: 0,
            step: 0.01,
            type: "slide",
          },
          {
            name: "strands",
            items: ["both", "plus"],
            disabled: "never",
            tooltip:
              "when comparing sequences with the cluster seed, check both strands (forward and reverse complementary) or the plus strand only",
            value: "both",
            type: "select",
          },
          {
            name: "minsize",
            value: 8,
            disabled: "never",
            tooltip: "minimum abundance of sequences for denoising",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
        ],
      },
      {
        tooltip:
          "assign taxonomy to zOTUs with BLAST against a selected database",
        scriptName: "taxonomy_BLAST.sh",
        imageName: "pipecraft/blast:2.14",
        serviceName: "assign taxonomy",
        manualLink:
          "https://pipecraft2-manual.readthedocs.io/en/latest/quicktools.html#blast-camacho-et-al-2009",
        selected: false,
        showExtra: false,
        extraInputs: [
          {
            name: "e_value",
            value: 10,
            disabled: "never",
            tooltip:
              "a parameter that describes the number of hits one can expect to see by chance when searching a database of a particular size. The lower the e-value the more 'significant' the match is",
            type: "numeric",
            default: 10,
            rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
          },
          {
            name: "word_size",
            value: 11,
            disabled: "never",
            tooltip:
              "the size of the initial word that must be matched between the database and the query sequence",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "reward",
            value: 2,
            disabled: "never",
            tooltip: "reward for a match",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
          },
          {
            name: "penalty",
            value: -3,
            disabled: "never",
            tooltip: "penalty for a mismatch",
            type: "numeric",
            rules: [(v) => v <= 0 || "ERROR: specify values <= 0"],
          },
          {
            name: "gap_open",
            value: 5,
            disabled: "never",
            tooltip: "cost to open a gap",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
          },
          {
            name: "gap_extend",
            value: 2,
            disabled: "never",
            tooltip: "cost to extend a gap",
            type: "numeric",
            // rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
          },
          {
            name: "cores",
            value: 4,
            disabled: "never",
            tooltip: "number of cores to use",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
        ],
        Inputs: [
          {
            name: "database_file",
            btnName: "select file",
            value: "undefined",
            disabled: "never",
            tooltip:
              "database file (may be fasta formated - automatically will convert to BLAST database format)",
            type: "file",
          },
          {
            name: "task",
            items: ["blastn", "megablast"],
            value: "blastn",
            disabled: "never",
            tooltip: "task (blastn or megablast)",
            type: "select",
          },
          {
            name: "strands",
            items: ["plus", "both"],
            value: "both",
            disabled: "never",
            tooltip:
              "query strand to search against database. Both = search also reverse complement",
            type: "select",
          },
        ],
      },
    ],
    // Metaworks_COI: [
    //   {
    //     tooltip:
    //       "MetaWorks v1.12.0 ASVs workflow for Illumina (paired-end) COI amplicons",
    //     scriptName: "metaworks_paired_end_ASV.sh",
    //     imageName: "pipecraft/metaworks:1.12.0",
    //     serviceName: "metaworks_COI",
    //     disabled: "never",
    //     selected: "always",
    //     showExtra: false,
    //     extraInputs: [
    //       {
    //         name: "quality_cutoff",
    //         value: 13,
    //         disabled: "never",
    //         tooltip:
    //           "Assemble paired-end reads setting. Phred score quality cutoff (default 20)",
    //         type: "numeric",
    //         rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
    //       },
    //       {
    //         name: "min_overlap",
    //         value: 25,
    //         disabled: "never",
    //         tooltip:
    //           "Assemble paired-end reads setting. Minimum overlap (bp) length between R1 and R2 reads when merging reads",
    //         type: "numeric",
    //         rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
    //       },
    //       {
    //         name: "match_fraction",
    //         value: 0.9,
    //         disabled: "never",
    //         tooltip:
    //           "Assemble paired-end reads setting. Minimum fraction of matching overlap (default 0.90) when merging reads",
    //         type: "numeric",
    //         rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
    //       },
    //       {
    //         name: "mismatch_fraction",
    //         value: [0.02],
    //         disabled: "never",
    //         tooltip:
    //           "Assemble paired-end reads setting. Maximum fraction of mismatches allowed in overlap (default 0.02) when merging reads",
    //         type: "numeric",
    //         rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
    //       },
    //       {
    //         name: "primer_mismatch",
    //         value: 1,
    //         disabled: "never",
    //         tooltip:
    //           "CUT PRIMERS setting. Maximum number of mismatches when searching and clipping primers",
    //         type: "numeric",
    //         rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
    //       },
    //       {
    //         name: "primer_overlap",
    //         value: 15,
    //         disabled: "never",
    //         tooltip:
    //           "CUT_PRIMERS setting. Minimum overlap to primer sequence when searching and clipping primers",
    //         type: "numeric",
    //         rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
    //       },
    //       {
    //         name: "min_seq_len",
    //         value: 150,
    //         disabled: "never",
    //         tooltip:
    //           "CUT_PRIMERS setting. Minimum sequence length (bp) to retain after trimming primers",
    //         type: "numeric",
    //         rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
    //       },
    //       {
    //         name: "qual_cutoff_3end",
    //         value: 20,
    //         disabled: "never",
    //         tooltip:
    //           "QUALITY FILT setting. Phred quality score cutoffs at the 3' end during quality filtering",
    //         type: "numeric",
    //         rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
    //       },
    //       {
    //         name: "qual_cutoff_5end",
    //         value: 20,
    //         disabled: "never",
    //         tooltip:
    //           "QUALITY FILT setting. Phred quality score cutoffs at the 5' end during quality filtering",
    //         type: "numeric",
    //         rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
    //       },
    //       {
    //         name: "maxNs",
    //         value: 0,
    //         disabled: "never",
    //         tooltip: "QUALITY FILT setting. Maximum number of Ns in the read",
    //         type: "numeric",
    //         rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
    //       },
    //       // denoise
    //       {
    //         name: "minsize",
    //         value: 8,
    //         disabled: "never",
    //         tooltip:
    //           "UNOISE denoising setting. minimum number of reads per cluster to retain (default 8)",
    //         type: "numeric",
    //         rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
    //       },
    //       {
    //         name: "pseudogene_filtering",
    //         value: true,
    //         disabled: "never",
    //         tooltip:
    //           "Filter out putative pseudogenes based on unusually short/long open reading frames (uses ORFfinder)",
    //         type: "bool",
    //       },
    //       {
    //         name: "genetic_code",
    //         items: ["1", "2", "5"],
    //         value: "5",
    //         disabled: "never",
    //         tooltip:
    //           "Pseudogene filtering setting. Genetic code translation table: 1 = standard code (use for rbcL); 2 = vertebrate mitochondrial (use for COI if targeting vertebrates); 5 = invertebrate mitochondrial (use for COI if targeting invertebrates)",
    //         type: "select",
    //       },
    //       {
    //         name: "orf_len",
    //         value: 75,
    //         disabled: "never",
    //         tooltip:
    //           "Pseudogene filtering setting. Minimum length of an open reading frame",
    //         type: "numeric",
    //         rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
    //       },
    //       {
    //         name: "cores",
    //         value: 4,
    //         disabled: "never",
    //         tooltip: "number of cores to use (for vsearch)",
    //         type: "numeric",
    //         rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
    //       },
    //     ],
    //     Inputs: [
    //       {
    //         name: "filename_structure",
    //         value: ["{sample}.R{read}"],
    //         disabled: "single_end",
    //         tooltip:
    //           "specify the sample filename structure. E.g. 'mysample1.R1.fastq' = {sample}.R{read}; 'mysample1_L001_R1_001.fastq' = {sample}_L001_R{read}_001",
    //         type: "chip",
    //         rules: [(v) => v.length <= 1 || "ADD ONLY ONE IDENTIFIER"],
    //       },
    //       // {
    //       //   name: "marker",
    //       //   items: [
    //       //     "16S",
    //       //     "18S_eukaryota",
    //       //     "18S_diatom",
    //       //     "12S_fish",
    //       //     "12S_vertebrate",
    //       //     "ITS_fungi",
    //       //     "28S_fungi",
    //       //     "rbcL_eukaryota",
    //       //     "rbcL_diatom",
    //       //     "rbcL_landPlant",
    //       //     "ITS_plants",
    //       //     "COI",
    //       //   ],
    //       //   value: "COI",
    //       //   disabled: "never",
    //       //   tooltip: "Which marker classifier will you be using?",
    //       //   type: "select",
    //       // },
    //       // {
    //       //   name: "ITS_region",
    //       //   items: ["ITS1", "ITS2"],
    //       //   value: "ITS2",
    //       //   disabled: "never",
    //       //   tooltip:
    //       //     "when marker = ITS, specify which region to extract with ITSx (if using other marker, then this setting is ignored)",
    //       //   type: "select",
    //       // },
    //       {
    //         name: "forward_primers",
    //         value: ["GGWACWGGWTGAACWGTWTAYCCYCC"],
    //         disabled: "never",
    //         tooltip: "specify forward primer (5'-3'); add up to 13 primers",
    //         type: "chip",
    //         iupac: true,
    //         rules: [(v) => v.length <= 13 || "TOO MANY PRIMERS"],
    //       },
    //       {
    //         name: "reverse_primers",
    //         value: ["TANACYTCNGGRTGNCCRAARAAYCA"],
    //         disabled: "never",
    //         tooltip: "specify reverse primer (3'-5'); add up to 13 primers",
    //         type: "chip",
    //         iupac: true,
    //         rules: [(v) => v.length <= 13 || "TOO MANY PRIMERS"],
    //       },
    //       {
    //         name: "database",
    //         btnName: "select RDP db",
    //         value: "undefined",
    //         disabled: "never",
    //         tooltip:
    //           "RDP-trained reference database for the RDP classifier. Click on the header to download trained reference databases the RDP classifier, link MetaWorks user guide: https://terrimporter.github.io/MetaWorksSite/#classifier_table",
    //         type: "file",
    //       },
    //     ],
    //   },
    // ],
    NextITS: [
      {
        tooltip:
          "Settings for STEP_1 (sequence filtering processes per sequencing run) in NextITS pipeline",
        scriptName: "",
        imageName: "vmikk/nextits:0.5.0",
        serviceName: "Step_1",
        manualLink: "https://next-its.github.io/parameters/#step-1",
        disabled: "never",
        selected: "always",
        showExtra: false,
        extraInputs: [
          {
            name: "qc_maxee",
            value: 1,
            disabled: "never",
            tooltip: "Maximum number of expected errors",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "qc_maxhomopolymerlen",
            value: 25,
            disabled: "never",
            tooltip:
              "Threshold for a homopolymer region length in a sequence (default, 25)",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "qc_maxn",
            value: 4,
            disabled: "never",
            tooltip:
              "Discard sequences with more than the specified number of N’s",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "ITSx_evalue",
            value: (0.01).toExponential(),
            disabled: "never",
            tooltip: "ITSx E-value cutoff threshold (default, 1e-1)",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify only values > 0"],
          },
          {
            name: "ITSx_partial",
            value: 0,
            disabled: "never",
            tooltip:
              "Keep partial ITS sequences (defalt, off), otherwise specify min length cutoff",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "chimera_database",
            active: false,
            btnName: "select file",
            value: "undefined",
            disabled: "never",
            tooltip: "Database for reference-based chimera removal",
            type: "boolfile",
          },
          {
            name: "ITSx_tax",
            items: [
              "all",
              "alveolata",
              "bryophyta",
              "bacillariophyta",
              "amoebozoa",
              "euglenozoa",
              "fungi",
              "chlorophyta",
              "rhodophyta",
              "phaeophyceae",
              "marchantiophyta",
              "metazoa",
              "oomycota",
              "haptophyceae",
              "raphidophyceae",
              "rhizaria",
              "synurophyceae",
              "tracheophyta",
              "eustigmatophyceae",
              "apusozoa",
              "parabasalia",
            ],
            value: ["all"],
            disabled: "never",
            tooltip: "ITSx taxonomy profile (default, 'all')",
            type: "combobox",
          },
          {
            name: "chimera_rescueoccurrence",
            value: 2,
            disabled: "never",
            tooltip:
              "Min occurrence of chimeric sequences required to rescue them (default, 2)",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "tj_f ",
            value: 0.01,
            disabled: "never",
            tooltip:
              "Tag-jump filtering, UNCROSS parameter `f` (default, 0.01)",
            max: 1,
            min: 0,
            step: 0.01,
            type: "slide",
          },
          {
            name: "tj_p",
            value: 1,
            disabled: "never",
            tooltip: "Tag-jump filtering parameter `p` (default, 1)",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "hp",
            value: true,
            disabled: "never",
            tooltip: "Homopolymer compression (default, true)",
            type: "bool",
          },
        ],
        Inputs: [
          {
            name: "primer_forward",
            value: ["GTACACACCGCCCGTCG"],
            disabled: "never",
            tooltip: "specify forward primer",
            type: "chip",
            iupac: true,
            rules: [(v) => v.length <= 1 || "TOO MANY PRIMERS"],
          },
          {
            name: "primer_reverse",
            value: ["CCTSCSCTTANTDATATGC"],
            disabled: "never",
            tooltip: "specify reverse primer",
            type: "chip",
            iupac: true,
            rules: [(v) => v.length <= 1 || "TOO MANY PRIMERS"],
          },
          {
            name: "primer_mismatches",
            value: 2,
            disabled: "never",
            tooltip: "Maximum number of mismatches when searching for primers",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "its_region",
            items: ["full", "ITS1", "ITS2"],
            value: "full",
            disabled: "never",
            tooltip: "sub-regions of the internal transcribed spacer",
            type: "select",
          },
        ],
      },
      {
        tooltip: "Settings for STEP_2 (clustering) in NextITS pipeline",
        scriptName: "",
        imageName: "vmikk/nextits:0.5.0",
        serviceName: "Step_2",
        manualLink: "https://next-its.github.io/parameters/#step-2",
        disabled: "never",
        selected: "always",
        showExtra: false,
        extraInputs: [
          {
            name: "otu_iddef",
            value: 2,
            disabled: "never",
            tooltip:
              "Sequence similarity definition for tag-jump removal step (default, 2)",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "otu_qmask",
            items: ["dust", "none"],
            value: "dust",
            disabled: "never",
            tooltip:
              'mask regions in sequences using the "dust" method, or do not mask ("none").',
            type: "select",
          },
          {
            name: "swarm_fastidious",
            value: true,
            disabled: "never",
            tooltip: "Link nearby low-abundance swarms (fastidious option)",
            type: "bool",
            depends_on:
              'state.NextITS[1].Inputs[0].value == "swarm" && state.NextITS[1].Inputs[2].value <= 1',
          },
          {
            name: "unoise_alpha",
            value: 2,
            disabled: "never",
            tooltip: "Alpha parameter of UNOISE",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
            depends_on: "state.NextITS[1].Inputs[4].value == true",
          },
          {
            name: "unoise_minsize",
            value: 8,
            disabled: "never",
            tooltip: "Minimum sequence abundance ",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
            depends_on: "state.NextITS[1].Inputs[4].value == true",
          },
          {
            name: "max_MEEP",
            value: 0.5,
            disabled: "never",
            tooltip: "Maximum allowed number of expected errors per 100 bp",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify values >= 1"],
          },
          {
            name: "max_ChimeraScore",
            value: 0.5,
            disabled: "never",
            tooltip: "Maximum allowed de novo chimera score",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify values >= 1"],
          },
          {
            name: "lulu_match",
            value: 95,
            disabled: "never",
            tooltip: "Minimum similarity threshold",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify values >= 1"],
            depends_on: "state.NextITS[1].Inputs[3].value == true",
          },
          {
            name: "lulu_ratio",
            value: 1,
            disabled: "never",
            tooltip: "Minimum abundance ratio",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify values >= 1"],
            depends_on: "state.NextITS[1].Inputs[3].value == true",
          },
          {
            name: "lulu_ratiotype",
            items: ["min", "avg"],
            value: "min",
            disabled: "never",
            tooltip: "Abundance ratio type - 'min' or 'avg'	",
            type: "select",
            depends_on: "state.NextITS[1].Inputs[3].value == true",
          },
          {
            name: "lulu_relcooc",
            value: 0.95,
            disabled: "never",
            tooltip: "Relative co-occurrence",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify values >= 1"],
            depends_on: "state.NextITS[1].Inputs[3].value == true",
          },
          {
            name: "lulu_maxhits",
            value: 0,
            disabled: "never",
            tooltip: "Maximum number of hits (0 = unlimited)",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify values >= 1"],
            depends_on: "state.NextITS[1].Inputs[3].value == true",
          },
        ],
        Inputs: [
          {
            name: "clustering_method",
            items: ["vsearch", "swarm", "unoise"],
            value: "vsearch",
            disabled: "never",
            tooltip: "Sequence clustering method",
            type: "select",
          },
          {
            name: "otu_id",
            value: 0.98,
            disabled: "never",
            tooltip: "Sequence similarity for OTU clustering (default, 0.98)",
            max: 1,
            min: 0,
            step: 0.01,
            type: "slide",
          },
          {
            name: "swarm_d",
            value: 1,
            disabled: "never",
            tooltip: "SWARM clustering resolution (d)",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify values >= 1"],
            depends_on: 'state.NextITS[1].Inputs[0].value == "swarm"',
          },
          {
            name: "lulu",
            value: true,
            disabled: "never",
            tooltip: "Run post-clustering curation with LULU",
            type: "bool",
          },
          {
            name: "unoise",
            value: false,
            disabled: "never",
            tooltip: "Perform denoising with UNOISE algorithm",
            type: "bool",
          },
        ],
      },
    ],
    DADA2_ASVs: [
      {
        tooltip: "remove primers sequences from the reads",
        scriptName: {
          FORWARD: "cut_primers_paired_end_reads.sh",
          MIXED: "cut_mixed_primers_paired_end_reads.sh",
          SINGLE_END: "cut_primers_single_end_reads.sh",
        },
        imageName: "pipecraft/cutadapt:4.4",
        serviceName: "cut primers",
        manualLink:
          "https://pipecraft2-manual.readthedocs.io/en/latest/quicktools.html#cut-primers",
        disabled: "never",
        selected: false,
        showExtra: false,
        extraInputs: [
          {
            name: "cores",
            value: 1,
            disabled: "never",
            tooltip:
              "number of cores to use. For paired-end data in fasta format, set to 1 [default]. For fastq formats you may set the value to 0 to use all cores",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "no_indels",
            value: true,
            disabled: "never",
            tooltip:
              "do not allow insertions or deletions is primer search. Mismatches are the only type of errors accounted in the error rate parameter. ",
            type: "bool",
          },
        ],
        Inputs: [
          {
            name: "forward_primers",
            value: [],
            disabled: "never",
            tooltip: "specify forward primer (5'-3'); add up to 13 primers",
            type: "chip",
            iupac: true,
            rules: [(v) => v.length <= 13 || "TOO MANY PRIMERS"],
          },
          {
            name: "reverse_primers",
            value: [],
            disabled: "never",
            tooltip: "specify reverse primer (3'-5'); add up to 13 primers",
            type: "chip",
            iupac: true,
            rules: [(v) => v.length <= 13 || "TOO MANY PRIMERS"],
          },
          {
            name: "mismatches",
            value: 1,
            disabled: "never",
            tooltip: "allowed mismatches in the primer search",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
          },
          {
            name: "min_overlap",
            value: 21,
            disabled: "never",
            tooltip:
              "number of overlap bases with the primer sequence. Partial matches are allowed, but short matches may occur by chance, leading to erroneously clipped bases. Specifying higher overlap than the length of primer sequnce will still clip the primer (e.g. primer length is 22 bp, but overlap is specified as 25 - this does not affect the identification and clipping of the primer as long as the match is in the specified mismatch error range)",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "seqs_to_keep",
            items: ["keep_all", "keep_only_linked"],
            value: "keep_all",
            disabled: "never",
            tooltip:
              "'keep_all' keeps the sequences where at least one primer was found (fwd or rev); recommended when cutting primers from paired-end data (unassembled), where individual R1 or R2 read lenghts are shorther than the expected amplicon length. 'keep_only_linked' = keep sequences if primers are found in both ends (fwd…rev); discards the read if both primers were not found in this read; maybe useful for single-end data",
            type: "select",
          },
          {
            name: "pair_filter",
            items: ["both", "any"],
            value: "both",
            disabled: "single_end",
            tooltip:
              "applies only for paired-end data. Option 'both' discards a read-pair when both reads do not contain a primer sequence. Option 'any' discards the read-pair when one of the reads does not contain a primer sequence",
            type: "select",
          },
        ],
      },
      {
        tooltip: "quality filtering with DADA2 'filterAndTrim' function",
        scriptName: {
          FORWARD: "quality_filtering_paired_end_dada2.sh",
          MIXED: "quality_filtering_paired_end_dada2_mixed.sh",
          SINGLE_END: "quality_filtering_single_end_dada2.sh",
        },
        imageName: "pipecraft/vsearch_dada2:2",
        serviceName: "quality filtering",
        manualLink:
          "https://pipecraft2-manual.readthedocs.io/en/latest/quicktools.html#quality-filtering",
        disabled: "never",
        selected: "always",
        showExtra: false,
        extraInputs: [],
        Inputs: [
          {
            name: "read_R1",
            value: ["_R1"],
            disabled: "single_end",
            tooltip:
              "[only for paired-end data] identifyer string that is common for all R1 reads. R1/R2 reads must contain R1/R2 strings in the file names; sample names must not contain R1/R2! (e.g. when all R1 files have '.R1' string, then enter '\\.R1'. Note that backslash is only needed to escape dot regex; e.g. when all R1 files have '_R1' string, then enter '_R1'.)",
            type: "chip",
            rules: [(v) => v.length <= 1 || "ADD ONLY ONE IDENTIFIER"],
          },
          {
            name: "read_R2",
            value: ["_R2"],
            disabled: "single_end",
            tooltip:
              "[only for paired-end data] identifyer string that is common for all R2 reads (e.g. when all R2 files have '.R2' string, then enter '\\.R2'. Note that backslash is only needed to escape dot regex; e.g. when all R2 files have '_R1' string, then enter '_R2'.)",
            type: "chip",
            rules: [(v) => v.length <= 1 || "ADD ONLY ONE IDENTIFIER"],
          },
          // {
          //   name: "samp_ID",
          //   value: ["_"],
          //   disabled: "single_end",
          //   tooltip:
          //     "[only for paired-end data] identifyer string that separates the sample name from redundant charachters (e.g. file name = sample1.R1.fastq, then '\\.' would be the 'identifier string' (sample name = sample1)); note that backslash is only needed to escape dot regex (e.g. when file name = sample1_R1.fastq then specify as '_')",
          //   type: "chip",
          //   rules: [(v) => v.length <= 1 || "ADD ONLY ONE IDENTIFIER"],
          // },
          {
            name: "maxEE",
            value: 2,
            disabled: "never",
            tooltip:
              "discard sequences with more than the specified number of expected errors",
            type: "numeric",
            rules: [(v) => v >= 0.1 || "ERROR: specify values >= 0.1"],
          },
          {
            name: "maxN",
            value: 0,
            disabled: "never",
            tooltip:
              "discard sequences with more than the specified number of Ns (ambiguous bases)",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
          },
          {
            name: "minLen",
            value: 20,
            disabled: "never",
            tooltip:
              "remove reads with length less than minLen. minLen is enforced after all other trimming and truncation",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "truncQ",
            value: 2,
            disabled: "never",
            tooltip:
              "truncate reads at the first instance of a quality score less than or equal to truncQ",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "truncLen",
            value: 0,
            disabled: "never",
            tooltip:
              "truncate reads after truncLen bases (applies to R1 reads when working with paired-end data). Reads shorter than this are discarded. Explore quality profiles (with QualityCheck module) see whether poor quality ends needs to truncated",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "truncLen_R2",
            value: 0,
            disabled: "single_end",
            tooltip:
              "truncate R2 reads after truncLen bases. Reads shorter than this are discarded. Explore quality profiles (with QualityCheck module) see whether poor quality ends needs to truncated",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "maxLen",
            value: 9999,
            disabled: "never",
            tooltip:
              "remove reads with length greater than maxLen. maxLen is enforced on the raw reads. In dada2, the default = Inf, but here set as 9999",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "minQ",
            value: 0,
            disabled: "never",
            tooltip:
              "after truncation, reads contain a quality score below minQ will be discarded",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
          },
          {
            name: "matchIDs",
            value: true,
            disabled: "single_end",
            tooltip:
              "applies only for paired-end data. If TRUE, then double-checking (with seqkit pair) that only paired reads that share ids are outputted",
            type: "bool",
          },
        ],
      },
      {
        tooltip: "select the denoising options for DADA2 'dada' function",
        scriptName: {
          FORWARD: "assemble_paired_end_data_dada2_wf.sh",
          MIXED: "assemble_paired_end_data_dada2_mixed_wf.sh",
          SINGLE_END: "denoise_single_end_data_dada2.sh",
        },
        imageName: "pipecraft/vsearch_dada2:2",
        serviceName: "denoise",
        manualLink:
          "https://pipecraft2-manual.readthedocs.io/en/latest/quicktools.html#dada2-denoise",
        selected: "always",
        disabled: "never",
        showExtra: false,
        extraInputs: [
          {
            name: "BAND_SIZE",
            value: 16,
            disabled: "never",
            tooltip:
              "default = 16. Banding for Needleman-Wunsch alignments. Default sets here = 16 for Illumina and 32 for PacBio",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
          },
          {
            name: "OMEGA_A",
            value: (0.0000000000000000000000000000000000000001).toExponential(),
            disabled: "never",
            tooltip:
              "default = 1e-40. Denoising setting; see DADA2 'setDadaOpt()' for detalis. Default value  is a conservative setting to avoid making false positive inferences, but comes at the cost of reducing the ability to identify some rare variants",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify values > 0"],
          },
          {
            name: "OMEGA_P",
            value: (0.0001).toExponential(),
            disabled: "never",
            tooltip:
              "default = 1e-4 (0.0001). Denoising setting; see DADA2 'setDadaOpt()' for detalis",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify values > 0"],
          },
          {
            name: "OMEGA_C",
            value: (0.0000000000000000000000000000000000000001).toExponential(),
            disabled: "never",
            tooltip:
              "default = 1e-40. Denoising setting; see DADA2 'setDadaOpt()' for detalis",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify values > 0"],
          },
          {
            name: "DETECT_SINGLETONS",
            disabled: "never",
            value: false,
            tooltip:
              "Default = FALSE. Denoising setting; see DADA2 'setDadaOpt()' for detalis. If set to TRUE, this removes the requirement for at least two reads with the same sequences to exist in order for a new ASV to be detected",
            type: "bool",
          },
        ],
        Inputs: [
          {
            name: "errorEstFun",
            items: ["PacBioErrfun", "loessErrfun"],
            value: "loessErrfun",
            disabled: "never",
            tooltip:
              "DADA2 errorEstimationFunction. 'loessErrfun' for Illumina data; 'PacBioErrfun' for PacBio data",
            type: "select",
          },
          {
            name: "pool",
            items: ["FALSE", "TRUE", "psuedo"],
            value: "FALSE",
            disabled: "never",
            tooltip:
              "if pool = TRUE, the algorithm will pool together all samples prior to sample inference. Pooling improves the detection of rare variants, but is computationally more expensive. If pool = 'pseudo', the algorithm will perform pseudo-pooling between individually processed samples. This argument has no effect if only 1 sample is provided, and pool does not affect error rates, which are always estimated from pooled observations across samples",
            type: "select",
          },
          // {
          //   name: "selfConsist",
          //   disabled: "never",
          //   value: false,
          //   tooltip:
          //     "if selfConsist = TRUE, the algorithm will alternate between sample inference and error rate estimation until convergence",
          //   type: "bool",
          // },
          {
            name: "qualityType",
            items: ["Auto", "FastqQuality"],
            value: "Auto",
            disabled: "never",
            tooltip:
              "Auto means to attempt to auto-detect the fastq quality encoding. This may fail for PacBio files with uniformly high quality scores, in which case use 'FastqQuality'",
            type: "select",
          },
        ],
      },
      {
        tooltip:
          "assemble paired-end reads (R1 and R2) with DADA2 'mergePairs' function",
        scriptName: {
          FORWARD: "assemble_paired_end_data_dada2_wf.sh",
          MIXED: "assemble_paired_end_data_dada2_mixed_wf.sh",
          SINGLE_END: "disabled.sh",
        },
        imageName: "pipecraft/vsearch_dada2:2",
        serviceName: "merge Pairs",
        manualLink:
          "https://pipecraft2-manual.readthedocs.io/en/latest/quicktools.html#assemble-paired-end-reads",
        selected: "always",
        disabled: "single_end",
        showExtra: false,
        extraInputs: [],
        Inputs: [
          {
            name: "minOverlap",
            value: 12,
            disabled: "never",
            tooltip:
              "the minimum length of the overlap required for merging the forward and reverse reads",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "maxMismatch",
            value: 0,
            disabled: "never",
            tooltip: "the maximum mismatches allowed in the overlap region",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
          },
          {
            name: "trimOverhang",
            value: false,
            disabled: "never",
            tooltip:
              "if TRUE, overhangs in the alignment between the forwards and reverse read are trimmed off. Overhangs are when the reverse read extends past the start of the forward read, and vice-versa, as can happen when reads are longer than the amplicon and read into the other-direction primer region",
            type: "bool",
          },
          {
            name: "justConcatenate",
            value: false,
            disabled: "never",
            tooltip:
              "if TRUE, the forward and reverse-complemented reverse read are concatenated rather than merged, with a NNNNNNNNNN (10 Ns) spacer inserted between them",
            type: "bool",
          },
        ],
      },
      {
        tooltip: "remove chimeras with DADA2 'removeBimeraDenovo' function",
        scriptName: {
          FORWARD: "chimera_filtering_dada2_wf.sh",
          MIXED: "chimera_filtering_dada2_mixed_wf.sh",
          SINGLE_END: "chimera_filtering_dada2_wf.sh",
        },
        imageName: "pipecraft/vsearch_dada2:2",
        serviceName: "chimera filtering",
        manualLink:
          "https://pipecraft2-manual.readthedocs.io/en/latest/quicktools.html#chimera-filtering",
        disabled: "never",
        selected: "always",
        showExtra: false,
        extraInputs: [],
        Inputs: [
          {
            name: "method",
            items: ["consensus", "pooled", "per-sample"],
            value: "consensus",
            disabled: "never",
            tooltip:
              "'consensus' - the samples are independently checked for chimeras, and a consensus decision on each sequence variant is made. If 'pooled', the samples are all pooled together for chimera identification. If 'per-sample', the samples are independently checked for chimeras",
            type: "select",
          },
        ],
      },
      {
        scriptName: "tag_jump_removal.sh",
        tooltip: "using UNCROSS2 to filter out putative tag-jumps",
        imageName: "vmikk/nextits:0.5.0",
        serviceName: "filter tag-jumps",
        manualLink:
          "https://pipecraft2-manual.readthedocs.io/en/latest/pre-defined_pipelines.html#tag-jump-correction",
        selected: false,
        showExtra: false,
        extraInputs: [],
        Inputs: [
          {
            name: "f_value",
            value: 0.01,
            disabled: "never",
            tooltip:
              "f-parameter of UNCROSS2, which defines the expected cross-talk rate. Default is 0.01 (equivalent to 1%). A higher value enforces stricter filtering",
            type: "numeric",
            rules: [(v) => v > 0 || "ERROR: specify values > 0"],
          },
          {
            name: "p_value",
            value: 1,
            disabled: "never",
            tooltip:
              "p-parameter, which controls the severity of tag-jump removal. It adjusts the exponent in the UNCROSS formula. Default is 1. Opt for 0.5 or 0.3 to steepen the curve",
            type: "numeric",
            rules: [(v) => v > 0 || "ERROR: specify values > 0"],
          },
        ],
      },
      {
        tooltip: "Collaplse identical ASVs or/and filter ASVs by length",
        scriptName: "table_filtering_dada2_wf.sh",
        imageName: "pipecraft/vsearch_dada2:2",
        serviceName: "filter ASV table",
        manualLink:
          "https://pipecraft2-manual.readthedocs.io/en/latest/quicktools.html#dada2-table-filtering",
        disabled: "never",
        selected: false,
        showExtra: false,
        extraInputs: [
          {
            name: "minOverlap",
            value: 20,
            disabled: "never",
            tooltip:
              "collapseNoMismatch setting. Default = 20. The minimum overlap of base pairs between ASV sequences required to collapse them together",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "vec",
            value: true,
            disabled: "never",
            tooltip:
              "collapseNoMismatch setting. Default = TRUE. Use the vectorized aligner. Should be turned off if sequences exceed 2kb in length",
            type: "bool",
          },
        ],
        Inputs: [
          {
            name: "collapseNoMismatch",
            value: true,
            disabled: "never",
            tooltip:
              "collapses ASVs that are identical up to shifts or length variation, i.e. that have no mismatches or internal indels (dada2 'collapseNoMismatch')",
            type: "bool",
          },
          {
            name: "by_length",
            value: 249,
            disabled: "never",
            tooltip:
              "discard ASVs from the ASV table that are shorter than specified value (in base pairs). Value 0 means OFF, no filtering by length",
            type: "numeric",
            rules: [(v) => v >= 0 || "ERROR: specify values >= 0"],
          },
        ],
      },
      {
        tooltip:
          "assign taxonomy with DADA2 'assignTaxonomy' function against a selected database. Untick the checkbox to skip this step",
        scriptName: "taxonomy_dada2.sh",
        imageName: "pipecraft/vsearch_dada2:2",
        serviceName: "assign Taxonomy",
        manualLink:
          "https://pipecraft2-manual.readthedocs.io/en/latest/quicktools.html#assign-taxonomy-dada2",
        disabled: "never",
        selected: false,
        showExtra: false,
        extraInputs: [],
        Inputs: [
          {
            name: "dada2_database",
            btnName: "select fasta",
            value: "undefined",
            disabled: "never",
            tooltip:
              "Select a reference database fasta file for taxonomy annotation. Click on the header to download DADA2-formatted reference databases https://benjjneb.github.io/dada2/training.html",
            type: "file",
          },
          {
            name: "minBoot",
            value: 50,
            disabled: "never",
            tooltip:
              "the minimum bootstrap confidence for assigning a taxonomic level",
            type: "numeric",
            rules: [(v) => v >= 1 || "ERROR: specify values >= 1"],
          },
          {
            name: "tryRC",
            value: false,
            disabled: "never",
            tooltip:
              "the reverse-complement of each sequences will be used for classification if it is a better match to the reference sequences than the forward sequence",
            type: "bool",
          },
        ],
      },
    ],
    customWorkflowInfo: {
      vsearch_OTUs: {
        info: "vsearch OTUs workflow for for demultiplexed Illumina data",
        link: "https://github.com/torognes/vsearch",
        title: "vsearch OTUs workflow",
      },
      DADA2_ASVs: {
        info: "DADA2 ASVs workflow for for demultiplexed Illumina or PacBio data. FORWARD: [for Illumina paired-end data] select this when all reads of interest are expected to be in 5-3 orient. MIXED: [for Illumina paired-end data] select this when reads of interest are expected to be in both 5-3 and 3-5 orient. PACBIO: [for PacBio data] select this when working with PacBio data",
        link: "https://benjjneb.github.io/dada2/index.html",
        title: "DADA2 ASVs workflow",
      },
      UNOISE_ASVs: {
        info: "vsearch ASVs (zOTUs) workflow for for demultiplexed Illumina data",
        link: "https://github.com/torognes/vsearch",
        title: "UNOISE3 ASVs workflow",
      },
      // Metaworks_COI: {
      //   info: "MetaWorks ASVs workflow for demultiplexed Illumina COI amplicons (paired-end)",
      //   link: "https://terrimporter.github.io/MetaWorksSite/quickstart/",
      //   title: "MetaWorks COI ASVs",
      // },
      NextITS: {
        info: "NextITS pipeline for demultiplexed PacBio ITS (single-end) amplicons",
        link: "https://github.com/vmikk/NextITS",
        title: "NextITS",
      },
    },
  },
  getters: {
    dada2modeIndex: (state) => {
      if (state.data.dada2mode == "FORWARD") {
        return 0;
      } else if (state.data.dada2mode == "MIXED") {
        return 1;
      } else if (state.data.dada2mode == "SINGLE_END") {
        return 2;
      }
    },
    check_depends_on: (state) => (input) => {
      if (input.depends_on && state) {
        return !eval(input.depends_on);
      } else {
        return false;
      }
    },
    linkify: () => (tooltip) => {
      var urlRegex = /(((https?:\/\/)|(www\.))[^\s]+)/g;
      return tooltip.match(urlRegex);
    },
    singleStepNames: (state) => state.steps.map(({ stepName }) => stepName),
    steps2Run: (state) => (id) => {
      let steps = state[id].filter(
        (el) => el.selected == true || el.selected == "always"
      );
      return steps.length;
    },
    selectedStepsReady: (state) => {
      let x = 0;
      let fileInputValues = [];
      for (let index of state.selectedSteps.entries()) {
        state.selectedSteps[index[0]].services.forEach((input) => {
          if (input.selected === true || input.selected == "always") {
            x = x + 1;
          }
        });
      }
      for (let index of state.selectedSteps.entries()) {
        state.selectedSteps[index[0]].services.forEach((input) => {
          if (input.selected === true || input.selected == "always") {
            input.Inputs.forEach((input) => {
              if (
                input.type == "file" ||
                input.type == "chip" ||
                (input.type == "boolfile" && input.active == true)
              ) {
                fileInputValues.push(input.value);
              }
            });
          }
        });
      }
      if (
        x == state.selectedSteps.length &&
        state.selectedSteps.length > 0 &&
        !fileInputValues.includes("undefined")
      ) {
        return true;
      } else {
        return false;
      }
    },
    customWorkflowReady: (state) => {
      let Ready = [];
      if (state.route.params.workflowName) {
        Ready.push(true);
        state[state.route.params.workflowName].forEach((step) => {
          if (step.selected == true || step.selected == "always") {
            step.Inputs.forEach((input) => {
              if (
                input.type == "file" ||
                input.type == "chip" ||
                (input.type == "boolfile" && input.active == true)
              ) {
                if (input.value == "undefined" || input.value.length < 1) {
                  Ready.push(false);
                }
              }
            });
          }
        });
      }
      if (Ready.includes(false)) {
        return false;
      } else {
        return true;
      }
    },
  },
  mutations: {
    activatePullLoader(state) {
      state.pullLoader.active = true;
    },
    deactivatePullLoader(state) {
      state.pullLoader.active = false;
    },
    // runInfo: { active: false, type: null, step: null, nrOfSteps: null },
    resetRunInfo(state) {
      state.runInfo = {
        active: false,
        type: null,
        step: null,
        nrOfSteps: null,
        containerID: null,
      };
    },
    addRunInfo(state, payload) {
      let result = Object.fromEntries(
        Object.keys(state.runInfo).map((_, i) => [
          Object.keys(state.runInfo)[i],
          payload[i],
        ])
      );
      state.runInfo = result;
    },
    toggle_PE_SE_scripts(state, payload) {
      for (const [key] of Object.entries(state.customWorkflowInfo)) {
        for (let i = 0; i < state[key].length; i++) {
          if (
            payload == "paired_end" &&
            typeof state[key][i].scriptName !== "object"
          ) {
            state[key][i].scriptName = state[key][i].scriptName.replace(
              "single_end",
              "paired_end"
            );
            if (state[key][i].disabled == "single_end") {
              state[key][i].selected = "always";
            }
          }
          if (
            payload == "single_end" &&
            typeof state[key][i].scriptName !== "object"
          ) {
            state[key][i].scriptName = state[key][i].scriptName.replace(
              "paired_end",
              "single_end"
            );
            if (state[key][i].disabled == "single_end") {
              state[key][i].selected = "never";
            }
          }
        }
      }
      for (let i = 0; i < state.steps.length; i++) {
        for (let j = 0; j < state.steps[i].services.length; j++) {
          if (payload == "paired_end") {
            state.steps[i].services[j].scriptName = state.steps[i].services[
              j
            ].scriptName.replace("single_end", "paired_end");
          }
          if (payload == "single_end") {
            state.steps[i].services[j].scriptName = state.steps[i].services[
              j
            ].scriptName.replace("paired_end", "single_end");
          }
        }
      }
      for (let i = 0; i < state.selectedSteps.length; i++) {
        for (let j = 0; j < state.selectedSteps[i].services.length; j++) {
          if (payload == "paired_end") {
            state.selectedSteps[i].services[j].scriptName = state.selectedSteps[
              i
            ].services[j].scriptName.replace("single_end", "paired_end");
          }
          if (payload == "single_end") {
            state.selectedSteps[i].services[j].scriptName = state.selectedSteps[
              i
            ].services[j].scriptName.replace("paired_end", "single_end");
          }
        }
      }
      for (let i = 0; i < state.selectedSteps.length; i++) {
        if (payload == "single_end") {
          state.selectedSteps = state.selectedSteps.filter(
            (item) => !(item.stepName == "assemble paired-end")
          );
          if (router.currentRoute != "/home") {
            router.push("/home").catch(() => {
              /* ignore */
            });
          }
        }
      }
    },
    toggle_demux_mux(state, payload) {
      for (const [key] of Object.entries(state.customWorkflowInfo)) {
        for (let i = 0; i < state[key].length; i++) {
          if (
            payload == "demultiplexed" &&
            state[key][i].disabled == "demultiplexed"
          ) {
            console.log(state[key][i].disabled);
            state[key][i].selected = false;
          }
          if (
            payload == "multiplexed" &&
            state[key][i].disabled == "demultiplexed"
          ) {
            state[key][i].selected = true;
          }
        }
      }
      for (let i = 0; i < state.selectedSteps.length; i++) {
        if (payload == "demultiplexed") {
          state.selectedSteps = state.selectedSteps.filter(
            (item) => !(item.stepName == "demultiplex")
          );
          if (router.currentRoute != "/home") {
            router.push("/home").catch(() => {
              /* ignore */
            });
          }
        }
      }
    },
    loadWorkflow(state, payload) {
      state.selectedSteps = payload;
    },
    loadCustomWorkflow(state, payload) {
      state[payload[1]] = payload[0];
    },
    toggleExtra(state, payload) {
      state.selectedSteps[payload.stepIndex].services[
        payload.serviceIndex
      ].showExtra =
        !state.selectedSteps[payload.stepIndex].services[payload.serviceIndex]
          .showExtra;
    },
    toggleExtraCustomWorkflow(state, payload) {
      state[payload.workflowName][payload.serviceIndex].showExtra =
        !state[payload.workflowName][payload.serviceIndex].showExtra;
    },
    addWorkingDir(state, filePath) {
      state.workingDir = filePath;
    },
    updateDockerStatus(state, payload) {
      state.dockerStatus = payload;
    },
    addInputDir(state, filePath) {
      state.inputDir = filePath;
    },
    addInputInfo(state, payload) {
      (state.data.fileFormat = payload.fileFormat),
        (state.data.readType = payload.readType);
    },
    removeStep(state, index) {
      state.selectedSteps.splice(index, 1);
    },
    addStep(state, payload) {
      state.selectedSteps = [];
      let step = _.cloneDeep(payload.step);
      state.selectedSteps.push(step);
    },
    DraggableUpdate(state, value) {
      state.selectedSteps = value;
    },
    blastSwitch(state, value) {
      if (value == "blastn") {
        state.vsearch_OTUs[6].extraInputs[1].value = 11;
        state.vsearch_OTUs[6].extraInputs[2].value = 2;
        state.vsearch_OTUs[6].extraInputs[3].value = -3;
        state.vsearch_OTUs[6].extraInputs[4].value = 5;
        state.vsearch_OTUs[6].extraInputs[5].value = 2;
        state.UNOISE_ASVs[5].extraInputs[1].value = 11;
        state.UNOISE_ASVs[5].extraInputs[2].value = 2;
        state.UNOISE_ASVs[5].extraInputs[3].value = -3;
        state.UNOISE_ASVs[5].extraInputs[4].value = 5;
        state.UNOISE_ASVs[5].extraInputs[5].value = 2;
      } else if (value == "megablast") {
        state.vsearch_OTUs[6].extraInputs[1].value = 28;
        state.vsearch_OTUs[6].extraInputs[2].value = 1;
        state.vsearch_OTUs[6].extraInputs[3].value = -2;
        state.vsearch_OTUs[6].extraInputs[4].value = 0;
        state.vsearch_OTUs[6].extraInputs[5].value = undefined;
        state.UNOISE_ASVs[5].extraInputs[1].value = 28;
        state.UNOISE_ASVs[5].extraInputs[2].value = 1;
        state.UNOISE_ASVs[5].extraInputs[3].value = -2;
        state.UNOISE_ASVs[5].extraInputs[4].value = 0;
        state.UNOISE_ASVs[5].extraInputs[5].value = undefined;
      }
    },
    blastSwitch2(state, payload) {
      if (payload.value == "blastn") {
        state.selectedSteps[payload.i1].services[
          payload.i2
        ].extraInputs[1].value = 11;
        state.selectedSteps[payload.i1].services[
          payload.i2
        ].extraInputs[2].value = 2;
        state.selectedSteps[payload.i1].services[
          payload.i2
        ].extraInputs[3].value = -3;
        state.selectedSteps[payload.i1].services[
          payload.i2
        ].extraInputs[4].value = 5;
        state.selectedSteps[payload.i1].services[
          payload.i2
        ].extraInputs[5].value = 2;
      } else if (payload.value == "megablast") {
        state.selectedSteps[payload.i1].services[
          payload.i2
        ].extraInputs[1].value = 28;
        state.selectedSteps[payload.i1].services[
          payload.i2
        ].extraInputs[2].value = 1;
        state.selectedSteps[payload.i1].services[
          payload.i2
        ].extraInputs[3].value = -2;
        state.selectedSteps[payload.i1].services[
          payload.i2
        ].extraInputs[4].value = 0;
        state.selectedSteps[payload.i1].services[
          payload.i2
        ].extraInputs[5].value = undefined;
      }
    },
    serviceInputUpdate(state, payload) {
      state.selectedSteps[payload.stepIndex].services[payload.serviceIndex] =
        payload.value;
    },
    setDADAmode(state, payload) {
      console.log(payload)
      if (payload == "SINGLE_END") {
        state.data.readType = "single_end";
        state.DADA2_ASVs[3].selected = false;
        state.DADA2_ASVs[2].Inputs[0].value = "PacBioErrfun";
        state.DADA2_ASVs[2].extraInputs[0].value = 32;
        state.DADA2_ASVs[2].Inputs[2].value = "FastqQuality";
      } else {
        state.DADA2_ASVs[3].selected = "always";
        state.DADA2_ASVs[2].Inputs[0].value = "loessErrfun";
        state.DADA2_ASVs[2].Inputs[2].value = "Auto";
        state.DADA2_ASVs[2].extraInputs[0].value = 16;
      }
      if (payload == "MIXED") {
        state.DADA2_ASVs[0].selected = "always";
      } else {
        state.DADA2_ASVs[0].selected = false;
      }
    },
    inputUpdate(state, payload) {
      state.selectedSteps[payload.stepIndex].services[payload.serviceIndex][
        payload.listName
      ][payload.inputIndex].value = payload.value;
    },
    premadeInputUpdate(state, payload) {
      state[payload.workflowName][payload.serviceIndex][payload.listName][
        payload.inputIndex
      ].value = payload.value;
    },
    toggleActive(state, payload) {
      state.selectedSteps[payload.stepIndex].services[payload.serviceIndex][
        payload.listName
      ][payload.inputIndex].active = payload.value;
      if (payload.value == false) {
        state.selectedSteps[payload.stepIndex].services[payload.serviceIndex][
          payload.listName
        ][payload.inputIndex].value = "undefined";
      }
    },
    premadeToggleActive(state, payload) {
      state[payload.workflowName][payload.serviceIndex][payload.listName][
        payload.inputIndex
      ].active = payload.value;
      if (payload.value == false) {
        state[payload.workflowName][payload.serviceIndex][payload.listName][
          payload.inputIndex
        ].value = "undefined";
      }
    },
    checkService(state, payload) {
      for (
        let index = 0;
        index < state.selectedSteps[payload.stepIndex].services.length;
        index++
      ) {
        if (index === payload.serviceIndex) {
          state.selectedSteps[payload.stepIndex].services[index].selected =
            payload.selected;
        } else {
          state.selectedSteps[payload.stepIndex].services[
            index
          ].selected = false;
        }
      }
    },
    checkCustomService(state, payload) {
      state[payload.name][payload.serviceIndex].selected = payload.selected;
    },
    toggleDebugger(state) {
      state.data.debugger = !state.data.debugger;
    },
  },
  actions: {},
  modules: {},
});
