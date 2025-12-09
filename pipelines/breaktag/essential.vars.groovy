// breaktag ESSENTIAL VARIABLES

// Define essential variables here.
// Further module-specific variables can be adjusted in the corresponding ".header" files for each module.
//

// General parameters
ESSENTIAL_PROJECT="/home/ljw/new_fold/projects/roukoslab/breaktag/my_test"
ESSENTIAL_SAMPLE_PREFIX="" 
ESSENTIAL_THREADS=16

// Mapping parameters
ESSENTIAL_BWA_REF="/home/ljw/sdc1/SRA_cache/human_dataset/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna"
ESSENTIAL_PAIRED="no"        // paired end design
ESSENTIAL_QUALITY=60          // min mapping quality of reads to be kept. Defaults to 60

// NOTE: you probably don't need to touch anything beyond this point
// further optional pipeline stages to include
RUN_IN_PAIRED_END_MODE=(ESSENTIAL_PAIRED == "yes")

// project folders
PROJECT=ESSENTIAL_PROJECT
LOGS=PROJECT + "/logs"
MAPPED=PROJECT + "/mapped"
QC=PROJECT + "/qc"
RAWDATA="/home/ljw/sdc1/roukos/sra"
REPORTS=PROJECT + "/reports"
RESULTS=PROJECT + "/results"
TMP=PROJECT + "/tmp"
TRACKS=PROJECT + "/tracks"
TARGETS=PROJECT + "/targets.txt"

