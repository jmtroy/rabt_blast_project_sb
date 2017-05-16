#!/bin/bash
#PBS -A simons
#PBS -l walltime=168:00:00
#PBS -l nodes=1:ppn=12
#PBS -N report
#PBS -o report.out
#PBS -e report.err
#PBS -m abe
#PBS -M jmtroy2@igb.illinois.edu

# IMPORTANT: use the below qsub command to run from the projects code folder
# qsub -v my_script_name=main_cufflinks_rabt_cuffcompare.sh -S /bin/bash main_cufflinks_rabt_cuffcompare.sh

# pool mouse bam files together and then run cufflinks

# change directory to torque working directory (the directory the script was run from)
cd $PBS_O_WORKDIR

# echo out script name and time
echo begin script "$my_script_name" at `date`

# Get project name, the project name is assumed to be same as the parent folder
# of the current folder.
# the current folder is assumed to be the code folder.
CURRENT_FOLDER=`pwd`
PARENT_FOLDER="$(dirname "$CURRENT_FOLDER")"
PROJECT_NAME="$(basename "$PARENT_FOLDER")"

# TO DO, set input_data_folder if needed
INPUT_DATA_FOLDER="/home/groups/simons/Joe/mm_blast_project/input_data"

# the PROJECT_FOLDER is in the simoms project foldes
PROJECT_FOLDER="$PARENT_FOLDER"
CODE_FOLDER="$CURRENT_FOLDER"
PROJECT_INPUT_DATA_FOLDER="$PROJECT_FOLDER"/project_input_data

## run id variable - output and intermediate files will go
## in the run id directory (for example /RUN_20130708_162650)
dt=`date +%Y%m%d`
tm=`date +%H%M%S`
RUN_ID=RUN_"$dt"_"$tm"

# TO DO, to help identify the contents of the outfolder, add something to the RUN_ID
RUN_ID="$RUN_ID"

## set a variable with the name of the directory the output (and interim data) will be placed, and then create the folders

# TO DO set the variable for the "interim data folder" if needed (remove the # at the beginning of line below)
# INTERIM_DATA_FOLDER="$PROJECT_FOLDER"/"$RUN_ID"/interim_data
# TO DO create the variable for the "interim data folder" if needed (remove the # at the beginning of line below)
# mkdir -p "$INTERIM_DATA_FOLDER" 

OUTPUT_DATA_FOLDER="$PROJECT_FOLDER"/output_code_nonpooled_cufflinks_rabt_cuffcompare_w_ensembl_gtf"_"$RUN_ID
SAVED_CODE="$OUTPUT_DATA_FOLDER"/saved_code
## (use -p option below) mkdir "$PROJECT_FOLDER"/"$RUN_ID"/output_data
mkdir -p "$OUTPUT_DATA_FOLDER" # the -p option will create any leading directories that do not already exist.
INTERIM_DATA_FOLDER="$OUTPUT_DATA_FOLDER"/interim_data
# not needed for this script # mkdir -p "$INTERIM_DATA_FOLDER"

# create a run log file in the output folder
# This file can be used to capture any log information you need in your script.
RUN_LOG_FILE="$OUTPUT_DATA_FOLDER"/"$RUN_ID"_LOG.txt
echo begin script "$my_script_name" at `date` >> "$RUN_LOG_FILE" 
echo "The qsub job name: $PBS_JOBNAME" at `date` >> "$RUN_LOG_FILE"
echo "The qsub job id: $PBS_JOBID" at `date` >> "$RUN_LOG_FILE"
echo "The project folder is: $PROJECT_FOLDER" >> "$RUN_LOG_FILE"
echo "The code folder is: $CODE_FOLDER" >> "$RUN_LOG_FILE"
echo "The output folder is: $OUTPUT_DATA_FOLDER" >> "$RUN_LOG_FILE"

#################################################################################################
### BEGIN THE REAL WORK NOW THAT THE INITIAL HOUSE KEEPING IS DONE ##############################
#################################################################################################

# define the reference GTF for cuffcompare
gtf="/home/groups/simons/Abbas/ensemble-genomes/fish/fish.gtf" 
# split chrome file into separate ones
CHR_FA="$OUTPUT_DATA_FOLDER"/fasta
mkdir -p "$CHR_FA"
# now make sure the fasta for each chromosome/segment are in the folder by calling split_multifasta.pl to split the fish.fa file
FASTA_IN_FILE="/home/groups/simons/Abbas/ensemble-genomes/fish/fish.fa"
perl split_multifasta.pl  --input_file="$FASTA_IN_FILE"  --output_dir="$CHR_FA" --compress_output=0

# load software modules
module load samtools/0.1.19
module load cufflinks/2.2.1

# 5D_TCCGGAGA-ATAGAGGC_L00M_R1_001.trimmed
# 12D_CGCTCATT-TAATCTTA_L00M_R1_001.trimmed
# 17D_GAGATTCC-GGCTCTGA_L00M_R1_001.trimmed
# 30D_GAATTCGT-GTACTGAC_L00M_R1_001.trimmed
# 33D_CTGAAGCT-GGCTCTGA_L00M_R1_001.trimmed

D30C_5D=/home/groups/simons/Joe/sb_blast_project/sb_D_30_C/aligned/5D_TCCGGAGA-ATAGAGGC_L00M_R1_001.trimmed.tophat/accepted_hits.bam
D30C_12D=/home/groups/simons/Joe/sb_blast_project/sb_D_30_C/aligned/12D_CGCTCATT-TAATCTTA_L00M_R1_001.trimmed.tophat/accepted_hits.bam
D30C_17D=/home/groups/simons/Joe/sb_blast_project/sb_D_30_C/aligned/17D_GAGATTCC-GGCTCTGA_L00M_R1_001.trimmed.tophat/accepted_hits.bam
D30C_30D=/home/groups/simons/Joe/sb_blast_project/sb_D_30_C/aligned/30D_GAATTCGT-GTACTGAC_L00M_R1_001.trimmed.tophat/accepted_hits.bam
D30C_33D=/home/groups/simons/Joe/sb_blast_project/sb_D_30_C/aligned/33D_CTGAAGCT-GGCTCTGA_L00M_R1_001.trimmed.tophat/accepted_hits.bam

CUFFLINKS_OUT="$OUTPUT_DATA_FOLDER"/5D_TCCGGAGA-ATAGAGGC_L00M_R1_001.trimmed.cufflinks_rabt
cufflinks -p 12 -o $CUFFLINKS_OUT --GTF-guide $gtf $D30C_5D
GTF1="$CUFFLINKS_OUT"/transcripts.gtf

CUFFLINKS_OUT="$OUTPUT_DATA_FOLDER"/12D_CGCTCATT-TAATCTTA_L00M_R1_001.trimmed.cufflinks_rabt
cufflinks -p 12 -o $CUFFLINKS_OUT  --GTF-guide $gtf $D30C_12D
GTF2="$CUFFLINKS_OUT"/transcripts.gtf

CUFFLINKS_OUT="$OUTPUT_DATA_FOLDER"/17D_GAGATTCC-GGCTCTGA_L00M_R1_001.trimmed.cufflinks_rabt
cufflinks -p 12 -o $CUFFLINKS_OUT  --GTF-guide $gtf $D30C_17D
GTF3="$CUFFLINKS_OUT"/transcripts.gtf

CUFFLINKS_OUT="$OUTPUT_DATA_FOLDER"/30D_GAATTCGT-GTACTGAC_L00M_R1_001.trimmed.cufflinks_rabt
cufflinks -p 12 -o $CUFFLINKS_OUT  --GTF-guide $gtf $D30C_30D
GTF4="$CUFFLINKS_OUT"/transcripts.gtf

CUFFLINKS_OUT="$OUTPUT_DATA_FOLDER"/33D_CTGAAGCT-GGCTCTGA_L00M_R1_001.trimmed.cufflinks_rabt
cufflinks -p 12 -o $CUFFLINKS_OUT  --GTF-guide $gtf $D30C_33D
GTF5="$CUFFLINKS_OUT"/transcripts.gtf

# do cuffcompare (see notes below)
CUFFCOMPARE_OUT="$OUTPUT_DATA_FOLDER"/D30C_cuffcompare_contained
cuffcompare -V -o "$CUFFCOMPARE_OUT" -C -s "$CHR_FA" -r "$gtf" -R "$GTF1" "$GTF2" "$GTF3" "$GTF4" "$GTF5" 2> "$OUTPUT_DATA_FOLDER"/D30c_contained.log.txt

#
# T
# 

T30C_5T=/home/groups/simons/Joe/sb_blast_project/sb_T_30_C/aligned/5T_TCCGGAGA-TATAGCCT_L00M_R1_001.trimmed.tophat/accepted_hits.bam
T30C_12T=/home/groups/simons/Joe/sb_blast_project/sb_T_30_C/aligned/12T_CGCTCATT-AGGCGAAG_L00M_R1_001.trimmed.tophat/accepted_hits.bam
T30C_17T=/home/groups/simons/Joe/sb_blast_project/sb_T_30_C/aligned/17T_GAGATTCC-CCTATCCT_L00M_R1_001.trimmed.tophat/accepted_hits.bam
T30C_30T=/home/groups/simons/Joe/sb_blast_project/sb_T_30_C/aligned/30T_GAATTCGT-CAGGACGT_L00M_R1_001.trimmed.tophat/accepted_hits.bam
T30C_33T=/home/groups/simons/Joe/sb_blast_project/sb_T_30_C/aligned/33T_CTGAAGCT-CCTATCCT_L00M_R1_001.trimmed.tophat/accepted_hits.bam

CUFFLINKS_OUT="$OUTPUT_DATA_FOLDER"/5T_TCCGGAGA-TATAGCCT_L00M_R1_001.trimmed.cufflinks_rabt
cufflinks -p 12 -o $CUFFLINKS_OUT  --GTF-guide $gtf $T30C_5T
GTF1="$CUFFLINKS_OUT"/transcripts.gtf

CUFFLINKS_OUT="$OUTPUT_DATA_FOLDER"/12T_CGCTCATT-AGGCGAAG_L00M_R1_001.trimmed.cufflinks_rabt
cufflinks -p 12 -o $CUFFLINKS_OUT  --GTF-guide $gtf $T30C_12T
GTF2="$CUFFLINKS_OUT"/transcripts.gtf

CUFFLINKS_OUT="$OUTPUT_DATA_FOLDER"/17T_GAGATTCC-CCTATCCT_L00M_R1_001.trimmed.cufflinks_rabt
cufflinks -p 12 -o $CUFFLINKS_OUT  --GTF-guide $gtf $T30C_17T
GTF3="$CUFFLINKS_OUT"/transcripts.gtf

CUFFLINKS_OUT="$OUTPUT_DATA_FOLDER"/30T_GAATTCGT-CAGGACGT_L00M_R1_001.trimmed.cufflinks_rabt
cufflinks -p 12 -o $CUFFLINKS_OUT  --GTF-guide $gtf $T30C_30T
GTF4="$CUFFLINKS_OUT"/transcripts.gtf

CUFFLINKS_OUT="$OUTPUT_DATA_FOLDER"/33T_CTGAAGCT-CCTATCCT_L00M_R1_001.trimmed.cufflinks_rabt
cufflinks -p 12 -o $CUFFLINKS_OUT  --GTF-guide $gtf $T30C_33T
GTF5="$CUFFLINKS_OUT"/transcripts.gtf

# do cuffcompare (see notes below)
CUFFCOMPARE_OUT="$OUTPUT_DATA_FOLDER"/T30C_cuffcompare_contained
cuffcompare -V -o "$CUFFCOMPARE_OUT" -C -s "$CHR_FA" -r "$gtf" -R "$GTF1" "$GTF2" "$GTF3" "$GTF4" "$GTF5"  2> "$OUTPUT_DATA_FOLDER"/T30c_contained.log.txt


# NOTES on cuffcompare
# from cuffcompare document at http://cole-trapnell-lab.github.io/cufflinks/cuffcompare/index.html
#
# -o <outprefix>
# 
# All output files created by Cuffcompare will have this prefix
# (e.g. .loci, .tracking, etc.). If this option is not provided
# the default output prefix being used is: "cuffcmp"
# 
# -r
# 
# An optional “reference” annotation GFF file. Each sample is matched
# against this file, and sample isoforms are tagged as overlapping, matching, or novel
# where appropriate. See the refmap and tmap output file descriptions below.
# 
# -R
# 
# If -r was specified, this option causes cuffcompare to ignore reference
# transcripts that are not overlapped by any transcript in one of cuff1.gtf,…,cuffN.gtf.
# Useful for ignoring annotated transcripts that are not present in your RNA-Seq samples
# and thus adjusting the “sensitivity” calculation in the accuracy report
# written in the file
# 
# -s <seq_dir>
# 
# Causes cuffcompare to look into for fasta files with the underlying genomic
# sequences (one file per contig) against which your reads were aligned for
# some optional classification functions. For example, Cufflinks transcripts consisting
# mostly of lower-case bases are classified as repeats. Note that must contain
# one fasta file per reference chromosome, and each file must be named after
# the chromosome, and have a .fa or .fasta extension.
# 
# -C
# 
# Enables the “contained” transcripts to be also written in the .combined.gtffile,
# with the attribute "contained_in" showing the first container transfrag found. By
# default, without this option, cuffcompare does not write in that file isoforms that
# were found to be fully contained/covered (with the same compatible intron structure)
# by other transfrags in the same locus.
#
# -V
# 
# Cuffcompare is a little more verbose about what it’s doing, printing messages 
# to stderr, and it will also show warning messages about any inconsistencies 
# or potential issues found while reading the given GFF file(s).


#####################################################################################
### END OF THE REAL WORK - DO the FINAL HOUSE KEEPING  ##############################
#####################################################################################

# copy the contents of the current folder (with this script and other code) to the saved code folder
cp -R "$CODE_FOLDER" "$SAVED_CODE"
echo end script "$my_script_name" at `date` >> "$RUN_LOG_FILE" 
