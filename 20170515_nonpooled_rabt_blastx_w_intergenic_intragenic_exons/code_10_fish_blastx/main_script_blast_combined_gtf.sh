#!/bin/bash
#PBS -A simons
#PBS -l walltime=168:00:00
#PBS -l nodes=1:ppn=1
#PBS -N blast
#PBS -o blast.out
#PBS -e blast.err
#PBS -m abe
#PBS -M jmtroy2@igb.illinois.edu

# IMPORTANT: use the below qsub command to run from the code directory
# qsub -v my_script_name=main_script_blast_combined_gtf.sh -S /bin/bash main_script_blast_combined_gtf.sh

#
#  This script...
#
#	Does a grep on the combined gtf file from cufflinks to find intergenic transcripts
#   Note that the combinded gtf file has been updated with exon RPKMs using homer software
#
#		We do a grep on 'class_code "u"' to get the intERgenic exons with no annotation (exons outside of a gene)
#		We do a grep on 'class_code "i"' to get the intRAgenic exons with no annotation (exons inside an intron of a gene)
#
#	use bedops to convert the intergenic .gtf file to an intergenic .bed file
#
#
#	filter out exons with RKPM less than 1
#
#	blastx the remaining exons against fruit fly proteins
#
#

# change directory to torque working directory (the directory "qsub -S /bin/bash main_script.sh" was run from)
cd $PBS_O_WORKDIR

# set a variable with this scripts name
echo begin script "$my_script_name" at `date`

# Get project name, the project name is assumed to be same as the parent folder
# of the current folder.
# the current folder is assumed to be the code folder.
CURRENT_FOLDER=`pwd`
PARENT_FOLDER="$(dirname "$CURRENT_FOLDER")"
PROJECT_NAME="$(basename "$PARENT_FOLDER")"

# TO DO, set input_data_folder if needed
# input data folder not needed for this script 
# INPUT_DATA_FOLDER=

# the PROJECT_FOLDER is in the simoms project foldes
PROJECT_FOLDER="$PARENT_FOLDER"
CODE_FOLDER="$CURRENT_FOLDER"
PROJECT_INPUT_DATA_FOLDER="$PROJECT_FOLDER"/project_input_data

## run id variable - run id variable -  has the timestamp
dt=`date +%Y%m%d`
tm=`date +%H%M%S`
RUN_ID=RUN_"$dt"_"$tm"

## set a variable with the name of the directory the output (and interim data) will be placed, and then create the folders
# TO DO set the variable for the "interim data folder" if needed (remove the # at the beginning of line below)
# INTERIM_DATA_FOLDER="$PROJECT_FOLDER"/"$RUN_ID"/interim_data
# TO DO create the variable for the "interim data folder" if needed (remove the # at the beginning of line below)
# mkdir -p "$INTERIM_DATA_FOLDER" 
OUTPUT_DATA_FOLDER="$PROJECT_FOLDER"/output_10_fish_blastx"_"$RUN_ID
# create a symbolic link to the output folder.
# SYMBL_LINK="$PROJECT_FOLDER"/output_10_fish_blastx_symbolic_link
# rm "$SYMBL_LINK"
# ln -s $OUTPUT_DATA_FOLDER "$SYMBL_LINK"
SAVED_CODE="$OUTPUT_DATA_FOLDER"/saved_code
## (use -p option below) mkdir "$PROJECT_FOLDER"/"$RUN_ID"/output_data
mkdir -p "$OUTPUT_DATA_FOLDER" # the -p option will create any leading directories that do not already exist.
# no interim data folder for this project # mkdir -p "$INTERIM_DATA_FOLDER"

# create a run log file in the output folder
# This file can be used to capture any log information you need in your script.
RUN_LOG_FILE="$OUTPUT_DATA_FOLDER"/"$RUN_ID"_LOG.txt
echo begin script "$my_script_name" at `date` >> "$RUN_LOG_FILE" 
echo "The qsub job name: $PBS_JOBNAME" at `date` >> "$RUN_LOG_FILE"
echo "The qsub job id: $PBS_JOBID" at `date` >> "$RUN_LOG_FILE"
echo "The project folder is: $PROJECT_FOLDER" >> "$RUN_LOG_FILE"
echo "The code folder is: $CODE_FOLDER" >> "$RUN_LOG_FILE"

#################################################################################################
### BEGIN THE REAL WORK NOW THAT THE INITIAL HOUSE KEEPING IS DONE ##############################
#################################################################################################

#
#  HEY, this is a bash function, please also see the code below the function
#
function DOBLAST {

mkdir $SUB_OUPUT_FLDR

# filter the combined gtf file from cuffcompare, using grep, to get all the rows with class_code "u"
EXONS_GTF="$SUB_OUPUT_FLDR"/"$FILE_PREFIX"_transcripts.gtf
grep "$GREP_VAR" $COMB_GTF > $EXONS_GTF

#	use bedops to convert the EXONS .gtf file to an EXONS .bed file
module load bedops/2.4.2
EXONS_BED="$SUB_OUPUT_FLDR"/"$FILE_PREFIX"_transcripts.bed
# the command below does:
# uses the gtf2bed command to convert the gtf file to a bed file
# uses awk to replace the contents of column 10 with just a "."
# uses sort to sort the results
# uses uniq to get only the unique rows
# uses sortBed to sort it as a properly sorted bed file
gtf2bed < $EXONS_GTF | awk -F $'\t'  'BEGIN {OFS=FS} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,".",$11}' | sort | uniq | sortBed > $EXONS_BED

# We have found that homer can calculate a slightly different rpkm value for the same
# exon (that is, chrom, start, end and strand  are the same) for different transcripts
# we consider those as "duplicates" and call the R code below to eliminate the duplicates
# and chose the first of the duplicate rows as the one to keep.  We do this without 
# consideration of the rpkm value because we have seen the values are only slightly different.
EXONS_BED_NO_DUPS="$SUB_OUPUT_FLDR"/"$FILE_PREFIX"_transcripts_no_dups.bed
Rscript keep_one_exon_with_multi_rpkm_values.R input_file="$EXONS_BED" output_file="$EXONS_BED_NO_DUPS"
rm "$EXONS_BED"
mv "$EXONS_BED_NO_DUPS" "$EXONS_BED"
echo "exons for $FILE_PREFIX, unduplicated"
wc -l "$EXONS_BED" >> "$RUN_LOG_FILE" 

#######  Start of Filter by repeats ####################################################
# now eliminate those cuffcompare exons that overlap with repeat masker
REPEAT_BED=/home/groups/simons/Joe/sb_blast_project/input_data/repeats/stickleback_repeatmsker_01_24_2017.bed
REPEAT_BED_GROUP_CHROMS=/home/groups/simons/Joe/sb_blast_project/input_data/repeats/stickleback_repeatmsker_group_chroms.bed
#
#
# Call the script below to change chrUn to scaffolds and 'chr' chromosome names to 'group' chromosome names
# $REPEAT_BED is the input file and $REPEAT_BED_GROUP_CHROMS is the output file, "$OUTPUT_DATA_FOLDER"/repeat_tmp_dir is just a temp directory to do some work.
sh stickleback_convert_chr_to_group_and_unchr_to_scaffolds.sh $REPEAT_BED $REPEAT_BED_GROUP_CHROMS "$OUTPUT_DATA_FOLDER"/repeat_tmp_dir
# now filter the exons for repeats
EXON_BED_REPEAT_FILTER="$SUB_OUPUT_FLDR"/"$FILE_PREFIX"_intergenic_transcripts_filtered_for_repeats.bed
bedtools intersect -v -a "$EXONS_BED" -b $REPEAT_BED_GROUP_CHROMS > $EXON_BED_REPEAT_FILTER
# move data back into EXONS_BED file
rm "$EXONS_BED"
mv "$EXON_BED_REPEAT_FILTER" "$EXONS_BED"
echo "exons for $FILE_PREFIX, unduplicated after repeat filter"
wc -l "$EXONS_BED" >> "$RUN_LOG_FILE" 

####### Now filter by the RPKM values supplied by HOMER ################################
EXONS_BED_RPKM_FILTER="$SUB_OUPUT_FLDR"/"$FILE_PREFIX"_transcripts_filtered_rpkm.bed
awk -F"\t" '($11 > 1.0) {print}' "$EXONS_BED" > $EXONS_BED_RPKM_FILTER
echo "exons for $FILE_PREFIX, unduplicated after repeat filter and rpkm > 1 filter"
wc -l "$EXONS_BED_RPKM_FILTER" >> "$RUN_LOG_FILE" 

# in column 4, add all the data the will be preserved in the fasta file
# change the delimeters from ":" to something more unique like ":+:"
TMP="$SUB_OUPUT_FLDR"/tmp
awk -F $'\t'  'BEGIN {OFS=FS} {$4 = $1":::"$2"-"$3":::"$4":::"$11":::"$6};1' $EXONS_BED_RPKM_FILTER  >$TMP && mv $TMP $EXONS_BED_RPKM_FILTER

# use bedtool getfasta to create a search fasta file of the EXONS regions seleted
GENOME_FA=/home/groups/simons/Abbas/ensemble-genomes/fish/fish.fa
SEARCH_FA="$SUB_OUPUT_FLDR"/blast_serch.fa
bedtools getfasta -name -fi $GENOME_FA -bed $EXONS_BED_RPKM_FILTER -fo $SEARCH_FA

# to get the ensemble protein files using the below commands ...
# cd /home/groups/simons/Joe/projects/20161004-mrsb-sb-D30C-cuffcompare-blast/project_input_data/ensembl_zebrafish_proteins
# wget ftp://ftp.ensembl.org/pub/release-86/fasta/danio_rerio/pep/Danio_rerio.GRCz10.pep.all.fa.gz
# wget ftp://ftp.ensembl.org/pub/release-86/fasta/danio_rerio/pep/Danio_rerio.GRCz10.pep.abinitio.fa.gz

# To get the ncbi protein files use the below commands  (This one looks to have more information in the header)
# cd /home/groups/simons/Joe/projects/20161004-mrsb-sb-D30C-cuffcompare-blast/project_input_data/ncbi_zebrafish_proteins
# wget ftp://ftp.ncbi.nih.gov/genomes/Danio_rerio/protein/protein.fa.gz
# gunzip protein.fa.gz


# make a blast DB of the protein sequences in the protein .faa file
module load blast+/2.3.0
FASTA_IN=/home/groups/simons/Joe/sb_blast_project/input_data/ncbi_zebrafish_proteins/protein.fa
DB_OUT=/home/groups/simons/Joe/sb_blast_project/input_data/ncbi_zebrafish_proteins/protein.db
# this has already been corrected # makeblastdb -in $FASTA_IN -input_type fasta -dbtype prot -out $DB_OUT -taxid 7955 -parse_seqids

# now use blast to get some results
BLASTX_OUT="$OUTPUT_DATA_FOLDER"/blast_results_fmt1.txt
# for now skip this guy (joe troy 9/20/2016) # blastx -db $DB_OUT -query $SEARCH_FA -out $BLASTX_OUT -outfmt 1


# now use blast to get some results
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
BLASTX_OUT="$SUB_OUPUT_FLDR"/blast_results_fmt6_ftp_db.txt
blastx -db $DB_OUT -query $SEARCH_FA -out $BLASTX_OUT -max_target_seqs 1 -max_hsps 1 -outfmt "6 qseqid sseqid evalue bitscore pident length mismatch gapopen qstart qend sstart send salltitles"

# split qseqid (which is the coordinates all as chrom:start-end) into 3 separate fields
# first split on the first occurance of ':' (write it to tmp, and them move tmp back to the original file name)
TMP="$SUB_OUPUT_FLDR"/tmp
awk '{ sub(/:::/, "\t"); print }' $BLASTX_OUT > $TMP && mv $TMP $BLASTX_OUT
# second split on the first occurance of "-" (write it to tmp, and them move tmp back to the original file name)
awk '{ sub(/-/, "\t"); print }' $BLASTX_OUT > $TMP && mv $TMP $BLASTX_OUT
# third splint on next occurance of :
awk '{ sub(/:::/, "\t"); print }' $BLASTX_OUT > $TMP && mv $TMP $BLASTX_OUT
# fourth splint on next occurance of :
awk '{ sub(/:::/, "\t"); print }' $BLASTX_OUT > $TMP && mv $TMP $BLASTX_OUT
# fifth splint on next occurance of :
awk '{ sub(/:::/, "\t"); print }' $BLASTX_OUT > $TMP && mv $TMP $BLASTX_OUT

# now add header

echo "# script name = $my_script_name" > $TMP
echo "# project folder = $PROJECT_FOLDER" >> $TMP
echo "# code folder = $CODE_FOLDER" >> $TMP
echo "# run id = $RUN_ID" >> $TMP
echo "# output folder = $SUB_OUPUT_FLDR" >> $TMP
VALUE=`wc -l $EXONS_BED`
echo "# $FILE_PREFIX Transfrags (exons) per cuffcompare $VALUE" >> TMP
VALUE=`wc -l $EXONS_BED_RPKM_FILTER`
echo "# Transfrags who belong to a transcript with an average across samples of  > 1 rpkm $VALUE" >> TMP
echo "#" >> $TMP
echo -e chrom"\t"start"\t"end"\t"name"\t"rpkm"\t"strand"\t"sseqid"\t"evalue"\t"bitscore"\t"pident"\t"length"\t"mismatch"\t"gapopen"\t"qstart"\t"qend"\t"sstart"\t"send"\t"salltitles >> $TMP

TMP2="$SUB_OUPUT_FLDR"/tmp2

cat $TMP $BLASTX_OUT > $TMP2 && mv $TMP2 $BLASTX_OUT

# THIS IS THE END OF THE DOBLAST function !!!!!!!!!!!!!
}

module load bedtools/2.25.0
module load R/3.2.3
# TODO make sure the below input files are correct
COMB_GTF=/home/groups/simons/Joe/rabt_blast_project_sb/input_data/rabt_nonpooled_ensembl_gtf_code/output/D30C_w_rpkm_counts.gtf
################# D30C_intERgenic ###########################
# set grep var to find intergenic exons
GREP_VAR='class_code "u"'
# set file prefix to name files as intergenic
FILE_PREFIX="D30C_intERgenic"
# not used # TRACKING="$PROJECT_INPUT_DATA_FOLDER"/contained.tracking
SUB_OUPUT_FLDR="$OUTPUT_DATA_FOLDER"/D30C_intERgenic
# call the blast function

DOBLAST

################# D30C_intRAgenic ###########################
# set grep var to find intRAgenic exons
GREP_VAR='class_code "i"'
# set file prefix to name files as intRAgenic
FILE_PREFIX="D30C_intRAgenic"
# not used # TRACKING="$PROJECT_INPUT_DATA_FOLDER"/contained.tracking
SUB_OUPUT_FLDR="$OUTPUT_DATA_FOLDER"/D30C_intRAgenic
# call the blast function
DOBLAST

# TODO make sure the below input files are correct
COMB_GTF=/home/groups/simons/Joe/rabt_blast_project_sb/input_data/rabt_nonpooled_ensembl_gtf_code/output/T30C_w_rpkm_counts.gtf
################# T30C_intERgenic ###########################
# set grep var to find intergenic exons
GREP_VAR='class_code "u"'
# set file prefix to name files as intergenic
FILE_PREFIX="T30C_intERgenic"
# not used # TRACKING="$PROJECT_INPUT_DATA_FOLDER"/contained.tracking
SUB_OUPUT_FLDR="$OUTPUT_DATA_FOLDER"/T30C_intERgenic
# call the blast function
DOBLAST

################# T30C_intRAgenic ###########################
# set grep var to find intRAgenic exons
GREP_VAR='class_code "i"'
# set file prefix to name files as intRAgenic
FILE_PREFIX="T30C_intRAgenic"
# not used # TRACKING="$PROJECT_INPUT_DATA_FOLDER"/contained.tracking
SUB_OUPUT_FLDR="$OUTPUT_DATA_FOLDER"/T30C_intRAgenic
# call the blast function
DOBLAST

#####################################################################################
### END OF THE REAL WORK - DO the FINAL HOUSE KEEPING  ##############################
#####################################################################################

# copy the contents of the current folder (with this script and other code) to the saved code folder
cp -R "$CODE_FOLDER" "$SAVED_CODE"
echo end script "$my_script_name" at `date` >> "$RUN_LOG_FILE" 
