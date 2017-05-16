#!/bin/bash

#  this routine takes a ".bed" file for stickleback fish with
#  the chromosomes named chrI, chrII ... chrXXI and
#  converts them to chromosomes named groupI, groupII ... groupXXI
#  also
#  the routine changes chrUN rows to the the proper scaffold id and the proper
#  start and end coordinates for the scaffold
#
#  The parameters passed to this routine are
#  1) the input file,
#  2) and the output file
#  3) a working folder where any temp files are created
#
#  This routine also uses the file /home/groups/simons/Joe/sb_blast_project/input_data/scaffolds_to_unchr_frm_abbas/ChrUn_to_scaffolds.txt

INPUT_BED=$1
OUTPUT_BED=$2
WRK_DIR=$3
CHRUN_TO_SCAFFOLDS_TXT=/home/groups/simons/Joe/sb_blast_project/input_data/scaffolds_to_unchr_frm_abbas/ChrUn_to_scaffolds.txt

# create the WRK_DIR if it does not exist
if ! [ -d "$WRK_DIR" ]
then
  mkdir -p "$WRK_DIR"
fi

# First step, convert the CHRUN_TO_SCAFFOLDS_TXT to a bed file, CHRUN_TO_SCAFFOLDS_BED
CHRUN_TO_SCAFFOLDS_BED="$WRK_DIR"/ChrUn_to_scaffolds.bed

# convert mapping file to a bed file that will be used to aid in the mapping of the chrUN repeats to scaffold locations
# 1) use tail -n +2 "$CHRUN_TO_SCAFFOLDS_TXT"  to remove the header from the original mapping file
# 2) then pipe the results to the awk command to reformat into a bed file.
# 2a) the    -F $'\t'    tells the awk command the file has a tab delimiter
# 2b) the    '{print $3"\t"$4-1"\t"$5"\t"$1"\t"$2}'   tells the awk command how to order the orginal file columns to create the bed file.
# 2c) note the  "\t"  between each column is to tell awk the resulting bed file is tab separated.
# 2d) also note    $4-1   to convert the start location from "1-based" to "0-based" as required for a bed file.  for more info on 1-based and 0-based see: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
tail -n +2 "$CHRUN_TO_SCAFFOLDS_TXT" | awk -F $'\t' '{print $3"\t"$4-1"\t"$5"\t"$1"\t"$2}'  > "$CHRUN_TO_SCAFFOLDS_BED"

# bedtools intersect aids in converting the chrUn positions to scaffold positions
module load bedtools/2.25.0

# get the chrUn (unknown {un-named?} chromosome) records from the input bed file
CHRUN_INPUT_BED="$WRK_DIR"/ChrUn_input.bed
NOT_CHRUN_INPUT_BED="$WRK_DIR"/not_ChrUn_input.bed
# grep for all rows in INPUT_BED beginning with chrUn
grep ^chrUn "$INPUT_BED" > "$CHRUN_INPUT_BED"
# grep all those that do not begin with chrUn;  use sed to change chr to group and write it to the non chrUn file
grep -v ^chrUn "$INPUT_BED" | sed 's/^chr/group/' > "$NOT_CHRUN_INPUT_BED"

TEMP="$WRK_DIR"/temp.bed
# find which mapping row (CHRUN_TO_SCAFFOLDS_BED) intersects with rows in chrUn bed ($CHRUN_INPUT_BED)
# use the -wo option to keep all the mapping data because we will use it in the next step 
bedtools intersect -wo -a "$CHRUN_INPUT_BED" -b "$CHRUN_TO_SCAFFOLDS_BED" > $TEMP

SCAFFOLD_OUTPUT_BED="$WRK_DIR"/scaffold.bed

# the below awk command re-orders the contents of the TEMP file (which was the results of the bedtools interset)
# 1) the    -F $'\t'    tells the awk command the file has a tab delimiter
# 2) we print $20  first because that is the name of the scaffold (which we want as the first column in our bed file).
# 3) next we print $2-$18  which is the starting location of repeat in terms of chrUn coordinates MINUS the starting position of the scaffold to covert it to scaffold coordinates
# 4) then we print $3-$18  which is the ending location of the repeat in terms of chrUn coordinates MINUX the starting location of the scaffold to convert it to scaffold coordinates
# awk -F $'\t' '{print $10"\t"$2-$8"\t"$3-$8"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' $TEMP > $SCAFFOLD_OUTPUT_BED
awk -F $'\t' '{print $10"\t"$2-$8"\t"$3-$8"\t"$4"\t"$5"\t"$6}' $TEMP > $SCAFFOLD_OUTPUT_BED

cat "$NOT_CHRUN_INPUT_BED" "$SCAFFOLD_OUTPUT_BED" > $OUTPUT_BED

# end of script