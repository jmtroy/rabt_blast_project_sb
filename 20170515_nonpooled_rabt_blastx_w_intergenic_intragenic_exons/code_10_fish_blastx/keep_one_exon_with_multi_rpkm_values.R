# keep_one_exon_with_multi_rpkm_values.R
# Joe Troy 01/17/2017

# remove rows that have duplicated exons

# STEP 1 read in paramaters
# process command lines
args <- commandArgs(trailingOnly = TRUE)
# set argument defaults
input_file = ''
output_file = ''
# rpkm_column is column 11


# loop through arguments and assign values
arg_count <- length(args)
for ( xxx in seq(along=args)) {
	# print(paste("xxx = ",xxx))
	re = regexpr("=",args[xxx],fixed=TRUE)
	argname = substr(args[xxx],1,re[1]-1)
	# print(paste("argname =",argname,sep=""))
	argvalue = substr(args[xxx],re[1]+1,nchar(args[xxx]))
	# print(paste("argvalue =",argvalue,sep=""))	
	
	if (argname == 'input_file') {
		input_file = argvalue
	} else if (argname == 'output_file') {
		output_file = argvalue
	} 
}
# print("end reading arguments")

D_IN <- read.delim(input_file, header = FALSE, sep = "\t", quote = "\"", stringsAsFactors=FALSE)

D_IN_COPY = D_IN

# rpkms are in col 11
summary(D_IN_COPY[,11])

# set the rpkm column to 0 so any columns that are different only because
# of the rpkm value will be identified as duplicates
D_IN_COPY[,11] = 0 
summary(D_IN_COPY[,11]) 

# D_OUT is everything in D_IN, except those identified as duplicates
D_OUT = D_IN[!duplicated(D_IN_COPY),]

write.table(D_OUT,file = output_file ,sep="\t", row.names=FALSE,col.names = FALSE,quote = FALSE)

