library(Rsubread)
library(limma)
library(edgeR)
library(tools)
options(digits=2)

command_args <- commandArgs(trailingOnly=TRUE)

if (length(command_args) != 6){
    print("ERROR: Submit the correct options!")
    print("Rscript <reference_path> <in_1.fq.gz> <in_2.fq.gz> <features.gtf> <out_prefix> <threads>")
    quit(save = "no", status = 1, runLast = FALSE)
}

referenceGenomeFastaFilePath = command_args[1]
inFilePath1 = command_args[2]
inFilePath2 = command_args[3] # NULL for single-end analyses or when a BAM file has been specified
gtfFilePath = command_args[4]
outFilePrefix = command_args[5]
nthreads = command_args[6]

input_format = "gzFASTQ"
if (file_ext(inFilePath1) == "bam")
  input_format = "BAM"

#Prep output files
outBamFilePath = paste(outFilePrefix, "bam", sep=".")
outCountsFilePath = paste(outFilePrefix, "featureCounts", sep=".")
outRpkmFilePath = paste(outFilePrefix, "rpkm", sep=".")
outRpkmLogFilePath = paste(outFilePrefix, "rpkmlog", sep=".")
outStatsFilePath = paste(outFilePrefix, "stats", sep=".")

#2nd file null if only fragment data
if (inFilePath2 == "NULL")
  inFilePath2 = NULL

#Create bam file if it does not exist
if (!file.exists(outBamFilePath))
  align(index=referenceGenomeFastaFilePath, readfile1=inFilePath1, readfile2=inFilePath2, output_file=outBamFilePath, nthreads=nthreads, input_format=input_format, tieBreakHamming=TRUE, unique=TRUE, indels=5)

#Create rpkm list
fCountsList = featureCounts(outBamFilePath, annot.ext=gtfFilePath, isGTFAnnotationFile=TRUE, nthreads=nthreads, isPairedEnd=!is.null(inFilePath2))
dgeList = DGEList(counts=fCountsList$counts, genes=fCountsList$annotation)
rpkm = rpkm(dgeList, dgeList$genes$Length)

#Write out stats
write.table(fCountsList$stat, outStatsFilePath, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

#Write out featureCounts file
featureCounts = cbind(fCountsList$annotation[,1], fCountsList$counts)
write.table(featureCounts, outCountsFilePath, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

#Write out rpkm file and rpkmlog file
write.table(cbind(fCountsList$annotation[,1], rpkm), outRpkmFilePath, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(cbind(fCountsList$annotation[,1], log2(rpkm + 1)), outRpkmLogFilePath, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

#Delete bam file and indel file
unlink(outBamFilePath)
unlink(paste(outBamFilePath, ".indel", sep=""))
