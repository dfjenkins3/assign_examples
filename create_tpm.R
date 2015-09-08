command_args <- commandArgs(trailingOnly=TRUE)

if (length(command_args) != 3){
  print("ERROR: Submit the correct options!")
  print("Rscript create_tpm.R <in.fpkm> <out.tpm> <out.tpmlog>")
  quit(save = "no", status = 1, runLast = FALSE)
}

rpkm <- read.table(command_args[1], sep="\t", header=FALSE)

tpm <- (rpkm$V2 / sum(rpkm$V2)) * 10^6

write.table(cbind(as.character(rpkm$V1), tpm), command_args[2], sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(cbind(as.character(rpkm$V1), log2(tpm+1)), command_args[3], sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
