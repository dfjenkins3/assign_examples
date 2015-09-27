library(biomaRt)

####
#Replaces pi3k_array rownames that are entrez ids to gene names. The output is in expr1
####

#Grab the hsapiens Mart object 
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

#get the rownames from the pi3k_array
entrezids <- rownames(pi3k_array)

#match the entrezgene ids to gene names
gene_symbols <- getBM(attributes=c('entrezgene', 'external_gene_name'), filters = 'entrezgene', values = entrezids, mart = ensembl)
sym <- gene_symbols[match(entrezids,gene_symbols[,1]),2] 

#find the non duplicate and not 'NA' values to keep in the final list
keep <- !duplicated(sym) & !is.na(sym)

#subset the pik_array to only the 'keep' values
expr1 <- pi3k_array[keep,]

#replace the rownames with the new gene ids
rownames(expr1)=sym[keep]
