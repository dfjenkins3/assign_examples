assign\_examples
================================================================================

The scripts in this repository are examples that show how to get started with
analyzing RNA-Seq data with the RSubread aligner and the ASSIGN package.

### Required R Packages

To run the software, make sure you have installed the following R packages:

```
[Rsubread](http://www.bioconductor.org/packages/release/bioc/html/Rsubread.html)
[limma](http://www.bioconductor.org/packages/release/bioc/html/limma.html)
[edgeR](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html)
[ASSIGN](http://www.bioconductor.org/packages/release/bioc/html/ASSIGN.html)
[sva](http://www.bioconductor.org/packages/release/bioc/html/sva.html)
```

### 1) Create RPKM Values using RSubread

The script ProcessRnaSeqFeatureCounts.R will create a list of RPKM values
based on the annotations in the given GTF file.

First, create an index of your reference genome using RSubread:

```
> library(Rsubread)
> ref <- "/path/to/reference.fa"
> buildindex(basename="reference_index",reference=ref)
```

Then, run ProcessRnaSeqFeatureCounts.R with the following options:

```
Rscript ProcessRnaSeqFeatureCounts.R <ref_index_prefix> <reads_1.fq.gz> \ 
                                     <reads_2.fq.gz> <annotations.gtf> \
                                     <out_prefix> <threads>
```

| Input              | Description |
|--------------------|-------------|
| <ref_index_prefix> | Creating the RSubread reference will create several files with the same prefix. Provide the path to the prefix of these files | 
| <reads_1.fq.gz>    | Input fastq.gz file |
| <reads_2.fq.gz>    | If the reads are paired end, provide the path the second read file, otherwise put 'NULL' |
| <annotations.gtf>  | Path to gene annotations in [GTF](http://www.ensembl.org/info/website/upload/gff.html) format |
| <out_prefix>       | The script will create several output files, provide a prefix where these files will be created |
| <threads>          | Number of threads to use while aligning the data |

### 2) Create a matrix of RPKM values

ASSIGN expects a matrix of input values, so you will need to merge the values 
created by running the ProcessRnaSeqFeatureCounts.R on each of your samples
into one large matrix.  Additionally, ASSIGN will fail if a gene has multiple
RPKM values of 0 across samples, so these will need to be removed. Another
potential source of ASSIGN failure would be an annotation name that is null.
You can merge these values manually or use the provided script: merge_rpkms.pl.
The script will automatically remove any annotations with more than one count
of zero and any annotations with an empty name:

```
./merge_rpkms.pl [file1.rpkmlog] [file2.rpkmlog] ... > merged.rpkmlog
```

### 3) Run ASSIGN

Running ASSIGN depends on your data, so you will need to edit the example script
ASSIGN.R to fit your use case. Comments in ASSIGN.R will help you get started.

##### Questions

Please email me with any questions <dfj@bu.edu> 
