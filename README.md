# Rbed2HGVS

### Description
appends an additional column of interval in HGVS nomenclature to a bedfile

### Usage
`Rbed2HGVS(bedfile, preferred_tx = NA, ncores = 1)`

### Arguments
```
bedfile	
granges object with seqnames=(1-22,X,Y,MT)

preferred_tx	
path tsv file where column1="gene symbol", column2="refseq transcript". File should be headerless / where there are multiple preferred transcripts from a gene, these should one-row per transcript.

ncores	
number of cores to use for hgvs calculation (default=1)
```

### Value
list object with [1] data.frame of original bed appended with HGVS [2] given preferred transcripts not available in db [3] preferred transcripts with different version to db

Note - exon values given in output refer to the closest CDS exon.


### Database

This package loads precomputed TxDb objects derived from an internal database of RefSeq transcripts. Stored data generated as follows:
```
# download latest RefSeq DB from UCSC
ucsc_hg19_ncbiRefSeq <- GenomicFeatures::makeTxDbFromUCSC(genome = "hg19", tablename = "ncbiRefSeq")

# tidy seqnames
ucsc_hg19_ncbiRefSeq <- AnnotationDbi::loadDb(db) %>%
    GenomeInfoDb::keepStandardChromosomes(.)

GenomeInfoDb::seqlevelsStyle(ucsc_hg19_ncbiRefSeq) <- "NCBI"

# extract UTR, CDS and Tx objects from database
tx_db        <- GenomicFeatures::transcripts(x = ucsc_hg19_ncbiRefSeq)
cds_db       <- GenomicFeatures::cdsBy(x = ucsc_hg19_ncbiRefSeq, by = "tx", use.name = T)
three_utr_db <- GenomicFeatures::threeUTRsByTranscript(x = ucsc_hg19_ncbiRefSeq, use.names=T)
five_utr_db  <- GenomicFeatures::fiveUTRsByTranscript(x = ucsc_hg19_ncbiRefSeq, use.names=T)

# save as internal data
usethis::use_data(tx_db, cds_db, three_utr_db, five_utr_db, internal = TRUE)
```

Current DB Object used:
```
# Db type: TxDb
# Supporting package: GenomicFeatures
# Data source: UCSC
# Genome: hg19
# Organism: Homo sapiens
# Taxonomy ID: 9606
# UCSC Table: ncbiRefSeq
# UCSC Track: NCBI RefSeq
# Resource URL: http://genome.ucsc.edu/
# Type of Gene ID: no gene ids
# Full dataset: yes
# miRBase build ID: NA
# transcript_nrow: 77810
# exon_nrow: 293707
# cds_nrow: 232165
# Db created by: GenomicFeatures package from Bioconductor
# Creation time: 2020-10-04 16:48:41 +0100 (Sun, 04 Oct 2020)
# GenomicFeatures version at creation time: 1.36.4
# RSQLite version at creation time: 2.2.0
# DBSCHEMAVERSION: 1.2
```
