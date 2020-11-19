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


This package loads precomputed TxDb objects derived from an internal database of RefSeq transcripts. Stored data generated as follows:
```

```

