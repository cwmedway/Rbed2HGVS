
#' appends an additional column of interval in HGVS nomenclature to a bedfile
#'
#' @param bedfile path to bedfile
#' @param db ucsc TxDb object
#' @param preferred_tx path tsv file where column1="gene symbol", column2="refseq transcript". File should be headerless /
#' where there are multiple preferred transcripts from a gene, these should one-row per transcript.
#'
#' @return list object with [1] data.frame of HGVS [2] preferred transcripts no available in db [3] preferred transcripts with different version to db
#' @export
#' @importFrom magrittr %>%
#' @import org.Hs.eg.db
Rbed2HGVS <- function(bedfile, db, preferred_tx = NA, ncores = NA) {

  # set number of cores if not given
  if (is.na(ncores)) {
    ncores <- parallel::detectCores() - 1
  }

  # load bedfile to GRanges and use NCBI chr convention
  bedfile <- rtracklayer::import.bed(bedfile)

  if (length(bedfile) > 0) {

    GenomeInfoDb::seqlevelsStyle(bedfile) <- "NCBI"

    # load refseq database
    ucsc_hg19_ncbiRefSeq <- AnnotationDbi::loadDb(db) %>%
      GenomeInfoDb::keepStandardChromosomes(.)

    # convert chr to NCBI style
    GenomeInfoDb::seqlevelsStyle(ucsc_hg19_ncbiRefSeq) <- "NCBI"

    # return list[3]. [1] = list of overlapping refseq transcripts indexed by bed entry
    # [2] = preferred transcripts missing in db
    # [3] = preferred transcripts different version to db
    # [2] & [3] both NA is preferred transcripts not given
    tx <- getTranscripts(preferred_tx = preferred_tx, db = ucsc_hg19_ncbiRefSeq, bedfile = bedfile)

    #
    tx_names <- unique(unlist(tx$model))

    # GRangesList indexed by RefSeq
    cds_by_tx <- GenomicFeatures::cdsBy(x = ucsc_hg19_ncbiRefSeq, by = "tx", use.name = T)[tx_names]
    three_utr <- GenomicFeatures::threeUTRsByTranscript(x = ucsc_hg19_ncbiRefSeq, use.names=T)[tx_names]
    five_utr  <- GenomicFeatures::fiveUTRsByTranscript(x = ucsc_hg19_ncbiRefSeq, use.names=T)[tx_names]

    # get cds info for start and end of bedfile
    hgvs <- getHgvs(
      bedfile = bedfile,
      cds_by_tx = cds_by_tx,
      cds_ol = tx[['model']],
      ncores = ncores,
      five_utr = five_utr,
      three_utr = three_utr
      ) %>% do.call(rbind, .)

    #append HGMD
    hgvs$gene <- getSymbolRefseq(refSeqId = hgvs$tx)

    list(
      "hgvs" = hgvs,
      "missing" = tx[['missing']],
      "version" = tx[['version']]
    ) %>% return()

  } else {
    # bedfile is empty
    warning("bedfile is empty")
    return(NA)
  }
}



getTranscripts <- function(preferred_tx, db, bedfile, flank_length=1000) {

  # GenomicRanges object of refseq transcripts
  refseq_tx <- GenomicFeatures::transcripts(x = db)

  # keep only 'NM_' transcripts
  refseq_tx <- stringr::str_detect(string = refseq_tx$tx_name, pattern = '^NM_') %>%
    refseq_tx[.]

  if (is.character(preferred_tx)) {
    # read preferred transcripts (ptx) if given
    df_ptx <- read.table(file = preferred_tx, header = F, stringsAsFactors = F)
    # make sure there are no duplicated rows
    df_ptx <- df_ptx[!duplicated(df_ptx),]
    # reformat refseq name without version suffix
    ptx_prefix <- stringr::str_remove(string = df_ptx[,2], pattern = '\\.\\d+')
    # get names of tx in RefSeq DB and get prefix
    cds_name   <- refseq_tx$tx_name
    cds_prefix <- stringr::str_remove(string = cds_name, pattern = "\\.\\d+")

    # are preferred tx in RedSeq DB
    refseq_ptx_db <- match(x = ptx_prefix, table = cds_prefix)

    # which preferred tx not in RefSeq DB - give warning
    ptx_missing <- df_ptx[is.na(refseq_ptx_db),]
    if (dim(ptx_missing)[1] > 0) {
      warning("Some preferred transcripts are missing from the RefSeq DB")
    }

    # get df containing only matched ptx
    ptx_not_missing <- df_ptx[!is.na(refseq_ptx_db),]

    # extract full tx names (inc. suffix) from RefSeq DB
    cds_not_missing <- cds_name[refseq_ptx_db[!is.na(refseq_ptx_db)]]

    # check if versions match - if not warn and make df of diferences
    not_same_version <- !(ptx_not_missing[,2] %in% cds_not_missing)
    if (sum(not_same_version) > 0) {
      warning("Some preferred transcripts are a different version from the RefSeq DB")
      ptx_version <- data.frame(
        "gene" = ptx_not_missing[not_same_version,1],
        "preferred_tx" = ptx_not_missing[not_same_version,2],
        "db_tx" = cds_not_missing[not_same_version] )
    } else { ptx_version <- NA }

    # overwite refseq object to only contain preferred transcripts
    refseq_tx <- refseq_tx[refseq_ptx_db[!is.na(refseq_ptx_db)]]

    # index of bed overlap with transcripts
    bed_ol_tx <- IRanges::findOverlaps(query = bedfile, subject = refseq_tx, maxgap = flank_length)

    # format output object as list containing transcripts indexed by bed line
    cds_ol <- lapply(
      seq(bedfile), function(x) {
        S4Vectors::queryHits(bed_ol_tx) %in% x %>%
          bed_ol_tx[.] %>%
          S4Vectors::subjectHits(.) %>%
          refseq_tx$tx_name[.]
      }
    )

    return(
      return(list("model" = cds_ol, "missing" = ptx_missing, "version" = ptx_version))
    )
  } else {
    # no preferred transcripts given
    # index of bed overlap with transcripts
    bed_ol_tx <- IRanges::findOverlaps(query = bedfile, subject = refseq_tx, maxgap = flank_length)

    # get transcript that overlap each bed entry
    # output should contain hgvs for these
    cds_ol <- lapply(
      seq(bedfile), function(x) {
        S4Vectors::queryHits(bed_ol_tx) %in% x %>%
          bed_ol_tx[.] %>%
          S4Vectors::subjectHits(.) %>%
          refseq_tx$tx_name[.]
      }
    )

    return(list("model" = cds_ol, "missing" = NA, "version" = NA))

  }
  }


getHgvs <- function(bedfile, cds_by_tx, cds_ol, ncores, five_utr, three_utr) {

  # loop over each element (bed entry)
  parallel::mclapply(seq(cds_ol), mc.cores = ncores, function(bedln) {

    # will be first three columns of output
    chr   <- GenomicRanges::seqnames(bedfile[bedln]) %>% as.vector()
    start <- GenomicRanges::start(bedfile[bedln]) - 1 # -1 because GRanges representation of BED different to UCSC (1 vs 0 based)
    end   <- GenomicRanges::end(bedfile[bedln])

    message(paste0(chr, ":", start, "-", end))

    # if bed entry overlaps at least one transcript
    if(!isEmpty(cds_ol)[bedln]) {

      # loop over each tx - get exon & hgvs
      cds_annot <- lapply(cds_ol[[bedln]], function(tx_i) {
        message(tx_i)
        mapCoordToCds(
          bedfile = bedfile[bedln],
          cds = cds_by_tx[[tx_i]],
          tx_i = tx_i,
          five_utr = five_utr[[tx_i]],
          three_utr = three_utr[[tx_i]]
          )
      })

      # parse fields to make df
      tx         <- cds_ol[[bedln]]
      hgvs_start <- lapply(cds_annot, function(x) {x$start$hgvs}) %>% unlist()
      hgvs_end   <- lapply(cds_annot, function(x) {x$end$hgvs}) %>% unlist()
      exon_start <- lapply(cds_annot, function(x) {x$start$exon_rank}) %>% unlist()
      exon_end   <- lapply(cds_annot, function(x) {x$end$exon_rank}) %>% unlist()
    } else {
      # null when no overlapping tx
      tx <- NA
      hgvs_start <- NA
      hgvs_end   <- NA
      exon_start <- NA
      exon_end   <- NA
    }

    data.frame(chr, start, end, tx, hgvs_start, exon_start, hgvs_end, exon_end) %>%
      return()
  })
}

mapCoordToCds <- function(bedfile, cds, tx_i, five_utr, three_utr) {

  # only perform if transcript has cds
  if (!is.null(cds)) {

    # make separate ranges object for start and end coordinates
    bedfile_start <- IRanges::narrow(x = bedfile, start = 1, width = 1)
    bedfile_end   <- IRanges::narrow(x = bedfile, start = -1, width = 1)

    hgvs_start <- getHgvs2(bedfile = bedfile_start, cds = cds, five_utr = five_utr, three_utr = three_utr)
    hgvs_end   <- getHgvs2(bedfile = bedfile_end, cds = cds, five_utr = five_utr, three_utr = three_utr)

    # check orientation of this transcript (samples from first cds)
    strand <- GenomicRanges::strand(cds)[1] %>% as.vector()

    # flip start and end HGVS for -ve strand
    if (strand == "+") {
      list("start" = hgvs_start, "end" = hgvs_end) %>% return()
    } else {
      list("start" = hgvs_end, "end" = hgvs_start) %>% return()
    }
  } else {
      # if transcript does not have cds (i.e. "NR_")
    list("start" = list("hgvs" = NA, "exon_rank" = NA), "end" = list("hgvs" = NA, "exon_rank" = NA)) %>% return()
    }
}



getHgvs2 <- function(bedfile, cds, five_utr, three_utr) {

  # distance to closest CDS (-/+)
  dist_dir <- dist_dir_to_nearest(bedfile = bedfile, ranges = cds)

  if (dist_dir$dist == 0) {

    # ----- WITHIN CDS -----
    if (dist_dir$strand == '+') {
      hgvs <- pos_within_cds_positive(dist_dir = dist_dir, cds = cds)
    } else if (dist_dir$strand == '-') {
      hgvs <- pos_within_cds_negative(dist_dir = dist_dir, cds = cds)
    }
  } else {

    # ----- OUTSIDE CDS -----
    if (dist_dir$strand == '+') {
      hgvs <- pos_outside_cds_positive(bedfile = bedfile, dist_dir = dist_dir, five_utr = five_utr, three_utr = three_utr)
    } else if (dist_dir$strand == "-") {
      hgvs <- pos_outside_cds_negative(bedfile = bedfile, dist_dir = dist_dir, five_utr = five_utr, three_utr = three_utr)
    }
  }

  return(list("hgvs" = hgvs, "exon_rank" = cds[dist_dir$index]$exon_rank))
}




pos_within_cds_positive <- function(dist_dir, cds) {
  # distance into cds hit
  cds_width_hit <- (dist_dir$event_bp + 1) - dist_dir$nearest_start
  # accum width of cds upstream of exon for + strand
  hgvs <- paste0(sum(cds_width_hit + dist_dir$upstream_size))

  return(hgvs)
}

pos_within_cds_negative <- function(dist_dir, cds) {
  # distance into cds hit
  cds_width_hit <- (dist_dir$nearest_end - (dist_dir$event_bp - 1))
  # accum width of cds upstream of exon for - strand
  hgvs <- paste0(sum(cds_width_hit + dist_dir$upstream_size))
  return(hgvs)
}

pos_outside_cds_positive <- function(bedfile, dist_dir, five_utr, three_utr) {

  if ( dist_dir$dir == '-' ) {
    # interval is upstream of exon. hgvs should not include nearest exon
    # +1 because in relation to first base of upstream exon
    entry_cds <- 1 + dist_dir$upstream_size

    if (entry_cds == 1) {
      # upstream of CDS#1 - could be 5'UTR
      hgvs <- pos_five_utr(bedfile = bedfile, five_utr = five_utr)
    } else {
      # not upstream of first CDS - UTR not a consideration
      hgvs <- paste0(entry_cds, dist_dir$dir, dist_dir$dist + 1)
    }
  } else if ( dist_dir$dir == '+' ) {
    # downstream, hgvs includes nearest exon
    entry_cds <- dist_dir$upstream_size + dist_dir$nearest_size
    if (entry_cds == dist_dir$total_size) {
      # downstream of last CDS = could be 3'UTR
      hgvs <- pos_three_utr(bedfile = bedfile, three_utr = three_utr)
    } else {
      # not dpwnstream of last CDS - UTR not consideration
      hgvs <- paste0(entry_cds, dist_dir$dir, dist_dir$dist + 1)
    }
  }

  return(hgvs)
}

pos_outside_cds_negative <- function(bedfile, dist_dir, five_utr, three_utr) {

  # upstream of CDS on negative strand
  if (dist_dir$dir == '+') {

    # amount of CDS upstream
    entry_cds <- dist_dir$upstream_size + 1

    # if no CDS upstream - 5'UTR?
    if ( entry_cds == 1 ) {
      hgvs <- neg_five_utr(bedfile = bedfile, five_utr = five_utr)
    } else {
      hgvs <- paste0(entry_cds, "-", (dist_dir$dist + 1))
    }

  } else if ( dist_dir$dir == '-' ) {
    # interval is downstream of CDS on negative strand
    # hgvs shouls not includes exon. +1 because in relation to first base of nearest exon
    entry_cds <- dist_dir$upstream_size + dist_dir$nearest_size #?
    if (entry_cds == dist_dir$total_size) {
      # upstream of CDS1 - 5' UTR?
      hgvs <- neg_three_utr(bedfile = bedfile, three_utr = three_utr)
    } else {
      hgvs <- paste0(entry_cds, "+", dist_dir$dist + 1)
    }
  }
  return(hgvs)
}


pos_five_utr <- function(bedfile, five_utr) {

  dist_dir_utr <- dist_dir_to_nearest(bedfile = bedfile, ranges = five_utr)

  if ( dist_dir_utr$dir == 'in' ) {
    # within UTR exon
    hgvs <- paste0(
      "-",
      abs(dist_dir_utr$event_bp - dist_dir_utr$nearest_end) + dist_dir_utr$downstream_size + 1
    )
  } else if ( dist_dir_utr$index == 1 && dist_dir_utr$dir == '-' ) {
    # upstream of utr
    # outside UTR
    hgvs <- paste0(
      "-",
      sum(
        dist_dir_utr$downstream_size,
        dist_dir_utr$nearest_size,
        dist_dir_utr$dist,
        1))
  } else {
    # between UTR exons
    if (dist_dir_utr$dir == '-') {
      hgvs <- paste0(
        "-",
        dist_dir_utr$downstream_size + dist_dir_utr$nearest_size,
        dist_dir_utr$dir,
        dist_dir_utr$dist + 1 )
      } else if (dist_dir_utr$dir == '+') {
        hgvs <- paste0(
          "-",
          dist_dir_utr$downstream_size,
          dist_dir_utr$dir,
          dist_dir_utr$dist + 1 )
      }
    }
  return(hgvs)
}


pos_three_utr <- function(bedfile, three_utr) {

  dist_dir_utr <- dist_dir_to_nearest(bedfile = bedfile, ranges = three_utr)

  if ( dist_dir_utr$dir == 'in' ) {
    # within UTR exon
    hgvs <- paste0(
      "*",
      abs(dist_dir_utr$event_bp - dist_dir_utr$nearest_start) + dist_dir_utr$downstream_size + 1
    )
  } else if ( dist_dir_utr$event_bp > IRanges::end(three_utr[length(three_utr)]))  {
    # event_bp is downstream of final UTR exon
    hgvs <- paste0(
      "*",
      sum(
        dist_dir_utr$downstream_size,
        dist_dir_utr$nearest_size,
        dist_dir_utr$dist,
        1))
  } else {
    # between UTR exons
    if (dist_dir_utr$dir == '-') {
      hgvs <- paste0(
        "*",
        dist_dir_utr$upstream_size + 1,
        dist_dir_utr$dir,
        dist_dir_utr$dist + 1)
    } else if (dist_dir_utr$dir == '+') {
      hgvs <- paste0(
        "*",
        dist_dir_utr$upstream_size + dist_dir_utr$nearest_size,
        dist_dir_utr$dir,
        dist_dir_utr$dist + 1)
    }
  }
  return(hgvs)
}



neg_five_utr <- function(bedfile, five_utr) {

  dist_dir_utr <- dist_dir_to_nearest(bedfile = bedfile, ranges = five_utr)

  if ( dist_dir_utr$dir == 'in' ) {
    # within UTR exon
    hgvs <- paste0(
      "-",
      abs(dist_dir_utr$event_bp - dist_dir_utr$nearest_start) + dist_dir_utr$downstream_size + 1
    )
  } else if (dist_dir_utr$event_bp > IRanges::end(five_utr[1])) {
    # upstream of first UTR exon
    hgvs <- paste0(
      "-",
      dist_dir_utr$downstream_size +
        dist_dir_utr$nearest_size +
        dist_dir_utr$dist )
    } else {
    # between UTR exons
    if (dist_dir_utr$dir == '-') {
      # downstream on -ve strand
      hgvs <- paste0(
        "-",
        dist_dir_utr$downstream_size + 1,
        '+',
        dist_dir_utr$dist + 1 )
    } else if (dist_dir_utr$dir == '+') {
      # upstream on -ve strand
      hgvs <- paste0(
        "-",
        dist_dir_utr$downstream_size + dist_dir_utr$nearest_size,
        "-",
        dist_dir_utr$dist + 1)
    }


  }
return(hgvs)
}

neg_three_utr <- function(bedfile, three_utr) {

  dist_dir_utr <- dist_dir_to_nearest(bedfile = bedfile, ranges = three_utr)

  if ( dist_dir_utr$dir == 'in' ) {
    # within UTR exon
    hgvs <- paste0(
      "*",
      abs(dist_dir_utr$event_bp - dist_dir_utr$nearest_end) + dist_dir_utr$downstream_size + 1
    )
  } else if ( dist_dir_utr$event_bp < IRanges::start(three_utr[length(three_utr)]))  {
    # event_bp is downstream of final UTR exon on -ve strand
    hgvs <- paste0(
      "*",
      sum(
        dist_dir_utr$total_size,
        dist_dir_utr$dist,
        1))
  } else {
    # between UTR exons
    if (dist_dir_utr$dir == '-') {
      # downstream on -ve strand
      hgvs <- paste0(
        "*",
        dist_dir_utr$upstream_size +
          dist_dir_utr$nearest_size,
        "+",
        dist_dir_utr$dist + 1)
    } else if (dist_dir_utr$dir == '+') {
      # upstream on -ve strand
      hgvs <- paste0(
        "*",
        dist_dir_utr$upstream_size + 1,
        '-',
        dist_dir_utr$dist + 1)
    }
  }
  return(hgvs)

}


dist_dir_to_nearest <- function(bedfile, ranges) {

  # total number of events (i.e. exons) in range object
  total_i <- length(ranges)
  total_size <- IRanges::width(ranges) %>% sum()

  strand <- GenomicRanges::strand(ranges) %>% as.character() %>% .[1]
  pos  <- IRanges::start(bedfile)
  dist <- IRanges::distance(x = bedfile, y = ranges)
  nearest_i    <- which.min(dist)
  nearest_dist <- dist[nearest_i]
  nearest_size <- IRanges::width(ranges[nearest_i])

  nearest_start <- IRanges::start(ranges[nearest_i])
  nearest_end <- IRanges::end(ranges[nearest_i])

  # pos within events
  if (nearest_dist == 0) {

    dir <- "in"  # in = pos within event

    # check if pos within first event
    if (nearest_i == 1) {
      # if within first event in range, there are no upstream
      upstream_size <- 0
    } else {
      upstream_size <- IRanges::width(ranges[1:(nearest_i-1)]) %>% sum()
    }

    # check if pos within last event
    if (nearest_i == total_i) {
      downstream_size <- 0
    } else {
      downstream_size <- IRanges::width(ranges[(nearest_i+1):total_i]) %>% sum()
    }

    # pos not within event
  } else if (nearest_dist != 0) {

    if (pos < nearest_start) {
      # upstream on positive strand
      dir <- "-"

      if (nearest_i == 1) {
        upstream_size <- 0
      } else {
        upstream_size   <- IRanges::width(ranges[1:(nearest_i - 1)]) %>% sum()
      }

      if (nearest_i == total_i) {
        downstream_size <- 0
      } else {
        downstream_size <- IRanges::width(ranges[(nearest_i+1):total_i]) %>% sum()
      }

    } else {
      # downstream of positive strand
      dir <- "+"

      if (nearest_i == 1) {
        upstream_size <- 0
      } else {
        upstream_size <- IRanges::width(ranges[1:(nearest_i - 1)]) %>% sum()
      }

      if (nearest_i == total_i) {
        downstream_size <- 0
      } else {
        # + 1 because size calculated from last cds of nearest exon
        downstream_size <- 1 + IRanges::width(ranges[(nearest_i+1):total_i]) %>% sum()
      }
    }
  }
  return(
    list(
      "event_bp"=pos,
      "strand"=strand,
      "index"=nearest_i,
      "total_i"=total_i,
      "total_size"=total_size,
      "upstream_size" = upstream_size,
      "downstream_size" = downstream_size,
      "nearest_start"=nearest_start,
      "nearest_end"=nearest_end,
      "nearest_size"=nearest_size,
      "dist"=nearest_dist,
      "dir"=dir)
  )
}




getSymbolRefseq <- function(refSeqId) {

  # remove suffix
  stringr::str_remove(string = refSeqId, pattern = "\\.\\d+") %>%
    AnnotationDbi::select(org.Hs.eg.db, keytype = "REFSEQ", keys = ., columns = "SYMBOL") %>%
    .[,"SYMBOL"] %>%
    return()
}



# -----------------------------
# TESTS
# BRCA2 (Positive Strand)



test <- function() {

test <- function(start, db, preferred_tx = 'BRCA2.preferredtx') {

  db <- 'inst/extdata/ucsc_hg19_ncbiRefSeq.sqlite'

  ucsc_hg19_ncbiRefSeq <- AnnotationDbi::loadDb(db) %>%
    GenomeInfoDb::keepStandardChromosomes(.)

  GenomeInfoDb::seqlevelsStyle(ucsc_hg19_ncbiRefSeq) <- "NCBI"

  gene <- c("BRCA2")
  tx   <- c("NM_000059.3")

  data.frame(gene, tx) %>%
    write.table(
      x = .,
      file = "BRCA2.preferredtx",
      quote = F,
      row.names = F,
      col.names = F,
      sep = "\t"
    )

  bedfile <- GenomicRanges::GRanges(seqnames = "13", ranges = IRanges::IRanges(start=start, width = 1), strand = '+')
  tx <- getTranscripts(preferred_tx = preferred_tx, db = ucsc_hg19_ncbiRefSeq, bedfile = bedfile)
  tx_names <- unique(unlist(tx$model))
  five_utr  <- GenomicFeatures::fiveUTRsByTranscript(x = ucsc_hg19_ncbiRefSeq, use.names=T)[tx_names]
  three_utr <- GenomicFeatures::threeUTRsByTranscript(x = ucsc_hg19_ncbiRefSeq, use.names=T)[tx_names]
  hgvs <- pos_five_utr(event_bp=GenomicRanges::start(bedfile), bedfile = bedfile, five_utr = five_utr[[1]])
  return(hgvs)
}

run_tests <- function() {
# 5' POSITIVE STRAND
# between UTR exons (upstream)
# expect -40+46
test(start = 32889850) # PASS

# between UTR exons (downstream)
# expect -39-59
test(start = 32890500) # PASS

# upstream of UTR
# expect -274
test(start = 32889570) # PASS

# within UTR exon
# expect -28
test(start = 32890570) # PASS
}


# 3'UTR POSITIVE STRAND
test <- function(start, db, preferred_tx = 'BRCA2.preferredtx') {

  db <- 'inst/extdata/ucsc_hg19_ncbiRefSeq.sqlite'

  ucsc_hg19_ncbiRefSeq <- AnnotationDbi::loadDb(db) %>%
    GenomeInfoDb::keepStandardChromosomes(.)

  GenomeInfoDb::seqlevelsStyle(ucsc_hg19_ncbiRefSeq) <- "NCBI"

  gene <- c("BRCA2")
  tx   <- c("NM_000059.3")

  data.frame(gene, tx) %>%
    write.table(
      x = .,
      file = "BRCA2.preferredtx",
      quote = F,
      row.names = F,
      col.names = F,
      sep = "\t"
    )

  bedfile <- GenomicRanges::GRanges(seqnames = "13", ranges = IRanges::IRanges(start=start, width = 1), strand = '+')
  tx <- getTranscripts(preferred_tx = './BRCA2.preferredtx', db = ucsc_hg19_ncbiRefSeq, bedfile = bedfile)
  tx_names <- unique(unlist(tx$model))
  five_utr  <- GenomicFeatures::fiveUTRsByTranscript(x = ucsc_hg19_ncbiRefSeq, use.names=T)[tx_names]
  three_utr <- GenomicFeatures::threeUTRsByTranscript(x = ucsc_hg19_ncbiRefSeq, use.names=T)[tx_names]
  hgvs <- pos_three_utr(event_bp=GenomicRanges::start(bedfile), bedfile = bedfile, three_utr = three_utr[[1]])
  return(hgvs)
}

# outside 3UTR
# expect *1093
test(start = 32974000) # PASS

#in 3'UTR exon
# expect *93
test(start = 32973000)





# 3'UTR POSITIVE STRAND
test <- function(start) {

  db <- 'inst/extdata/ucsc_hg19_ncbiRefSeq.sqlite'

  ucsc_hg19_ncbiRefSeq <- AnnotationDbi::loadDb(db) %>%
    GenomeInfoDb::keepStandardChromosomes(.)

  GenomeInfoDb::seqlevelsStyle(ucsc_hg19_ncbiRefSeq) <- "NCBI"

  gene <- c("VWA1")
  tx   <- c("NM_199121.2")

  data.frame(gene, tx) %>%
    write.table(
      x = .,
      file = "VWA1.preferredtx",
      quote = F,
      row.names = F,
      col.names = F,
      sep = "\t"
    )

  bedfile <- GenomicRanges::GRanges(seqnames = "1", ranges = IRanges::IRanges(start=start, width = 1), strand = '+')
  tx <- getTranscripts(preferred_tx = './VWA1.preferredtx', db = ucsc_hg19_ncbiRefSeq, bedfile = bedfile)
  tx_names <- unique(unlist(tx$model))
  three_utr <- GenomicFeatures::threeUTRsByTranscript(x = ucsc_hg19_ncbiRefSeq, use.names=T)[tx_names]
  hgvs <- pos_three_utr(event_bp=GenomicRanges::start(bedfile), bedfile = bedfile, three_utr = three_utr[[1]])
  return(hgvs)
}

#in 3'UTR exon
# expect *42-461
test(start = 1374000)

#in 3'UTR exon
# expect *41+136
test(start = 1373000)



# NEGATIVE STRAND

test <- function(start) {

  db <- 'inst/extdata/ucsc_hg19_ncbiRefSeq.sqlite'

  ucsc_hg19_ncbiRefSeq <- AnnotationDbi::loadDb(db) %>%
    GenomeInfoDb::keepStandardChromosomes(.)

  GenomeInfoDb::seqlevelsStyle(ucsc_hg19_ncbiRefSeq) <- "NCBI"

  gene <- c("BRCA1")
  tx   <- c("NM_007300.4")

  data.frame(gene, tx) %>%
    write.table(
      x = .,
      file = "BRCA1.preferredtx",
      quote = F,
      row.names = F,
      col.names = F,
      sep = "\t"
    )

  bedfile <- GenomicRanges::GRanges(seqnames = "17", ranges = IRanges::IRanges(start=start, width = 1), strand = '-')
  tx <- getTranscripts(preferred_tx = 'BRCA1.preferredtx', db = ucsc_hg19_ncbiRefSeq, bedfile = bedfile, flank_length = 1000)
  tx_names <- unique(unlist(tx$model))
  five_utr  <- GenomicFeatures::fiveUTRsByTranscript(x = ucsc_hg19_ncbiRefSeq, use.names=T)[tx_names]
  hgvs <- neg_five_utr(bedfile = bedfile, five_utr = five_utr[[1]])
  return(hgvs)
}

# inside 5'UTR
#-7
test(start = 41276120) # PASS

#-32
test(start = 41277300) # PASS

# BETWEEN UTR
#-20+288
test(start = 41277000) # PASS

#-19-368
test(start = 41276500) # PASS

# UPSTREAM of UTR
#-232
test(start = 41277500) # PASS






test <- function(start) {

  db <- 'inst/extdata/ucsc_hg19_ncbiRefSeq.sqlite'

  ucsc_hg19_ncbiRefSeq <- AnnotationDbi::loadDb(db) %>%
    GenomeInfoDb::keepStandardChromosomes(.)

  GenomeInfoDb::seqlevelsStyle(ucsc_hg19_ncbiRefSeq) <- "NCBI"

  gene <- c("BRCA1")
  tx   <- c("NM_007300.4")

  data.frame(gene, tx) %>%
    write.table(
      x = .,
      file = "BRCA1.preferredtx",
      quote = F,
      row.names = F,
      col.names = F,
      sep = "\t"
    )

  bedfile <- GenomicRanges::GRanges(seqnames = "17", ranges = IRanges::IRanges(start=start, width = 1), strand = '-')
  tx <- getTranscripts(preferred_tx = 'BRCA1.preferredtx', db = ucsc_hg19_ncbiRefSeq, bedfile = bedfile, flank_length = 1000)
  tx_names <- unique(unlist(tx$model))
  three_utr  <- GenomicFeatures::threeUTRsByTranscript(x = ucsc_hg19_ncbiRefSeq, use.names=T)[tx_names]
  hgvs <- neg_three_utr(bedfile = bedfile, three_utr = three_utr[[1]])
  return(hgvs)
}

# 3UTR NEGATIVE STRAND

# downstream of 3UTR on -ve strand
# expect *1395
test(start = 41196300) # PASS

# in 3UTR
# *1295
test(start = 41196400) # PASS

# GENE ON -VE STRAND with 2 x 3UTR EXONS
# GNG5
test <- function(start) {

  db <- 'inst/extdata/ucsc_hg19_ncbiRefSeq.sqlite'

  ucsc_hg19_ncbiRefSeq <- AnnotationDbi::loadDb(db) %>%
    GenomeInfoDb::keepStandardChromosomes(.)

  GenomeInfoDb::seqlevelsStyle(ucsc_hg19_ncbiRefSeq) <- "NCBI"

  gene <- c("GNG5")
  tx   <- c("NM_005274.3")

  data.frame(gene, tx) %>%
    write.table(
      x = .,
      file = "GNG5.preferredtx",
      quote = F,
      row.names = F,
      col.names = F,
      sep = "\t"
    )

  bedfile <- GenomicRanges::GRanges(seqnames = "1", ranges = IRanges::IRanges(start=start, width = 1), strand = '-')
  tx <- getTranscripts(preferred_tx = 'GNG5.preferredtx', db = ucsc_hg19_ncbiRefSeq, bedfile = bedfile, flank_length = 1000)
  tx_names <- unique(unlist(tx$model))
  three_utr  <- GenomicFeatures::threeUTRsByTranscript(x = ucsc_hg19_ncbiRefSeq, use.names=T)[tx_names]
  hgvs <- neg_three_utr(bedfile = bedfile, three_utr = three_utr[[1]])
  return(hgvs)
}

# Upstream og UTR
#*251
test(start = 84964000) # PASS

#*20-769
test(start = 84965000)

#*19+509
test(start = 84967000)


}
