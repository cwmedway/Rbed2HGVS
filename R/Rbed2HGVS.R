
# parseCmdArgs <- function() {
#
#   parser <- argparse::ArgumentParser(description="Add Transcript, CDS & Exon annotations to bedfile")
#
#   parser$add_argument("-b", "--bedfile", help="path to BED file", required=T, type="character")
#   parser$add_argument("-o", "--outname", help="output file name", default="hgvs", type="character")
#   parser$add_argument("-O", "--outdir", help="output file directory", default="./", type="character")
#
#   parser$parse_args() %>%
#     return()
# }

# ucsc_hg19_ncbiRefSeq <- GenomicFeatures::makeTxDbFromUCSC(genome = "hg19", tablename = "ncbiRefSeq") %>%
#   GenomeInfoDb::keepStandardChromosomes(.)
# GenomeInfoDb::seqlevelsStyle(ucsc_hg19_ncbiRefSeq) <- "NCBI"
# AnnotationDbi::saveDb(x = ucsc_hg19_ncbiRefSeq, file = './data/ucsc_hg19_ncbiRefSeq.sqlite')



Rbed2HGVS <- function(bedfile, db) {

  # load bedfile
  bedfile <- rtracklayer::import.bed(bedfile)

  # loadDb
  ucsc_hg19_ncbiRefSeq <- AnnotationDbi::loadDb(db)
  GenomeInfoDb::seqlevelsStyle(ucsc_hg19_ncbiRefSeq) <- "NCBI"

  # get cds fron db
  cds_by_tx <- GenomicFeatures::cdsBy(x = ucsc_hg19_ncbiRefSeq, by = "tx", use.name = T)

  # make separate ranges object for start and end
  bedfile_start <- IRanges::narrow(x = bedfile, start = 1, width = 1)
  bedfile_end   <- IRanges::narrow(x = bedfile, start = -1, width = 1)

  # get cds info for start and end of bedfile
  s <- getHgvs(bedfile = bedfile_start, cds = cds_by_tx)
  e <- getHgvs(bedfile = bedfile_end, cds = cds_by_tx)

  # combine start and end for each bed entry
  lapply(seq(bedfile), function(bed_i) {

    chr   <- GenomicRanges::seqnames(bedfile[bed_i]) %>% as.vector()
    start <- GenomicRanges::start(bedfile[bed_i])
    end   <- GenomicRanges::end(bedfile[bed_i])

    # check if tx/start has a corresponding tx/end - ignore NA
    s_tx <- s[[bed_i]]$tx %>% as.vector() %>% .[!is.na(.)]
    e_tx <- e[[bed_i]]$tx %>% as.vector() %>% .[!is.na(.)]

    # tx missing cds for end
    miss_e <- s_tx[!(s_tx %in% e_tx)]

    # tx missing cds for end
    miss_s <-e_tx[!(e_tx %in% s_tx)]

    if (!isEmpty(miss_s)) {
      df_miss_s <- getUsDs(missing_tx = miss_s, bedfile = bedfile_start[bed_i], cds = cds_by_tx)
      s[[bed_i]] <- rbind(
        s[[bed_i]][complete.cases(s[[bed_i]]),],
        df_miss_s
        )
    }

    if (!isEmpty(miss_e)) {
      df_miss_e <- getUsDs(missing_tx = miss_e, bedfile = bedfile_end[bed_i], cds = cds_by_tx)
      e[[bed_i]] <- rbind(
        e[[bed_i]][complete.cases(e[[bed_i]]),],
        df_miss_e
        )
    }

    df <- merge.data.frame(
      x = s[[bed_i]],
      y = e[[bed_i]],
      by = "tx",
      all = TRUE, suffixes = c(".start",".end"))

    gene <- getSymbolRefseq(refSeqId = df$tx)[,"SYMBOL"]

    cbind(chr, start, end, gene, df) %>% return()
    }) %>%
    do.call(rbind, .)
}


getUsDs <- function(missing_tx, bedfile, cds) {

  out <- lapply(missing_tx, function(tx) {

    out <- GenomicRanges::distance(x = bedfile, y = cds[[tx]])
    dist <- min(out)
    side <- GenomicRanges::follow(subject = bedfile, x = cds[[tx]])[which.min(out)]
    exon <- cds[[tx]]$exon_rank[which.min(out)]

    if (is.na(side)) {
      # interval is downstream of exon
      dist <- (0 - dist)
      # cds should not include this exon
      entry_cds <- 1 + GenomicRanges::width(cds[[tx]][1:which.min(out)-1]) %>% sum
      hgvs <- paste0("c.", entry_cds, "-", dist)
    } else {
      # interval is upstream of exon
      # cds includes exon
      entry_cds <- GenomicRanges::width(cds[[tx]][1:which.min(out)]) %>% sum
      hgvs <- paste0("c.", entry_cds, "+", dist)
    }

    data.frame(tx, exon, hgvs) %>% return()

    })

  do.call(rbind, out) %>%
    return()
  }

getHgvs <- function(bedfile, cds) {

  # overlap with cds
  ol <- IRanges::findOverlaps(query = bedfile, subject = cds)

  # list transcripts indexed by bed entry
  cds_ol <- lapply(seq(bedfile), function(x) {
    S4Vectors::queryHits(ol) %in% x %>%
      ol[.] %>%
      S4Vectors::subjectHits(.) %>%
      cds[.]
  })

  # loop over bed entry
  lapply(seq(cds_ol), function(bedln) {

    if(!isEmpty(cds_ol)[bedln]) {

      # loop over tx entry - get cds
      cds_annot <- lapply(names(cds_ol[[bedln]]), function(tx) {
        mapCoordToCds(range = bedfile[bedln], cds = cds_ol[[bedln]][[tx]])
      })

      tx   <- names(cds_ol[[bedln]])
      #gene <- getSymbolRefseq(refSeqId = tx)[,"SYMBOL"]
      exon <- lapply(cds_annot, function(x) {x$exon_rank}) %>% unlist()
      hgvs  <- lapply(cds_annot, function(x) {x$cds}) %>% unlist() %>% paste0("c.", .)
    } else {
      tx   <- NA
      exon <- NA
      hgvs <- NA
    }

    data.frame(tx, exon, hgvs)
  })

}





mapCoordToCds <- function(range, cds) {

  # find which exon overlaps
  hits <- IRanges::findOverlaps(query = range, subject = cds)

  # get that exon
  cds_hit_ln <-cds[S4Vectors::subjectHits(hits)]

  # how far into exon hit is bed entry
  cds_width_hit <- IRanges::end(cds_hit_ln) - IRanges::start(range)

  # accum width of cds upstream of exon
  cds_width_us <- IRanges::width(cds)[1:(S4Vectors::subjectHits(hits)) - 1] %>% sum()

  # total cds
  cds_width_final <- sum(cds_width_us, cds_width_hit)

  cds_hit_ln$cds <- cds_width_final

  return(cds_hit_ln)
}

getSymbolRefseq <- function(refSeqId) {

  # remove suffix
  stringr::str_remove(string = refSeqId, pattern = "\\.\\d+") %>%
    AnnotationDbi::select(x = org.Hs.eg.db, keytype = "REFSEQ", keys = ., columns = "SYMBOL") %>%
    return()

}

out <- Rbed2HGVS(bedfile = './data/test.bed', db = './data/ucsc_hg19_ncbiRefSeq.sqlite')

