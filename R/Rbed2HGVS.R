
makeCdsDb <- function() {

  ucsc_hg19_ncbiRefSeq <- GenomicFeatures::makeTxDbFromUCSC(genome = "hg19", tablename = "ncbiRefSeq") %>%
    GenomeInfoDb::keepStandardChromosomes(.)
  GenomeInfoDb::seqlevelsStyle(ucsc_hg19_ncbiRefSeq) <- "NCBI"
  AnnotationDbi::saveDb(x = ucsc_hg19_ncbiRefSeq, file = './data/ucsc_hg19_ncbiRefSeq.sqlite')
}

#' appends an additional column of interval in HGVS nomenclature to a bedfile
#'
#' @param bedfile path to bedfile
#'
#' @return GRanges object with HGVS appended
#' @export
#' @importFrom magrittr %>%
Rbed2HGVS <- function(bedfile, db = './data/ucsc_hg19_ncbiRefSeq.sqlite') {

  # load bedfile
  bedfile <- rtracklayer::import.bed(bedfile)

  # load refseq database
  ucsc_hg19_ncbiRefSeq <- AnnotationDbi::loadDb(db) %>%
    GenomeInfoDb::keepStandardChromosomes(.)

  # convert chr to NCBI style
  GenomeInfoDb::seqlevelsStyle(ucsc_hg19_ncbiRefSeq) <- "NCBI"

  # get cds indexed by refseq fron db
  cds_by_tx <- GenomicFeatures::cdsBy(x = ucsc_hg19_ncbiRefSeq, by = "tx", use.name = T)

  # make separate ranges object for start and end coordinates
  bedfile_start <- IRanges::narrow(x = bedfile, start = 1, width = 1)
  bedfile_end   <- IRanges::narrow(x = bedfile, start = -1, width = 1)

  # get cds info for start and end of bedfile
  hgvs_start <- getHgvs(bedfile = bedfile_start, cds = cds_by_tx)
  hgvs_end   <- getHgvs(bedfile = bedfile_end, cds = cds_by_tx)

  # combine start and end for each bed entry
  lapply(seq(bedfile), function(bed_i) {

    # get bed line info
    chr   <- GenomicRanges::seqnames(bedfile[bed_i]) %>% as.vector()
    start <- GenomicRanges::start(bedfile[bed_i])
    end   <- GenomicRanges::end(bedfile[bed_i])

    # check if tx/start has a corresponding tx/end - ignore NA
    s_tx <- hgvs_start[[bed_i]]$tx %>% as.vector() %>% .[!is.na(.)]
    e_tx <- hgvs_end[[bed_i]]$tx %>% as.vector() %>% .[!is.na(.)]

    # get tx missing cds for end
    miss_e <- s_tx[!(s_tx %in% e_tx)]

    # get tx missing cds for end
    miss_s <-e_tx[!(e_tx %in% s_tx)]

    if (!isEmpty(miss_s)) {
      df_miss_s <- getUsDs(missing_tx = miss_s, bedfile = bedfile_start[bed_i], cds = cds_by_tx)
      hgvs_start[[bed_i]] <- rbind(
        hgvs_start[[bed_i]][stats::complete.cases(hgvs_start[[bed_i]]),],
        df_miss_s
        )
    }

    if (!isEmpty(miss_e)) {
      df_miss_e <- getUsDs(missing_tx = miss_e, bedfile = bedfile_end[bed_i], cds = cds_by_tx)
      hgvs_end[[bed_i]] <- rbind(
        hgvs_end[[bed_i]][complete.cases(hgvs_end[[bed_i]]),],
        df_miss_e
        )
    }

    df <- merge.data.frame(
      x = hgvs_start[[bed_i]],
      y = hgvs_end[[bed_i]],
      by = "tx",
      all = TRUE, suffixes = c(".start",".end"))

    gene <- getSymbolRefseq(refSeqId = df$tx)[,"SYMBOL"]

    cbind(chr, start, end, gene, df) %>% return()
    }) %>%
    do.call(rbind, .)
}


#' getUsDs
#'
#' @param missing_tx character vector of transcript names
#' @param bedfile GRanges object containing a single bedfile
#' @param cds
#'
#' @return
getUsDs <- function(missing_tx, bedfile, cds) {

  out <- lapply(missing_tx, function(tx) {

    out <- GenomicRanges::distance(x = bedfile, y = cds[[tx]])
    dist <- min(out)
    side <- GenomicRanges::follow(subject = bedfile, x = cds[[tx]])[which.min(out)]
    exon <- cds[[tx]]$exon_rank[which.min(out)]

    # if positive strand
    if (GenomicRanges::strand(cds[[tx]])[1] %>% as.vector() == '+') {
      # event occuring upstream of cds
      if (is.na(side)) {
        # interval is downstream of exon
        # hgvs includes nearest exon
        entry_cds <- GenomicRanges::width(cds[[tx]][1:which.min(out)]) %>% sum
        hgvs <- paste0("c.", entry_cds, "+", dist)
      } else {
        # hgvs should not include nearest exon. +1 because in relation to first base of
        # nearest exon
        entry_cds <- 1 + GenomicRanges::width(cds[[tx]][1:which.min(out)-1]) %>% sum
        hgvs <- paste0("c.", entry_cds, "-", dist+1)
      }
    } else {
      # if negative strand
      if (is.na(side)) {
        # interval is downstream of exon
        # hgvs should include this exon
        entry_cds <- GenomicRanges::width(cds[[tx]][1:which.min(out)]) %>% sum
        hgvs <- paste0("c.", entry_cds, "+", dist + 1)
      } else {
        # interval is upstream of exon
        # hgvs shouls not includes exon. +1 because in relation to first base of nearest exon
        entry_cds <- 1 + GenomicRanges::width(cds[[tx]][1:which.min(out) - 1]) %>% sum
        hgvs <- paste0("c.", entry_cds, "-", dist + 1) #
      }

    }

    data.frame(tx, exon, hgvs) %>% return()
    })

  do.call(rbind, out) %>%
    return()
  }


getHgvs <- function(bedfile, cds_by_tx) {

  # index of bed overlap with cds
  bed_ol_tx <- IRanges::findOverlaps(query = bedfile, subject = cds_by_tx)

  # extract transcripts and index by bed entry
  cds_ol <- lapply(seq(bedfile), function(x) {
    S4Vectors::queryHits(bed_ol_tx) %in% x %>%
      bed_ol_tx[.] %>%
      S4Vectors::subjectHits(.) %>%
      cds_by_tx[.]
  })

  # loop over each element (bed entry)
  lapply(seq(cds_ol), function(bedln) {

    # if not empty
    if(!isEmpty(cds_ol)[bedln]) {

      # all transcripts for bed entry
      tx   <- names(cds_ol[[bedln]])

      # loop over tx entry - get exon & hgvs
      cds_annot <- lapply(names(cds_ol[[bedln]]), function(tx) {
        mapCoordToCds(range = bedfile[bedln], cds = cds_ol[[bedln]][[tx]])
      })

      exon  <- lapply(cds_annot, function(x) {x$exon_rank}) %>% unlist()
      hgvs  <- lapply(cds_annot, function(x) {x$cds}) %>% unlist() %>% paste0("c.", .)

    } else {
      tx   <- NA
      exon <- NA
      hgvs <- NA
    }

    data.frame(tx, exon, hgvs)
  })

}



#' mapCoordToCds
#'
#' @param range
#' @param cds
#'
#' @return
mapCoordToCds <- function(range, cds) {

  # find which exon overlaps
  hits <- IRanges::findOverlaps(query = range, subject = cds)

  # get that exon
  cds_hit_ln <-cds[S4Vectors::subjectHits(hits)]

  # how far into exon hit is bed entry
  # for positive strand
  if (GenomicRanges::strand(cds_hit_ln) %>% as.vector() == '+') {
    cds_width_hit <- IRanges::start(range) - IRanges::start(cds_hit_ln)
    # accum width of cds upstream of exon for + strand
    cds_width_us <- IRanges::width(cds)[1:(S4Vectors::subjectHits(hits)) - 1] %>% sum()
  } else {
    # negative strand
    cds_width_hit <- IRanges::end(cds_hit_ln) - IRanges::start(range)
    cds_width_us <- IRanges::width(cds)[1:(S4Vectors::subjectHits(hits)) - 1] %>% sum()
  }

  # total cds
  cds_width_final <- sum(cds_width_us, cds_width_hit) + 1

  cds_hit_ln$cds <- cds_width_final

  return(cds_hit_ln)
}


#' getSymbolRefSeq
#'
#' @param refSeqId
#'
#' @return
#' @import org.Hs.eg.db
getSymbolRefseq <- function(refSeqId) {

  # remove suffix
  stringr::str_remove(string = refSeqId, pattern = "\\.\\d+") %>%
    AnnotationDbi::select(org.Hs.eg.db, keytype = "REFSEQ", keys = ., columns = "SYMBOL") %>%
    return()

}

