
# makeCdsDb <- function() {
#
#   ucsc_hg19_ncbiRefSeq <- GenomicFeatures::makeTxDbFromUCSC(genome = "hg19", tablename = "ncbiRefSeq") %>%
#     GenomeInfoDb::keepStandardChromosomes(.)
#   GenomeInfoDb::seqlevelsStyle(ucsc_hg19_ncbiRefSeq) <- "NCBI"
#   AnnotationDbi::saveDb(x = ucsc_hg19_ncbiRefSeq, file = './data/ucsc_hg19_ncbiRefSeq.sqlite')
# }

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

  # cds GRangesList indexed by RefSeq
  cds_by_tx <- GenomicFeatures::cdsBy(x = ucsc_hg19_ncbiRefSeq, by = "tx", use.name = T)

  # trim to only include necessary transcripts
  cds_by_tx <- cds_by_tx[names(cds_by_tx) %in% unique(unlist(tx$model))]

  # get cds info for start and end of bedfile
  hgvs <- getHgvs(bedfile = bedfile, cds_by_tx = cds_by_tx, cds_ol = tx[['model']], ncores = ncores) %>% do.call(rbind, .)

  #append HGMD
  hgvs$gene <- getSymbolRefseq(refSeqId = hgvs$tx)

  list(
    "hgvs" = hgvs,
    "missing" = tx[['missing']],
    "version" = tx[['version']]
    ) %>% return()
}


getTranscripts <- function(preferred_tx, db, bedfile) {

  refseq_tx <- GenomicFeatures::transcripts(x = db)

  if (is.character(preferred_tx)) {
    # read preferred transcripts (ptx)
    df_ptx <- read.table(file = preferred_tx, header = F, stringsAsFactors = F)
    # ptx name without version suffix
    ptx_prefix <- stringr::str_remove(string = df_ptx[,2], pattern = '\\.\\d+')
    # names of tx in RefSeq DB
    cds_name   <- refseq_tx$tx_name
    cds_prefix <- stringr::str_remove(string = cds_name, pattern = "\\.\\d+")

    # find preferred tx in RedSeq DB
    m <- match(x = ptx_prefix, table = cds_prefix)

    # preferred tx not in RefSeq DB
    ptx_missing <- df_ptx[is.na(m),]

    if (dim(ptx_missing)[1] > 0) {
      warning("Some preferred transcripts are missing from the RefSeq DB")
    }

    # df containing only matched ptx
    ptx_not_missing <- df_ptx[!is.na(m),]
    # extract full tx names from RefSeq DB (inc version)
    cds_not_missing <- cds_name[m[!is.na(m)]]
    # check if versions match
    not_same_version <- !(ptx_not_missing[,2] %in% cds_not_missing)

    if (sum(not_same_version) > 0) {
      warning("Some preferred transcripts are a different version from the RefSeq DB")
      ptx_version <- data.frame(
        "gene" = ptx_not_missing[not_same_version,1],
        "preferred_tx" = ptx_not_missing[not_same_version,2],
        "db_tx" = cds_not_missing[not_same_version] )
    } else { ptx_version <- NA }

    refseq_tx <- refseq_tx[m[!is.na(m)]]
    # index of bed overlap with transcripts
    bed_ol_tx <- IRanges::findOverlaps(query = bedfile, subject = refseq_tx)

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
    bed_ol_tx <- IRanges::findOverlaps(query = bedfile, subject = refseq_tx)

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


getHgvs <- function(bedfile, cds_by_tx, cds_ol, ncores) {

  # loop over each element (bed entry)
  parallel::mclapply(seq(cds_ol), mc.cores = ncores, function(bedln) {

    # will be first three columns of output
    chr   <- GenomicRanges::seqnames(bedfile[bedln]) %>% as.vector()
    start <- GenomicRanges::start(bedfile[bedln]) - 1 # -1 because GRanges representation of BED different to UCSC (1 vs 0 based)
    end   <- GenomicRanges::end(bedfile[bedln])

    # if bed entry overlaps at least one transcript
    if(!isEmpty(cds_ol)[bedln]) {

      # loop over each tx - get exon & hgvs
      cds_annot <- lapply(cds_ol[[bedln]], function(tx_i) {
        mapCoordToCds(bedfile = bedfile[bedln], cds = cds_by_tx[[tx_i]])
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

mapCoordToCds <- function(bedfile, cds) {

  # only perform if transcript has cds
  if (!is.null(cds)) {

    # make separate ranges object for start and end coordinates
    bedfile_start <- IRanges::narrow(x = bedfile, start = 1, width = 1)
    bedfile_end   <- IRanges::narrow(x = bedfile, start = -1, width = 1)

    hgvs_start <- getHgvs2(bedfile = bedfile_start, cds = cds)
    hgvs_end   <- getHgvs2(bedfile = bedfile_end, cds = cds)

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


getHgvs2 <- function(bedfile, cds) {

  # find which exon overlaps
  dist_to_cds <- IRanges::distance(x = bedfile, y = cds)

  # which exon
  near_ex_i <- which.min(dist_to_cds)

  # get dist to event (0-based - if using in HGVS need to add 1bp)
  dist <- dist_to_cds[near_ex_i]

  # get that exon
  near_ex <- cds[near_ex_i]

  # exon number
  exon_rank <- near_ex$exon_rank

  # which strand
  strand  <- GenomicRanges::strand(near_ex) %>% as.vector()

  if ( dist == 0 ) {
    # coordinate within cds
    if ( strand == "+") {
      # on positive strand
      # distance into cds hit
      cds_width_hit <- (1 + IRanges::start(bedfile) - IRanges::start(near_ex))
      # accum width of cds upstream of exon for + strand
      cds_width_us <- IRanges::width(cds)[1:near_ex_i - 1] %>% sum()
      hgvs <- paste0(sum(cds_width_hit + cds_width_us))
    } else {
      # on negative strand
      # distance into cds hit
      cds_width_hit <- (1 + IRanges::end(near_ex) - IRanges::start(bedfile))
      # accum width of cds upstream of exon for - strand
      cds_width_us <- IRanges::width(cds)[1:near_ex_i - 1] %>% sum()
      hgvs <- paste0(sum(cds_width_hit + cds_width_us))
    }
  } else {
    # coordinate is flanking cds
    # if positive strand
    if (strand == '+') {
      # event occuring upstream of cds
      if (IRanges::start(bedfile) < IRanges::start(near_ex)) {
        # interval is upstream of exon. hgvs should not include nearest exon
        # +1 because in relation to first base of upstream exon
        entry_cds <- 1 + IRanges::width(cds[1:near_ex_i - 1]) %>% sum
        hgvs <- paste0(entry_cds, "-", dist + 1)
      } else {
        # downstream, hgvs includes nearest exon
        entry_cds <- IRanges::width(cds[1:(near_ex_i)]) %>% sum
        hgvs <- paste0(entry_cds, "+", dist + 1)
      }
    } else {
      # if negative strand
      if (IRanges::start(bedfile) < IRanges::start(near_ex)) {
        # interval is downstream of exon
        # hgvs should include this exon
        entry_cds <- IRanges::width(cds[1:near_ex_i]) %>% sum
        hgvs <- paste0(entry_cds, "+", dist + 1)
      } else {
        # interval is upstream of exon
        # hgvs shouls not includes exon. +1 because in relation to first base of nearest exon
        entry_cds <- 1 + IRanges::width(cds[1:near_ex_i - 1]) %>% sum
        hgvs <- paste0(entry_cds, "-", dist + 1)
      }
    }
  }
  return(list("hgvs" = hgvs, "exon_rank" = exon_rank))
}


getSymbolRefseq <- function(refSeqId) {

  # remove suffix
  stringr::str_remove(string = refSeqId, pattern = "\\.\\d+") %>%
    AnnotationDbi::select(org.Hs.eg.db, keytype = "REFSEQ", keys = ., columns = "SYMBOL") %>%
    .[,"SYMBOL"] %>%
    return()
}

