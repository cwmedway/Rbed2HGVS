
parseCmdArgs <- function() {

  parser <- argparse::ArgumentParser(description="Add Transcript, CDS & Exon annotations to bedfile")

  parser$add_argument("-b", "--bedfile", help="path to BED file", required=T, type="character")
  parser$add_argument("-o", "--outname", help="output file name", default="hgvs", type="character")
  parser$add_argument("-O", "--outdir", help="output file directory", default="./", type="character")

  parser$parse_args() %>%
    return()
}


#' load a bedfile
#'
#' @param bedfile path to bedfile
#'
#' @return GRanges object
loadBedFile <- function(bedfile) {

  # load bedfile as granges object
  rtracklayer::import.bed(bedfile) %>%
    return()
}


#' create ensembldb from GTF
#'
#' @param gtf
#'
#' @return
#'
makeDbFromGtf <- function(gtf) {

  db <- ensembldb::ensDbFromGtf(gtf)

    edb <- ensembldb::EnsDb("./Homo_sapiens.GRCh37.87.sqlite")

}


#' load db
#'
#' @param db_path path
#'
#' @return IRanges object
loadDbHGVS <- function(db_path) {

  ensembldb::EnsDb(db_path) %>%
    return()
}


#' annotate bedfile with transcript information
#'
#' @param bedfile GRanges
#' @param db ensembl db object
#'
#' @return IRanges object
bedToTx <- function(bedfile, db) {

  results <- ensembldb::genomeToTranscript(x = bedfile, db = db)
  return(results)
}


#' annotate bedfile with CDS information
#'
#' @param results IRanges
#' @param db ebsembl db object
#'
#' @return IRanges object
txToCds <- function(results, db) {

  lapply(seq(results), function(i) {

    # db query
    tx2cds <- ensembldb::transcriptToCds(x = results[[i]], db = db)

    # parse db results
    cds <- tx2cds %>% IRanges::start()
    tx_id     <- tx2cds %>% GenomicRanges::mcols() %>% .['tx_id']
    exon_rank <- tx2cds %>% GenomicRanges::mcols() %>% .['exon_rank']

    # format report
    data.frame(tx_id, cds, exon_rank) %>% return()
  })
}


#' annotates a bedfile with Transcript, exon & CDS information
#'
#' @param bedfile path to bedfile
#' @param outname output file name
#' @param outdir output file location
#' @param db ensembl db object
#'
#' @return writes report
#' @importFrom magrittr %>%
#' @export
bedToCds <- function(bedfile, outname, outdir, db_path) {

  bedfile <- loadBedFile(bedfile = args$bedfile)

  edb <- loadDbHGVS(db_path)

  bedfile_start <- IRanges::narrow(x = bedfile, start = 1, width = 1)
  bedfile_end   <- IRanges::narrow(x = bedfile, start = -1, width = 1)

  results_start <- bedToTx(bedfile = bedfile_start, db = db) %>% txToCds(results = ., db = db)
  results_end   <- bedToTx(bedfile = bedfile_end, db = db) %>% txToCds(results = ., db = db)

  combineEndsBed(bed_a = results_start, bed_b = results_end, bedfile = bedfile) %>%
    do.call(rbind,.) %>%
    readr::write_tsv(x = ., path = paste0(outdir,"/", outname, ".tsv"))
}


#' Combines start & end coordinates
#'
#' @param bed_a
#' @param bed_b
#' @param bedfile
#'
#' @return data.frame
combineEndsBed <- function(bed_a, bed_b, bedfile) {

  lapply(seq(bed_a), function(i) {

    # parse bedfile
    chr   <- bedfile[i] %>% GenomicRanges::seqnames()
    start <- bedfile[i] %>% IRanges::start()
    end   <- bedfile[i] %>% IRanges::end()

    df <- merge.data.frame(x = bed_a[[i]], y = bed_b[[i]], by = "tx_id", suffixes = c('_start','_end'), all = T)

    cbind(chr, start, end, df) %>% return()

  })
}

# args <- parseCmdArgs()
#
# o <- bedToCds(
#   bedfile = bedfile,
#   outname = args$outname,
#   outdir = args$outdir,
#   db_path = '/home/chris/Apps/Rbed2HGVS/data/Homo_sapiens.GRCh37.87.sqlite'
#   )


