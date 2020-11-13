# make bed file
make_test_bed <- function(chr,bp) {
   GenomicRanges::GRanges(seqnames = chr, ranges = IRanges::IRanges(start = bp, width = 1)) %>%
      return()
}

make_test_ptx <- function(gene, tx) {
   data.frame(gene, tx) %>%
      write.table(
         .,
         file = 'ptx.txt',
         quote = F, sep = '\t',
         col.names = F,
         row.names = F
      )
}

run_test <- function(chr, bp, gene=NA, tx=NA) {

   if (is.na(gene)) {
      ptx <- NA
   } else {
      ptx <- 'ptx.txt'
      make_test_ptx(gene,tx)
   }

   Rbed2HGVS::Rbed2HGVS(
      bedfile = make_test_bed(chr = chr, bp = bp),
      preferred_tx = ptx
   ) %>%
      return()
}

get_hgvs <- function(out) {
   out$hgvs$hgvs_start %>% as.character() %>%
      return()
}

get_exon <- function(out) {
   out$hgvs$exon_start %>% as.character() %>%
      return()
}

check_test <- function(test, exp_hgvs, exp_exon) {
   hgvs <- get_hgvs(test)
   exon <- get_exon(test)
   testthat::expect_equal(object = hgvs, expected = exp_hgvs)
   testthat::expect_equal(object = exon, expected = exp_exon)
   }



testthat::test_that("position downstream of cds exon on negative strand", {

   out <- run_test(chr = 17, bp = 7573912, gene = "TP53", tx = "NM_000546.5")
   check_test(test = out, exp_hgvs = '1100+15', exp_exon = "10")
})



# test_that("cds gap that finishes downstream on negative strand", {
#
#   data.frame("17",	"7573912",	"7573940") %>%
#     write.table(
#       x = .,
#       file = "test-neg-ds.bed",
#       quote = F,
#       row.names = F,
#       col.names = F,
#       sep = "\t"
#     )
#
#   data.frame("TP53", "NM_000546.5") %>%
#     write.table(
#       x = .,
#       file = "test-neg-ds.preferredtx",
#       quote = F,
#       row.names = F,
#       col.names = F,
#       sep = "\t"
#     )
#
#   out <- Rbed2HGVS::Rbed2HGVS(
#     bedfile = 'test-neg-ds.bed',
#     db = '../../inst/extdata/ucsc_hg19_ncbiRefSeq.sqlite',
#     preferred_tx = 'test-neg-ds.preferredtx'
#   )
#
#   testthat::expect_equal(
#     out$hgvs$hgvs_end %>% as.vector(),
#     c("1100+14")
#     )
#
#   testthat::expect_equal(out$hgvs$exon_start, c(10))
#
#   testthat::expect_equal(
#     out$hgvs$hgvs_start %>% as.vector(),
#     c("1087")
#     )
#
#   testthat::expect_equal(out$hgvs$exon_end, c(10))
#
#   unlink(c('test-neg-ds.bed', 'test-neg-ds.preferredtx'))
# })
#
#
#
# test_that("cds gap that finishes upstream on positive strand", {
#
#   data.frame("13", "32910371", "32910431") %>%
#     write.table(
#       x = .,
#       file = "test-pos-us.bed",
#       quote = F,
#       row.names = F,
#       col.names = F, sep = "\t"
#     )
#
#   out <- Rbed2HGVS(
#     bedfile = 'test-pos-us.bed',
#     db = '../../inst/extdata/ucsc_hg19_ncbiRefSeq.sqlite'
#   )
#
#   testthat::expect_equal(out$hgvs$hgvs_start %>% as.vector(), "1910-30")
#   testthat::expect_equal(out$hgvs$exon_start, 11)
#   testthat::expect_equal(out$hgvs$hgvs_end %>% as.vector(), "1939")
#   testthat::expect_equal(out$hgvs$exon_end, 11)
#
#   unlink('test-pos-us.bed')
# })
#
#
# test_that("cds gap that finishes upstream on negative strand", {
#
#   data.frame("17",	"7574017",	"7574052") %>%
#     write.table(
#       x = .,
#       file = "test-neg-us.bed",
#       quote = F,
#       row.names = F,
#       col.names = F,
#       sep = "\t"
#     )
#
#   data.frame("TP53", "NM_000546.5") %>%
#     write.table(
#       x = .,
#       file = "test-neg-us.preferredtx",
#       quote = F,
#       row.names = F,
#       col.names = F,
#       sep = "\t"
#     )
#
#   out <- Rbed2HGVS(
#     bedfile = 'test-neg-us.bed',
#     db = '../../inst/extdata/ucsc_hg19_ncbiRefSeq.sqlite',
#     preferred_tx = 'test-neg-us.preferredtx'
#   )
#
#   testthat::expect_equal(
#     out$hgvs$hgvs_end %>% as.vector(),
#     c("1009")
#     )
#
#   testthat::expect_equal(out$hgvs$exon_start, c(10))
#
#   testthat::expect_equal(
#     out$hgvs$hgvs_start %>% as.vector(),
#     c("994-19")
#     )
#
#   testthat::expect_equal(out$hgvs$exon_end, c(10))
#
#   unlink(c('test-neg-us.bed', 'test-neg-us.preferredtx'))
# })
#
# # ---------------------------------------------------------
# # Issue 1: gap falls outside transcript (promotor i.e TERT)
#
#  test_that("cds gap that finishes upstream on negative strand", {
#
#    chr   <- c(5)
#    start <- c(1295226, 1295248)
#    end   <- c(1295228, 1295250)
#    meta  <- c("TERT(NM_198253.2):c.-124_-124", "TERT(NM_198253.2):c.-146_-146")
#
#    data.frame(chr, start, end, meta) %>%
#      write.table(
#        x = .,
#        file = "issue1.bed",
#        quote = F,
#        row.names = F,
#        col.names = F,
#        sep = "\t"
#      )
#
#    data.frame("TERT", "NM_198253.2") %>%
#      write.table(
#        x = .,
#        file = "issue1.preferredtx",
#        quote = F,
#        row.names = F,
#         col.names = F,
#        sep = "\t"
#      )
#
#    out <- Rbed2HGVS(
#      bedfile = 'issue1.bed',
#      db = '../../inst/extdata/ucsc_hg19_ncbiRefSeq.sqlite',
#      preferred_tx = 'issue1.preferredtx'
#    )
#
#    testthat::expect_equal(
#      out$hgvs$hgvs_end %>% as.vector(),
#      c("1-123","1-145")
#    )
#
#    testthat::expect_equal(out$hgvs$exon_start, c(1,1))
#
#    testthat::expect_equal(
#      out$hgvs$hgvs_start %>% as.vector(),
#      c("1-124", "1-146")
#    )
#
#    testthat::expect_equal(out$hgvs$exon_end, c(1,1))
#
#    unlink(c('issue1.bed', 'issue1.preferredtx'))
#  })
#
# # # ----------------------------------------------------
# # # ISSUE 2 - gap between UTR exons
#
#  chr   <- c("13")
#  start <- c(32889850)
#  end   <- c(32890500)
#  meta  <- c(
#    "BRCA2(NM_000059.3):c.-40+46_-39-59"
#   )
#
# data.frame(chr, start, end, meta) %>%
#   write.table(
#     x = .,
#     file = "issue2.bed",
#     quote = F,
#     row.names = F,
#     col.names = F,
#     sep = "\t"
#   )
#
# gene <- c("BRCA2")
# tx   <- c("NM_000059.3")
#
# data.frame(gene, tx) %>%
#   write.table(
#     x = .,
#     file = "issue2.preferredtx",
#     quote = F,
#     row.names = F,
#     col.names = F,
#     sep = "\t"
#   )
#
# out <- Rbed2HGVS(
#   bedfile = 'issue2.bed',
#   db = '../../inst/extdata/ucsc_hg19_ncbiRefSeq.sqlite',
#   preferred_tx = 'issue2.preferredtx'
# )
#
#  testthat::expect_equal(
#    out$hgvs$hgvs_start %>% as.vector(),
#    c("-40+46")
#  )
#
#  testthat::expect_equal(
#    out$hgvs$hgvs_end %>% as.vector(),
#    c("-39-59")
#  )
#
#
#  unlink(c("issue2.bed","issue2.preferredtx"))
#
# # # ----------------------
# # Upstream of 5UTR to inside UTR
# chr   <- c("13")
# start <- c(32889570)
# end   <- c(32890570)
# meta  <- c(
#    "BRCA2(NM_000059.3):c.-274_-28"
# )
#
#  data.frame(chr, start, end, meta) %>%
#    write.table(
#      x = .,
#      file = "issue3.bed",
#      quote = F,
#      row.names = F,
#      col.names = F,
#      sep = "\t"
#    )
#
#  gene <- c("BRCA2")
#  tx   <- c("NM_000059.3")
#
#  data.frame(gene, tx) %>%
#    write.table(
#      x = .,
#      file = "issue3.preferredtx",
#      quote = F,
#      row.names = F,
#      col.names = F,
#      sep = "\t"
#    )
#
#  out <- Rbed2HGVS(
#    bedfile = 'issue3.bed',
#    db = '../../inst/extdata/ucsc_hg19_ncbiRefSeq.sqlite',
#    preferred_tx = 'issue3.preferredtx'
#  )
#
#  testthat::expect_equal(
#    out$hgvs$hgvs_start %>% as.vector(),
#    c("-274")
#  )
#
#  testthat::expect_equal(
#    out$hgvs$hgvs_end %>% as.vector(),
#    c("-28")
#  )
#
#  unlink(c("issue3.bed","issue3.preferredtx"))
#
# # ----------------------------
# #5' UTR on negative strand
#  chr   <- c("17")
#  start <- c(41276000, 41276050, 41277000)
#  end   <- c(41276120, 41276500, 41278000)
#  meta  <- c(
#    "BRCA1(NM_007300.4):c.80+34_-7",
#    "BRCA1(NM_007300.4):c.64_-19-368",
#    "BRCA1(NM_007300.4):c.-20+288_-732"
#  )
#
#  data.frame(chr, start, end, meta) %>%
#    write.table(
#      x = .,
#      file = "issue4.bed",
#      quote = F,
#      row.names = F,
#      col.names = F,
#      sep = "\t"
#    )
#
#  gene <- c("BRCA1")
#  tx   <- c("NM_007300.4")
#
#  data.frame(gene, tx) %>%
#    write.table(
#      x = .,
#      file = "issue4.preferredtx",
#      quote = F,
#      row.names = F,
#      col.names = F,
#      sep = "\t"
#    )
#
#  out <- Rbed2HGVS(
#    bedfile = 'issue4.bed',
#    db = '../../inst/extdata/ucsc_hg19_ncbiRefSeq.sqlite',
#    preferred_tx = 'issue4.preferredtx'
#  )
#
#
#
# # # -----------------
# # # 3'UTR
# # # Positive Strand
# #
# # chr   <- c("1")
# # start <- c(1373000)
# # end   <- c(1375000)
# # meta  <- c(
# #   "VWA1(NM_199121.2):c.*41+136_*581"
# # )
# #
# # data.frame(chr, start, end, meta) %>%
# #   write.table(
# #     x = .,
# #     file = "issue4.bed",
# #     quote = F,
# #     row.names = F,
# #     col.names = F,
# #     sep = "\t"
# #   )
# #
# # gene <- c("VWA1")
# # tx   <- c("NM_199121.2")
# #
# # data.frame(gene, tx) %>%
# #   write.table(
# #     x = .,
# #     file = "issue4.preferredtx",
# #     quote = F,
# #     row.names = F,
# #     col.names = F,
# #     sep = "\t"
# #   )
# #
# # out <- Rbed2HGVS(
# #   bedfile = 'issue4.bed',
# #   db = '../../inst/extdata/ucsc_hg19_ncbiRefSeq.sqlite',
# #   preferred_tx = 'issue4.preferredtx'
# # )
#
#
#
#  # -----------------------------
#  # TESTS
#  # BRCA2 (Positive Strand)
#
#
#
#  test <- function() {
#
#    test <- function(start, db, preferred_tx = 'BRCA2.preferredtx') {
#
#      db <- 'inst/extdata/ucsc_hg19_ncbiRefSeq.sqlite'
#
#      ucsc_hg19_ncbiRefSeq <- AnnotationDbi::loadDb(db) %>%
#        GenomeInfoDb::keepStandardChromosomes(.)
#
#      GenomeInfoDb::seqlevelsStyle(ucsc_hg19_ncbiRefSeq) <- "NCBI"
#
#      gene <- c("BRCA2")
#      tx   <- c("NM_000059.3")
#
#      data.frame(gene, tx) %>%
#        write.table(
#          x = .,
#          file = "BRCA2.preferredtx",
#          quote = F,
#          row.names = F,
#          col.names = F,
#          sep = "\t"
#        )
#
#      bedfile <- GenomicRanges::GRanges(seqnames = "13", ranges = IRanges::IRanges(start=start, width = 1), strand = '+')
#      tx <- getTranscripts(preferred_tx = preferred_tx, db = ucsc_hg19_ncbiRefSeq, bedfile = bedfile)
#      tx_names <- unique(unlist(tx$model))
#      five_utr  <- GenomicFeatures::fiveUTRsByTranscript(x = ucsc_hg19_ncbiRefSeq, use.names=T)[tx_names]
#      three_utr <- GenomicFeatures::threeUTRsByTranscript(x = ucsc_hg19_ncbiRefSeq, use.names=T)[tx_names]
#      hgvs <- pos_five_utr(event_bp=GenomicRanges::start(bedfile), bedfile = bedfile, five_utr = five_utr[[1]])
#      return(hgvs)
#    }
#
#    run_tests <- function() {
#      # 5' POSITIVE STRAND
#      # between UTR exons (upstream)
#      # expect -40+46
#      test(start = 32889850) # PASS
#
#      # between UTR exons (downstream)
#      # expect -39-59
#      test(start = 32890500) # PASS
#
#      # upstream of UTR
#      # expect -274
#      test(start = 32889570) # PASS
#
#      # within UTR exon
#      # expect -28
#      test(start = 32890570) # PASS
#    }
#
#
#    # 3'UTR POSITIVE STRAND
#    test <- function(start, db, preferred_tx = 'BRCA2.preferredtx') {
#
#      db <- 'inst/extdata/ucsc_hg19_ncbiRefSeq.sqlite'
#
#      ucsc_hg19_ncbiRefSeq <- AnnotationDbi::loadDb(db) %>%
#        GenomeInfoDb::keepStandardChromosomes(.)
#
#      GenomeInfoDb::seqlevelsStyle(ucsc_hg19_ncbiRefSeq) <- "NCBI"
#
#      gene <- c("BRCA2")
#      tx   <- c("NM_000059.3")
#
#      data.frame(gene, tx) %>%
#        write.table(
#          x = .,
#          file = "BRCA2.preferredtx",
#          quote = F,
#          row.names = F,
#          col.names = F,
#          sep = "\t"
#        )
#
#      bedfile <- GenomicRanges::GRanges(seqnames = "13", ranges = IRanges::IRanges(start=start, width = 1), strand = '+')
#      tx <- getTranscripts(preferred_tx = './BRCA2.preferredtx', db = ucsc_hg19_ncbiRefSeq, bedfile = bedfile)
#      tx_names <- unique(unlist(tx$model))
#      five_utr  <- GenomicFeatures::fiveUTRsByTranscript(x = ucsc_hg19_ncbiRefSeq, use.names=T)[tx_names]
#      three_utr <- GenomicFeatures::threeUTRsByTranscript(x = ucsc_hg19_ncbiRefSeq, use.names=T)[tx_names]
#      hgvs <- pos_three_utr(event_bp=GenomicRanges::start(bedfile), bedfile = bedfile, three_utr = three_utr[[1]])
#      return(hgvs)
#    }
#
#    # outside 3UTR
#    # expect *1093
#    test(start = 32974000) # PASS
#
#    #in 3'UTR exon
#    # expect *93
#    test(start = 32973000)
#
#
#
#
#
#    # 3'UTR POSITIVE STRAND
#    test <- function(start) {
#
#      db <- 'inst/extdata/ucsc_hg19_ncbiRefSeq.sqlite'
#
#      ucsc_hg19_ncbiRefSeq <- AnnotationDbi::loadDb(db) %>%
#        GenomeInfoDb::keepStandardChromosomes(.)
#
#      GenomeInfoDb::seqlevelsStyle(ucsc_hg19_ncbiRefSeq) <- "NCBI"
#
#      gene <- c("VWA1")
#      tx   <- c("NM_199121.2")
#
#      data.frame(gene, tx) %>%
#        write.table(
#          x = .,
#          file = "VWA1.preferredtx",
#          quote = F,
#          row.names = F,
#          col.names = F,
#          sep = "\t"
#        )
#
#      bedfile <- GenomicRanges::GRanges(seqnames = "1", ranges = IRanges::IRanges(start=start, width = 1), strand = '+')
#      tx <- getTranscripts(preferred_tx = './VWA1.preferredtx', db = ucsc_hg19_ncbiRefSeq, bedfile = bedfile)
#      tx_names <- unique(unlist(tx$model))
#      three_utr <- GenomicFeatures::threeUTRsByTranscript(x = ucsc_hg19_ncbiRefSeq, use.names=T)[tx_names]
#      hgvs <- pos_three_utr(event_bp=GenomicRanges::start(bedfile), bedfile = bedfile, three_utr = three_utr[[1]])
#      return(hgvs)
#    }
#
#    #in 3'UTR exon
#    # expect *42-461
#    test(start = 1374000)
#
#    #in 3'UTR exon
#    # expect *41+136
#    test(start = 1373000)
#
#
#
#    # NEGATIVE STRAND
#
#    test <- function(start) {
#
#      db <- 'inst/extdata/ucsc_hg19_ncbiRefSeq.sqlite'
#
#      ucsc_hg19_ncbiRefSeq <- AnnotationDbi::loadDb(db) %>%
#        GenomeInfoDb::keepStandardChromosomes(.)
#
#      GenomeInfoDb::seqlevelsStyle(ucsc_hg19_ncbiRefSeq) <- "NCBI"
#
#      gene <- c("BRCA1")
#      tx   <- c("NM_007300.4")
#
#      data.frame(gene, tx) %>%
#        write.table(
#          x = .,
#          file = "BRCA1.preferredtx",
#          quote = F,
#          row.names = F,
#          col.names = F,
#          sep = "\t"
#        )
#
#      bedfile <- GenomicRanges::GRanges(seqnames = "17", ranges = IRanges::IRanges(start=start, width = 1), strand = '-')
#      tx <- getTranscripts(preferred_tx = 'BRCA1.preferredtx', db = ucsc_hg19_ncbiRefSeq, bedfile = bedfile, flank_length = 1000)
#      tx_names <- unique(unlist(tx$model))
#      five_utr  <- GenomicFeatures::fiveUTRsByTranscript(x = ucsc_hg19_ncbiRefSeq, use.names=T)[tx_names]
#      hgvs <- neg_five_utr(bedfile = bedfile, five_utr = five_utr[[1]])
#      return(hgvs)
#    }
#
#    # inside 5'UTR
#    #-7
#    test(start = 41276120) # PASS
#
#    #-32
#    test(start = 41277300) # PASS
#
#    # BETWEEN UTR
#    #-20+288
#    test(start = 41277000) # PASS
#
#    #-19-368
#    test(start = 41276500) # PASS
#
#    # UPSTREAM of UTR
#    #-232
#    test(start = 41277500) # PASS
#
#
#
#
#
#
#    test <- function(start) {
#
#      db <- 'inst/extdata/ucsc_hg19_ncbiRefSeq.sqlite'
#
#      ucsc_hg19_ncbiRefSeq <- AnnotationDbi::loadDb(db) %>%
#        GenomeInfoDb::keepStandardChromosomes(.)
#
#      GenomeInfoDb::seqlevelsStyle(ucsc_hg19_ncbiRefSeq) <- "NCBI"
#
#      gene <- c("BRCA1")
#      tx   <- c("NM_007300.4")
#
#      data.frame(gene, tx) %>%
#        write.table(
#          x = .,
#          file = "BRCA1.preferredtx",
#          quote = F,
#          row.names = F,
#          col.names = F,
#          sep = "\t"
#        )
#
#      bedfile <- GenomicRanges::GRanges(seqnames = "17", ranges = IRanges::IRanges(start=start, width = 1), strand = '-')
#      tx <- getTranscripts(preferred_tx = 'BRCA1.preferredtx', db = ucsc_hg19_ncbiRefSeq, bedfile = bedfile, flank_length = 1000)
#      tx_names <- unique(unlist(tx$model))
#      three_utr  <- GenomicFeatures::threeUTRsByTranscript(x = ucsc_hg19_ncbiRefSeq, use.names=T)[tx_names]
#      hgvs <- neg_three_utr(bedfile = bedfile, three_utr = three_utr[[1]])
#      return(hgvs)
#    }
#
#    # 3UTR NEGATIVE STRAND
#
#    # downstream of 3UTR on -ve strand
#    # expect *1395
#    test(start = 41196300) # PASS
#
#    # in 3UTR
#    # *1295
#    test(start = 41196400) # PASS
#
#    # GENE ON -VE STRAND with 2 x 3UTR EXONS
#    # GNG5
#    test <- function(start) {
#
#      db <- 'inst/extdata/ucsc_hg19_ncbiRefSeq.sqlite'
#
#      ucsc_hg19_ncbiRefSeq <- AnnotationDbi::loadDb(db) %>%
#        GenomeInfoDb::keepStandardChromosomes(.)
#
#      GenomeInfoDb::seqlevelsStyle(ucsc_hg19_ncbiRefSeq) <- "NCBI"
#
#      gene <- c("GNG5")
#      tx   <- c("NM_005274.3")
#
#      data.frame(gene, tx) %>%
#        write.table(
#          x = .,
#          file = "GNG5.preferredtx",
#          quote = F,
#          row.names = F,
#          col.names = F,
#          sep = "\t"
#        )
#
#      bedfile <- GenomicRanges::GRanges(seqnames = "1", ranges = IRanges::IRanges(start=start, width = 1), strand = '-')
#      tx <- getTranscripts(preferred_tx = 'GNG5.preferredtx', db = ucsc_hg19_ncbiRefSeq, bedfile = bedfile, flank_length = 1000)
#      tx_names <- unique(unlist(tx$model))
#      three_utr  <- GenomicFeatures::threeUTRsByTranscript(x = ucsc_hg19_ncbiRefSeq, use.names=T)[tx_names]
#      hgvs <- neg_three_utr(bedfile = bedfile, three_utr = three_utr[[1]])
#      return(hgvs)
#    }
#
#    # Upstream og UTR
#    #*251
#    test(start = 84964000) # PASS
#
#    #*20-769
#    test(start = 84965000)
#
#    #*19+509
#    test(start = 84967000)
#
#
#  }
#
#
