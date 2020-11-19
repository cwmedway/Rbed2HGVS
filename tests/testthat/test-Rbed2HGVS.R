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


# ----- Inside CDS -----
testthat::test_that("position inside cds exon on positive strand", {

   out <- run_test(chr = 13, bp = 32915300, gene = "BRCA2", tx = "NM_000059.3")
   check_test(test = out, exp_hgvs = "6808", exp_exon = "11")
})

testthat::test_that("position inside cds exon on negative strand", {

   out <- run_test(chr = 17, bp = 7577550, gene = "TP53", tx = "NM_000546.5")
   check_test(test = out, exp_hgvs = '731', exp_exon = "7")
})


# ----- Upstream / Downstream of CDS Exons -----
testthat::test_that("position downstream of cds exon on positive strand", {

   out <- run_test(chr = 13, bp = 32915400, gene = "BRCA2", tx = "NM_000059.3")
   check_test(test = out, exp_hgvs = '6841+67', exp_exon = "11")
})

testthat::test_that("position downstream of cds exon on negative strand", {

   out <- run_test(chr = 17, bp = 7573912, gene = "TP53", tx = "NM_000546.5")
   check_test(test = out, exp_hgvs = '1100+15', exp_exon = "10")
})

testthat::test_that("cds gap that finishes upstream on positive strand", {

   out <- run_test(chr = 13, bp = 32910371, gene = 'BRCA2', tx = 'NM_000059.3')
   check_test(test = out, exp_hgvs = '1910-31', exp_exon = "11")
})

testthat::test_that("cds gap that finishes upstream on negative strand", {

   out <- run_test(chr = 17, bp = 7574052, gene = 'TP53', tx = 'NM_000546.5')
   check_test(test = out, exp_hgvs = '994-19', exp_exon = "10")
})


# ----- Upstream / Downstream of UTRs -----
# ----- Positive Strand -----
# ----- 5'UTR -----
testthat::test_that("inside 5' UTR exon positive strand", {

   out <- run_test(chr = 13, bp = 32889700, gene = 'BRCA2', tx = 'NM_000059.3')
   check_test(test = out, exp_hgvs = '-144', exp_exon = "2")
})

testthat::test_that("upstream of 1st 5' UTR exon positive strand", {

   out <- run_test(chr = 13, bp = 32889570, gene = 'BRCA2', tx = 'NM_000059.3')
   check_test(test = out, exp_hgvs = '-274', exp_exon = "2")
})

testthat::test_that("upstream of 5' UTR exon positive strand", {

   out <- run_test(chr = 13, bp = 32889850, gene = 'BRCA2', tx = 'NM_000059.3')
   check_test(test = out, exp_hgvs = '-40+46', exp_exon = "2")
})

testthat::test_that("downstream of 5' UTR positive strand", {

   out <- run_test(chr = 13, bp = 32890500, gene = 'BRCA2', tx = 'NM_000059.3')
   check_test(test = out, exp_hgvs = '-39-59', exp_exon = "2")
})

# ----- 3'UTR -----
testthat::test_that("inside 3'UTR exon positive strand", {

   out <- run_test(chr = 1, bp = 1375000, gene = 'VWA1', tx = 'NM_199121.2')
   check_test(test = out, exp_hgvs = '*581', exp_exon = "2")
})

testthat::test_that("downstream of final 3' UTR exon positive strand", {

   out <- run_test(chr = 1, bp = 1378300, gene = 'VWA1', tx = 'NM_199121.2')
   check_test(test = out, exp_hgvs = '*3881', exp_exon = "2")
})

testthat::test_that("upstream of 3' UTR exon positive strand", {

   out <- run_test(chr = 1, bp = 1374000, gene = 'VWA1', tx = 'NM_199121.2')
   check_test(test = out, exp_hgvs = '*42-461', exp_exon = "2")
})

testthat::test_that("downstream of 3' UTR exon positive strand", {

   out <- run_test(chr = 1, bp = 1373000, gene = 'VWA1', tx = 'NM_199121.2')
   check_test(test = out, exp_hgvs = '*41+136', exp_exon = "2")
})

# ----- Negative Strand -----
# ----- 5'UTR -----
testthat::test_that("within 5' UTR exon nagative strand", {

   out <- run_test(chr = 17, bp = 41276120, gene = 'BRCA1', tx = 'NM_007300.4')
   check_test(test = out, exp_hgvs = '-7', exp_exon = "2")
})

testthat::test_that("upstream of 1st 5' UTR exon negative strand", {

   out <- run_test(chr = 17, bp = 41278000, gene = 'BRCA1', tx = 'NM_007300.4')
   check_test(test = out, exp_hgvs = '-732', exp_exon = "2")
})

testthat::test_that("upstream of 5' UTR exon negative strand", {

   out <- run_test(chr = 17, bp = 41276500, gene = 'BRCA1', tx = 'NM_007300.4')
   check_test(test = out, exp_hgvs = '-19-368', exp_exon = "2")
})

testthat::test_that("downstream of 5' UTR exon negative strand", {

   out <- run_test(chr = 17, bp = 41277000, gene = 'BRCA1', tx = 'NM_007300.4')
   check_test(test = out, exp_hgvs = '-20+288', exp_exon = "2")
})

# ----- 3'UTR -----
testthat::test_that("within 3' UTR exon nagative strand", {

   out <- run_test(chr = 1, bp = 84964020, gene = 'GNG5', tx = 'NM_005274.3')
   check_test(test = out, exp_hgvs = '*212', exp_exon = "3")
})

testthat::test_that("downstream of last 3' UTR exon negative strand", {

   out <- run_test(chr = 1, bp = 84964000, gene = 'GNG5', tx = 'NM_005274.3')
   check_test(test = out, exp_hgvs = '*251', exp_exon = "3")
})

testthat::test_that("upstream of 3' UTR exon negative strand", {

   out <- run_test(chr = 1, bp = 84965000, gene = 'GNG5', tx = 'NM_005274.3')
   check_test(test = out, exp_hgvs = '*20-769', exp_exon = "3")
})

testthat::test_that("downstream of 3' UTR exon negative strand", {

   out <- run_test(chr = 1, bp = 84967000, gene = 'GNG5', tx = 'NM_005274.3')
   check_test(test = out, exp_hgvs = '*19+509', exp_exon = "3")
})

# ----- Edge Cases -----
testthat::test_that("no interval near gene", {

   out <- run_test(chr = 1, bp = 86700000)
   testthat::expect_equal(object = out$hgvs$hgvs_start, expected = NA)
   testthat::expect_equal(object = out$hgvs$exon_start, expected = NA)
})


testthat::test_that("no preferred tx found for bed entry", {

   out <- run_test(chr = 17, bp = 41277000, gene = "BRCA2", 'NM_000059.3')
   testthat::expect_equal(object = out$hgvs$hgvs_start, expected = NA)
   testthat::expect_equal(object = out$hgvs$exon_start, expected = NA)
})




