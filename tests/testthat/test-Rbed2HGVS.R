test_that("cds gap that finishes downstream on positive strand", {

  out <- Rbed2HGVS(
    bedfile = '../../inst/extdata/positive-strand-downstream.bed',
    db = '../../inst/extdata/ucsc_hg19_ncbiRefSeq.sqlite'
    )

  testthat::expect_equal(out$hgvs_start %>% as.vector(), "6829")
  testthat::expect_equal(out$exon_start, 11)
  testthat::expect_equal(out$hgvs_end %>% as.vector(), "6841+19")
  testthat::expect_equal(out$exon_end, 11)
})

test_that("cds gap that finishes downstream on negative strand", {

  out <- Rbed2HGVS(
    bedfile = '../../inst/extdata/negative-strand-downstream.bed',
    db = '../../inst/extdata/ucsc_hg19_ncbiRefSeq.sqlite'
  )

  testthat::expect_equal(
    out$hgvs_start %>% as.vector(),
    c("623+14",  "704+14",  "983+14",  "983+14",  "1100+14", "1100+14", "983+14")
    )

  testthat::expect_equal(out$exon_start, c(6, 6, 10, 10, 10, 10, 9))

  testthat::expect_equal(
    out$hgvs_end %>% as.vector(),
    c("610", "691", "970", "970", "1087", "1087", "970")
    )

  testthat::expect_equal(out$exon_end, c(6, 6, 10, 10, 10, 10, 9))
})


test_that("cds gap that finishes upstream on positive strand", {

  out <- Rbed2HGVS(
    bedfile = '../../inst/extdata/positive-strand-upstream.bed',
    db = '../../inst/extdata/ucsc_hg19_ncbiRefSeq.sqlite'
  )

  testthat::expect_equal(out$hgvs_start %>% as.vector(), "1910-30")
  testthat::expect_equal(out$exon_start, 11)
  testthat::expect_equal(out$hgvs_end %>% as.vector(), "1939")
  testthat::expect_equal(out$exon_end, 11)
})


test_that("cds gap that finishes upstream on negative strand", {

  out <- Rbed2HGVS(
    bedfile = '../../inst/extdata/negative-strand-upstream.bed',
    db = '../../inst/extdata/ucsc_hg19_ncbiRefSeq.sqlite'
  )

  testthat::expect_equal(
    out$hgvs_start %>% as.vector(),
    c("532", "613", "892", "892", "1009", "1009", "892")
    )

  testthat::expect_equal(out$exon_start, c(6, 6, 10, 10, 10, 10, 9))

  testthat::expect_equal(
    out$hgvs_end %>% as.vector(),
    c("517-19", "598-19", "877-19", "877-19", "994-19", "994-19", "877-19")
    )

  testthat::expect_equal(out$exon_end, c(6, 6, 10, 10, 10, 10, 9))
})


test_that("cds gap that encompasses whole exon", {

  out <- Rbed2HGVS(
    bedfile = '../../inst/extdata/whole-exon-gap.bed',
    db = '../../inst/extdata/ucsc_hg19_ncbiRefSeq.sqlite'
  )

  testthat::expect_equal(
    out$hgvs_start %>% as.vector(),
    c("623+11", "704+11", "983+11", "983+11", "1100+11", "1100+11", "983+11"))

  testthat::expect_equal(out$exon_start, c(6, 6, 10, 10, 10, 10, 9))

  testthat::expect_equal(
    out$hgvs_end %>% as.vector(),
    c("517-16", "598-16", "877-16", "877-16", "994-16", "994-16", "877-16"))

  testthat::expect_equal(out$exon_end, c(6, 6, 10, 10, 10, 10, 9))
})

