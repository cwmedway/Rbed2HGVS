test_that("cds gap that finishes downstream on positive strand", {

  out <- Rbed2HGVS(
    bedfile = '/home/chris/Apps/Rbed2HGVS/data/positive-strand-downstream.bed',
    db = '/home/chris/Apps/Rbed2HGVS/data/ucsc_hg19_ncbiRefSeq.sqlite'
    )

  testthat::expect_equal(out$hgvs.start %>% as.vector(), "c.6829")
  testthat::expect_equal(out$exon.start, 11)
  testthat::expect_equal(out$hgvs.end %>% as.vector(), "c.6841+18")
  testthat::expect_equal(out$exon.end, 11)
})

test_that("cds gap that finishes downstream on negative strand", {

  out <- Rbed2HGVS(
    bedfile = '/home/chris/Apps/Rbed2HGVS/data/negative-strand-downstream.bed',
    db = '/home/chris/Apps/Rbed2HGVS/data/ucsc_hg19_ncbiRefSeq.sqlite'
  )

  testthat::expect_equal(
    out$hgvs.start %>% as.vector(),
    c("c.1100+14","c.1100+14","c.704+14","c.983+14","c.623+14","c.983+14","c.983+14")
    )

  testthat::expect_equal(out$exon.start, c(10,10,6,9,6,10,10))

  testthat::expect_equal(
    out$hgvs.end %>% as.vector(),
    c("c.1087","c.1087","c.691","c.970","c.610","c.970","c.970")
    )

  testthat::expect_equal(out$exon.end, c(10,10,6,9,6,10,10))
})


test_that("cds gap that finishes upstream on positive strand", {

  out <- Rbed2HGVS(
    bedfile = '/home/chris/Apps/Rbed2HGVS/data/positive-strand-upstream.bed',
    db = '/home/chris/Apps/Rbed2HGVS/data/ucsc_hg19_ncbiRefSeq.sqlite'
  )

  testthat::expect_equal(out$hgvs.start %>% as.vector(), "c.1910-30")
  testthat::expect_equal(out$exon.start, 11)
  testthat::expect_equal(out$hgvs.end %>% as.vector(), "c.1939")
  testthat::expect_equal(out$exon.end, 11)
})


test_that("cds gap that finishes upstream on negative strand", {

  out <- Rbed2HGVS(
    bedfile = '/home/chris/Apps/Rbed2HGVS/data/negative-strand-upstream.bed',
    db = '/home/chris/Apps/Rbed2HGVS/data/ucsc_hg19_ncbiRefSeq.sqlite'
  )

  testthat::expect_equal(
    out$hgvs.start %>% as.vector(),
    c("c.1009","c.1009","c.613","c.892","c.532","c.892","c.892")
    )

  testthat::expect_equal(out$exon.start, c(10, 10, 6, 9, 6, 10, 10))

  testthat::expect_equal(
    out$hgvs.end %>% as.vector(),
    c("c.994-19", "c.994-19", "c.598-19", "c.877-19", "c.517-19", "c.877-19", "c.877-19")
    )

  testthat::expect_equal(out$exon.end, c(10, 10, 6, 9, 6, 10, 10))
})

