test_that("cds gap that finishes downstream on positive strand", {

  # make test bed
  # 13	32915320	32915352
  data.frame("13", "32915320", "32915352") %>%
    write.table(
      x = .,
      file = "test-pos-ds.bed",
      quote = F,
      row.names = F,
      col.names = F, sep = "\t"
      )

  out <- Rbed2HGVS::Rbed2HGVS(
    bedfile = 'test-pos-ds.bed',
    db = '../../inst/extdata/ucsc_hg19_ncbiRefSeq.sqlite'
    )

  testthat::expect_equal(out$hgvs$hgvs_start %>% as.vector(), "6829")
  testthat::expect_equal(out$hgvs$exon_start, 11)
  testthat::expect_equal(out$hgvs$hgvs_end %>% as.vector(), "6841+19")
  testthat::expect_equal(out$hgvs$exon_end, 11)

  unlink("test-pos-ds.bed")
})



test_that("cds gap that finishes downstream on negative strand", {

  data.frame("17",	"7573912",	"7573940") %>%
    write.table(
      x = .,
      file = "test-neg-ds.bed",
      quote = F,
      row.names = F,
      col.names = F,
      sep = "\t"
    )

  data.frame("TP53", "NM_000546.5") %>%
    write.table(
      x = .,
      file = "test-neg-ds.preferredtx",
      quote = F,
      row.names = F,
      col.names = F,
      sep = "\t"
    )

  out <- Rbed2HGVS::Rbed2HGVS(
    bedfile = 'test-neg-ds.bed',
    db = '../../inst/extdata/ucsc_hg19_ncbiRefSeq.sqlite',
    preferred_tx = 'test-neg-ds.preferredtx'
  )

  testthat::expect_equal(
    out$hgvs$hgvs_end %>% as.vector(),
    c("1100+14")
    )

  testthat::expect_equal(out$hgvs$exon_start, c(10))

  testthat::expect_equal(
    out$hgvs$hgvs_start %>% as.vector(),
    c("1087")
    )

  testthat::expect_equal(out$hgvs$exon_end, c(10))

  unlink(c('test-neg-ds.bed', 'test-neg-ds.preferredtx'))
})



test_that("cds gap that finishes upstream on positive strand", {

  data.frame("13", "32910371", "32910431") %>%
    write.table(
      x = .,
      file = "test-pos-us.bed",
      quote = F,
      row.names = F,
      col.names = F, sep = "\t"
    )

  out <- Rbed2HGVS(
    bedfile = 'test-pos-us.bed',
    db = '../../inst/extdata/ucsc_hg19_ncbiRefSeq.sqlite'
  )

  testthat::expect_equal(out$hgvs$hgvs_start %>% as.vector(), "1910-30")
  testthat::expect_equal(out$hgvs$exon_start, 11)
  testthat::expect_equal(out$hgvs$hgvs_end %>% as.vector(), "1939")
  testthat::expect_equal(out$hgvs$exon_end, 11)

  unlink('test-pos-us.bed')
})


test_that("cds gap that finishes upstream on negative strand", {

  data.frame("17",	"7574017",	"7574052") %>%
    write.table(
      x = .,
      file = "test-neg-us.bed",
      quote = F,
      row.names = F,
      col.names = F,
      sep = "\t"
    )

  data.frame("TP53", "NM_000546.5") %>%
    write.table(
      x = .,
      file = "test-neg-us.preferredtx",
      quote = F,
      row.names = F,
      col.names = F,
      sep = "\t"
    )

  out <- Rbed2HGVS(
    bedfile = 'test-neg-us.bed',
    db = '../../inst/extdata/ucsc_hg19_ncbiRefSeq.sqlite',
    preferred_tx = 'test-neg-us.preferredtx'
  )

  testthat::expect_equal(
    out$hgvs$hgvs_end %>% as.vector(),
    c("1009")
    )

  testthat::expect_equal(out$hgvs$exon_start, c(10))

  testthat::expect_equal(
    out$hgvs$hgvs_start %>% as.vector(),
    c("994-19")
    )

  testthat::expect_equal(out$hgvs$exon_end, c(10))
})



