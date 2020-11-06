test_that("cds gap that finishes downstream on positive strand", {

  # make test bed
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

  unlink(c('test-neg-us.bed', 'test-neg-us.preferredtx'))
})

# ---------------------------------------------------------
# Issue 1: gap falls outside transcript (promotor i.e TERT)

 test_that("cds gap that finishes upstream on negative strand", {

   chr   <- c(5)
   start <- c(1295226, 1295248)
   end   <- c(1295228, 1295250)
   meta  <- c("TERT(NM_198253.2):c.-124_-124", "TERT(NM_198253.2):c.-146_-146")

   data.frame(chr, start, end, meta) %>%
     write.table(
       x = .,
       file = "issue1.bed",
       quote = F,
       row.names = F,
       col.names = F,
       sep = "\t"
     )

   data.frame("TERT", "NM_198253.2") %>%
     write.table(
       x = .,
       file = "issue1.preferredtx",
       quote = F,
       row.names = F,
        col.names = F,
       sep = "\t"
     )

   out <- Rbed2HGVS(
     bedfile = 'issue1.bed',
     db = '../../inst/extdata/ucsc_hg19_ncbiRefSeq.sqlite',
     preferred_tx = 'issue1.preferredtx'
   )

   testthat::expect_equal(
     out$hgvs$hgvs_end %>% as.vector(),
     c("1-123","1-145")
   )

   testthat::expect_equal(out$hgvs$exon_start, c(1,1))

   testthat::expect_equal(
     out$hgvs$hgvs_start %>% as.vector(),
     c("1-124", "1-146")
   )

   testthat::expect_equal(out$hgvs$exon_end, c(1,1))

   unlink(c('issue1.bed', 'issue1.preferredtx'))
 })

# # ----------------------------------------------------
# # ISSUE 2 - gap between UTR exons

 chr   <- c("13")
 start <- c(32889850)
 end   <- c(32890500)
 meta  <- c(
   "BRCA2(NM_000059.3):c.-40+46_-39-59"
  )

data.frame(chr, start, end, meta) %>%
  write.table(
    x = .,
    file = "issue2.bed",
    quote = F,
    row.names = F,
    col.names = F,
    sep = "\t"
  )

gene <- c("BRCA2")
tx   <- c("NM_000059.3")

data.frame(gene, tx) %>%
  write.table(
    x = .,
    file = "issue2.preferredtx",
    quote = F,
    row.names = F,
    col.names = F,
    sep = "\t"
  )

out <- Rbed2HGVS(
  bedfile = 'issue2.bed',
  db = '../../inst/extdata/ucsc_hg19_ncbiRefSeq.sqlite',
  preferred_tx = 'issue2.preferredtx'
)

 testthat::expect_equal(
   out$hgvs$hgvs_start %>% as.vector(),
   c("-40+46")
 )

 testthat::expect_equal(
   out$hgvs$hgvs_end %>% as.vector(),
   c("-39-59")
 )


 unlink(c("issue2.bed","issue2.preferredtx"))

# # ----------------------
# Upstream of 5UTR to inside UTR
chr   <- c("13")
start <- c(32889570)
end   <- c(32890570)
meta  <- c(
   "BRCA2(NM_000059.3):c.-274_-28"
)

 data.frame(chr, start, end, meta) %>%
   write.table(
     x = .,
     file = "issue3.bed",
     quote = F,
     row.names = F,
     col.names = F,
     sep = "\t"
   )

 gene <- c("BRCA2")
 tx   <- c("NM_000059.3")

 data.frame(gene, tx) %>%
   write.table(
     x = .,
     file = "issue3.preferredtx",
     quote = F,
     row.names = F,
     col.names = F,
     sep = "\t"
   )

 out <- Rbed2HGVS(
   bedfile = 'issue3.bed',
   db = '../../inst/extdata/ucsc_hg19_ncbiRefSeq.sqlite',
   preferred_tx = 'issue3.preferredtx'
 )

 testthat::expect_equal(
   out$hgvs$hgvs_start %>% as.vector(),
   c("-274")
 )

 testthat::expect_equal(
   out$hgvs$hgvs_end %>% as.vector(),
   c("-28")
 )

 unlink(c("issue3.bed","issue3.preferredtx"))

# ----------------------------
#5' UTR on negative strand
 chr   <- c("17")
 start <- c(41276000, 41276050, 41277000)
 end   <- c(41276120, 41276500, 41278000)
 meta  <- c(
   "BRCA1(NM_007300.4):c.80+34_-7",
   "BRCA1(NM_007300.4):c.64_-19-368",
   "BRCA1(NM_007300.4):c.-20+288_-732"
 )

 data.frame(chr, start, end, meta) %>%
   write.table(
     x = .,
     file = "issue4.bed",
     quote = F,
     row.names = F,
     col.names = F,
     sep = "\t"
   )

 gene <- c("BRCA1")
 tx   <- c("NM_007300.4")

 data.frame(gene, tx) %>%
   write.table(
     x = .,
     file = "issue4.preferredtx",
     quote = F,
     row.names = F,
     col.names = F,
     sep = "\t"
   )

 out <- Rbed2HGVS(
   bedfile = 'issue4.bed',
   db = '../../inst/extdata/ucsc_hg19_ncbiRefSeq.sqlite',
   preferred_tx = 'issue4.preferredtx'
 )



# # -----------------
# # 3'UTR
# # Positive Strand
#
# chr   <- c("1")
# start <- c(1373000)
# end   <- c(1375000)
# meta  <- c(
#   "VWA1(NM_199121.2):c.*41+136_*581"
# )
#
# data.frame(chr, start, end, meta) %>%
#   write.table(
#     x = .,
#     file = "issue4.bed",
#     quote = F,
#     row.names = F,
#     col.names = F,
#     sep = "\t"
#   )
#
# gene <- c("VWA1")
# tx   <- c("NM_199121.2")
#
# data.frame(gene, tx) %>%
#   write.table(
#     x = .,
#     file = "issue4.preferredtx",
#     quote = F,
#     row.names = F,
#     col.names = F,
#     sep = "\t"
#   )
#
# out <- Rbed2HGVS(
#   bedfile = 'issue4.bed',
#   db = '../../inst/extdata/ucsc_hg19_ncbiRefSeq.sqlite',
#   preferred_tx = 'issue4.preferredtx'
# )

