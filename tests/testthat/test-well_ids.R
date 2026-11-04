fake_clones <- function(wells, n_alpha = 10, n_beta = 10) {
  mk <- function(n) data.table::data.table(readCount = rep(1, n), readFraction = rep(1 / n, n))
  a <- lapply(wells, function(w) mk(n_alpha))
  b <- lapply(wells, function(w) mk(n_beta))
  names(a) <- paste0(wells, "_TCRa_ID10_clones_TRA_tsv")
  names(b) <- paste0(wells, "_TCRb_ID10_clones_TRB_tsv")
  list(a = a, b = b)
}

test_that("normalize_well_id strips leading zeros and upper-cases the row", {
  expect_equal(normalize_well_id(c("A01", "A1", "a01", "P24", "H009", "B10")),
               c("A1", "A1", "A1", "P24", "H9", "B10"))
  expect_equal(normalize_well_id(c("TCRa", "ID10", "A1B")), c("TCRa", "ID10", "A1B"))
})

test_that("get_good_wells_sub matches zero-padded file names to get_well_subset()", {
  padded <- sprintf("%s%02d", rep(LETTERS[1:2], each = 12), rep(1:12, times = 2))
  cl <- fake_clones(padded)
  qc <- TIRTL:::get_good_wells_sub(cl$a, cl$b, thres = 5, pos = 1,
                                   wellset = get_well_subset(1:2, 1:12))
  expect_true(all(qc$a))
  expect_true(all(qc$b))
  expect_setequal(qc$well_ids, get_well_subset(1:2, 1:12))
})

test_that("padded and unpadded file names select the same wells", {
  cols <- c(1:3, 10:12)
  padded <- sprintf("A%02d", cols)
  unpadded <- sprintf("A%d", cols)
  ws <- get_well_subset(1, c(1, 2, 10))
  qp <- TIRTL:::get_good_wells_sub(fake_clones(padded)$a, fake_clones(padded)$b, 5, pos = 1, wellset = ws)
  qu <- TIRTL:::get_good_wells_sub(fake_clones(unpadded)$a, fake_clones(unpadded)$b, 5, pos = 1, wellset = ws)
  expect_equal(qp$a, qu$a)
  expect_equal(qp$b, qu$b)
  expect_equal(qp$well_ids, qu$well_ids)
  expect_equal(qp$a, c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE))
})

test_that("a padded wellset matches unpadded file names", {
  cl <- fake_clones(c("A1", "A2", "A10"))
  qc <- TIRTL:::get_good_wells_sub(cl$a, cl$b, 5, pos = 1, wellset = c("A01", "A02"))
  expect_equal(qc$a, c(TRUE, TRUE, FALSE))
})

test_that("the clone threshold still filters wells after normalization", {
  cl <- fake_clones(c("A01", "A02"))
  cl$a[[2]] <- cl$a[[2]][1:2]
  qc <- TIRTL:::get_good_wells_sub(cl$a, cl$b, 5, pos = 1, wellset = get_well_subset(1, 1:2))
  expect_equal(qc$well_ids, "A1")
  expect_equal(qc$a, c(TRUE, FALSE))
  expect_equal(qc$b, c(TRUE, FALSE))
})

test_that("A1 and A01 in the same folder is an error", {
  cl <- fake_clones(c("A1", "A01"))
  expect_error(TIRTL:::get_good_wells_sub(cl$a, cl$b, 5, pos = 1),
               "map to the same well.*A1/A01")
})

test_that("run_single_point_analysis_sub_gpu stops when no requested wells are present", {
  dir <- tempfile("tirtl_wells_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  for (w in c("A01", "A02")) {
    for (ch in c("TRA", "TRB")) {
      data.table::fwrite(data.table::data.table(readCount = 1, readFraction = 1),
                         file.path(dir, sprintf("%s.TCR_ID10.clones_%s.tsv", w, ch)), sep = "\t")
    }
  }
  expect_error(
    capture.output(run_single_point_analysis_sub_gpu(dir, outdir = dir, well_pos = 1,
                                                      wellset1 = get_well_subset(2, 1:2),
                                                      compute = FALSE)),
    "None of the requested wells were found"
  )
})
