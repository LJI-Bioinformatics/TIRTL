test_that("get_well_subset works", {
    expected_wellset1 <- readRDS(test_path('testdata/expected_wellset1.rds'))
    expect_equal(get_well_subset(1:16,1:24), expected_wellset1)
})
