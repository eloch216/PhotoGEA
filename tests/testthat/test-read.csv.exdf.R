test_that('data agrees with original exdf object', {
    exdf_from_excel <- read_gasex_file(
        system.file('extdata', 'ball_berry_1.xlsx', package = 'PhotoGEA', mustWork = TRUE)
    )

    exdf_from_csv <- read.csv.exdf(
        system.file('extdata', 'ball_berry_1.csv', package = 'PhotoGEA', mustWork = TRUE)
    )

    expect_equal(colnames(exdf_from_csv), colnames(exdf_from_excel))
    expect_equal(exdf_from_csv$units, exdf_from_excel$units)
    expect_equal(exdf_from_csv$categories, exdf_from_excel$categories)
    expect_equal(exdf_from_csv[, 'A'], exdf_from_excel[, 'A'])
})
