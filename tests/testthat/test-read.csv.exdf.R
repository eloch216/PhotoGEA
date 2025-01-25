test_that('data agrees with original exdf object', {
    exdf_from_excel <- read_gasex_file(
        PhotoGEA_example_file_path('ball_berry_1.xlsx')
    )

    exdf_from_csv <- read.csv.exdf(
        PhotoGEA_example_file_path('ball_berry_1.csv')
    )

    expect_equal(colnames(exdf_from_csv),        colnames(exdf_from_excel))
    expect_equal(exdf_from_csv$units,            exdf_from_excel$units)
    expect_equal(exdf_from_csv$categories,       exdf_from_excel$categories)
    expect_equal(exdf_from_csv[, 'A'],           exdf_from_excel[, 'A'])
    expect_equal(exdf_from_csv[, 'user_remark'], exdf_from_excel[, 'user_remark'])
})
