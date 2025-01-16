test_that('plaintext Licor files that have been reopened are read properly', {
    v1 <- expect_silent(read_gasex_file(PhotoGEA_example_file_path('plaintext_licor_file')))
    v2 <- expect_silent(read_gasex_file(PhotoGEA_example_file_path('plaintext_licor_file_v2')))

    expect_equal(nrow(v2), 96)

    expect_equal(colnames(v2), colnames(v1))
    expect_equal(v2$units, v1$units)
    expect_equal(v2$categories, v1$categories)
    expect_equal(v2[, 'A'], v1[, 'A'])
    expect_equal(v2[, 'oxygen'], v1[, 'oxygen'])
})
