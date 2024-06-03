test_that('plaintext Licor files that have been reopened are read properly', {
    v1 <- expect_silent(read_gasex_file(PhotoGEA_example_file_path('plaintext_licor_file')))
    v2 <- expect_silent(read_gasex_file(PhotoGEA_example_file_path('plaintext_licor_file_v2')))

    expect_equal(colnames(v1), colnames(v2))
    expect_equal(v1$units, v2$units)
    expect_equal(v1$categories, v2$categories)
    expect_equal(v1[, 'A'], v2[, 'A'])
    expect_equal(v1[, 'oxygen'], v2[, 'oxygen'])
})
