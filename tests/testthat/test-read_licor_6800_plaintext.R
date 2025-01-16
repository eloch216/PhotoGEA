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

test_that('plaintext Licor files with user remarks are read properly', {
    licor_file <- expect_silent(read_gasex_file(PhotoGEA_example_file_path('plaintext_licor_file_v2')))

    expect_true('user_remarks' %in% names(licor_file))

    expect_equal(
        licor_file$user_remarks$remark_time,
        c('09:22:24', '11:04:20', '13:10:42')
    )

    expect_equal(
        licor_file$user_remarks$remark_value,
        c('Stability Definition: gsw (GasEx): Slp<0.2 Std<0.02 Per=30 A (GasEx): Slp<1 Std<0.2 Per=30', 'a user remark', 'another user remark')
    )
})
