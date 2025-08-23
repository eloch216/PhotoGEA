test_that('Excel Licor files with user remarks are read properly', {
    licor_file <- expect_silent(read_gasex_file(PhotoGEA_example_file_path('c3_aci_1.xlsx')))

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

test_that('Excel Licor files are not required to have a Remarks sheet', {
    licor_file <- expect_silent(read_gasex_file(PhotoGEA_example_file_path('c3_aci_2.xlsx')))

    expect_true('user_remarks' %in% names(licor_file))

    expect_equal(
        nrow(licor_file$user_remarks),
        0
    )

    expect_true(
        all(is.na(licor_file[, 'user_remark']))
    )
})

test_that('Changes to oxygen percentage are properly read', {
    licor_file <- expect_silent(read_gasex_file(PhotoGEA_example_file_path('c3_aci_2.xlsx')))

    expect_true('Oxygen' %in% colnames(licor_file))

    expect_equal(
        licor_file[, 'Oxygen'],
        c(rep_len(21, 10), rep_len(22, 38))
    )
})
