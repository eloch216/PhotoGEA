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
