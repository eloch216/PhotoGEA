test_that('Ball-Berry index is calculated for exdf and data frame objects', {
    # Load an example file and calculate additional properties
    licor_file <- read_gasex_file(
        PhotoGEA_example_file_path('ball_berry_1.xlsx')
    )

    licor_file <- calculate_total_pressure(licor_file)

    licor_file <- calculate_gas_properties(licor_file)

    # Calculate BB index from exdf
    bb_exdf <- expect_silent(
        calculate_ball_berry_index(licor_file)
    )

    # Calculate BB index from data frame
    bb_df <- expect_silent(
        calculate_ball_berry_index(licor_file$main_data)
    )

    # Check to make sure results have not changed
    expect_equal(
        as.numeric(bb_exdf[seq_len(4), 'bb_index']),
        c(0.064874428, 0.042172990, 0.047174397, 0.042301287)
    )

    # Check to make sure results match
    expect_equal(
        as.numeric(bb_exdf[, 'bb_index']),
        as.numeric(bb_df[, 'bb_index'])
    )
})
