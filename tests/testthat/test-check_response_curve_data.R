# Read an example Licor file included in the PhotoGEA package
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('c3_aci_1.xlsx')
)

# Define a new column that uniquely identifies each curve
licor_file[, 'curve_identifier'] <-
  paste(licor_file[, 'species'], licor_file[, 'plot'], sep = ' - ')

test_that('check_response_curve_data produces messages only when expected', {
    # Open a tempfile to redirect the printed output
    sink(tempfile())

    # Run tests
    expect_silent(
        check_response_curve_data(licor_file, 'curve_identifier', 16, 'CO2_r_sp')
    )

    expect_error(
        check_response_curve_data(licor_file, 'curve_identifier', 15, 'CO2_r_sp'),
        'One or more curves does not have the expected number of points.'
    )

    expect_error(
        check_response_curve_data(licor_file, 'curve_identifier', 16, 'Ci'),
        'The curves do not all follow the same sequence of the driving variable.'
    )

    expect_warning(
        check_response_curve_data(licor_file, 'curve_identifier', 15, 'CO2_r_sp', error_on_failure = FALSE),
        'One or more curves does not have the expected number of points.'
    )

    expect_warning(
        check_response_curve_data(licor_file, 'curve_identifier', 16, 'Ci', error_on_failure = FALSE),
        'The curves do not all follow the same sequence of the driving variable.'
    )

    # Close the tempfile
    sink()
})
