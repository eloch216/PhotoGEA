# Read an example Licor file included in the PhotoGEA package
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('c3_aci_1.xlsx')
)

# Define a new column that uniquely identifies each curve
licor_file[, 'species_plot'] <-
  paste(licor_file[, 'species'], '-', licor_file[, 'plot'] )

# Organize the data
licor_file <- organize_response_curve_data(
    licor_file,
    'species_plot',
    c(9, 10, 16),
    'CO2_r_sp'
)

# Apply the Laisk method. Note: this is a bad example because these curves were
# measured at the same light intensity, but from different species. Because of
# this, the results are not meaningful.
laisk_results <- fit_laisk(
  licor_file, 20, 150,
  ppfd_column_name = 'species_plot'
)

# Choose test tolerance
TOLERANCE <- 1e-4

test_that('fit results have not changed', {
    expect_equal(
        as.numeric(laisk_results$first_fit_parameters[1, c('laisk_intercept', 'laisk_intercept_err', 'laisk_slope', 'laisk_slope_err', 'r_squared', 'p_value')]),
        c(-1.074544e+01, 5.639684e-01, 1.828644e-01, 6.673785e-03, 9.947005e-01, 1.055060e-05),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(laisk_results$second_fit_parameters[1, c('RL', 'RL_err', 'Ci_star', 'Ci_star_err', 'r_squared', 'p_value')]),
        c(-1.61745253 , 1.84213648, 69.37656570, 9.89203385, 0.98007470, 0.09016446),
        tolerance = TOLERANCE
    )
})
