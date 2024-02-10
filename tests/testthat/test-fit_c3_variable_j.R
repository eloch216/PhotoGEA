# Read an example Licor file included in the PhotoGEA package
licor_file <- read_gasex_file(
  system.file('extdata', 'c3_aci_1.xlsx', package = 'PhotoGEA', mustWork = TRUE)
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

# Calculate the total pressure in the Licor chamber
licor_file <- calculate_total_pressure(licor_file)

# Calculate temperature-dependent values of C3 photosynthetic parameters
licor_file <- calculate_arrhenius(licor_file, c3_arrhenius_bernacchi)

# Get just one curve
one_curve <- licor_file[licor_file[, 'species_plot'] == 'tobacco - 1', , TRUE]

# Purposely introduce negative Ci values
one_curve_bad <- one_curve
one_curve_bad[one_curve_bad[, 'CO2_r_sp'] == 20, 'Ci'] <- -5

test_that('fit failures are handled properly', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res_bad <- expect_silent(
        fit_c3_variable_j(one_curve_bad, Ca_atmospheric = 420)
    )

    expect_equal(unique(fit_res_bad$fits[, 'c3_assimilation_msg']), 'Cc must be >= 0')
    expect_equal(unique(fit_res_bad$fits[, 'c3_variable_j_msg']), 'Ci must be >= 0')
    expect_equal(fit_res_bad$parameters[, 'c3_assimilation_msg'], 'Cc must be >= 0')
    expect_equal(fit_res_bad$parameters[, 'c3_variable_j_msg'], 'Ci must be >= 0')
    expect_true(all(is.na(fit_res_bad$fits[, 'A_fit'])))
    expect_true(all(is.na(fit_res_bad$fits[, 'gmc'])))
    expect_true(all(is.na(fit_res_bad$parameters[, c('Vcmax_at_25', 'J_at_25', 'Rd_at_25', 'TPU', 'tau')])))
})

test_that('fit results have not changed', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- fit_c3_variable_j(one_curve, Ca_atmospheric = 420)

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25', 'J_at_25', 'Rd_at_25', 'TPU', 'tau')]),
        c(240.5297, 253.8052, 1.8880, 31.3243, 0.4048),
        tolerance = 1e-6
    )
})
