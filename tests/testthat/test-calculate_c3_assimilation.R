# Set up some inputs
inputs <- exdf(data.frame(
  Cc = seq(1, 601, by = 60),
  oxygen = 21,
  Tleaf = 30,
  total_pressure = 1
))

inputs <- document_variables(
  inputs,
  c('', 'Cc',             'micromol mol^(-1)'),
  c('', 'oxygen',         'percent'),
  c('', 'Tleaf',          'degrees C'),
  c('', 'total_pressure', 'bar')
)

inputs <- calculate_temperature_response(inputs, c3_temperature_param_sharkey, 'Tleaf')

inputs2 <- set_variable(inputs, 'alpha_g', 'dimensionless', '', 0.5)

inputs3 <- inputs
inputs3[, 'Cc'] <- inputs[, 'Cc'] - 50

test_that('c3 assimilation works for numeric values of flexible inputs', {
    expect_silent(
        calculate_c3_assimilation(inputs, 0, 0, 0, 0, '', 150, '', '', 1, 12, 120)
    )

    res <- expect_silent(
        calculate_c3_assimilation(inputs2, 0, 0, 0, 0, '', 150, '', '', 1, 12, 120)
    )

    # The alpha_g value in inputs2 should be overwritten by the supplied numeric
    # value
    expect_equal(res[1, 'alpha_g'], 0)
})

test_that('c3 assimilation works for non-numeric values of flexible inputs', {
    expect_error(
        calculate_c3_assimilation(inputs, '', 0, 0, 0, '', 150, '', '', 1, 12, 120),
        'The following required columns are not present: alpha_g'
    )

    res <- expect_silent(
        calculate_c3_assimilation(inputs2, '', 0, 0, 0, '', 150, '', '', 1, 12, 120)
    )

    # The alpha_g value should be as in inputs2
    expect_equal(res[1, 'alpha_g'], 0.5)
})

test_that('alpha restrictions are enforced', {
    # All alphas nonzero
    expect_error(
        calculate_c3_assimilation(inputs, 0.1, 0.1, 0.1, 0.1, '', 150, '', '', 1, 12, 120),
        'Cannot specify nonzero alpha_old and nonzero alpha_g / alpha_s / alpha_t',
        fixed = TRUE
    )

    # alpha_old and alpha_s nonzero
    expect_error(
        calculate_c3_assimilation(inputs, 0.0, 0.1, 0.1, 0.0, '', 150, '', '', 1, 12, 120),
        'Cannot specify nonzero alpha_old and nonzero alpha_g / alpha_s / alpha_t',
        fixed = TRUE
    )

    # alpha_g and alpha_old nonzero
    expect_error(
        calculate_c3_assimilation(inputs, 0.1, 0.1, 0.0, 0.0, '', 150, '', '', 1, 12, 120),
        'Cannot specify nonzero alpha_old and nonzero alpha_g / alpha_s / alpha_t',
        fixed = TRUE
    )

    # alphas_s too high for supplied value of alpha_g (only when hard
    # constraints are applied)
    expect_error(
        calculate_c3_assimilation(inputs, 0.5, 0.0, 0.5, 0.5, '', 150, '', '', 1, 12, 120, hard_constraints = 2),
        'alpha_g + 2 * alpha_t + 4 * alpha_s / 3 must be <= 1',
        fixed = TRUE
    )

    expect_silent(
        calculate_c3_assimilation(inputs, 0.5, 0.0, 0.5, 0.5, '', 150, '', '', 1, 12, 120)
    )

    # alpha_g and alpha_s nonzero and Wj_coef_C not 4
    expect_error(
        calculate_c3_assimilation(inputs, 0.5, 0.0, 0.1, 0.0, '', 150, '', '', 1, 12, 120, Wj_coef_C = 4.5),
        'Wj_coef_C must be 4 and Wj_coef_Gamma_star must be 8 when alpha_g / alpha_s / alpha_t are nonzero',
        fixed = TRUE
    )


    # alpha_g and alpha_s nonzero and Wj_coef_Gamma_star not 8
    expect_error(
        calculate_c3_assimilation(inputs, 0.5, 0.0, 0.1, 0.0, '', 150, '', '', 1, 12, 120, Wj_coef_Gamma_star = 10.5),
        'Wj_coef_C must be 4 and Wj_coef_Gamma_star must be 8 when alpha_g / alpha_s / alpha_t are nonzero',
        fixed = TRUE
    )

    # alpha_g and alpha_s nonzero, Wj_coef_C not 4, and Wj_coef_Gamma_star not 8
    expect_error(
        calculate_c3_assimilation(inputs, 0.5, 0.0, 0.1, 0.0, '', 150, '', '', 1, 12, 120, Wj_coef_C = 4.5, Wj_coef_Gamma_star = 10.5),
        'Wj_coef_C must be 4 and Wj_coef_Gamma_star must be 8 when alpha_g / alpha_s / alpha_t are nonzero',
        fixed = TRUE
    )
})

test_that('fitting parameter limits can be bypassed', {
    expect_error(
        calculate_c3_assimilation(inputs, -1, -1, -1, 0, '', 150, '', '', 1, 12, 120, hard_constraints = 2),
        'alpha_g must be >= 0 and <= 1. alpha_old must be >= 0 and <= 1. alpha_s must be >= 0 and <= 1',
        fixed = TRUE
    )

    expect_silent(
        calculate_c3_assimilation(inputs, -1, -1, -1, 0, '', -150, -400, -275, -1, -12, -120, hard_constraints = 1)
    )

    expect_error(
        calculate_c3_assimilation(inputs3, -1, -1, -1, 0, '', 150, '', '', 1, 12, 120, hard_constraints = 1),
        'Cc must be >= 0'
    )
})

test_that('Cc limits can be bypassed', {
    expect_error(
        calculate_c3_assimilation(inputs3, 0, 0, 0, 0, '', 150, '', '', 1, 12, 120, hard_constraints = 2),
        'Cc must be >= 0'
    )

    expect_error(
        calculate_c3_assimilation(inputs3, 0, 0, 0, 0, '', 150, '', '', 1, 12, 120, hard_constraints = 1),
        'Cc must be >= 0'
    )

    expect_silent(
        calculate_c3_assimilation(inputs3, 0, 0, 0, 0, '', 150, '', '', 1, 12, 120, hard_constraints = 0)
    )
})

test_that('c3 assimilation works for exdf and data frame objects', {
    res_exdf <- expect_silent(
        calculate_c3_assimilation(inputs, 0, 0, 0, 0, '', 150, '', '', 1, 12, 120)
    )

    res_df <- expect_silent(
        calculate_c3_assimilation(inputs$main_data, 0, 0, 0, 0, '', 150, '', '', 1, 12, 120)
    )

    expect_equal(
        as.numeric(res_df[, 'An']),
        as.numeric(res_exdf[, 'An'])
    )
})
