# Set up some inputs
inputs <- exdf(data.frame(
  PCm = seq(1, 601, by = 60),
  Qin = 1800,
  Tleaf = 30,
  total_pressure = 1,
  oxygen = 21
))

inputs <- document_variables(
  inputs,
  c('', 'PCm',            'microbar'),
  c('', 'Qin',            'micromol m^(-2) s^(-1)'),
  c('', 'Tleaf',          'degrees C'),
  c('', 'total_pressure', 'bar'),
  c('', 'oxygen',         'percent')
)

inputs <- calculate_temperature_response(inputs, c4_temperature_param_vc, 'Tleaf')

inputs2 <- set_variable(inputs, 'alpha_psii', 'dimensionless', '', 0.5)

inputs3 <- inputs
inputs3[, 'PCm'] <- inputs[, 'PCm'] - 50

test_that('c4 assimilation works for numeric values of flexible inputs', {
    expect_silent(
        calculate_c4_assimilation(inputs, 0, 0.003, 400, 1, 0.5, 80, 120, 400)
    )

    res <- expect_silent(
        calculate_c4_assimilation(inputs2, 0, 0.003, 400, 1, 0.5, 80, 120, 400)
    )

    # The alpha_psii value in inputs2 should be overwritten by the supplied
    # numeric value
    expect_equal(res[1, 'alpha_psii'], 0)
})

test_that('c4 assimilation works for non-numeric values of flexible inputs', {
    expect_error(
        calculate_c4_assimilation(inputs, '', 0.003, 400, 1, 0.5, 80, 120, 400),
        'The following required columns are not present: alpha_psii'
    )

    res <- expect_silent(
        calculate_c4_assimilation(inputs2, '', 0.003, 400, 1, 0.5, 80, 120, 400)
    )

    # The alpha_psii value should be as in inputs2
    expect_equal(res[1, 'alpha_psii'], 0.5)
})

test_that('fitting parameter limits can be bypassed', {
    expect_error(
        calculate_c4_assimilation(inputs, 0, -0.003, -400, 1, 0.5, 80, 120, 400, hard_constraints = 2),
        'gbs must be >= 0. Jmax_at_25 must be >= 0',
        fixed = TRUE
    )

    expect_silent(
        calculate_c4_assimilation(inputs, 0, -0.003, -400, 1, 0.5, 80, 120, 400, hard_constraints = 1)
    )

    expect_error(
        calculate_c4_assimilation(inputs3, 0, -0.003, -400, 1, 0.5, 80, 120, 400, hard_constraints = 1),
        'PCm must be >= 0'
    )
})

test_that('PCm limits can be bypassed', {
    expect_error(
        calculate_c4_assimilation(inputs3, 0, 0.003, 400, 1, 0.5, 80, 120, 400, hard_constraints = 2),
        'PCm must be >= 0'
    )

    expect_error(
        calculate_c4_assimilation(inputs3, 0, 0.003, 400, 1, 0.5, 80, 120, 400, hard_constraints = 1),
        'PCm must be >= 0'
    )

    expect_silent(
        calculate_c4_assimilation(inputs3, 0, 0.003, 400, 1, 0.5, 80, 120, 400, hard_constraints = 0)
    )
})
