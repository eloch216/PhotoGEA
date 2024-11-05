# Set up some inputs
inputs <- exdf(data.frame(
  Ci = seq(1, 601, by = 60),
  Qin = 1800
))

inputs <- document_variables(
  inputs,
  c('', 'Ci',  'micromol mol^(-1)'),
  c('', 'Qin', 'micromol m^(-2) s^(-1)')
)

inputs2 <- set_variable(inputs, 'c4_curvature', 'dimensionless', '', 0.5)

inputs3 <- inputs
inputs3[, 'Ci'] <- inputs[, 'Ci'] - 50

test_that('c4 assimilation works for numeric values of flexible inputs', {
    expect_silent(
        calculate_c4_assimilation_hyperbola(inputs, 0.8, 0.5, 1.0, 65)
    )

    res <- expect_silent(
        calculate_c4_assimilation_hyperbola(inputs2, 0.8, 0.5, 1.0, 65)
    )

    # The c4_curvature value in inputs2 should be overwritten by the supplied
    # numeric value
    expect_equal(res[1, 'c4_curvature'], 0.8)
})

test_that('c4 assimilation works for non-numeric values of flexible inputs', {
    expect_error(
        calculate_c4_assimilation_hyperbola(inputs, '', 0.5, 1.0, 65),
        'The following required columns are not present: c4_curvature'
    )

    res <- expect_silent(
        calculate_c4_assimilation_hyperbola(inputs2, '', 0.5, 1.0, 65)
    )

    # The c4_curvature value should be as in inputs2
    expect_equal(res[1, 'c4_curvature'], 0.5)
})

test_that('fitting parameter limits can be bypassed', {
    expect_error(
        calculate_c4_assimilation_hyperbola(inputs, -0.8, -0.5, 1.0, 65, hard_constraints = 2),
        'c4_curvature must be >= 0 and <= 1. c4_slope must be >= 0',
        fixed = TRUE
    )

    expect_silent(
        calculate_c4_assimilation_hyperbola(inputs, -0.8, -0.5, 1.0, 65, hard_constraints = 1)
    )

    expect_error(
        calculate_c4_assimilation_hyperbola(inputs3, -0.8, -0.5, 1.0, 65, hard_constraints = 1),
        'Ci must be >= 0'
    )
})

test_that('Ci limits can be bypassed', {
    expect_error(
        calculate_c4_assimilation_hyperbola(inputs3, 0.8, 0.5, 1.0, 65, hard_constraints = 2),
        'Ci must be >= 0'
    )

    expect_error(
        calculate_c4_assimilation_hyperbola(inputs3, 0.8, 0.5, 1.0, 65, hard_constraints = 1),
        'Ci must be >= 0'
    )

    expect_silent(
        calculate_c4_assimilation_hyperbola(inputs3, 0.8, 0.5, 1.0, 65, hard_constraints = 0)
    )
})
