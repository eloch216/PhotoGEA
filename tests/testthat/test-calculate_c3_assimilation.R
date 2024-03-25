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

inputs <- calculate_arrhenius(inputs, c3_arrhenius_sharkey, 'Tleaf')

inputs2 <- set_variable(inputs, 'alpha_g', 'dimensionless', '', 0.5)

test_that('c3 assimilation works for numeric values of flexible inputs', {
    expect_silent(
        calculate_c3_assimilation(inputs, 0, 40, 150, 1, 12, 120)
    )

    res <- expect_silent(
        calculate_c3_assimilation(inputs2, 0, 40, 150, 1, 12, 120)
    )

    # The alpha_g value in inputs2 should be overwritten by the supplied numeric
    # value
    expect_equal(res[1, 'alpha_g'], 0)

    # The Gamma_star value in inputs2 should be overwritten by the supplied
    # numeric value
    expect_equal(res[1, 'Gamma_star'], 40)
})

test_that('c3 assimilation works for non-numeric values of flexible inputs', {
    expect_error(
        calculate_c3_assimilation(inputs, '', '', 150, 1, 12, 120),
        'The following columns are undefined: alpha_g'
    )

    res <- expect_silent(
        calculate_c3_assimilation(inputs2, '', '', 150, 1, 12, 120)
    )

    # The alpha_g value should be as in inputs2
    expect_equal(res[1, 'alpha_g'], 0.5)
})
