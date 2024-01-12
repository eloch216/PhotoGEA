# Set up some inputs
inputs <- exdf(data.frame(
  Cc = seq(1, 601, by = 60),
  Tleaf = 30,
  total_pressure = 1
))

inputs <- document_variables(
  inputs,
  c('', 'Cc',             'micromol mol^(-1)'),
  c('', 'Tleaf',          'degrees C'),
  c('', 'total_pressure', 'bar')
)

inputs <- calculate_arrhenius(inputs, c3_arrhenius_sharkey, 'Tleaf')

inputs2 <- set_variable(inputs, 'alpha', 'dimensionless', '', 0.5)

test_that('c3 assimilation works for numeric values of flexible inputs', {
    expect_silent(
        calculate_c3_assimilation(inputs, 0, 40, 150, 1, 12, 120)
    )

    res <- expect_silent(
        calculate_c3_assimilation(inputs2, 0, 40, 150, 1, 12, 120)
    )

    # The alpha value in inputs2 should be overwritten by the supplied numeric
    # value
    expect_equal(res[1, 'alpha'], 0)

    # The Gamma_star value in inputs2 should be overwritten by the supplied
    # numeric value
    expect_equal(res[1, 'Gamma_star'], 40)
})

test_that('c3 assimilation works for non-numeric values of flexible inputs', {
    expect_error(
        calculate_c3_assimilation(inputs, '', '', 150, 1, 12, 120),
        'The following columns are undefined: alpha'
    )

    res <- expect_silent(
        calculate_c3_assimilation(inputs2, '', '', 150, 1, 12, 120)
    )

    # The alpha value should be as in inputs2
    expect_equal(res[1, 'alpha'], 0.5)
})
