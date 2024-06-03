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
        calculate_c3_assimilation(inputs, 0, 0, 0, 40, 150, 1, 12, 120)
    )

    res <- expect_silent(
        calculate_c3_assimilation(inputs2, 0, 0, 0, 40, 150, 1, 12, 120)
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
        calculate_c3_assimilation(inputs, '', 0, 0, '', 150, 1, 12, 120),
        'The following columns are undefined: alpha_g'
    )

    res <- expect_silent(
        calculate_c3_assimilation(inputs2, '', 0, 0, '', 150, 1, 12, 120)
    )

    # The alpha_g value should be as in inputs2
    expect_equal(res[1, 'alpha_g'], 0.5)
})

test_that('alpha restrictions are enforced', {
    # All alphas nonzero
    expect_error(
        calculate_c3_assimilation(inputs, 0.1, 0.1, 0.1, '', 150, 1, 12, 120),
        'Cannot specify nonzero alpha_old and nonzero alpha_s or alpha_g'
    )

    # alpha_old and alpha_s nonzero
    expect_error(
        calculate_c3_assimilation(inputs, 0.0, 0.1, 0.1, '', 150, 1, 12, 120),
        'Cannot specify nonzero alpha_old and nonzero alpha_s or alpha_g'
    )

    # alpha_g and alpha_old nonzero
    expect_error(
        calculate_c3_assimilation(inputs, 0.1, 0.1, 0.0, '', 150, 1, 12, 120),
        'Cannot specify nonzero alpha_old and nonzero alpha_s or alpha_g'
    )

    # alphas_s too high for supplied value of alpha_g
    expect_error(
        calculate_c3_assimilation(inputs, 0.5, 0.0, 0.5, '', 150, 1, 12, 120),
        'alpha_s must be >= 0 and <= 0.75 * (1 - alpha_g)',
        fixed = TRUE
    )

    # alpha_g and alpha_s nonzero and atp_use not 4
    expect_error(
        calculate_c3_assimilation(inputs, 0.5, 0.0, 0.1, '', 150, 1, 12, 120, atp_use = 4.5),
        'atp_use must be 4 and nadph_use must be 8 when alpha_s or alpha_s are nonzero'
    )


    # alpha_g and alpha_s nonzero and nadph_use not 8
    expect_error(
        calculate_c3_assimilation(inputs, 0.5, 0.0, 0.1, '', 150, 1, 12, 120, nadph_use = 10.5),
        'atp_use must be 4 and nadph_use must be 8 when alpha_s or alpha_s are nonzero'
    )

    # alpha_g and alpha_s nonzero, atp_use not 4, and nadph_use not 8
    expect_error(
        calculate_c3_assimilation(inputs, 0.5, 0.0, 0.1, '', 150, 1, 12, 120, atp_use = 4.5, nadph_use = 10.5),
        'atp_use must be 4 and nadph_use must be 8 when alpha_s or alpha_s are nonzero'
    )
})
