example_exdf <- exdf(
    data.frame(
        Cc             = seq(10, 400, length.out = 3),
        Gamma_star     = 40,
        gmc            = 0.5,
        gsc            = 0.8,
        J_tl           = 200,
        Kc             = 270,
        Ko             = 160,
        oxygen         = 21,
        total_pressure = 1,
        Vcmax_tl       = 150

    ),
    units = data.frame(
        Cc             = 'micromol mol^(-1)',
        Gamma_star     = 'micromol mol^(-1)',
        gmc            = 'mol m^(-2) s^(-1) bar^(-1)',
        gsc            = 'mol m^(-2) s^(-1)',
        J_tl           = 'micromol m^(-2) s^(-1)',
        Kc             = 'micromol mol^(-1)',
        Ko             = 'mmol mol^(-1)',
        oxygen         = 'percent',
        total_pressure = 'bar',
        Vcmax_tl       = 'micromol m^(-2) s^(-1)',
        stringsAsFactors = FALSE
    )
)

test_that('j limitations are not calculated by default', {
    limit_res <- calculate_c3_limitations_grassi(example_exdf)
    expect_true(is.exdf(limit_res))
    expect_true('dAdC_rubisco' %in% colnames(limit_res))
    expect_true('ls_rubisco_grassi' %in% colnames(limit_res))
    expect_true('lm_rubisco_grassi' %in% colnames(limit_res))
    expect_true('lb_rubisco_grassi' %in% colnames(limit_res))
    expect_false('dAdC_j' %in% colnames(limit_res))
    expect_false('ls_j_grassi' %in% colnames(limit_res))
    expect_false('lm_j_grassi' %in% colnames(limit_res))
    expect_false('lb_j_grassi' %in% colnames(limit_res))
})

test_that('limitations add to 1', {
    limit_res <- calculate_c3_limitations_grassi(example_exdf, j_column_name = 'J_tl')

    limit_res[, 'l_rubisco'] <- limit_res[, 'ls_rubisco_grassi'] + limit_res[, 'lm_rubisco_grassi'] + limit_res[, 'lb_rubisco_grassi']
    limit_res[, 'l_j'] <- limit_res[, 'ls_j_grassi'] + limit_res[, 'lm_j_grassi'] + limit_res[, 'lb_j_grassi']

    expect_equal(limit_res[, 'l_rubisco'], c(1, 1, 1))
    expect_equal(limit_res[, 'l_j'], c(1, 1, 1))
})

# Get test curve to use
source('one_curve_c3_aci.R')

# Specify mesophyll conductance
one_curve <- set_variable(
  one_curve,
  'gmc', 'mol m^(-2) s^(-1) bar^(-1)', value = 1.0
)

one_curve_bad <- set_variable(
  one_curve_bad,
  'gmc', 'mol m^(-2) s^(-1) bar^(-1)', value = 1.0
)

# Calculate Cc
one_curve <- apply_gm(one_curve)
one_curve_bad <- apply_gm(one_curve_bad)

# Calculate additional gas properties
one_curve <- calculate_gas_properties(one_curve)
one_curve_bad <- calculate_gas_properties(one_curve_bad)


test_that('fit failures are handled properly', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res_bad <- fit_c3_aci(one_curve_bad, Ca_atmospheric = 420)

    limit_res_bad <- expect_no_error(
        calculate_c3_limitations_grassi(fit_res_bad$fits)
    )

    expect_true(all(is.na(limit_res_bad[, 'lm_rubisco_grassi'])))
    expect_true(all(is.na(limit_res_bad[, 'dAdC_rubisco'])))
})

test_that('fit results have not changed', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- fit_c3_aci(one_curve, Ca_atmospheric = 420)

    limit_res <- expect_silent(
        calculate_c3_limitations_grassi(fit_res$fits)
    )

    expect_equal(
        as.numeric(limit_res[1, c('dAdC_rubisco', 'ls_rubisco_grassi', 'lm_rubisco_grassi', 'lb_rubisco_grassi')]),
        c(0.39900, 0.43480, 0.16148, 0.40372),
        tolerance = 1e-5
    )
})
