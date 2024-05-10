# Get test curves to use
source('one_curve_c4_aci.R')

# Calculate PCm
one_curve <- apply_gm(one_curve, 'C4')
one_curve_bad <- apply_gm(one_curve_bad, 'C4')

test_that('fit failures are handled properly', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res_bad <- expect_silent(
        fit_c4_aci(
            one_curve_bad,
            Ca_atmospheric = 420,
            calculate_confidence_intervals = TRUE,
            remove_unreliable_param = TRUE
        )
    )

    expect_equal(unique(fit_res_bad$fits[, 'c4_assimilation_msg']), 'PCm must be >= 0')
    expect_equal(fit_res_bad$parameters[, 'c4_assimilation_msg'], 'PCm must be >= 0')
    expect_true(all(is.na(fit_res_bad$fits[, c('A_fit', 'Apr', 'Apc', 'Ar')])))
    expect_true(all(is.na(fit_res_bad$fits_interpolated[, c('An', 'Apr', 'Apc', 'Ar')])))
    expect_true(all(is.na(fit_res_bad$parameters[, c('Vcmax_at_25', 'Vpmax_at_25', 'Rd_at_25', 'AIC')])))
    expect_true(all(is.na(fit_res_bad$parameters[, c('Vcmax_at_25_upper', 'Vpmax_at_25_upper', 'Rd_at_25_upper')])))
})

test_that('some fit settings cause warnings', {
    fit_res_warn <- expect_warning(
        fit_c4_aci(
            one_curve_bad,
            Ca_atmospheric = 420,
            fit_options = list(Jmax_at_opt = 'fit')
        ),
        'The following C4 fit options have been selected: Jmax_at_opt: fit, Vcmax_at_25: fit, Vpr: 1000. It is not recommended to fit more than one of Jmax_at_opt, Vcmax_at_25, and Vpr.'
    )
})

test_that('fit results have not changed (Vcmax)', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- fit_c4_aci(
        one_curve,
        Ca_atmospheric = 420,
        fit_options = list(Vcmax_at_25 = 'fit', Vpr = 1000, Jmax_at_opt = 1000),
        calculate_confidence_intervals = TRUE,
        remove_unreliable_param = TRUE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25', 'Vpmax_at_25', 'Rd_at_25', 'AIC')]),
        c(3.630116e+01, 1.804791e+02, 1.069116e-08, 8.026640e+01),
        tolerance = 1e-5
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25_upper', 'Vpmax_at_25_upper', 'Rd_at_25_upper')]),
        c(38.434695, 214.046523, 1.568026),
        tolerance = 1e-5
    )
})

test_that('fit results have not changed (Vpr)', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- fit_c4_aci(
        one_curve,
        Ca_atmospheric = 420,
        fit_options = list(Vcmax_at_25 = 1000, Vpr = 'fit', Jmax_at_opt = 1000),
        calculate_confidence_intervals = TRUE,
        remove_unreliable_param = TRUE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vpr', 'Vpmax_at_25', 'Rd_at_25', 'AIC')]),
        c(NA, 202.95187, 18.47438, 103.06962),
        tolerance = 1e-5
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vpr_upper', 'Vpmax_at_25_upper', 'Rd_at_25_upper')]),
        c(Inf, 246.89935, 23.17948),
        tolerance = 1e-5
    )
})

test_that('fit results have not changed (Jmax)', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- fit_c4_aci(
        one_curve,
        Ca_atmospheric = 420,
        fit_options = list(Vcmax_at_25 = 1000, Vpr = 1000, Jmax_at_opt = 'fit'),
        calculate_confidence_intervals = TRUE,
        remove_unreliable_param = TRUE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Jmax_at_opt', 'Vpmax_at_25', 'Rd_at_25', 'AIC')]),
        c(5.215746e+02, 1.338467e+02, 1.475187e-08, 8.675720e+01),
        tolerance = 1e-5
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Jmax_at_opt_upper', 'Vpmax_at_25_upper', 'Rd_at_25_upper')]),
        c(573.15632, 157.30750, 2.24672),
        tolerance = 1e-5
    )
})
