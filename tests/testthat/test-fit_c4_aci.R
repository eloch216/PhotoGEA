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
    expect_true(all(is.na(fit_res_bad$parameters[, c('Vcmax_at_25', 'Vpmax_at_25', 'Rd_at_25')])))
    expect_true(all(is.na(fit_res_bad$parameters[, c('Vcmax_at_25_upper', 'Vpmax_at_25_upper', 'Rd_at_25_upper')])))
})

test_that('fit results have not changed', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- fit_c4_aci(
        one_curve,
        Ca_atmospheric = 420,
        calculate_confidence_intervals = TRUE,
        remove_unreliable_param = TRUE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25', 'Vpmax_at_25', 'Rd_at_25')]),
        c(3.630116e+01, 1.804791e+02, 1.069116e-08),
        tolerance = 1e-5
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25_upper', 'Vpmax_at_25_upper', 'Rd_at_25_upper')]),
        c(38.434695, 214.046523, 1.568026),
        tolerance = 1e-5
    )
})
