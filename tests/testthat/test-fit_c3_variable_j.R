# Get test curves to use
source('one_curve_c3_aci.R')

test_that('fit failures are handled properly', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res_bad <- expect_silent(
        fit_c3_variable_j(
            one_curve_bad,
            Ca_atmospheric = 420,
            calculate_confidence_intervals = TRUE,
            remove_unreliable_param = TRUE
        )
    )

    expect_equal(unique(fit_res_bad$fits[, 'c3_assimilation_msg']), 'Cc must be >= 0')
    expect_equal(unique(fit_res_bad$fits[, 'c3_variable_j_msg']), 'Ci must be >= 0')
    expect_equal(fit_res_bad$parameters[, 'c3_assimilation_msg'], 'Cc must be >= 0')
    expect_equal(fit_res_bad$parameters[, 'c3_variable_j_msg'], 'Ci must be >= 0')
    expect_true(all(is.na(fit_res_bad$fits[, 'A_fit'])))
    expect_true(all(is.na(fit_res_bad$fits[, 'gmc'])))
    expect_true(all(is.na(fit_res_bad$parameters[, c('Vcmax_at_25', 'J_at_25', 'Rd_at_25', 'Tp', 'tau')])))
    expect_true(all(is.na(fit_res_bad$parameters[, c('Vcmax_at_25_upper', 'J_at_25_upper', 'Rd_at_25_upper', 'Tp', 'tau_upper')])))
})

test_that('fit results have not changed', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- fit_c3_variable_j(
        one_curve,
        Ca_atmospheric = 420,
        calculate_confidence_intervals = TRUE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25', 'J_at_25', 'Rd_at_25', 'tau', 'Tp')]),
        c(242.4517893, 255.3335038, 1.8934618, 0.4072432, 31.7388206),
        tolerance = 1e-5
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25_upper', 'J_at_25_upper', 'Rd_at_25_upper', 'tau_upper', 'Tp_upper')]),
        c(249.2159655, 257.8865971, 1.9041862, 0.4108754, Inf),
        tolerance = 1e-5
    )
})
