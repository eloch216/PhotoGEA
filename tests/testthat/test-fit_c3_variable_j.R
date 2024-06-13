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
            OPTIM_FUN = optimizer_deoptim(200),
            calculate_confidence_intervals = TRUE,
            remove_unreliable_param = TRUE
        )
    )

    expect_equal(unique(fit_res_bad$fits[, 'c3_assimilation_msg']), 'Cc must be >= 0. Rd_at_25 must be >= 0. Vcmax_at_25 must be >= 0')
    expect_equal(unique(fit_res_bad$fits[, 'c3_variable_j_msg']), 'Ci must be >= 0. Rd_at_25 must be >= 0')
    expect_equal(fit_res_bad$parameters[, 'c3_assimilation_msg'], 'Cc must be >= 0. Rd_at_25 must be >= 0. Vcmax_at_25 must be >= 0')
    expect_equal(fit_res_bad$parameters[, 'c3_variable_j_msg'], 'Ci must be >= 0. Rd_at_25 must be >= 0')
    expect_true(all(is.na(fit_res_bad$fits[, c('A_fit', 'Ac', 'Aj', 'Ap', 'gmc', 'Cc')])))
    expect_true(all(is.na(fit_res_bad$fits_interpolated[, c('An', 'Ac', 'Aj', 'Ap', 'gmc', 'Cc')])))
    expect_true(all(is.na(fit_res_bad$parameters[, c('Vcmax_at_25', 'J_at_25', 'Rd_at_25', 'Tp', 'tau', 'AIC')])))
    expect_true(all(is.na(fit_res_bad$parameters[, c('Vcmax_at_25_upper', 'J_at_25_upper', 'Rd_at_25_upper', 'Tp', 'tau_upper')])))
})

test_that('Ci and Cc limits can be bypassed', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- expect_silent(
        fit_c3_variable_j(
            one_curve_bad,
            Ca_atmospheric = 420,
            OPTIM_FUN = optimizer_deoptim(200),
            hard_constraints = 0,
            calculate_confidence_intervals = TRUE,
            remove_unreliable_param = TRUE
        )
    )

    expect_equal(unique(fit_res$fits[, 'c3_assimilation_msg']), '')
    expect_equal(unique(fit_res$fits[, 'c3_variable_j_msg']), '')
    expect_equal(fit_res$parameters[, 'c3_assimilation_msg'], '')
    expect_equal(fit_res$parameters[, 'c3_variable_j_msg'], '')
    expect_true(all(!is.na(fit_res$fits[, c('A_fit', 'Ac', 'Aj', 'gmc', 'Cc')])))
})

test_that('fit results have not changed (no alpha)', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- fit_c3_variable_j(
        one_curve,
        Ca_atmospheric = 420,
        fit_options = list(alpha_old = 0, alpha_g = 0, alpha_s = 0),
        OPTIM_FUN = optimizer_deoptim(200),
        calculate_confidence_intervals = TRUE,
        remove_unreliable_param = TRUE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25', 'J_at_25', 'Rd_at_25', 'tau', 'Tp', 'AIC')]),
        c(240.718, 254.101, 1.885, 0.405, NA, 38.416),
        tolerance = 1e-5
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25_upper', 'J_at_25_upper', 'Rd_at_25_upper', 'tau_upper', 'Tp_upper')]),
        c(247.455, 256.611, 1.892, 0.409, Inf),
        tolerance = 1e-5
    )
})

test_that('fit results have not changed (alpha_old)', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- fit_c3_variable_j(
        one_curve,
        Ca_atmospheric = 420,
        fit_options = list(alpha_old = 'fit', alpha_g = 0, alpha_s = 0),
        OPTIM_FUN = optimizer_deoptim(200),
        calculate_confidence_intervals = TRUE,
        remove_unreliable_param = TRUE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25', 'J_at_25', 'Rd_at_25', 'tau', 'Tp', 'AIC')]),
        c(240.8281, 253.9907, 1.8771, 0.4049, NA, 40.4394),
        tolerance = 1e-5
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25_upper', 'J_at_25_upper', 'Rd_at_25_upper', 'tau_upper', 'Tp_upper')]),
        c(247.2329, 256.4144, 1.8906, 0.4088, Inf),
        tolerance = 1e-5
    )
})

test_that('fit results have not changed (alpha_g and alpha_s)', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- fit_c3_variable_j(
        one_curve,
        Ca_atmospheric = 420,
        fit_options = list(alpha_old = 0, alpha_g = 'fit', alpha_s = 'fit'),
        OPTIM_FUN = optimizer_deoptim(200),
        calculate_confidence_intervals = TRUE,
        remove_unreliable_param = TRUE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25', 'J_at_25', 'Rd_at_25', 'tau', 'Tp', 'AIC')]),
        c(216.444, 268.318, 2.360, 0.428, NA, 43.524),
        tolerance = 1e-5
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25_upper', 'J_at_25_upper', 'Rd_at_25_upper', 'tau_upper', 'Tp_upper')]),
        c(223.693, 271.054, 2.457, 0.432, Inf),
        tolerance = 1e-5
    )
})
