# Get test curves to use
source('one_curve_c3_aci.R')

test_that('fit failures are handled properly', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res_bad <- expect_silent(
        fit_c3_variable_j(one_curve_bad, Ca_atmospheric = 420)
    )

    expect_equal(unique(fit_res_bad$fits[, 'c3_assimilation_msg']), 'Cc must be >= 0')
    expect_equal(unique(fit_res_bad$fits[, 'c3_variable_j_msg']), 'Ci must be >= 0')
    expect_equal(fit_res_bad$parameters[, 'c3_assimilation_msg'], 'Cc must be >= 0')
    expect_equal(fit_res_bad$parameters[, 'c3_variable_j_msg'], 'Ci must be >= 0')
    expect_true(all(is.na(fit_res_bad$fits[, 'A_fit'])))
    expect_true(all(is.na(fit_res_bad$fits[, 'gmc'])))
    expect_true(all(is.na(fit_res_bad$parameters[, c('Vcmax_at_25', 'J_at_25', 'Rd_at_25', 'Tp', 'tau')])))
})

test_that('fit results have not changed', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- fit_c3_variable_j(one_curve, Ca_atmospheric = 420)

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25', 'J_at_25', 'Rd_at_25', 'Tp', 'tau')]),
        c(240.5297, 253.8052, 1.8880, 31.3243, 0.4048),
        tolerance = 1e-6
    )
})
