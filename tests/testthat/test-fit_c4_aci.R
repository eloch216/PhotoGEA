# Get test curves to use
source('one_curve_c4_aci.R')

# Load helping function
source('get_duplicated_colnames.R')

# Choose test tolerance
TOLERANCE <- 1e-4

test_that('fit failures are handled properly', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res_bad <- expect_silent(
        fit_c4_aci(
            one_curve_bad,
            Ca_atmospheric = 420,
            optim_fun = optimizer_nmkb(1e-7),
            hard_constraints = 2,
            calculate_confidence_intervals = TRUE,
            remove_unreliable_param = 2
        )
    )

    expect_equal(unique(fit_res_bad$fits[, 'c4_assimilation_msg']), 'PCm must be >= 0')
    expect_equal(fit_res_bad$parameters[, 'c4_assimilation_msg'], 'PCm must be >= 0')
    expect_true(all(is.na(fit_res_bad$fits[, c('A_fit', 'Apr', 'Apc', 'Ar')])))
    expect_true(all(is.na(fit_res_bad$fits_interpolated[, c('An', 'Apr', 'Apc', 'Ar')])))
    expect_true(all(is.na(fit_res_bad$parameters[, c('Vcmax_at_25', 'Vpmax_at_25', 'RL_at_25', 'AIC')])))
    expect_true(all(is.na(fit_res_bad$parameters[, c('Vcmax_at_25_upper', 'Vpmax_at_25_upper', 'RL_at_25_upper')])))
})

test_that('PCm limits can be bypassed', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- expect_silent(
        fit_c4_aci(
            one_curve_bad,
            Ca_atmospheric = 420,
            optim_fun = optimizer_nmkb(1e-7),
            hard_constraints = 0,
            calculate_confidence_intervals = TRUE,
            remove_unreliable_param = 2
        )
    )

    expect_equal(unique(fit_res$fits[, 'c4_assimilation_msg']), '')
    expect_equal(fit_res$parameters[, 'c4_assimilation_msg'], '')
    expect_true(all(!is.na(fit_res$fits[, c('A_fit')])))
})

test_that('fit results have not changed (Vcmax)', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- fit_c4_aci(
        one_curve,
        Ca_atmospheric = 420,
        fit_options = list(Vcmax_at_25 = 'fit', Vpr = 1000, J_at_25 = 1000),
        optim_fun = optimizer_nmkb(1e-7),
        hard_constraints = 2,
        calculate_confidence_intervals = TRUE,
        remove_unreliable_param = 2
    )

    expect_equal(
        get_duplicated_colnames(fit_res$fits),
        character(0)
    )

    expect_equal(
        get_duplicated_colnames(fit_res$parameters),
        character(0)
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25', 'Vpmax_at_25', 'RL_at_25', 'AIC', 'TleafCnd_avg')]),
        c(3.630116e+01, 1.804791e+02, 1.069116e-08, 8.226640e+01, 3.030823e+01),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25_upper', 'Vpmax_at_25_upper', 'RL_at_25_upper')]),
        c(38.434695, 214.046523, 1.568026),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('npts', 'nparam', 'dof')]),
        c(13, 3, 10)
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vpmax_trust', 'Vcmax_trust', 'Vpr_trust', 'J_trust')]),
        c(2, 2, 0, 0)
    )
})

test_that('fit results have not changed (Vpr)', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- fit_c4_aci(
        one_curve,
        Ca_atmospheric = 420,
        fit_options = list(Vcmax_at_25 = 1000, Vpr = 'fit', J_at_25 = 1000),
        optim_fun = optimizer_nmkb(1e-7),
        hard_constraints = 2,
        calculate_confidence_intervals = TRUE,
        remove_unreliable_param = 2
    )

    expect_equal(
        get_duplicated_colnames(fit_res$fits),
        character(0)
    )

    expect_equal(
        get_duplicated_colnames(fit_res$parameters),
        character(0)
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vpr', 'Vpmax_at_25', 'RL_at_25', 'AIC')]),
        c(58.1503, 133.8474, 0.0000, 88.3427),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vpr_upper', 'Vpmax_at_25_upper', 'RL_at_25_upper')]),
        c(62.43, 156.94, 2.76),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('npts', 'nparam', 'dof')]),
        c(13, 3, 10)
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vpmax_trust', 'Vcmax_trust', 'Vpr_trust', 'J_trust')]),
        c(2, 1, 2, 0)
    )
})

test_that('fit results have not changed (J)', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- fit_c4_aci(
        one_curve,
        Ca_atmospheric = 420,
        fit_options = list(Vcmax_at_25 = 1000, Vpr = 1000, J_at_25 = 'fit'),
        optim_fun = optimizer_nmkb(1e-7),
        hard_constraints = 2,
        calculate_confidence_intervals = TRUE,
        remove_unreliable_param = 2
    )

    fit_res$parameters <- calculate_temperature_response(
        fit_res$parameters,
        jmax_temperature_param_bernacchi,
        'TleafCnd_avg'
    )

    fit_res$parameters <- calculate_jmax(
        fit_res$parameters,
        0.6895,
        0.97875
    )

    expect_equal(
        get_duplicated_colnames(fit_res$fits),
        character(0)
    )

    expect_equal(
        get_duplicated_colnames(fit_res$parameters),
        character(0)
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('J_at_25', 'Vpmax_at_25', 'RL_at_25', 'AIC', 'Jmax_at_25')]),
        c(258.1464, 135.7058, 0.0000, 88.6061, 259.4098),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('J_at_25_upper', 'Vpmax_at_25_upper', 'RL_at_25_upper')]),
        c(275.66, 157.25, 2.35),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('npts', 'nparam', 'dof')]),
        c(13, 3, 10)
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vpmax_trust', 'Vcmax_trust', 'Vpr_trust', 'J_trust')]),
        c(2, 1, 0, 2)
    )
})

test_that('fit results have not changed (gmc with temperature dependence)', {
    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res <- fit_c4_aci(
        one_curve,
        Ca_atmospheric = 420,
        fit_options = list(gmc_at_25 = 'fit'),
        optim_fun = optimizer_nmkb(1e-7),
        hard_constraints = 2,
        calculate_confidence_intervals = TRUE,
        remove_unreliable_param = 2
    )

    expect_equal(
        get_duplicated_colnames(fit_res$fits),
        character(0)
    )

    expect_equal(
        get_duplicated_colnames(fit_res$parameters),
        character(0)
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25', 'Vpmax_at_25', 'RL_at_25', 'gmc_at_25', 'AIC')]),
        c(42.177755, 114.264158, 3.029476, 9.999997, 75.003441),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_at_25_upper', 'Vpmax_at_25_upper', 'RL_at_25_upper', 'gmc_at_25_upper')]),
        c(43.925521, 124.343968, 4.295661, Inf),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_tl_avg', 'Vpmax_tl_avg', 'RL_tl_avg', 'gmc_tl_avg')]),
        c(73.155628, 162.740193, 4.841121, 14.212310),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vcmax_tl_avg_lower', 'Vpmax_tl_avg_lower', 'RL_tl_avg_lower', 'gmc_tl_avg_lower')]),
        c(70.179640, 149.671916, 2.812345, 4.772084),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('npts', 'nparam', 'dof')]),
        c(13, 4, 9)
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('Vpmax_trust', 'Vcmax_trust', 'Vpr_trust', 'J_trust')]),
        c(2, 2, 0, 0)
    )
})

test_that('removing and excluding points produce the same fit results', {
    pts_to_remove <- c(3, 5, 13)

    one_curve_remove <- remove_points(
        one_curve,
        list(seq_num = pts_to_remove),
        method = 'remove'
    )

    one_curve_exclude <- remove_points(
        one_curve,
        list(seq_num = pts_to_remove),
        method = 'exclude'
    )

    expect_equal(nrow(one_curve_remove), 10)
    expect_equal(nrow(one_curve_exclude), 13)

    # Set a seed before fitting since there is randomness involved with the
    # default optimizer
    set.seed(1234)

    fit_res_remove <- fit_c4_aci(
        one_curve_remove,
        Ca_atmospheric = 420,
        optim_fun = optimizer_nmkb(1e-7)
    )

    set.seed(1234)

    fit_res_exclude <- fit_c4_aci(
        one_curve_exclude,
        Ca_atmospheric = 420,
        optim_fun = optimizer_nmkb(1e-7)
    )

    # Check that results haven't changed
    expect_equal(
        as.numeric(fit_res_remove$parameters[1, c('Vcmax_at_25', 'Vpmax_at_25', 'RL_at_25', 'AIC')]),
        c(37.01, 211.08, 0.70, 59.67),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res_remove$parameters[1, c('npts', 'nparam', 'dof')]),
        c(10, 3, 7)
    )

    expect_equal(
        as.numeric(fit_res_remove$parameters[1, c('RSS', 'RMSE')]),
        c(102.677, 3.204),
        tolerance = TOLERANCE
    )

    # Check that remove/exclude results are the same
    expect_equal(
        as.numeric(fit_res_remove$parameters[1, c('Vcmax_at_25', 'Vpmax_at_25', 'RL_at_25', 'AIC')]),
        as.numeric(fit_res_exclude$parameters[1, c('Vcmax_at_25', 'Vpmax_at_25', 'RL_at_25', 'AIC')]),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res_remove$parameters[1, c('npts', 'nparam', 'dof')]),
        as.numeric(fit_res_exclude$parameters[1, c('npts', 'nparam', 'dof')])
    )

    expect_equal(
        as.numeric(fit_res_remove$parameters[1, c('RSS', 'RMSE')]),
        as.numeric(fit_res_exclude$parameters[1, c('RSS', 'RMSE')]),
        tolerance = TOLERANCE
    )
})

