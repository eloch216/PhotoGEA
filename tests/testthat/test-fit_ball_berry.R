# Get test curves to use
source('one_curve_ball_berry.R')

# Choose test tolerance
TOLERANCE <- 1e-4

test_that('fit results have not changed', {
    fit_res <- fit_ball_berry(one_curve)

    expect_equal(
        as.numeric(fit_res$parameters[1, c('bb_intercept', 'bb_intercept_err', 'bb_slope', 'bb_slope_err')]),
        c(0.1033, 0.0308, 3.9951, 0.7595),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('r_squared', 'p_value')]),
        c(0.8469, 0.0033),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res$parameters[1, c('npts', 'nparam', 'dof')]),
        c(7, 2, 5)
    )
})

test_that('removing and excluding points produce the same fit results', {
    pts_to_remove <- c(3, 5)

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

    expect_equal(nrow(one_curve_remove), 5)
    expect_equal(nrow(one_curve_exclude), 7)

    fit_res_remove <- fit_ball_berry(one_curve_remove)
    fit_res_exclude <- fit_ball_berry(one_curve_exclude)

    # Check that results haven't changed
    expect_equal(
        as.numeric(fit_res_remove$parameters[1, c('bb_intercept', 'bb_intercept_err', 'bb_slope', 'bb_slope_err')]),
        c(0.1066, 0.0405, 4.0265, 0.9909),
        tolerance = TOLERANCE
    )

    expect_equal(
        as.numeric(fit_res_remove$parameters[1, c('npts', 'nparam', 'dof')]),
        c(5, 2, 3)
    )

    expect_equal(
        as.numeric(fit_res_remove$parameters[1, c('RSS', 'RMSE')]),
        c(0.005213733, 0.032291588),
        tolerance = TOLERANCE
    )

    # Check that remove/exclude results are the same
    expect_equal(
        as.numeric(fit_res_remove$parameters[1, c('bb_intercept', 'bb_intercept_err', 'bb_slope', 'bb_slope_err')]),
        as.numeric(fit_res_exclude$parameters[1, c('bb_intercept', 'bb_intercept_err', 'bb_slope', 'bb_slope_err')]),
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
