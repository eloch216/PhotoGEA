# Choose test tolerance
TOLERANCE <- 1e-4

test_that('all stats are calculated properly', {
    # Generate some residuals
    residuals <- seq_len(10)

    # Calculate residual stats as if these values had units of `kg` and were
    # related to a model with 3 free parameters
    result <- residual_stats(residuals, 'kg', 3)

    expect_equal(
        as.numeric(result[1, c('npts', 'nparam', 'dof', 'RSS', 'MSE', 'RMSE', 'RSE', 'AIC')]),
        c(10, 3, 7, 385, 38.5, 6.204837, 7.416198, 72.885353),
        tolerance = TOLERANCE
    )
})
