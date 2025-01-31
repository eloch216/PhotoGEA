test_that('Residuals are calculated for exdf and data frame objects', {
    testing_df <- data.frame(
        A = c(1, 2, 3, 4, 5),
        A_fit = c(1.1, 2.2, 3.3, 4.4, 5.5)
    )

    testing_exdf <- exdf(testing_df)

    resid_df <- expect_silent(
        PhotoGEA:::calculate_residuals(testing_df, 'A')
    )

    resid_exdf <- expect_silent(
        PhotoGEA:::calculate_residuals(testing_exdf, 'A')
    )

    expect_equal(
        as.numeric(resid_df[, 'A_residuals']),
        as.numeric(resid_exdf[, 'A_residuals'])
    )
})
