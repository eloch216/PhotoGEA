# Define a data frame where there is a `y_val` outlier when we split the entries
# by `x_val`
test_df <- data.frame(
    y_val = c(1,   2,   NA,  NA,  3,   100, 6,   NA,  7,   8,   NA,   10),
    x_val = c('A', 'A', 'A', 'A', 'A', 'A', 'B', 'B', 'B', 'B', 'B',  'B')
)

test_that('outliers can be removed when making barcharts', {
    keep_outliers_result <- barchart_with_errorbars(
        test_df[, 'y_val'],
        test_df[, 'x_val']
    )

    expect_equal(
        keep_outliers_result[['panel.args']][[1]][['y']],
        c(26.5, 7.75)
    )

    exclude_outliers_result <- barchart_with_errorbars(
        test_df[, 'y_val'],
        test_df[, 'x_val'],
        remove_outliers = TRUE
    )

    expect_equal(
        exclude_outliers_result[['panel.args']][[1]][['y']],
        c(2, 7.75)
    )
})
