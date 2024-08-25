# Define a data frame where there is a `y_val` outlier when we split the entries
# by `x_val`
test_df <- data.frame(
    y_val = c(1,   2,   NA,  NA,  2,   100, 6,   NA,  7,   8,   NA,   10),
    x_val = c('A', 'A', 'A', 'A', 'A', 'A', 'B', 'B', 'B', 'B', 'B',  'B')
)

# Define a data frame where the outlier has been removed
test_df_removed <- data.frame(
    y_val = c(1,   2,   NA,  NA,  2,   6,   NA,  7,   8,   NA,   10),
    x_val = c('A', 'A', 'A', 'A', 'A', 'B', 'B', 'B', 'B', 'B',  'B')
)

# Define a data frame where the outlier has been replaced by NA
test_df_excluded <- data.frame(
    y_val = c(1,   2,   NA,  NA,  2,   NA,  6,   NA,  7,   8,   NA,   10),
    x_val = c('A', 'A', 'A', 'A', 'A', 'A', 'B', 'B', 'B', 'B', 'B',  'B')
)

test_that('outliers are removed and NAs are retained (data frame)', {
    # Exclude outliers
    res <- exclude_outliers(test_df, 'y_val', test_df$x_val, method = 'remove')

    # Ensure row names are the same because we don't care about them
    expected <- test_df_removed
    rownames(expected) <- NULL
    rownames(res) <- NULL

    expect_true(is.data.frame(res))
    expect_identical(expected, res)
})

test_that('outliers are removed and NAs are retained (exdf)', {
    # Exclude outliers
    res <- exclude_outliers(exdf(test_df), 'y_val', test_df$x_val, method = 'remove')

    # Ensure row names are the same because we don't care about them
    expected <- exdf(test_df_removed)
    rownames(expected$main_data) <- NULL
    rownames(res$main_data) <- NULL

    expect_true(is.exdf(res))
    expect_identical(expected, res)
})

test_that('outliers are excluded and NAs are retained (data frame)', {
    # Exclude outliers
    res <- exclude_outliers(test_df, 'y_val', test_df$x_val, method = 'exclude')

    # Ensure row names are the same because we don't care about them
    expected <- test_df_excluded
    rownames(expected) <- NULL
    rownames(res) <- NULL

    expect_true(is.data.frame(res))
    expect_identical(expected, res)
})

test_that('outliers are excluded and NAs are retained (exdf)', {
    # Exclude outliers
    res <- exclude_outliers(exdf(test_df), 'y_val', test_df$x_val, method = 'exclude')

    # Ensure row names are the same because we don't care about them
    expected <- exdf(test_df_excluded)
    rownames(expected$main_data) <- NULL
    rownames(res$main_data) <- NULL

    expect_true(is.exdf(res))
    expect_identical(expected, res)
})

