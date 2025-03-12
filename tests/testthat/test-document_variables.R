test_that('Variables can be documented in exdf and data frame objects', {
    testing_exdf <- exdf(data.frame(A = c(3, 2, 7, 9), B = c(4, 5, 1, 8)))

    # Specify units and categories for the `A` and `B` columns, and add a new `C`
    # column.
    res_exdf <- document_variables(
        testing_exdf,
        c('cat1', 'A', 'm'), # The category of `A` is `cat1` and its units are `m`
        c('cat2', 'B', 's'), # The category of `B` is `cat2` and its units are `s`
        c('cat3', 'C', 'g')  # The category of `C` is `cat3` and its units are `g`
    )

    # Do the same but for a data frame
    res_df <- document_variables(
          testing_exdf$main_data,
          c('cat1', 'A', 'm'), # The category of `A` is `cat1` and its units are `m`
          c('cat2', 'B', 's'), # The category of `B` is `cat2` and its units are `s`
          c('cat3', 'C', 'g')  # The category of `C` is `cat3` and its units are `g`
    )

    expect_equal(
        res_df,
        res_exdf$main_data
    )
})
