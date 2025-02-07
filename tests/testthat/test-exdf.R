test_that('an empty exdf can be created and printed', {
    empty_exdf <- expect_silent(
        exdf()
    )

    expect_equal(
        length(PhotoGEA:::fancy_column_names(empty_exdf)),
        ncol(empty_exdf)
    )
})

test_that('an exdf with one row can be created directly or by subsetting', {
    # Direct creation
    v1 <- expect_silent(
        exdf(data.frame(A = seq_len(5)))
    )

    # Taking a subset of an exdf
    big_exdf <- exdf(
        data.frame(
            A = seq_len(5),
            B = seq_len(5) * 2
        )
    )

    v2 <- expect_silent(
        big_exdf[, 'A', TRUE]
    )

    expect_equal(v1, v2)
})

test_that('exdf units and categories cannot have multiple rows', {
    expect_error(
        exdf(
            data.frame(
                A = seq_len(5),
                B = seq_len(5) * 2
            ),
            units = data.frame(
                A = c('A1 units', 'A2 units'),
                B = c('B1 units', 'B2 units')
            ),
            categories = data.frame(
                A = c('A1 category', 'A2 category'),
                B = c('B1 category', 'B2 category')
            )
        ),
        '`units` must have exactly one row\n  `categories` must have exactly one row'
    )
})

test_that('units and categories cannot have duplicated column names', {
    bad_units <- data.frame(
        A = 'A1 units',
        A = 'A2 units'
    )
    colnames(bad_units) <- c('A', 'A')

    bad_categories <- data.frame(
        B = 'B1 category',
        B = 'B2 category'
    )
    colnames(bad_categories) <- c('B', 'B')

    expect_error(
        exdf(
            data.frame(
                A = seq_len(5),
                B = seq_len(5) * 2
            ),
            units = bad_units,
            categories = bad_categories
        ),
        'All columns of `units` must have unique names, but the following names are duplicated: A\n  All columns of `categories` must have unique names, but the following names are duplicated: B'
    )
})

test_that('Not all units and categories are required for exdf creation', {
    test_exdf <- expect_silent(
        exdf(
            data.frame(
                A = seq_len(5),
                B = seq_len(5) * 2,
                C = seq_len(5) * 3 - 6
            ),
            units = data.frame(
                A = 'A units',
                stringsAsFactors = FALSE
            ),
            categories = data.frame(
                B = 'B category',
                stringsAsFactors = FALSE
            )
        )
    )

    expect_true(is.exdf(test_exdf, consistency_check = TRUE))
})
