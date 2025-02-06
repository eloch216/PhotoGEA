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
        "'units' must have exactly one row\n  'categories' must have exactly one row"
    )
})
