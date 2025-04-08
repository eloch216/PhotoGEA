simple_exdf <- exdf(
  data.frame(A = c(3, 2, 7, 9), B = c(4, 5, 1, 8)),
  data.frame(A = 'm', B = 's'),
  data.frame(A = 'Cat1', B = 'Cat2')
)

test_that('column names must match', {
    simple_exdf_diff_colname <- simple_exdf
    simple_exdf_diff_colname[, 'C'] <- 3

    expect_message(
        expect_error(
            rbind(simple_exdf, simple_exdf_diff_colname),
            'exdf objects must all have the same column names when using rbind'
        ),
        'colnames from first exdf object:\nA, B\ncolnames from current exdf object:\nA, B, C'
    )
})

test_that('categories must match', {
    simple_exdf_diff_category <- simple_exdf
    simple_exdf_diff_category$categories$A = 'Cat3'

    expect_message(
        expect_error(
            rbind(simple_exdf, simple_exdf_diff_category),
            'exdf objects must all have the same categories when using rbind'
        ),
        'categories from first exdf object:\nA = Cat1, B = Cat2\ncategories from current exdf object:\nA = Cat3, B = Cat2'
    )
})

test_that('units must match', {
    simple_exdf_diff_units <- simple_exdf
    simple_exdf_diff_units$units$A = 'kg'

    expect_message(
        expect_error(
            rbind(simple_exdf, simple_exdf_diff_units),
            'exdf objects must all have the same units when using rbind'
        ),
        'units from first exdf object:\nA = m, B = s\nunits from current exdf object:\nA = kg, B = s'
    )
})
