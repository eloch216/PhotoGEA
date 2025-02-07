test_that('genotype (WT) identifiers are factorized properly', {
    ids <- c('4', 'WT', '2', 'Wt', '8')

    factorized_ids <- factorize_id_column(ids)

    expect_equal(
        levels(factorized_ids),
        c('WT', '2',  '4',  '8')
    )

    expect_equal(
        as.character(factorized_ids),
        c('4', 'WT', '2', 'WT', '8')
    )
})

test_that('genotype - replicate (WT) identifiers are factorized properly', {
    ids <- c('4 - 4', 'wT - 2', 'a - 2', 'WT - 1', '4 - 8', 'wt - 9')

    factorized_ids <- factorize_id_column(ids)

    expect_equal(
        levels(factorized_ids),
        c("WT - 1", "WT - 2", "WT - 9", "4 - 4", "4 - 8", "a - 2")
    )

    expect_equal(
        as.character(factorized_ids),
        c("4 - 4", "WT - 2", "a - 2", "WT - 1", "4 - 8", "WT - 9")
    )
})

test_that('table columns are factorized', {
    ids <- c('4 - 4', 'wT - 2', 'a - 2', 'WT - 1', '4 - 8', 'wt - 9')

    dat <- data.frame(replicate_id = ids, val = seq_along(ids))

    exdf_obj <- exdf(dat, units = data.frame(replicate_id = '', val = 'm / s'))

    dat <- expect_silent(
        factorize_id_column(dat, 'replicate_id')
    )

    exdf_obj <- expect_silent(
        factorize_id_column(exdf_obj, 'replicate_id')
    )

    expect_equal(
        levels(dat[, 'replicate_id']),
        levels(exdf_obj[, 'replicate_id'])
    )
})
